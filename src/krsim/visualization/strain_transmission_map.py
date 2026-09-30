import os
from collections import defaultdict

import pandas as pd


def _build_viz_data(lineage, node_map):
    """
    Turn a lineage table (pre-simplify) + node ID map into the layout data
    the strain-transmission-map page renders: one entry per strain (generation,
    host slot, position within its host, recomb ID, post-simplify ID) plus
    the parent-child edges between them.

    Parameters
    ----------
    lineage : pd.DataFrame
        As written by output.lineage.save_lineage — columns gen,
        child_node_ID, child_ind_ID, recomb_ID, parent_node_ID, parent_ind_ID
        (plus left, right, which this function ignores).
    node_map : pd.DataFrame
        As written by output.node_id_map.save_node_id_map — columns
        pre_simplify_node_ID, post_simplify_node_ID.

    Returns
    -------
    dict
        JSON-serialisable payload for the page's embedded data script.
    """
    T = int(lineage["gen"].max()) + 1
    F = int(lineage["child_node_ID"].min())

    child_info = lineage.drop_duplicates("child_node_ID")[["child_node_ID", "child_ind_ID", "gen", "recomb_ID"]]
    node_to_ind = dict(zip(child_info["child_node_ID"], child_info["child_ind_ID"]))
    node_to_gen = dict(zip(child_info["child_node_ID"], child_info["gen"]))
    node_to_recomb = dict(zip(child_info["child_node_ID"], child_info["recomb_ID"]))

    # founders never appear as a child, so their host has to be recovered
    # from the parent side of gen (T-1) edges; any founder that was never
    # sampled as a parent at all has no recorded host and is dropped into
    # a fallback "unknown" lane
    founder_rows = lineage[lineage["parent_node_ID"] < F][["parent_node_ID", "parent_ind_ID"]].drop_duplicates()
    founder_known = dict(zip(founder_rows["parent_node_ID"], founder_rows["parent_ind_ID"]))
    never_used_founders = F - len(founder_known)
    ids_sorted = sorted(founder_known)
    for m in [i for i in range(F) if i not in founder_known]:
        before = max([i for i in ids_sorted if i < m], default=None)
        after = min([i for i in ids_sorted if i > m], default=None)
        bi, ai = founder_known.get(before), founder_known.get(after)
        founder_known[m] = bi if bi == ai and bi is not None else None

    for k, v in founder_known.items():
        node_to_ind[k] = v
        node_to_gen[k] = float(T)
        node_to_recomb[k] = None

    post_map = dict(zip(node_map["pre_simplify_node_ID"], node_map["post_simplify_node_ID"]))

    gens = sorted(set(node_to_gen.values()), reverse=True)
    gen_to_inds = {}
    for g in gens:
        inds = sorted({node_to_ind[n] for n in node_to_ind if node_to_gen[n] == g and node_to_ind[n] is not None})
        gen_to_inds[g] = inds
    Nh = max((len(v) for v in gen_to_inds.values()), default=0)

    cell_nodes = defaultdict(list)
    for n, g in node_to_gen.items():
        cell_nodes[(g, node_to_ind[n])].append(n)
    for k in cell_nodes:
        cell_nodes[k].sort()

    nodes_out = []
    for n in sorted(node_to_gen):
        g = node_to_gen[n]
        ind = node_to_ind[n]
        lane = gen_to_inds[g].index(ind) if ind is not None else None
        cell = cell_nodes[(g, ind)]
        slot = cell.index(n)
        slot_count = len(cell)
        post = post_map.get(n)
        post = None if pd.isna(post) else int(post)
        recomb = node_to_recomb.get(n)
        nodes_out.append({
            "id": int(n),
            "gen": float(g),
            "lane": lane,
            "slot": slot,
            "slotCount": slot_count,
            "recomb": None if pd.isna(recomb) else int(recomb),
            "post": post,
            "isFounder": n < F,
            "isSample": g == 0,
        })

    edges_out = [{"child": int(c), "parent": int(p)} for c, p in
                 lineage[["child_node_ID", "parent_node_ID"]].drop_duplicates().itertuples(index=False)]

    return {
        "T": T,
        "Nh": Nh,
        "neverUsedFounders": never_used_founders,
        "nodes": nodes_out,
        "edges": edges_out,
        "genOrder": gens,
    }


def _build_genome_data(ts, lineage, node_map, nodes):
    """
    Place every mutation of the mutated tree sequence on the pedigree, and list
    the trees along each chromosome, for the map's genome view.

    A mutation sits on a post-simplify node at a genome position and a time.
    Starting from that node's pre-simplify ID, we follow the pedigree upwards
    along the parent that passed on that genome position, until reaching the
    parent-child link whose generations bracket the mutation's time.

    Parameters
    ----------
    ts : tskit.TreeSequence
        Mutated tree sequence from `krsim mutate` (ts_mut_<seed>.trees).
    lineage : pd.DataFrame
        As written by output.lineage.save_lineage (with left, right).
    node_map : pd.DataFrame
        As written by output.node_id_map.save_node_id_map.
    nodes : list of dict
        The "nodes" entry of _build_viz_data's output (for each node's generation).

    Returns
    -------
    dict
        chroms, trees, segments, mutations and unplacedMutations entries for
        the page's embedded data script.
    """
    import numpy as np
    from krsim.sim.params import COMP_BREAKPOINTS, GENOME_LENGTH

    starts = [0] + [int(b) for b in COMP_BREAKPOINTS]
    ends = [int(b) for b in COMP_BREAKPOINTS] + [int(GENOME_LENGTH)]
    chroms = [{"chrom": i + 1, "start": s, "end": e} for i, (s, e) in enumerate(zip(starts, ends))]

    # a tree can span a chromosome boundary (independent assortment keeps the
    # same parent half the time), so it is listed once per chromosome it covers
    trees = []
    for tree in ts.trees():
        left, right = tree.interval
        for c in chroms:
            if left < c["end"] and right > c["start"]:
                trees.append({"index": tree.index, "chrom": c["chrom"],
                              "left": max(left, c["start"]), "right": min(right, c["end"])})

    seg_by_child = defaultdict(list)
    segments = []
    for c, p, left, right in lineage[["child_node_ID", "parent_node_ID", "left", "right"]].itertuples(index=False):
        seg_by_child[int(c)].append((float(left), float(right), int(p)))
        segments.append([int(c), int(p), float(left), float(right)])

    gen_of = {n["id"]: n["gen"] for n in nodes}
    mapped = node_map.dropna(subset=["post_simplify_node_ID"])
    pre_of_post = dict(zip(mapped["post_simplify_node_ID"].astype(int), mapped["pre_simplify_node_ID"].astype(int)))

    breakpoints = np.asarray(ts.breakpoints(as_array=True))
    site_pos = ts.tables.sites.position
    mutations = []
    unplaced = 0
    for m in ts.mutations():
        x = float(site_pos[m.site])
        t = float(m.time)
        cur = pre_of_post.get(int(m.node))
        placed = None
        while cur is not None and not np.isnan(t):
            parent = next((p for left, right, p in seg_by_child.get(cur, []) if left <= x < right), None)
            if parent is None:
                break
            g, gp = gen_of[cur], gen_of[parent]
            if g <= t < gp:
                placed = (cur, parent, (t - g) / (gp - g))
                break
            cur = parent
        if placed is None:
            unplaced += 1
            continue
        ci = int(np.searchsorted(ends, x, side="right"))
        mutations.append({
            "pos": round(x, 1),
            "chrom": ci + 1,
            "local": round(x - starts[ci], 1),
            "time": round(t, 3),
            "tree": int(np.searchsorted(breakpoints, x, side="right") - 1),
            "child": placed[0],
            "parent": placed[1],
            "frac": round(placed[2], 4),
        })

    return {
        "chroms": chroms,
        "trees": trees,
        "numTrees": ts.num_trees,
        "segments": segments,
        "mutations": mutations,
        "unplacedMutations": unplaced,
    }


_HTML_TEMPLATE = r"""<meta charset="utf-8">
<title>Strain transmission map — sim __SEED__</title>
<style>
  .viz-root {
    color-scheme: light;
    --surface:        #fcfcfb;
    --page:           #f9f9f7;
    --ink:            #0b0b0b;
    --ink-2:          #52514e;
    --ink-muted:      #898781;
    --gridline:       #e1e0d9;
    --baseline:       #c3c2b7;
    --border:         rgba(11,11,11,0.10);
    --lane-1: #2a78d6; --lane-2: #eb6834; --lane-3: #1baf7a; --lane-4: #eda100;
    --lane-5: #e87ba4; --lane-6: #008300; --lane-7: #4a3aa7; --lane-8: #e34948;
    --lane-overflow:  #898781;
    --card:           #ffffff;
    --edge:           #c9c8bf;
    --warn-bg:        #fdf3e4;
    --mut: #E23137; --seg-a: #DEE9F6; --seg-b: #C8E3F3; --warn-ink:       #7a4a00;
  }
  @media (prefers-color-scheme: dark) {
    :root:where(:not([data-theme="light"])) .viz-root {
      color-scheme: dark;
      --surface:        #1a1a19;
      --page:           #0d0d0d;
      --ink:            #ffffff;
      --ink-2:          #c3c2b7;
      --ink-muted:      #898781;
      --gridline:       #2c2c2a;
      --baseline:       #383835;
      --border:         rgba(255,255,255,0.10);
      --lane-1: #3987e5; --lane-2: #d95926; --lane-3: #199e70; --lane-4: #c98500;
      --lane-5: #d55181; --lane-6: #008300; --lane-7: #9085e9; --lane-8: #e66767;
      --lane-overflow:  #898781;
      --card:           #222221;
      --edge:           #3c3c39;
      --warn-bg:        #2c2410;
      --mut: #E23137; --seg-a: #22475F; --seg-b: #152F40; --warn-ink:       #e0b25c;
    }
  }
  :root[data-theme="dark"] .viz-root {
    color-scheme: dark;
    --surface: #1a1a19; --page: #0d0d0d; --ink: #ffffff; --ink-2: #c3c2b7; --ink-muted: #898781;
    --gridline: #2c2c2a; --baseline: #383835; --border: rgba(255,255,255,0.10);
    --lane-1: #3987e5; --lane-2: #d95926; --lane-3: #199e70; --lane-4: #c98500;
    --lane-5: #d55181; --lane-6: #008300; --lane-7: #9085e9; --lane-8: #e66767;
    --lane-overflow: #898781; --card: #222221; --edge: #3c3c39;
    --warn-bg: #2c2410; --mut: #E23137; --seg-a: #22475F; --seg-b: #152F40; --warn-ink: #e0b25c;
  }
  :root[data-theme="light"] .viz-root {
    color-scheme: light;
    --surface: #fcfcfb; --page: #f9f9f7; --ink: #0b0b0b; --ink-2: #52514e; --ink-muted: #898781;
    --gridline: #e1e0d9; --baseline: #c3c2b7; --border: rgba(11,11,11,0.10);
    --lane-1: #2a78d6; --lane-2: #eb6834; --lane-3: #1baf7a; --lane-4: #eda100;
    --lane-5: #e87ba4; --lane-6: #008300; --lane-7: #4a3aa7; --lane-8: #e34948;
    --lane-overflow: #898781; --card: #ffffff; --edge: #c9c8bf;
    --warn-bg: #fdf3e4; --mut: #E23137; --seg-a: #DEE9F6; --seg-b: #C8E3F3; --warn-ink: #7a4a00;
  }

  * { box-sizing: border-box; }

  .viz-root {
    background: var(--page);
    color: var(--ink);
    font-family: system-ui, -apple-system, "Segoe UI", sans-serif;
    min-height: 100vh;
    padding: 28px 20px 60px;
  }

  .mono {
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    font-variant-numeric: tabular-nums;
  }

  .wrap { max-width: 1180px; margin: 0 auto; }

  header.page { margin-bottom: 22px; }
  .eyebrow {
    font-size: 11px; font-weight: 600; letter-spacing: 0.08em; text-transform: uppercase;
    color: var(--ink-muted); margin: 0 0 6px;
  }
  h1 { font-size: 22px; font-weight: 650; margin: 0 0 8px; text-wrap: balance; }
  .dek { font-size: 14px; line-height: 1.55; color: var(--ink-2); max-width: 66ch; margin: 0; }
  .dek strong { color: var(--ink); font-weight: 600; }

  .meta-row { display: flex; flex-wrap: wrap; gap: 8px 10px; margin-top: 14px; }
  .chip {
    font-size: 12px; color: var(--ink-2); background: var(--card);
    border: 1px solid var(--border); border-radius: 6px; padding: 4px 9px;
  }
  .chip b { color: var(--ink); font-weight: 600; }

  .warn-banner {
    display: none;
    margin-top: 14px;
    font-size: 12.5px;
    line-height: 1.5;
    background: var(--warn-bg);
    color: var(--warn-ink);
    border-radius: 8px;
    padding: 10px 14px;
  }
  .warn-banner.show { display: block; }

  .panel { background: var(--surface); border: 1px solid var(--border); border-radius: 12px; }

  .toolbar {
    display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap;
    gap: 10px; padding: 14px 16px; border-bottom: 1px solid var(--border);
  }
  .legend { display: flex; align-items: center; flex-wrap: wrap; gap: 14px; font-size: 12px; color: var(--ink-2); }
  .legend-group { display: flex; align-items: center; gap: 7px; }
  .swatch { width: 11px; height: 11px; border-radius: 50%; flex: none; display: inline-block; }
  .swatch.founder { border: 1.5px solid var(--ink-muted); background: transparent; }
  .swatch.sample { border: 2px solid var(--ink); background: var(--lane-1); }

  #clearBtn {
    font: inherit; font-size: 12px; font-weight: 600; color: var(--ink); background: var(--card);
    border: 1px solid var(--border); border-radius: 6px; padding: 6px 12px; cursor: pointer;
  }
  #clearBtn:hover { border-color: var(--ink-muted); }
  #clearBtn[disabled] { opacity: 0.4; cursor: default; }

  .stage { display: flex; flex-direction: column; }

  .canvas-scroll { position: relative; overflow: auto; max-height: 78vh; padding: 18px 10px 22px; }

  .zoom-toolbar {
    position: sticky; top: 8px; left: 8px; z-index: 5; width: fit-content;
    display: inline-flex; align-items: center; gap: 4px; margin-bottom: 8px;
    background: var(--card); border: 1px solid var(--border); border-radius: 8px;
    padding: 4px 6px; box-shadow: 0 1px 4px rgba(0,0,0,0.12);
  }
  .zoom-toolbar button {
    font: inherit; font-size: 13px; font-weight: 600; color: var(--ink); background: var(--surface);
    border: 1px solid var(--border); border-radius: 5px; width: 24px; height: 24px; line-height: 1;
    cursor: pointer; display: inline-flex; align-items: center; justify-content: center;
  }
  .zoom-toolbar button:hover { border-color: var(--ink-muted); }
  .zoom-toolbar #zoomResetBtn { width: auto; padding: 0 8px; font-size: 11px; }
  #zoomLevel { font-size: 11.5px; color: var(--ink-2); min-width: 38px; text-align: center; }
  #zoomWrap { display: inline-block; transform-origin: 0 0; }

  .side {
    border-bottom: 1px solid var(--border); padding: 14px 16px;
    display: flex; flex-direction: row; flex-wrap: wrap; align-items: flex-start; gap: 14px 24px;
  }
  @media (max-width: 700px) { .side { flex-direction: column; } }

  .side h2 {
    font-size: 11px; font-weight: 600; letter-spacing: 0.06em; text-transform: uppercase;
    color: var(--ink-muted); margin: 0;
  }

  .detail-empty { font-size: 13px; color: var(--ink-muted); line-height: 1.6; }

  #detail { display: flex; flex-direction: row; flex-wrap: wrap; gap: 14px 24px; align-items: flex-start; flex: 1; }

  .detail-card { display: flex; flex-direction: row; flex-wrap: wrap; gap: 14px 24px; align-items: flex-start; }
  .detail-section { flex: 1 1 170px; min-width: 150px; }
  .id-row { display: flex; align-items: baseline; gap: 8px; flex-wrap: wrap; }
  .id-big { font-size: 26px; font-weight: 650; }
  .id-tag {
    font-size: 11px; font-weight: 600; padding: 2px 7px; border-radius: 5px;
    background: var(--card); border: 1px solid var(--border); color: var(--ink-2);
  }
  .id-tag.post { color: var(--ink); border-color: var(--ink); }

  .kv { display: grid; grid-template-columns: auto 1fr; gap: 5px 12px; font-size: 12.5px; margin-top: 8px; }
  .kv dt { color: var(--ink-muted); }
  .kv dd { margin: 0; color: var(--ink); font-weight: 500; }

  .rel-list { display: flex; flex-direction: column; gap: 4px; font-size: 12.5px; }
  .rel-list .none { color: var(--ink-muted); font-style: italic; }
  .rel-chip { display: inline-flex; align-items: center; gap: 5px; padding: 1px 0; }
  .dot { width: 7px; height: 7px; border-radius: 50%; flex: none; }

  .stat-row { display: flex; gap: 10px; min-width: 170px; flex: 0 0 auto; }
  .stat { flex: 1; background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 8px 10px; }
  .stat .n { font-size: 18px; font-weight: 650; }
  .stat .l { font-size: 10.5px; color: var(--ink-muted); text-transform: uppercase; letter-spacing: 0.05em; }

  svg { display: block; }
  .row-label { font-size: 11px; fill: var(--ink-muted); }
  .row-label.founder-label { font-weight: 600; fill: var(--ink-2); }
  .lane-label {
    font-size: 10.5px; font-weight: 600; letter-spacing: 0.04em; fill: var(--ink-muted); text-transform: uppercase;
  }
  .gridline { stroke: var(--gridline); stroke-width: 1; }

  .edge-path { fill: none; stroke: var(--edge); stroke-width: 1.4; transition: stroke 0.12s ease, opacity 0.12s ease; }
  .node-group { cursor: pointer; transition: opacity 0.12s ease; }
  .node-dot { stroke-width: 1.6; }
  .node-dot.founder-ring { stroke-dasharray: 2.2 2.2; }
  .node-hit { fill: transparent; cursor: pointer; }
  .pre-label {
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    font-variant-numeric: tabular-nums;
    font-size: 9.5px; font-weight: 700; text-anchor: middle; pointer-events: none;
  }
  .pre-label.on-fill { fill: #ffffff; }
  .pre-label.on-ring { fill: var(--ink); }
  .post-badge {
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    font-variant-numeric: tabular-nums;
    font-size: 8.5px; font-weight: 700; fill: var(--ink); text-anchor: middle; pointer-events: none;
  }
  .post-badge-bg { fill: var(--surface); stroke: var(--ink); stroke-width: 1; }

  .dimmed .edge-path { opacity: 0.12; }
  .node-group.dimmed { opacity: 0.16; }
  .hi .edge-path { opacity: 1; stroke: var(--ink-2); stroke-width: 1.8; }
  .node-group.hi { opacity: 1; }

  #tooltip {
    position: fixed; pointer-events: none; background: var(--ink); color: var(--surface);
    font-size: 12px; line-height: 1.5; padding: 7px 10px; border-radius: 7px; max-width: 240px;
    z-index: 50; opacity: 0; transform: translate(-9999px, -9999px); transition: opacity 0.08s ease;
  }
  #tooltip b { font-weight: 700; }
  #tooltip .tt-muted { opacity: 0.7; }

  details.table-toggle { margin-top: 18px; }
  details.table-toggle > summary {
    cursor: pointer; font-size: 13px; font-weight: 600; padding: 12px 16px; list-style: none;
    display: flex; align-items: center; gap: 8px;
  }
  details.table-toggle > summary::-webkit-details-marker { display: none; }
  details.table-toggle > summary::before {
    content: "▸"; font-size: 11px; color: var(--ink-muted); transition: transform 0.12s ease;
  }
  details.table-toggle[open] > summary::before { transform: rotate(90deg); }
  .table-scroll { overflow-x: auto; border-top: 1px solid var(--border); max-height: 420px; overflow-y: auto; }
  table.raw { border-collapse: collapse; width: 100%; font-size: 12px; }
  table.raw th {
    position: sticky; top: 0; background: var(--surface); text-align: left; font-weight: 600;
    color: var(--ink-muted); text-transform: uppercase; font-size: 10.5px; letter-spacing: 0.04em;
    padding: 8px 12px; border-bottom: 1px solid var(--border); white-space: nowrap;
  }
  table.raw td { padding: 6px 12px; border-bottom: 1px solid var(--gridline); white-space: nowrap; color: var(--ink-2); }
  table.raw tr.sample-row td { color: var(--ink); font-weight: 500; }
  table.raw tbody tr:hover { background: var(--card); }

  /* genome view: chromosome / tree pickers and mutation crosses */
  .genome-bar {
    display: flex; align-items: center; flex-wrap: wrap; gap: 8px 16px;
    padding: 12px 16px; border-bottom: 1px solid var(--border); font-size: 12.5px; color: var(--ink-2);
  }
  .genome-bar label { display: inline-flex; align-items: center; gap: 6px; font-weight: 600; color: var(--ink); }
  .genome-bar select {
    font: inherit; font-size: 12px; color: var(--ink); background: var(--card);
    border: 1px solid var(--border); border-radius: 6px; padding: 4px 8px; max-width: 320px;
  }
  #genomeInfo { color: var(--ink-muted); }
  /* bar of trees (or chromosomes), each segment as wide as the stretch of genome it covers */
  .tree-bar {
    flex: 1 1 100%; display: flex; height: 30px; border: 1px solid var(--border);
    border-radius: 6px; overflow: hidden; background: var(--card);
  }
  .tree-seg {
    flex: none; height: 100%; display: flex; align-items: center; justify-content: center;
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace; font-size: 10.5px; font-weight: 600;
    color: var(--ink-2); white-space: nowrap; overflow: hidden; cursor: pointer;
    background: var(--seg-a); border-right: 1px solid var(--surface);
  }
  .tree-seg:nth-child(even) { background: var(--seg-b); }
  .tree-seg:last-child { border-right: none; }
  .tree-seg:hover { filter: brightness(0.93); }
  .tree-seg.selected { background: var(--ink); color: var(--surface); }
  .swatch.mut { border-radius: 0; background: none; position: relative; }
  .swatch.mut::before, .swatch.mut::after {
    content: ""; position: absolute; left: 5px; top: -1px; width: 2px; height: 13px; background: var(--mut);
  }
  .swatch.mut::before { transform: rotate(45deg); }
  .swatch.mut::after { transform: rotate(-45deg); }
  .mut-x { stroke: var(--mut); stroke-width: 2; stroke-linecap: round; pointer-events: none; }
  .mut-hit { fill: transparent; cursor: help; }
  .mut-group.hidden { display: none; }
  .edge-path.off { opacity: 0.06; }
  .node-group.off { opacity: 0.18; }

  footer.note { margin-top: 16px; font-size: 12px; color: var(--ink-muted); line-height: 1.6; max-width: 74ch; }
</style>

<div class="viz-root">
  <div class="wrap">
    <header class="page">
      <p class="eyebrow">KR malaria simulator — pedigree, pre-simplify</p>
      <h1>Strain transmission map — seed __SEED__</h1>
      <p class="dek">Every strain, every generation, drawn from the raw (pre-<code>simplify()</code>) tables — so nothing is compressed away. The number in each circle is its <strong>pre-simplify node ID</strong>; the small corner badge, where present, is the <strong>post-simplify ID</strong> from the final tree sequence. <strong>Click a strain</strong> to trace its full ancestry and descendants; hover for detail. Each <strong>&times;</strong> is a mutation, drawn on the link and at the time it happened. <strong>Pick a chromosome, then a tree in the bar</strong> (each tree is as wide as the stretch of genome it covers) to see just the genealogy of that stretch and its mutations.</p>
      <div class="meta-row" id="metaRow"></div>
      <div class="warn-banner" id="warnBanner"></div>
    </header>

    <div class="panel">
      <div class="toolbar">
        <div class="legend" id="legend"></div>
        <button id="clearBtn" disabled>Clear selection</button>
      </div>

      <div class="genome-bar">
        <label>Chromosome <select id="chromSel"></select></label>
        <span id="genomeInfo"></span>
        <div class="tree-bar" id="treeBar"></div>
      </div>

      <div class="stage">
        <div class="side">
          <h2>Selection</h2>
          <div id="detail"><p class="detail-empty">Click any strain in the map to see its identity, host, and full lineage — ancestors traced back to the founders, descendants traced forward to today.</p></div>
        </div>
        <div class="canvas-scroll">
          <div class="zoom-toolbar">
            <button id="zoomOutBtn" type="button" aria-label="Zoom out">&minus;</button>
            <span id="zoomLevel" class="mono">100%</span>
            <button id="zoomInBtn" type="button" aria-label="Zoom in">+</button>
            <button id="zoomResetBtn" type="button">Reset</button>
          </div>
          <div id="zoomWrap">
            <svg id="canvas"></svg>
          </div>
        </div>
      </div>

      <details class="table-toggle">
        <summary>Raw lineage rows (<span id="rowCount"></span> transmission edges)</summary>
        <div class="table-scroll">
          <table class="raw">
            <thead>
              <tr><th>gen</th><th>child</th><th>child post-ID</th><th>recomb ID</th><th>parent</th><th>parent post-ID</th></tr>
            </thead>
            <tbody id="rawBody"></tbody>
          </table>
        </div>
      </details>
    </div>

    <footer class="note" id="footerNote"></footer>
  </div>
</div>

<div id="tooltip"></div>

<script id="lineage-data" type="application/json">
__DATA_JSON__
</script>

<script>
(function () {
  const data = JSON.parse(document.getElementById("lineage-data").textContent);
  const { T, Nh, neverUsedFounders, nodes, edges, genOrder,
          chroms, trees, numTrees, segments, mutations, unplacedMutations } = data;

  const LANE_PALETTE = ["var(--lane-1)", "var(--lane-2)", "var(--lane-3)", "var(--lane-4)",
                         "var(--lane-5)", "var(--lane-6)", "var(--lane-7)", "var(--lane-8)"];
  const OVERFLOW_COLOR = "var(--lane-overflow)";
  function laneColor(i) { return i === null || i === undefined ? OVERFLOW_COLOR : (i < 8 ? LANE_PALETTE[i] : OVERFLOW_COLOR); }

  document.getElementById("rowCount").textContent = edges.length;
  document.getElementById("metaRow").innerHTML = `
    <span class="chip"><b>${nodes.length}</b> strains</span>
    <span class="chip"><b>${edges.length}</b> transmission edges</span>
    <span class="chip"><b>${T + 1}</b> generations</span>
    <span class="chip"><b>${Nh}</b> hosts</span>
    <span class="chip"><b>${nodes.filter(n => n.isSample).length}</b> current samples</span>
    <span class="chip"><b>${mutations.length}</b> mutations</span>
    <span class="chip"><b>${numTrees}</b> trees</span>
  `;

  const legend = document.getElementById("legend");
  const visibleLanes = Math.min(Nh, 8);
  let legendHTML = "";
  for (let i = 0; i < visibleLanes; i++) {
    legendHTML += `<div class="legend-group"><span class="swatch" style="background:${LANE_PALETTE[i]}"></span>Host slot ${i + 1}</div>`;
  }
  if (Nh > 8) {
    legendHTML += `<div class="legend-group"><span class="swatch" style="background:${OVERFLOW_COLOR}"></span>+${Nh - 8} more (shared color)</div>`;
  }
  legendHTML += `<div class="legend-group"><span class="swatch founder"></span>Founder</div>`;
  legendHTML += `<div class="legend-group"><span class="swatch sample"></span>Current sample</div>`;
  legendHTML += `<div class="legend-group"><span class="swatch mut"></span>Mutation</div>`;
  legendHTML += `<div class="legend-group">Number = pre-simplify ID &middot; badge = post-simplify ID</div>`;
  legend.innerHTML = legendHTML;

  const warnBanner = document.getElementById("warnBanner");
  const warnings = [];
  if (nodes.length > 400) {
    warnings.push(`This run has ${nodes.length} strains across ${T + 1} generations — this diagram is built for small debug runs (a handful of hosts, a few dozen generations) and gets hard to read at this scale. Consider a smaller --Nh/--T for a run you specifically want to visualize.`);
  }
  if (Nh > 8) {
    warnings.push(`${Nh} hosts exceeds the 8-color host-slot palette — slots beyond 8 share a single neutral color and are only distinguishable by hovering.`);
  }
  if (unplacedMutations > 0) {
    warnings.push(`${unplacedMutations} mutation(s) could not be placed on the pedigree and are not drawn.`);
  }
  if (warnings.length) {
    warnBanner.classList.add("show");
    warnBanner.innerHTML = warnings.join("<br>");
  }

  document.getElementById("footerNote").innerHTML = `
    Built from <span class="mono">lineage___SEED__.csv</span> + <span class="mono">node_id_map___SEED__.csv</span> + <span class="mono">ts_mut___SEED__.trees</span> (seed __SEED__).
    "Host slot" is this generation's draw order, not a persistent host identity — the simulator creates a fresh individual
    record for every host every generation, so slot 1 at one generation and slot 1 at another are different tskit
    individuals that happen to occupy the same position in the transmission loop.
    ${neverUsedFounders > 0 ? `${neverUsedFounders} founder strain(s) were never sampled as a parent in this run — that's expected sampling variance, not a data gap.` : ""}
  `;

  // ---- layout constants ----
  const MARGIN_LEFT = 96;
  const MARGIN_TOP = 46;
  const MARGIN_RIGHT = 30;
  const MARGIN_BOTTOM = 20;
  const ROW_H = 70;
  const LANE_W = 216;
  const SLOT_SPACING = 32;
  const R = 10.5;
  const R_SAMPLE = 11.5;

  const laneCount = Math.max(1, Nh);
  const width = MARGIN_LEFT + laneCount * LANE_W + MARGIN_RIGHT;
  const height = MARGIN_TOP + genOrder.length * ROW_H + MARGIN_BOTTOM;

  function rowIndex(gen) { return T - gen; }
  function nodeX(n) {
    const lane = n.lane === null ? laneCount - 1 : n.lane;
    const laneCenter = MARGIN_LEFT + lane * LANE_W + LANE_W / 2;
    const offset = (n.slot - (n.slotCount - 1) / 2) * SLOT_SPACING;
    return laneCenter + offset;
  }
  function nodeY(n) { return MARGIN_TOP + rowIndex(n.gen) * ROW_H; }

  const byId = new Map(nodes.map(n => [n.id, n]));
  const parentsOf = new Map();
  const childrenOf = new Map();
  nodes.forEach(n => { parentsOf.set(n.id, []); childrenOf.set(n.id, []); });
  edges.forEach(e => {
    parentsOf.get(e.child).push(e.parent);
    childrenOf.get(e.parent).push(e.child);
  });

  const svg = document.getElementById("canvas");
  svg.setAttribute("width", width);
  svg.setAttribute("height", height);
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

  // ---- zoom ----
  const canvasScroll = document.querySelector(".canvas-scroll");
  const zoomWrap = document.getElementById("zoomWrap");
  const zoomLevelEl = document.getElementById("zoomLevel");
  let zoom = 1;
  const ZOOM_MIN = 0.2, ZOOM_MAX = 3, ZOOM_STEP = 0.2;

  function zoomAt(newZoom, clientX, clientY) {
    newZoom = Math.min(ZOOM_MAX, Math.max(ZOOM_MIN, newZoom));
    if (newZoom === zoom) return;
    const rect = canvasScroll.getBoundingClientRect();
    const px = clientX - rect.left, py = clientY - rect.top;
    const ratio = newZoom / zoom;
    const newScrollLeft = (canvasScroll.scrollLeft + px) * ratio - px;
    const newScrollTop  = (canvasScroll.scrollTop  + py) * ratio - py;
    zoom = newZoom;
    zoomWrap.style.transform = `scale(${zoom})`;
    zoomLevelEl.textContent = Math.round(zoom * 100) + "%";
    canvasScroll.scrollLeft = newScrollLeft;
    canvasScroll.scrollTop  = newScrollTop;
  }
  function zoomAtCenter(newZoom) {
    const rect = canvasScroll.getBoundingClientRect();
    zoomAt(newZoom, rect.left + rect.width / 2, rect.top + rect.height / 2);
  }
  document.getElementById("zoomInBtn").addEventListener("click", () => zoomAtCenter(zoom + ZOOM_STEP));
  document.getElementById("zoomOutBtn").addEventListener("click", () => zoomAtCenter(zoom - ZOOM_STEP));
  document.getElementById("zoomResetBtn").addEventListener("click", () => zoomAtCenter(1));
  canvasScroll.addEventListener("wheel", (ev) => {
    if (!ev.ctrlKey && !ev.metaKey) return;
    ev.preventDefault();
    const factor = ev.deltaY < 0 ? 1.1 : 1 / 1.1;
    zoomAt(zoom * factor, ev.clientX, ev.clientY);
  }, { passive: false });

  const NS = "http://www.w3.org/2000/svg";
  function el(tag, attrs) {
    const e = document.createElementNS(NS, tag);
    for (const k in attrs) e.setAttribute(k, attrs[k]);
    return e;
  }

  const gridLayer = el("g");
  genOrder.forEach(g => {
    const y = MARGIN_TOP + rowIndex(g) * ROW_H;
    gridLayer.appendChild(el("line", { class: "gridline", x1: MARGIN_LEFT - 14, x2: width - MARGIN_RIGHT, y1: y, y2: y }));
    const label = el("text", { class: "row-label" + (g === T ? " founder-label" : ""), x: MARGIN_LEFT - 20, y: y + 4, "text-anchor": "end" });
    label.textContent = g === T ? "founders" : `gen ${g}`;
    gridLayer.appendChild(label);
  });
  for (let i = 0; i < laneCount; i++) {
    const label = el("text", { class: "lane-label", x: MARGIN_LEFT + i * LANE_W + LANE_W / 2, y: MARGIN_TOP - 20, "text-anchor": "middle" });
    label.textContent = "slot " + (i + 1);
    gridLayer.appendChild(label);
  }
  svg.appendChild(gridLayer);

  const edgeLayer = el("g");
  const edgeEls = [];
  edges.forEach(e => {
    const c = byId.get(e.child), p = byId.get(e.parent);
    const x1 = nodeX(p), y1 = nodeY(p), x2 = nodeX(c), y2 = nodeY(c);
    const my = (y1 + y2) / 2;
    const d = `M ${x1} ${y1} C ${x1} ${my}, ${x2} ${my}, ${x2} ${y2}`;
    const path = el("path", { class: "edge-path", d });
    path.dataset.child = e.child;
    path.dataset.parent = e.parent;
    edgeLayer.appendChild(path);
    edgeEls.push(path);
  });
  svg.appendChild(edgeLayer);

  // ---- mutations: a cross on the parent-child link, at the height of its time ----
  const edgeKey = (c, p) => c + "-" + p;
  function pointOnLink(c, p, frac) {
    // the link is a cubic curve from parent (s = 0) to child (s = 1); find the
    // point whose height matches the mutation's time, frac of the way up from the child
    const x1 = nodeX(p), y1 = nodeY(p), x2 = nodeX(c), y2 = nodeY(c), my = (y1 + y2) / 2;
    const bez = (s, a, b, cc, d) => (1 - s) ** 3 * a + 3 * (1 - s) ** 2 * s * b + 3 * (1 - s) * s ** 2 * cc + s ** 3 * d;
    const yt = y2 + frac * (y1 - y2);
    let lo = 0, hi = 1;
    for (let i = 0; i < 30; i++) {
      const mid = (lo + hi) / 2;
      if (bez(mid, y1, my, my, y2) < yt) lo = mid; else hi = mid;
    }
    const sp = (lo + hi) / 2;
    return [bez(sp, x1, x1, x2, x2), bez(sp, y1, my, my, y2)];
  }
  const mutLayer = el("g");
  const mutEls = [];
  const ARM = 4.5;
  mutations.forEach(m => {
    const [x, y] = pointOnLink(byId.get(m.child), byId.get(m.parent), m.frac);
    const g = el("g", { class: "mut-group" });
    g.appendChild(el("line", { class: "mut-x", x1: x - ARM, y1: y - ARM, x2: x + ARM, y2: y + ARM }));
    g.appendChild(el("line", { class: "mut-x", x1: x - ARM, y1: y + ARM, x2: x + ARM, y2: y - ARM }));
    const hit = el("circle", { class: "mut-hit", cx: x, cy: y, r: ARM + 3 });
    g.appendChild(hit);
    hit.addEventListener("mouseenter", (ev) => showMutTooltip(m, ev));
    hit.addEventListener("mousemove", (ev) => moveTooltip(ev));
    hit.addEventListener("mouseleave", hideTooltip);
    mutLayer.appendChild(g);
    mutEls.push({ m, g });
  });
  svg.appendChild(mutLayer);

  const nodeLayer = el("g");
  const nodeEls = new Map();
  nodes.forEach(n => {
    const g = el("g", { class: "node-group" });
    const x = nodeX(n), y = nodeY(n);
    const r = n.isSample ? R_SAMPLE : R;
    const color = laneColor(n.lane);

    const dot = el("circle", {
      class: "node-dot" + (n.isFounder ? " founder-ring" : ""),
      cx: x, cy: y, r,
      fill: n.isFounder ? "none" : color,
      stroke: n.isSample ? "var(--ink)" : color,
    });
    if (n.isSample) dot.setAttribute("stroke-width", "2.5");
    g.appendChild(dot);

    const preLabel = el("text", { class: "pre-label " + (n.isFounder ? "on-ring" : "on-fill"), x, y: y + 3.3 });
    preLabel.textContent = n.id;
    g.appendChild(preLabel);

    if (n.post !== null) {
      const bx = x + r * 0.78, by = y - r * 0.78;
      const bg = el("circle", { class: "post-badge-bg", cx: bx, cy: by, r: 7.5 });
      const txt = el("text", { class: "post-badge", x: bx, y: by + 2.9 });
      txt.textContent = n.post;
      g.appendChild(bg);
      g.appendChild(txt);
    }

    const hit = el("circle", { class: "node-hit", cx: x, cy: y, r: r + 4 });
    g.appendChild(hit);

    nodeLayer.appendChild(g);
    nodeEls.set(n.id, { g, dot, x, y });

    hit.addEventListener("mouseenter", (ev) => showTooltip(n, ev));
    hit.addEventListener("mousemove", (ev) => moveTooltip(ev));
    hit.addEventListener("mouseleave", hideTooltip);
    hit.addEventListener("click", (ev) => { ev.stopPropagation(); selectNode(n.id); });
  });
  svg.appendChild(nodeLayer);

  const tooltip = document.getElementById("tooltip");
  function showTooltip(n, ev) {
    const lane = n.lane === null ? "unknown" : n.lane + 1;
    tooltip.innerHTML = `
      <b>node ${n.id}</b>${n.post !== null ? ` <span class="tt-muted">(post-simplify ${n.post})</span>` : ""}<br>
      ${n.isFounder ? "founder" : "gen " + n.gen} &middot; host slot ${lane}${n.recomb !== null ? " &middot; recomb " + n.recomb : ""}
    `;
    tooltip.style.opacity = 1;
    moveTooltip(ev);
  }
  function moveTooltip(ev) {
    const pad = 14;
    tooltip.style.transform = `translate(${ev.clientX + pad}px, ${ev.clientY + pad}px)`;
  }
  function hideTooltip() { tooltip.style.opacity = 0; }
  function showMutTooltip(m, ev) {
    tooltip.innerHTML = `
      <b>mutation</b> <span class="tt-muted">(tree ${m.tree})</span><br>
      chr ${m.chrom}: ${Math.round(m.local).toLocaleString()} bp<br>
      ${m.time.toFixed(2)} generations ago<br>
      on the link node ${m.child} &rarr; node ${m.parent}
    `;
    tooltip.style.opacity = 1;
    moveTooltip(ev);
  }

  // ---- genome view: whole genome, one chromosome, or one tree ----
  const segByChild = new Map();
  segments.forEach(([c, p, left, right]) => {
    if (!segByChild.has(c)) segByChild.set(c, []);
    segByChild.get(c).push({ p, left, right });
  });
  const sampleIds = nodes.filter(n => n.isSample).map(n => n.id);

  // pedigree links that carry the sampled strains' ancestry for [L, R)
  function linksCarrying(L, R) {
    const cuts = new Set([L, R]);
    segments.forEach(([, , left, right]) => {
      if (left > L && left < R) cuts.add(left);
      if (right > L && right < R) cuts.add(right);
    });
    const sorted = [...cuts].sort((a, b) => a - b);
    const keys = new Set(), onPath = new Set(sampleIds);
    for (let i = 0; i + 1 < sorted.length; i++) {
      const x = (sorted[i] + sorted[i + 1]) / 2;
      sampleIds.forEach(sid => {
        let cur = sid;
        for (;;) {
          const seg = (segByChild.get(cur) || []).find(s => s.left <= x && x < s.right);
          if (!seg) break;
          keys.add(edgeKey(cur, seg.p));
          onPath.add(seg.p);
          cur = seg.p;
        }
      });
    }
    return { keys, onPath };
  }

  const chromSel = document.getElementById("chromSel");
  const treeBar = document.getElementById("treeBar");
  const genomeInfo = document.getElementById("genomeInfo");
  const kb = x => (x / 1000).toLocaleString(undefined, { maximumFractionDigits: 1 });
  const mutCount = (L, R) => mutations.filter(m => m.pos >= L && m.pos < R).length;
  const genomeLength = chroms[chroms.length - 1].end;
  let selectedTree = null; // index into the current chromosome's trees, or null for the whole chromosome

  chromSel.innerHTML = `<option value="">Whole genome</option>` + chroms.map(c =>
    `<option value="${c.chrom}">Chromosome ${c.chrom} (${mutCount(c.start, c.end)} mutations)</option>`).join("");

  // one segment per tree on the chosen chromosome (or per chromosome for the whole genome),
  // its width the share of the bar that its stretch of genome covers
  function drawBar() {
    const c = chroms.find(ch => String(ch.chrom) === chromSel.value);
    const items = c
      ? trees.filter(t => t.chrom === c.chrom).map((t, i) => ({
          label: "t" + (i + 1), left: t.left, right: t.right,
          title: `t${i + 1} (tree ${t.index}): ${kb(t.left - c.start)}–${kb(t.right - c.start)} kb, ${mutCount(t.left, t.right)} mutations`,
        }))
      : chroms.map(ch => ({
          label: "chr" + ch.chrom, left: ch.start, right: ch.end,
          title: `Chromosome ${ch.chrom}: ${kb(ch.end - ch.start)} kb, ${mutCount(ch.start, ch.end)} mutations`,
        }));
    const span = c ? c.end - c.start : genomeLength;
    treeBar.innerHTML = "";
    items.forEach((it, i) => {
      const pct = (it.right - it.left) / span * 100;
      const seg = document.createElement("div");
      seg.className = "tree-seg" + (c && selectedTree === i ? " selected" : "");
      seg.style.width = pct + "%";
      seg.title = it.title;
      if (pct >= 2.5) seg.textContent = it.label;
      seg.addEventListener("click", () => {
        if (c) {
          selectedTree = selectedTree === i ? null : i;
        } else {
          chromSel.value = String(chroms[i].chrom);
          selectedTree = null;
        }
        drawBar();
        applyGenomeView();
      });
      treeBar.appendChild(seg);
    });
  }

  function applyGenomeView() {
    const c = chroms.find(ch => String(ch.chrom) === chromSel.value);
    if (!c) {
      edgeEls.forEach(pth => pth.classList.remove("off"));
      nodeEls.forEach(({ g }) => g.classList.remove("off"));
      mutEls.forEach(({ g }) => g.classList.remove("hidden"));
      genomeInfo.textContent = `Showing the whole pedigree and all ${mutations.length} mutations across ${numTrees} trees. Click a chromosome below.`;
      return;
    }
    const onChrom = trees.filter(t => t.chrom === c.chrom);
    const t = selectedTree === null ? null : onChrom[selectedTree];
    const L = t ? t.left : c.start, R = t ? t.right : c.end;
    const { keys, onPath } = linksCarrying(L, R);
    edgeEls.forEach(pth => pth.classList.toggle("off", !keys.has(edgeKey(pth.dataset.child, pth.dataset.parent))));
    nodeEls.forEach(({ g }, id) => g.classList.toggle("off", !onPath.has(id)));
    let shown = 0;
    mutEls.forEach(({ m, g }) => {
      const on = m.pos >= L && m.pos < R;
      g.classList.toggle("hidden", !on);
      if (on) shown++;
    });
    genomeInfo.textContent = t
      ? `t${selectedTree + 1} (tree ${t.index}) covers ${kb(L - c.start)}–${kb(R - c.start)} kb of chromosome ${c.chrom}: ${shown} mutations. Click it again for the whole chromosome.`
      : `Chromosome ${c.chrom}: ${onChrom.length} trees, ${shown} mutations. Click a tree below.`;
  }

  chromSel.addEventListener("change", () => { selectedTree = null; drawBar(); applyGenomeView(); });
  drawBar();
  applyGenomeView();

  const clearBtn = document.getElementById("clearBtn");
  const detail = document.getElementById("detail");
  let selectedId = null;

  function ancestorsOf(id) {
    const seen = new Set([id]);
    let frontier = [id];
    while (frontier.length) {
      const next = [];
      frontier.forEach(f => (parentsOf.get(f) || []).forEach(p => { if (!seen.has(p)) { seen.add(p); next.push(p); } }));
      frontier = next;
    }
    seen.delete(id);
    return seen;
  }
  function descendantsOf(id) {
    const seen = new Set([id]);
    let frontier = [id];
    while (frontier.length) {
      const next = [];
      frontier.forEach(f => (childrenOf.get(f) || []).forEach(c => { if (!seen.has(c)) { seen.add(c); next.push(c); } }));
      frontier = next;
    }
    seen.delete(id);
    return seen;
  }

  function clearSelection() {
    selectedId = null;
    nodeEls.forEach(({ g }) => g.classList.remove("hi", "dimmed"));
    edgeEls.forEach(p => p.classList.remove("hi", "dimmed"));
    clearBtn.disabled = true;
    detail.innerHTML = `<p class="detail-empty">Click any strain in the map to see its identity, host, and full lineage — ancestors traced back to the founders, descendants traced forward to today.</p>`;
  }

  function selectNode(id) {
    if (selectedId === id) { clearSelection(); return; }
    selectedId = id;
    const n = byId.get(id);
    const ancestors = ancestorsOf(id);
    const descendants = descendantsOf(id);
    const highlighted = new Set([id, ...ancestors, ...descendants]);

    nodeEls.forEach(({ g }, nid) => {
      g.classList.toggle("hi", highlighted.has(nid));
      g.classList.toggle("dimmed", !highlighted.has(nid));
    });
    edgeEls.forEach(p => {
      const c = Number(p.dataset.child), pa = Number(p.dataset.parent);
      const on = highlighted.has(c) && highlighted.has(pa);
      p.classList.toggle("hi", on);
      p.classList.toggle("dimmed", !on);
    });
    clearBtn.disabled = false;

    const lane = n.lane === null ? "unknown" : n.lane + 1;
    const directParents = parentsOf.get(id) || [];
    const directChildren = childrenOf.get(id) || [];

    function chipList(ids, emptyText) {
      if (!ids.length) return `<span class="none">${emptyText}</span>`;
      return ids.map(pid => {
        const pn = byId.get(pid);
        const c = laneColor(pn.lane);
        return `<span class="rel-chip"><span class="dot" style="background:${c}"></span>node ${pid}${pn.post !== null ? ` <span class="tt-muted">(→${pn.post})</span>` : ""}</span>`;
      }).join("");
    }

    detail.innerHTML = `
      <div class="detail-card">
        <div class="detail-section">
          <div class="id-row">
            <span class="id-big mono">${id}</span>
            ${n.post !== null ? `<span class="id-tag post">post-simplify: ${n.post}</span>` : `<span class="id-tag">not in final ts</span>`}
          </div>
          <dl class="kv">
            <dt>Generation</dt><dd>${n.isFounder ? "founder (gen " + T + ")" : "gen " + n.gen}</dd>
            <dt>Host slot</dt><dd>${lane}</dd>
            ${n.recomb !== null ? `<dt>Recomb ID</dt><dd>${n.recomb}</dd>` : ""}
          </dl>
        </div>
        <div class="stat-row">
          <div class="stat"><div class="n mono">${ancestors.size}</div><div class="l">Ancestors</div></div>
          <div class="stat"><div class="n mono">${descendants.size}</div><div class="l">Descendants</div></div>
        </div>
        <div class="detail-section">
          <h2 style="margin-bottom:6px">Direct parent(s)</h2>
          <div class="rel-list">${chipList(directParents, "none — this is a founder")}</div>
        </div>
        <div class="detail-section">
          <h2 style="margin-bottom:6px">Direct child(ren)</h2>
          <div class="rel-list">${chipList(directChildren, "none — never sampled as a parent")}</div>
        </div>
      </div>
    `;
  }

  clearBtn.addEventListener("click", clearSelection);
  svg.addEventListener("click", clearSelection);

  const rawBody = document.getElementById("rawBody");
  const rowsByChild = new Map();
  edges.forEach(e => {
    const c = byId.get(e.child), p = byId.get(e.parent);
    if (!rowsByChild.has(e.child)) rowsByChild.set(e.child, []);
    rowsByChild.get(e.child).push({ c, p });
  });
  [...nodes].sort((a, b) => b.gen - a.gen || a.id - b.id).forEach(n => {
    const rows = rowsByChild.get(n.id);
    if (!rows) return;
    rows.forEach(({ c, p }) => {
      const tr = document.createElement("tr");
      if (c.isSample) tr.className = "sample-row";
      tr.innerHTML = `
        <td class="mono">${c.gen}</td>
        <td class="mono">${c.id}</td>
        <td class="mono">${c.post !== null ? c.post : "—"}</td>
        <td class="mono">${c.recomb !== null ? c.recomb : "—"}</td>
        <td class="mono">${p.id}</td>
        <td class="mono">${p.post !== null ? p.post : "—"}</td>
      `;
      rawBody.appendChild(tr);
    });
  });
})();
</script>
"""


def save_strain_transmission_map(lineage, node_map, ts, sim_id, output_dir="."):
    """
    Render an interactive, self-contained HTML transmission map for one
    simulation run and write it to disk.

    Parameters
    ----------
    lineage : pd.DataFrame
        As written by output.lineage.save_lineage for this sim_id.
    node_map : pd.DataFrame
        As written by output.node_id_map.save_node_id_map for this sim_id.
    ts : tskit.TreeSequence
        Mutated tree sequence for this sim_id (ts_mut_<seed>.trees).
    sim_id : int
        Simulation / random seed identifier.
    output_dir : str
        Directory to write the HTML file to.
    """
    import json

    data = _build_viz_data(lineage, node_map)
    data.update(_build_genome_data(ts, lineage, node_map, data["nodes"]))
    seed = str(sim_id)
    html = (_HTML_TEMPLATE
            .replace("__DATA_JSON__", json.dumps(data))
            .replace("__SEED__", seed))

    output_path = os.path.join(output_dir, f"strain_transmission_map_{sim_id}.html")
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(html)

    return output_path
