import os
import logging

from output.transmission_map import _build_viz_data


_HTML_TEMPLATE = r"""<title>Pairwise coalescence — sim __SEED__</title>
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
    --card:           #ffffff;
    --edge:           #c9c8bf;
    --a-color:        #2a78d6;
    --b-color:        #eb6834;
    --common-color:   #9a5fd6;
    --mrca-color:     #c9a800;
    --warn-bg:        #fdf3e4;
    --warn-ink:       #7a4a00;
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
      --card:           #222221;
      --edge:           #3c3c39;
      --a-color:        #3987e5;
      --b-color:        #d95926;
      --common-color:   #b087e6;
      --mrca-color:     #e0c030;
      --warn-bg:        #2c2410;
      --warn-ink:       #e0b25c;
    }
  }
  :root[data-theme="dark"] .viz-root {
    color-scheme: dark;
    --surface: #1a1a19; --page: #0d0d0d; --ink: #ffffff; --ink-2: #c3c2b7; --ink-muted: #898781;
    --gridline: #2c2c2a; --baseline: #383835; --border: rgba(255,255,255,0.10);
    --card: #222221; --edge: #3c3c39;
    --a-color: #3987e5; --b-color: #d95926; --common-color: #b087e6; --mrca-color: #e0c030;
    --warn-bg: #2c2410; --warn-ink: #e0b25c;
  }
  :root[data-theme="light"] .viz-root {
    color-scheme: light;
    --surface: #fcfcfb; --page: #f9f9f7; --ink: #0b0b0b; --ink-2: #52514e; --ink-muted: #898781;
    --gridline: #e1e0d9; --baseline: #c3c2b7; --border: rgba(11,11,11,0.10);
    --card: #ffffff; --edge: #c9c8bf;
    --a-color: #2a78d6; --b-color: #eb6834; --common-color: #9a5fd6; --mrca-color: #c9a800;
    --warn-bg: #fdf3e4; --warn-ink: #7a4a00;
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
  .dek { font-size: 14px; line-height: 1.55; color: var(--ink-2); max-width: 70ch; margin: 0; }
  .dek strong { color: var(--ink); font-weight: 600; }

  .caveat {
    margin-top: 14px; font-size: 12.5px; line-height: 1.55; color: var(--ink-2);
    background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 10px 14px;
    max-width: 74ch;
  }
  .caveat b { color: var(--ink); }

  .warn-banner {
    display: none; margin-top: 14px; font-size: 12.5px; line-height: 1.5;
    background: var(--warn-bg); color: var(--warn-ink); border-radius: 8px; padding: 10px 14px;
  }
  .warn-banner.show { display: block; }

  .panel { background: var(--surface); border: 1px solid var(--border); border-radius: 12px; }

  .picker {
    display: flex; align-items: center; flex-wrap: wrap; gap: 14px 22px;
    padding: 14px 16px; border-bottom: 1px solid var(--border);
  }
  .picker-field { display: flex; align-items: center; gap: 8px; font-size: 12.5px; color: var(--ink-2); }
  .picker-field label { font-weight: 600; display: flex; align-items: center; gap: 6px; }
  .swatch { width: 11px; height: 11px; border-radius: 50%; flex: none; display: inline-block; }
  .swatch.a { background: var(--a-color); }
  .swatch.b { background: var(--b-color); }
  .swatch.common { background: var(--common-color); }
  .swatch.mrca { background: var(--mrca-color); }
  select {
    font: inherit; font-size: 12.5px; color: var(--ink); background: var(--card);
    border: 1px solid var(--border); border-radius: 6px; padding: 5px 8px;
  }
  #traceBtn, #clearBtn {
    font: inherit; font-size: 12px; font-weight: 600; color: var(--ink); background: var(--card);
    border: 1px solid var(--border); border-radius: 6px; padding: 6px 14px; cursor: pointer;
  }
  #traceBtn { color: #fff; background: var(--ink); border-color: var(--ink); }
  #traceBtn:hover, #clearBtn:hover { opacity: 0.85; }
  #traceBtn[disabled] { opacity: 0.4; cursor: default; }
  .picker-hint { font-size: 11.5px; color: var(--ink-muted); }

  .stage { display: grid; grid-template-columns: minmax(0, 1fr) 300px; }
  @media (max-width: 900px) { .stage { grid-template-columns: 1fr; } }

  .canvas-scroll { overflow-x: auto; padding: 18px 10px 22px; }

  .side {
    border-left: 1px solid var(--border); padding: 16px;
    display: flex; flex-direction: column; gap: 14px; min-height: 320px;
  }
  @media (max-width: 900px) { .side { border-left: none; border-top: 1px solid var(--border); } }

  .side h2 {
    font-size: 11px; font-weight: 600; letter-spacing: 0.06em; text-transform: uppercase;
    color: var(--ink-muted); margin: 0;
  }

  .detail-empty { font-size: 13px; color: var(--ink-muted); line-height: 1.6; }

  .sample-card { display: flex; flex-direction: column; gap: 4px; font-size: 12.5px; padding: 8px 10px; border-radius: 8px; background: var(--card); border: 1px solid var(--border); }
  .sample-card .id-big { font-size: 18px; font-weight: 650; }
  .sample-card .id-tag { font-size: 10.5px; color: var(--ink-muted); }

  .mrca-list { display: flex; flex-direction: column; gap: 8px; }
  .mrca-card { font-size: 12.5px; padding: 8px 10px; border-radius: 8px; background: var(--card); border: 1px solid var(--border); }
  .mrca-card .id-big { font-size: 15px; font-weight: 650; color: var(--mrca-color); }
  .mrca-card dl { display: grid; grid-template-columns: auto 1fr; gap: 3px 10px; margin: 6px 0 0; font-size: 12px; }
  .mrca-card dt { color: var(--ink-muted); }
  .mrca-card dd { margin: 0; color: var(--ink); font-weight: 500; }

  .stat-row { display: flex; gap: 10px; }
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

  .edge-path { fill: none; stroke: var(--edge); stroke-width: 1.4; opacity: 0.1; transition: stroke 0.12s ease, opacity 0.12s ease, stroke-width 0.12s ease; }
  .edge-path.a { stroke: var(--a-color); opacity: 0.85; }
  .edge-path.b { stroke: var(--b-color); opacity: 0.85; }
  .edge-path.common { stroke: var(--common-color); opacity: 1; stroke-width: 2.2; }
  .node-group { cursor: pointer; }
  .node-dot { stroke-width: 1.6; opacity: 0.16; transition: opacity 0.12s ease, stroke-width 0.12s ease; }
  .node-dot.founder-ring { stroke-dasharray: 2.2 2.2; }
  .node-hit { fill: transparent; cursor: pointer; }
  .node-group.a .node-dot { fill: var(--a-color); stroke: var(--a-color); opacity: 0.9; }
  .node-group.b .node-dot { fill: var(--b-color); stroke: var(--b-color); opacity: 0.9; }
  .node-group.common .node-dot { fill: var(--common-color); stroke: var(--common-color); opacity: 1; }
  .node-group.mrca .node-dot { stroke: var(--mrca-color); stroke-width: 3.5; opacity: 1; }
  .node-group.endpoint .node-dot { stroke: var(--ink); stroke-width: 3; }
  .pre-label {
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    font-variant-numeric: tabular-nums;
    font-size: 9.5px; font-weight: 700; text-anchor: middle; pointer-events: none; fill: #fff;
    opacity: 0.16; transition: opacity 0.12s ease;
  }
  .node-group.a .pre-label, .node-group.b .pre-label, .node-group.common .pre-label { opacity: 1; }
  .post-badge {
    font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace;
    font-variant-numeric: tabular-nums;
    font-size: 8.5px; font-weight: 700; fill: var(--ink); text-anchor: middle; pointer-events: none;
  }
  .post-badge-bg { fill: var(--surface); stroke: var(--ink); stroke-width: 1; }

  #tooltip {
    position: fixed; pointer-events: none; background: var(--ink); color: var(--surface);
    font-size: 12px; line-height: 1.5; padding: 7px 10px; border-radius: 7px; max-width: 240px;
    z-index: 50; opacity: 0; transform: translate(-9999px, -9999px); transition: opacity 0.08s ease;
  }
  #tooltip b { font-weight: 700; }
  #tooltip .tt-muted { opacity: 0.7; }

  footer.note { margin-top: 16px; font-size: 12px; color: var(--ink-muted); line-height: 1.6; max-width: 74ch; }
</style>

<div class="viz-root">
  <div class="wrap">
    <header class="page">
      <p class="eyebrow">KR malaria simulator — pedigree, pre-simplify</p>
      <h1>Pairwise coalescence — seed __SEED__</h1>
      <p class="dek">Pick two current-generation parasite strains. Each strain's full ancestry back to the founders is traced and drawn in its own color; where those two ancestries <strong>overlap</strong> is a shared ancestor, and the most recent shared ancestor(s) are the pedigree's <strong>coalescence point(s)</strong> for this pair.</p>
      <div class="caveat">
        <b>Reading this correctly:</b> this is the transmission <em>pedigree</em> — it shows when the two strains' host lineages last shared a parasite, generation by generation, regardless of genome position. A shared ancestor here does <em>not</em> guarantee any specific genomic region was actually inherited from it — recombination can replace a segment with DNA from the other parent at every step. For the true, recombination-aware, per-genomic-interval coalescence time, use <span class="mono">coal_T_{seed}.csv</span> / <span class="mono">pedigree_relate_strain_genome_{seed}.csv</span> instead. Also note a strain can have up to <strong>two</strong> immediate parents (one per gamete in the oocyst) — so the ancestor count can up to double each generation back, fewer if the same strain is drawn twice (a clonal self-cross) or if lineages merge.
      </div>
      <div class="warn-banner" id="warnBanner"></div>
    </header>

    <div class="panel">
      <div class="picker">
        <div class="picker-field">
          <label><span class="swatch a"></span>Strain 1</label>
          <select id="selA"></select>
        </div>
        <div class="picker-field">
          <label><span class="swatch b"></span>Strain 2</label>
          <select id="selB"></select>
        </div>
        <button id="traceBtn">Trace coalescence</button>
        <button id="clearBtn" disabled>Clear</button>
        <span class="picker-hint">or click any two current-sample strains in the diagram</span>
      </div>

      <div class="stage">
        <div class="canvas-scroll">
          <svg id="canvas"></svg>
        </div>
        <div class="side">
          <h2>Pair</h2>
          <div id="detail"><p class="detail-empty">Choose two strains above (or click two current-generation strains in the diagram) and hit "Trace coalescence".</p></div>
        </div>
      </div>
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
  const { T, Nh, neverUsedFounders, nodes, edges, genOrder } = data;

  document.getElementById("footerNote").innerHTML = `
    Built from <span class="mono">lineage___SEED__.csv</span> + <span class="mono">node_id_map___SEED__.csv</span> (seed __SEED__), same source as the transmission map.
    "Host slot" is this generation's draw order, not a persistent host identity.
    ${neverUsedFounders > 0 ? `${neverUsedFounders} founder strain(s) were never sampled as a parent in this run — expected sampling variance.` : ""}
  `;

  const warnBanner = document.getElementById("warnBanner");
  if (nodes.length > 400) {
    warnBanner.classList.add("show");
    warnBanner.innerHTML = `This run has ${nodes.length} strains across ${T + 1} generations — built for small debug runs, gets hard to read at this scale.`;
  }

  // ---- layout constants (matches transmission_map.py) ----
  const MARGIN_LEFT = 96;
  const MARGIN_TOP = 46;
  const MARGIN_RIGHT = 30;
  const MARGIN_BOTTOM = 20;
  const ROW_H = 100;
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

  function ancestorsOf(id) {
    const seen = new Set([id]);
    let frontier = [id];
    while (frontier.length) {
      const next = [];
      frontier.forEach(f => (parentsOf.get(f) || []).forEach(p => { if (!seen.has(p)) { seen.add(p); next.push(p); } }));
      frontier = next;
    }
    return seen; // includes id itself
  }

  const svg = document.getElementById("canvas");
  svg.setAttribute("width", width);
  svg.setAttribute("height", height);
  svg.setAttribute("viewBox", `0 0 ${width} ${height}`);

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

  const nodeLayer = el("g");
  const nodeEls = new Map();
  nodes.forEach(n => {
    const g = el("g", { class: "node-group" });
    const x = nodeX(n), y = nodeY(n);
    const r = n.isSample ? R_SAMPLE : R;

    const dot = el("circle", {
      class: "node-dot" + (n.isFounder ? " founder-ring" : ""),
      cx: x, cy: y, r,
      fill: "var(--edge)",
      stroke: "var(--edge)",
    });
    g.appendChild(dot);

    const preLabel = el("text", { class: "pre-label", x, y: y + 3.3 });
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
    nodeEls.set(n.id, { g, dot });

    hit.addEventListener("mouseenter", (ev) => showTooltip(n, ev));
    hit.addEventListener("mousemove", (ev) => moveTooltip(ev));
    hit.addEventListener("mouseleave", hideTooltip);
    hit.addEventListener("click", (ev) => { ev.stopPropagation(); clickNode(n); });
  });
  svg.appendChild(nodeLayer);

  const tooltip = document.getElementById("tooltip");
  function showTooltip(n, ev) {
    const lane = n.lane === null ? "unknown" : n.lane + 1;
    tooltip.innerHTML = `
      <b>node ${n.id}</b>${n.post !== null ? ` <span class="tt-muted">(post-simplify ${n.post})</span>` : ""}<br>
      ${n.isFounder ? "founder" : "gen " + n.gen} &middot; host slot ${lane}
    `;
    tooltip.style.opacity = 1;
    moveTooltip(ev);
  }
  function moveTooltip(ev) {
    const pad = 14;
    tooltip.style.transform = `translate(${ev.clientX + pad}px, ${ev.clientY + pad}px)`;
  }
  function hideTooltip() { tooltip.style.opacity = 0; }

  // ---- strain pickers ----
  const samples = nodes.filter(n => n.isSample).sort((a, b) => a.id - b.id);
  const selA = document.getElementById("selA");
  const selB = document.getElementById("selB");
  function fillSelect(sel) {
    sel.innerHTML = samples.map(n =>
      `<option value="${n.id}">pre ${n.id}${n.post !== null ? ` (post ${n.post})` : ""} — host slot ${n.lane + 1}</option>`
    ).join("");
  }
  fillSelect(selA);
  fillSelect(selB);
  if (samples.length > 1) selB.selectedIndex = 1;

  const traceBtn = document.getElementById("traceBtn");
  const clearBtn = document.getElementById("clearBtn");
  const detail = document.getElementById("detail");

  let pendingClicks = [];
  function clickNode(n) {
    if (!n.isSample) return; // only current-generation strains are valid endpoints
    pendingClicks.push(n.id);
    if (pendingClicks.length === 1) {
      selA.value = String(n.id);
      clearHighlight();
      detail.innerHTML = `<p class="detail-empty">Strain 1 set to node ${n.id}. Click one more current-generation strain to trace coalescence.</p>`;
    } else {
      selB.value = String(n.id);
      pendingClicks = [];
      trace();
    }
  }

  function clearHighlight() {
    nodeEls.forEach(({ g }) => { g.setAttribute("class", "node-group"); });
    edgeEls.forEach(p => { p.setAttribute("class", "edge-path"); });
  }

  function nGen(id) { return byId.get(id).gen; }

  function trace() {
    pendingClicks = [];
    const idA = Number(selA.value), idB = Number(selB.value);
    clearHighlight();
    clearBtn.disabled = false;

    if (samples.length < 2) {
      detail.innerHTML = `<p class="detail-empty">This run only has one current-generation strain — nothing to compare.</p>`;
      return;
    }
    if (idA === idB) {
      detail.innerHTML = `<p class="detail-empty">Pick two <em>different</em> strains to trace a coalescence.</p>`;
      return;
    }

    const ancA = ancestorsOf(idA); // includes idA
    const ancB = ancestorsOf(idB); // includes idB
    const common = new Set([...ancA].filter(x => ancB.has(x)));

    // minimal common ancestors (MRCAs): process common set most-recent-first (smallest gen),
    // an ancestor is dominated (not an MRCA) if a more-recent common ancestor already covers it
    const commonSorted = [...common].sort((x, y) => nGen(x) - nGen(y));
    const covered = new Set();
    const mrcas = [];
    commonSorted.forEach(x => {
      if (covered.has(x)) return;
      mrcas.push(x);
      ancestorsOf(x).forEach(a => { if (a !== x && common.has(a)) covered.add(a); });
    });

    // paint nodes
    nodeEls.forEach(({ g }, id) => {
      const classes = ["node-group"];
      const inA = ancA.has(id), inB = ancB.has(id);
      if (common.has(id)) classes.push("common");
      else if (inA) classes.push("a");
      else if (inB) classes.push("b");
      if (mrcas.includes(id)) classes.push("mrca");
      if (id === idA || id === idB) classes.push("endpoint");
      g.setAttribute("class", classes.join(" "));
    });

    // paint edges: color by which trace(s) both endpoints belong to
    edgeEls.forEach(p => {
      const c = Number(p.dataset.child), pa = Number(p.dataset.parent);
      const classes = ["edge-path"];
      const bothCommon = common.has(c) && common.has(pa);
      const bothA = ancA.has(c) && ancA.has(pa);
      const bothB = ancB.has(c) && ancB.has(pa);
      if (bothCommon) classes.push("common");
      else if (bothA) classes.push("a");
      else if (bothB) classes.push("b");
      p.setAttribute("class", classes.join(" "));
    });

    const nA = byId.get(idA), nB = byId.get(idB);
    function sampleCard(n, colorClass, label) {
      return `
        <div class="sample-card">
          <div><span class="swatch ${colorClass}"></span> <b>${label}</b></div>
          <div class="id-big mono">${n.id}</div>
          <div class="id-tag">${n.post !== null ? `post-simplify ${n.post} &middot; ` : ""}host slot ${n.lane + 1}</div>
        </div>`;
    }

    let mrcaHTML;
    if (mrcas.length === 0) {
      mrcaHTML = `<p class="detail-empty">No shared ancestor found anywhere in the recorded pedigree — these two lineages remain fully independent all the way back to the founder generation.</p>`;
    } else {
      mrcaHTML = `<div class="mrca-list">` + mrcas.map(id => {
        const n = byId.get(id);
        const genBackA = nGen(id) - nGen(idA);
        const genBackB = nGen(id) - nGen(idB);
        return `
          <div class="mrca-card">
            <div class="id-big mono">${n.isFounder ? "founder " : ""}${n.id}</div>
            <dl>
              <dt>Generation</dt><dd>${n.isFounder ? "founder (gen " + T + ")" : "gen " + n.gen}</dd>
              <dt>Post-simplify</dt><dd>${n.post !== null ? n.post : "not in final ts"}</dd>
              <dt>Gens back (Strain 1)</dt><dd>${genBackA}</dd>
              <dt>Gens back (Strain 2)</dt><dd>${genBackB}</dd>
            </dl>
          </div>`;
      }).join("") + `</div>`;
    }

    detail.innerHTML = `
      <div class="sample-card-row" style="display:flex; flex-direction:column; gap:8px;">
        ${sampleCard(nA, "a", "Strain 1")}
        ${sampleCard(nB, "b", "Strain 2")}
      </div>
      <div class="stat-row">
        <div class="stat"><div class="n mono">${common.size}</div><div class="l">Shared ancestors</div></div>
        <div class="stat"><div class="n mono">${mrcas.length}</div><div class="l">Coalescence pt${mrcas.length === 1 ? "" : "s"}</div></div>
      </div>
      <div>
        <h2 style="margin-bottom:6px">Most recent common ancestor${mrcas.length === 1 ? "" : "s"}</h2>
        ${mrcaHTML}
      </div>
      ${mrcas.length > 1 ? `<p class="detail-empty">More than one coalescence point means the two strains' ancestries merge along more than one independent path (e.g. via different recombined parents) — consistent with each strain having up to two immediate parents.</p>` : ""}
    `;
  }

  function clearSelection() {
    pendingClicks = [];
    clearHighlight();
    clearBtn.disabled = true;
    detail.innerHTML = `<p class="detail-empty">Choose two strains above (or click two current-generation strains in the diagram) and hit "Trace coalescence".</p>`;
  }

  traceBtn.addEventListener("click", trace);
  clearBtn.addEventListener("click", clearSelection);
  svg.addEventListener("click", (ev) => { if (ev.target === svg) clearSelection(); });

  if (samples.length >= 2) trace();
})();
</script>
"""


def save_pairwise_coalescence(lineage, node_map, sim_id, output_dir="."):
    """
    Render an interactive, self-contained HTML page for tracing the
    pedigree coalescence of any two current-generation parasite strains
    and write it to disk.

    Parameters
    ----------
    lineage : pd.DataFrame
        As written by output.lineage.save_lineage for this sim_id.
    node_map : pd.DataFrame
        As written by output.node_id_map.save_node_id_map for this sim_id.
    sim_id : int
        Simulation / random seed identifier.
    output_dir : str
        Directory to write the HTML file to.
    """
    import json

    data = _build_viz_data(lineage, node_map)
    seed = str(sim_id)
    html = (_HTML_TEMPLATE
            .replace("__DATA_JSON__", json.dumps(data))
            .replace("__SEED__", seed))

    output_path = os.path.join(output_dir, f"pairwise_coalescence_{sim_id}.html")
    with open(output_path, "w") as f:
        f.write(html)

    logging.info(f"sim_id={sim_id} | wrote pairwise coalescence map to {output_path}")
    return output_path
