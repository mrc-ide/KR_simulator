# krsim — detailed notes

Material moved out of the main README. Not linked from anywhere yet.

## Commands in detail

### Stage 1 — Simulate genealogy

```bash
krsim simulate \
    --Nh 100 --T 50 --K 2 --chi 0.18 --C 1 \
    --recombination_rate 7.4e-7 \
    --random_seed 42 \
    --output_dir results/ \
    --save_pairs
```

Always writes the tree sequence `ts_{seed}.trees`, which every later stage starts from, and the individual–node mapping `sim_ind_node_map_{seed}.csv` (which strains each host carries, with each node's `time` in generations before present; `0` = most recent/current generation), needed by `vcf-host` and `relate-host`. The `--save_*` flags below add extra outputs.

| Flag | Description |
|---|---|
| `--Nh` | Number of human hosts |
| `--T` | Transmission generations |
| `--K` | Parasites transmitted per mosquito bite |
| `--chi` | Superinfection probability (0–1) |
| `--C` | Oocysts per mosquito |
| `--recombination_rate` | Per-bp recombination rate — uses `7.4e-7` (*P. falciparum*) if not given |
| `--random_seed` | Seed for reproducibility |
| `--output_dir` | Output directory |
| `--save_pairs` | Save recombinant node pairs to `node_pairs_{seed}.csv` |
| `--save_lineage` | Save per-generation child→parent strain lineage to `lineage_{seed}.csv`, plus a `node_id_map_{seed}.csv` pre/post-`simplify()` ID mapping — turn this on if you want true generation-by-generation parentage (for `krsim strain-transmission-map`, `krsim pairwise-strain-pedigree-map`, or your own pedigree analysis), since the final tree sequence alone can skip generations |

### Strain transmission map — visualize the pedigree

```bash
krsim strain-transmission-map \
    --lineage_csv results/lineage_42.csv \
    --node_id_map_csv results/node_id_map_42.csv \
    --random_seed 42 \
    --output_dir results/
```

Needs `--save_lineage` from Stage 1. Renders an interactive, self-contained `strain_transmission_map_{seed}.html` — a generation-by-generation pedigree diagram (founders at top, current samples at bottom, one lane per host) built from the pre-simplify `lineage_{seed}.csv`, so it shows the true parentage that `simplify()` compresses away. Each strain is labeled with its pre-simplify ID (with its post-simplify ID as a small corner badge, where one exists); clicking a strain traces and highlights its full ancestry and descendants. Open the file directly in a browser — no server or dependencies needed. Intended for small debug runs (a handful of hosts, a few dozen generations); it warns in-page and degrades gracefully past that, but isn't meant for production-scale runs (`--Nh 100`, etc.). It only needs Stage 1, so you can inspect the raw transmission history before running mutations or VCFs.

### Pairwise strain pedigree map — trace two strains back to their shared ancestor

```bash
krsim pairwise-strain-pedigree-map \
    --lineage_csv results/lineage_42.csv \
    --node_id_map_csv results/node_id_map_42.csv \
    --random_seed 42 \
    --output_dir results/
```

Also only needs `--save_lineage` from Stage 1 — same inputs as the transmission map. Renders `pairwise_strain_pedigree_map_{seed}.html`: pick any two current-generation strains (dropdowns, or click two strains in the diagram) and it traces each one's full ancestry back through the pedigree, drawing strain 1's lineage in one color and strain 2's in another. Where the two colors overlap is a shared ancestor; the *most recent* shared ancestor(s) are highlighted as the coalescence point(s), with the number of generations back reported for each strain. Because a strain can have up to two immediate parents (one per gamete in the oocyst that produced it), there can be more than one coalescence point — the two ancestries can independently merge via more than one path.

This is a **pedigree** trace, not a genomic one: a shared ancestor here means the two strains' host lineages shared a parasite at that generation, not that any specific genomic region was actually inherited from it — recombination can replace a segment with DNA from the other parent at every step. Same small-debug-run scope as the transmission map.

### Stage 2 — Add mutations

```bash
krsim mutate \
    --ts results/ts_42.trees \
    --mutation_rate 1e-6 \
    --random_seed 42 \
    --output_dir results/
```

If `--mutation_rate` is not given, `1e-6` per bp is used. Outputs `mutated_{seed}.trees`. All stages from here onwards use the mutated tree sequence.

### Stage 3 — Coalescence & mutation times

```bash
krsim coal \
    --ts results/mutated_42.trees \
    --random_seed 42 \
    --output_dir results/
```

Outputs `coal_T_{seed}.csv` (pairwise MRCA times per tree interval) and `mut_T_{seed}.csv` (position and age of every mutation). Both tables include a `chromosome` column (1-indexed, from `krsim/sim/params.py`'s `COMP_BREAKPOINTS`), with `position`/`left`/`right` given as chromosome-local coordinates rather than genome-wide ones. `coal_T` only contains rows where the sample pair has a common ancestor in that tree interval — pairs with no shared ancestry anywhere in the ARG are simply absent from the table.

### Stage 4 — Strain-level VCF

```bash
krsim vcf-strain \
    --ts results/mutated_42.trees \
    --random_seed 42 \
    --output_dir results/
```

Outputs `output_{seed}_strain.vcf` — one sample column per parasite strain (node), genotypes `0`/`1`.

### Stage 5 — Host-level VCF

```bash
krsim vcf-host \
    --strain_vcf results/output_42_strain.vcf \
    --ind_node_map_csv results/sim_ind_node_map_42.csv \
    --random_seed 42 \
    --depth 10 \
    --output_dir results/ \
    --save_proportions
```

Aggregates the strain VCF to individual (host) level using a **Dirichlet-Multinomial** model:
1. Within-host strain proportions sampled once per host from Dirichlet(α=1 per strain)
2. `--depth` sequencing reads drawn once per host via Multinomial from those proportions
3. At each variant site, genotype called as `0/0` (all reads carry ref), `1/1` (all reads carry alt), or `0/1` (mixed) from the fixed per-strain read counts

Outputs `output_{seed}_host.vcf`. With `--save_proportions`, also writes `host_proportions_{seed}.csv`.

---

## Relatedness statistics

Optional ground-truth relatedness, computed from Stage 3's `coal_T` table. Not part of the tutorial — run these when you need true IBD.

```bash
krsim relate-strain \
    --coal_csv results/coal_T_42.csv \
    --random_seed 42 \
    --time_boundary 2 \
    --output_dir results/

krsim relate-host \
    --relate_csv results/relate_strain_42.csv \
    --ind_node_map_csv results/sim_ind_node_map_42.csv \
    --random_seed 42 \
    --output_dir results/
```

`krsim relate-strain`: for each strain pair in `coal_T_{seed}.csv`, sums the length of tree intervals where `coal_time <= --time_boundary` and divides by `GENOME_LENGTH` to give the fraction of the genome that is IBD below that time threshold. Pairs with no recorded coalescence anywhere in the ARG are not written as a row (rather than a `0` row). Outputs `relate_strain_{seed}.csv` (`sample_1`, `sample_2`, `frac_IBD`).

`krsim relate-host`: for each currently-infected host (from `sim_ind_node_map_{seed}.csv`, restricted to `time == 0`, i.e. the most recent generation), takes every pairwise combination of that host's strain nodes and averages their `frac_IBD` from `relate_strain_{seed}.csv`. Pairs missing from that table (no recorded coalescence) are treated as `frac_IBD = 0`. Hosts with fewer than two strains get an empty average. Outputs `relate_host_{seed}.csv` (`host`, `avg_frac_IBD`).

---

## Output files (full table)

| File | Stage | Description |
|---|---|---|
| `ts_{seed}.trees` | 1 | Tree sequence (no mutations) — always written |
| `sim_ind_node_map_{seed}.csv` | 1 | Host → strain node mapping, with each node's `time` (generations before present) — always written |
| `node_pairs_{seed}.csv` | 1 | Recombinant node pairs from same meiosis |
| `lineage_{seed}.csv` | 1 | Per-generation child→parent strain lineage, pre-simplify node IDs (`--save_lineage`) |
| `node_id_map_{seed}.csv` | 1 | Pre-simplify → post-simplify node ID mapping (`--save_lineage`) |
| `mutated_{seed}.trees` | 2 | Tree sequence with mutations |
| `coal_T_{seed}.csv` | 3 | Pairwise coalescence times per tree interval (chromosome-local coordinates) |
| `mut_T_{seed}.csv` | 3 | Mutation chromosome, position, and age |
| `output_{seed}_strain.vcf` | 4 | Per-strain genotypes (0/1) |
| `output_{seed}_host.vcf` | 5 | Per-host genotypes (0/0, 0/1, 1/1) |
| `host_proportions_{seed}.csv` | 5 | Dirichlet proportions and read counts per host |
| `strain_transmission_map_{seed}.html` | additional | Interactive strain-transmission-map visualization |
| `pairwise_strain_pedigree_map_{seed}.html` | additional | Interactive pedigree coalescence trace for any two strains |
| `relate_strain_{seed}.csv` | stats | Pairwise strain IBD fraction (`krsim relate-strain`) |
| `relate_host_{seed}.csv` | stats | Average pairwise IBD fraction per host, most recent generation only (`krsim relate-host`) |
| `KR_sim_{seed}.log` | all | Run log — every command adds a block with the time, krsim version, exact command line, files written and run time (and the full error details if it fails) |

---

## Project structure

```
KR_SIM/
├── pyproject.toml        # Package definition: dependencies + the `krsim` command
├── src/krsim/            # The installable package
│   ├── cli/              #   The `krsim` command: one module per subcommand (flags, logging, file saving)
│   ├── sim/              #   Core simulation engine
│   │   ├── params.py         #   SimParams dataclass + genome constants
│   │   ├── simulator.py      #   KRSimulator class
│   │   ├── recombination.py  #   Compulsory and Poisson recombination
│   │   └── transmission.py   #   Human-mosquito-human transmission
│   ├── output/           #   Standard outputs built from the simulation
│   │   ├── sim_ind_node_id.py
│   │   ├── node_pairs.py
│   │   ├── lineage.py        #   Builds lineage_{seed}.csv (pre-simplify)
│   │   ├── node_id_map.py    #   Builds node_id_map_{seed}.csv
│   │   ├── coal_mut_times.py #   coal_T / mut_T tables
│   │   ├── vcf_strain_single.py
│   │   └── vcf_host_single.py
│   ├── stats/            #   Ground-truth relatedness
│   │   ├── relate_strain.py  #   Pairwise strain IBD fraction (from coal_T)
│   │   └── relate_host.py    #   Per-host average of relate_strain
│   ├── visualization/    #   Interactive pedigree HTML rendering
│   │   ├── strain_transmission_map.py
│   │   └── pairwise_strain_pedigree_map.py
│   └── config.py         #   Load SimParams from a YAML file
│
├── tutorial/             # Tutorial
│   ├── tutorial.ipynb    #   Runs every command on a small simulation, with outputs
│   └── README.md         #   How to open the tutorial
│
├── tests/                # Testing and regression tools
│   ├── test_sim_mut.sh   #   Run full pipeline + regression check
│   ├── compare_outputs.py#   Diff outputs against baseline/
│   ├── baseline/         #   Committed reference outputs
│   └── output/           #   Latest test run (not committed)
│
├── environment.yml       # Conda environment
└── requirements.txt      # Exact pinned versions (for reproducing results)
```

---

## Install options

### Option 1 — uv (fastest)
```bash
uv venv
source .venv/bin/activate
uv pip install -e .
```

### Option 2 — conda
```bash
conda env create -f environment.yml   # also runs `pip install -e .`
conda activate kr_malaria
```

### Option 3 — pip
```bash
python -m venv .venv
source .venv/bin/activate   # Mac / Linux
.venv\Scripts\activate      # Windows
pip install -e .
```

Check it worked with `krsim --help`, which lists every command.

`requirements.txt` pins the exact package versions the simulator was developed and tested with; install from it (`pip install -r requirements.txt`, then `pip install -e .`) if you need to reproduce results exactly.
