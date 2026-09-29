# krsim — Malaria Transmission Simulator

A forward simulator that brings together malaria modelling and population genetics to model parasite transmission. Whole-genome information is tracked through every recombination event as an ancestral recombination graph, stored efficiently with [tskit](https://tskit.dev). From the recorded genealogy (a tree sequence), krsim overlays mutations and outputs coalescence times across the genome and VCF files at the strain level and at the (observed) host level.

## Installation

Requires Python 3.9 or newer.

```bash
git clone https://github.com/mrc-ide/KR_simulator.git
cd KR_simulator
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -e .
```

Upgrading pip matters: older pip versions (such as the 21.x that ships with macOS's built-in Python) fail with `Directory cannot be installed in editable mode`.

Or with conda: `conda env create -f environment.yml && conda activate kr_malaria`.

To reproduce results with the exact package versions used in development, run `pip install -r requirements.txt` before `pip install -e .`.

### Check the install

Still in the `KR_simulator` folder with the environment active, run:

```bash
krsim --help
bash tests/test_sim_mut.sh
```

The test runs every command on a small simulation and compares the results in `tests/output/` with the saved results in `tests/baseline/` (decimal numbers only need to match to 12 significant digits, since the last digit can vary between machines). If everything is installed correctly, it ends with `All 13 check(s) passed.` Answer `N` (or press Enter) at the final prompt, which asks whether to replace the baseline.

## Quick start

From the `KR_simulator` folder, run the pipeline on a small example. Results go to `krsim_output/`, a folder next to the repo (krsim creates it), so they stay out of git:

```bash
OUT=../krsim_output

krsim simulate   --Nh 3 --T 10 --K 2 --chi 0.18 --C 1 --random_seed 123 --output_dir $OUT
krsim mutate     --ts $OUT/ts_123.trees --random_seed 123 --output_dir $OUT
krsim coal       --ts $OUT/mutated_123.trees --random_seed 123 --output_dir $OUT
krsim vcf-strain --ts $OUT/mutated_123.trees --random_seed 123 --output_dir $OUT
krsim vcf-host   --strain_vcf $OUT/output_123_strain.vcf --ind_node_map_csv $OUT/sim_ind_node_map_123.csv \
                 --depth 10 --random_seed 123 --output_dir $OUT
```

`krsim <command> --help` lists every option for a command. The [tutorial](tutorial/tutorial.ipynb) runs each command and shows what it writes.

## Pipeline

```mermaid
%%{init: {"theme": "base", "themeVariables": {"fontFamily": "ui-sans-serif, system-ui, -apple-system, Segoe UI, Helvetica, Arial, sans-serif", "fontSize": "14px", "lineColor": "#475569"}, "flowchart": {"curve": "basis", "nodeSpacing": 30, "rankSpacing": 50, "padding": 16}}}%%
flowchart TD
    subgraph SIM["Stage 1 · krsim simulate"]
        ts("ts_{seed}.trees"):::file
        sinm("sim_ind_node_map_{seed}.csv"):::file
        lin("lineage_{seed}.csv<br/>node_id_map_{seed}.csv<br/><i>with --save_lineage</i>"):::optfile
        pairs("node_pairs_{seed}.csv<br/><i>with --save_pairs</i>"):::optfile
    end

    subgraph MUT["Stage 2 · krsim mutate"]
        mts("mutated_{seed}.trees"):::file
    end

    subgraph COAL["Stage 3 · krsim coal"]
        coalT("coal_T_{seed}.csv"):::file
        mutT("mut_T_{seed}.csv"):::file
    end

    subgraph VCFS["Stage 4 · krsim vcf-strain"]
        svcf("output_{seed}_strain.vcf"):::file
    end

    subgraph VCFH["Stage 5 · krsim vcf-host"]
        hvcf("output_{seed}_host.vcf"):::file
        props("host_proportions_{seed}.csv<br/><i>with --save_proportions</i>"):::optfile
    end

    subgraph VIS["Pedigree visualizations · optional"]
        subgraph TMAP["krsim strain-transmission-map"]
            tmapH("strain_transmission_map_{seed}.html"):::file
        end
        subgraph PW["krsim pairwise-strain-pedigree-map"]
            pwH("pairwise_strain_pedigree_map_{seed}.html"):::file
        end
    end

    %% arrows point at the files inside each stage (not the stage box itself),
    %% which keeps mermaid from stretching the boxes
    ts --> mts
    mts --> coalT
    mts --> svcf
    svcf --> hvcf
    sinm --> hvcf
    lin -.-> tmapH
    lin -.-> pwH

    classDef file fill:#ffffff,stroke:#94a3b8,stroke-width:1px,color:#0f172a
    classDef optfile fill:#ffffff,stroke:#94a3b8,stroke-width:1px,stroke-dasharray:4 3,color:#334155
    classDef stage fill:#eff6ff,stroke:#3b82f6,stroke-width:1.5px,color:#1e3a8a
    class SIM,MUT,COAL,VCFS,VCFH,TMAP,PW stage
    style VIS fill:#f8fafc,stroke:#94a3b8,stroke-width:1px,stroke-dasharray:4 3,color:#475569
    linkStyle default stroke:#475569,stroke-width:1.6px
    linkStyle 5,6 stroke:#94a3b8,stroke-width:1.6px
```

Each blue box is a `krsim` command and the files it writes; dashed files need the flag shown. Arrows show which files each command reads. Every command also records what it did in `KR_sim_{seed}.log` in its `--output_dir`.

## Output files

| File | Command | Contents |
|---|---|---|
| `ts_{seed}.trees` | simulate | Tree sequence (no mutations) |
| `sim_ind_node_map_{seed}.csv` | simulate | Which strains each host carries |
| `lineage_{seed}.csv`, `node_id_map_{seed}.csv` | simulate `--save_lineage` | True generation-by-generation pedigree |
| `node_pairs_{seed}.csv` | simulate `--save_pairs` | Strain pairs from the same meiosis |
| `mutated_{seed}.trees` | mutate | Tree sequence with mutations |
| `coal_T_{seed}.csv` | coal | Pairwise coalescence times per tree interval |
| `mut_T_{seed}.csv` | coal | Position and age of every mutation |
| `output_{seed}_strain.vcf` | vcf-strain | Per-strain genotypes (0/1) |
| `output_{seed}_host.vcf` | vcf-host | Per-host genotypes (0/0, 0/1, 1/1) |
| `host_proportions_{seed}.csv` | vcf-host `--save_proportions` | Strain proportions and reads per host |
| `strain_transmission_map_{seed}.html` | strain-transmission-map | Interactive pedigree of the whole run |
| `pairwise_strain_pedigree_map_{seed}.html` | pairwise-strain-pedigree-map | Traces two strains back to shared ancestors |
| `KR_sim_{seed}.log` | all | What each command ran and wrote |
