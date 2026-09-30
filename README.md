# krsim — Malaria Transmission Simulator

A forward simulator that brings together malaria modelling and population genetics to model parasite transmission. Whole-genome information is tracked through every recombination event as an ancestral recombination graph, stored efficiently with [tskit](https://tskit.dev). From the recorded genealogy (a tree sequence), krsim overlays mutations and outputs coalescence times across the genome and VCF files at the strain level and at the (observed) host level.

## Installation

You need [conda](https://conda-forge.org/download/) (Miniforge, Miniconda or Anaconda). Download the environment file and create the environment:

```bash
curl -O https://raw.githubusercontent.com/mrc-ide/KR_simulator/main/environment-user.yml
conda env create -f environment-user.yml
conda activate krsim
```

Check it worked with `krsim --help`, which lists every command. 

To update krsim later, run `pip install --upgrade --force-reinstall --no-deps git+https://github.com/mrc-ide/KR_simulator.git` with the environment active.

## Tutorial

The tutorial is a Jupyter notebook, [tutorial.ipynb](tutorial/tutorial.ipynb), that goes through each stage of the pipeline and visualises the outputs. To run it, activate the `krsim` environment, then run the commands below. 

```bash
curl -O https://raw.githubusercontent.com/mrc-ide/KR_simulator/main/tutorial/tutorial.ipynb
jupyter notebook tutorial.ipynb
```

It writes its files to `output/` next to the notebook.

## Quick start

With the `krsim` environment active, run the pipeline on a small example. Results go to `krsim_quickstart_output/` in the current folder (krsim creates it):

```bash
OUT=krsim_quickstart_output

krsim simulate   --Nh 3 --T 10 --K 2 --chi 0.18 --C 1 --random_seed 123 --output_dir $OUT
krsim mutate     --ts $OUT/ts_123.trees --random_seed 123 --output_dir $OUT
krsim coal       --ts $OUT/mutated_123.trees --random_seed 123 --output_dir $OUT
krsim vcf-strain --ts $OUT/mutated_123.trees --random_seed 123 --output_dir $OUT
krsim vcf-host   --strain_vcf $OUT/output_123_strain.vcf --ind_node_map_csv $OUT/sim_ind_node_map_123.csv \
                 --depth 10 --random_seed 123 --output_dir $OUT
```

`krsim <command> --help` lists every option for a command.

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

## Development

To work on the code itself, see [DEVELOPMENT.md](DEVELOPMENT.md).
