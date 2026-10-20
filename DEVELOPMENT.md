# krsim — development

Notes for working on the krsim code. To just use krsim, follow the [README](README.md).

## Install from a clone

Clone the repo and install krsim in editable mode, so changes to `src/krsim/` take effect without reinstalling. Requires Python 3.9 or newer.

```bash
git clone https://github.com/mrc-ide/KR_simulator.git
cd KR_simulator
```

Then pick one:

### conda
```bash
conda env create -f environment.yml   # also runs `pip install -e .`
conda activate kr_malaria
```

### pip
```bash
python3 -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -e ".[tutorial]"
```

Upgrading pip matters: older pip versions (such as the 21.x that ships with macOS's built-in Python) fail with `Directory cannot be installed in editable mode`.

### uv (fastest)
```bash
uv venv
source .venv/bin/activate
uv pip install -e ".[tutorial]"
```

Check it worked with `krsim --help`, which lists every command.

`requirements.txt` pins the exact package versions the simulator was developed and tested with; install from it (`pip install -r requirements.txt`, then `pip install -e .`) if you need to reproduce results exactly.

## Run the tests

From the `KR_simulator` folder with the environment active, run:

```bash
bash tests/test_sim_mut.sh
```

The test runs every command on a small simulation and compares the results in `tests/output/` with the saved results in `tests/baseline/` (decimal numbers only need to match to 12 significant digits, since the last digit can vary between machines). If everything is installed correctly, it ends with `All 16 check(s) passed.` Answer `N` (or press Enter) at the final prompt, which asks whether to replace the baseline.

## Website

The website source is in `docs/`, built with [Quarto](https://quarto.org/docs/get-started/). To preview it locally while editing:

```bash
quarto preview docs
```

Pushing to `main` rebuilds and publishes it through `.github/workflows/docs.yml`. The tutorial pages (`docs/tutorial-*.qmd`) link to the downloadable notebook in `notebooks/` rather than showing its outputs.

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
│
├── docs/                 # Quarto website (published to GitHub Pages)
│   ├── _quarto.yml       #   Site config and navigation
│   ├── *.qmd             #   Pages: welcome, background, installation, tutorial (3 pages), command line, output files
│   └── maps/             #   Interactive pedigree maps embedded in the tutorial
│
├── notebooks/            # Notebooks for users to download and run
│   └── tutorial.ipynb    #   Runs every command on a small simulation, with outputs
│
├── tests/                # Testing and regression tools
│   ├── test_sim_mut.sh   #   Run full pipeline + regression check
│   ├── compare_outputs.py#   Diff outputs against baseline/
│   ├── baseline/         #   Committed reference outputs
│   └── output/           #   Latest test run (not committed)
│
├── environment.yml       # Conda environment for development (editable install)
├── environment-user.yml  # Conda environment for users (installs krsim from GitHub)
└── requirements.txt      # Exact pinned versions (for reproducing results)
```
