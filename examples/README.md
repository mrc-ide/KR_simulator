# Examples

A worked example running the full KR Malaria Transmission Simulator pipeline with seed 123.

Install the package first (`pip install -e .` from the repo root) — see the main [README](../README.md#installation).

## How to run

**Option 1 — Jupyter notebook** (recommended for exploration):

```bash
# from the repo root
jupyter notebook examples/example.ipynb
```

**Option 2 — script** (runs all five stages in one go):

```bash
# from the repo root
python examples/run_example.py
```

**Option 3 — stage by stage** (run each command independently):

```bash
# Stage 1: simulate genealogy
krsim simulate \
    --Nh 3 --T 10 --K 2 --chi 0.18 --C 1 \
    --recombination_rate 7.4e-7 --random_seed 123 \
    --output_dir examples/output \
    --save_pairs

# Stage 2: add mutations
krsim mutate \
    --ts examples/output/ts_123.trees \
    --mutation_rate 1e-6 --random_seed 123 \
    --output_dir examples/output

# Stage 3: coalescence & mutation times
krsim coal \
    --ts examples/output/mutated_123.trees \
    --random_seed 123 --output_dir examples/output

# Stage 4: per-strain VCF
krsim vcf-strain \
    --ts examples/output/mutated_123.trees \
    --random_seed 123 --output_dir examples/output

# Stage 5: per-host VCF + strain proportions
krsim vcf-host \
    --strain_vcf examples/output/output_123_strain.vcf \
    --random_seed 123 --depth 10 --output_dir examples/output \
    --save_proportions
```

### Example parameters

| Flag | Value | Meaning |
|---|---|---|
| `--Nh` | 3 | Number of human hosts |
| `--T` | 10 | Transmission generations |
| `--K` | 2 | Parasite strains per human host (2K if superinfected) |
| `--chi` | 0.18 | Superinfection probability |
| `--C` | 1 | Oocysts per mosquito |
| `--recombination_rate` | 7.4e-7 | Per-bp recombination rate |
| `--mutation_rate` | 1e-6 | Per-bp mutation rate |
| `--random_seed` | 123 | Seed for reproducibility |

---

## Output files

All outputs land in `output/` named with the random seed (here `123`).
Reference copies live in `expected/` — diff against them to verify a clean run.
Log entries from all stages are appended to `KR_sim_123.log`.

### `sim_ind_node_map_123.csv`

Maps each individual (host) to its parasite strain node IDs.

| Column | Description |
|---|---|
| `sim_id` | Random seed identifying the simulation run |
| `ind_id` | Host individual ID |
| `node_id` | Haploid parasite genome (strain) carried by this host |

A host with 2 rows carries K strains; a host with 4 rows carries 2K strains (superinfected).

### `node_pairs_123.csv`

Pairs of nodes that are the two haploid products of the same meiosis event.

| Column | Description |
|---|---|
| `sim_id` | Simulation run ID |
| `node_id` | First node of the pair (smaller ID) |
| `paired_node_id` | Second node of the pair |

Each pair appears once. Comparing their genotypes across the VCF reveals recombination breakpoints.

### `coal_T_123.csv`

Per-interval coalescence times for every pair of sample nodes.

| Column | Description |
|---|---|
| `sample_1`, `sample_2` | The two sample node IDs |
| `ancestor` | MRCA node ID in this tree |
| `chromosome` | Chromosome number (1–14) |
| `left`, `right` | Start and end of the tree interval (bp, local to chromosome) |
| `coal_time` | Time to most recent common ancestor (generations) |

### `mut_T_123.csv`

Position and time of every mutation on the mutated tree sequence.

| Column | Description |
|---|---|
| `node_id` | Node the mutation was placed on |
| `chromosome` | Chromosome number (1–14) |
| `position` | Position within the chromosome (bp) |
| `time` | Age of the mutation (generations) |

### `output_123_strain.vcf` — per-strain genotypes

Standard VCF with one sample column per parasite strain (node).

- Genotype field: `0` = ref allele, `1` = alt allele
- Use with `sim_ind_node_map_123.csv` to group strains by host

### `output_123_host.vcf` — per-host genotypes

VCF aggregated to the individual (host) level using a **Dirichlet-Multinomial** model:

1. Strain proportions within each host are sampled from Dirichlet(α=1 per strain)
2. 10 sequencing reads are drawn from those proportions via Multinomial
3. Genotype is called from the read counts:
   - `0/0` — all reads carry the ref allele
   - `1/1` — all reads carry the alt allele
   - `0/1` — reads carry both alleles (within-host mixed infection detected)

### `host_proportions_123.csv`

Records the Dirichlet-sampled strain proportions and Multinomial read counts used to generate `output_123_host.vcf`.

| Column | Description |
|---|---|
| `sim_id` | Simulation run ID |
| `ind_id` | Host individual ID |
| `node_ids` | List of strain node IDs in this host |
| `proportions` | Dirichlet-sampled within-host proportion of each strain |
| `reads` | Multinomial read count assigned to each strain (sums to 10) |

---

## Verifying your run

```bash
diff examples/output/sim_ind_node_map_123.csv  examples/expected/sim_ind_node_map_123.csv
diff examples/output/node_pairs_123.csv         examples/expected/node_pairs_123.csv
diff examples/output/coal_T_123.csv             examples/expected/coal_T_123.csv
diff examples/output/mut_T_123.csv              examples/expected/mut_T_123.csv
diff examples/output/output_123_strain.vcf      examples/expected/output_123_strain.vcf
diff examples/output/output_123_host.vcf        examples/expected/output_123_host.vcf
diff examples/output/host_proportions_123.csv   examples/expected/host_proportions_123.csv
```

All diffs should be empty for a clean install with the same seed.
