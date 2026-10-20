"""krsim vcf-host-depth — Stage 5 (alternative): per-host VCF with read depth from real Pf8 samples."""

from krsim.cli._common import KrsimError, Run

HELP = "Stage 5 (alternative): per-host VCF, read depth from real Pf8 samples"

EPILOG = """\
Each host's depth follows one real MalariaGEN Pf8 sample (one --samples ID per
host, in order). A sample's depth profile is made the first time it is used and
kept in <pf8_dir>/output/; making one needs mosdepth (conda install -c bioconda
mosdepth) and downloads its CRAM (0.5-11 GB, deleted afterwards).

Outputs:
  always  output_<seed>_host_depth.vcf
          host_depth_profiles_<seed>.csv     (which Pf8 sample each host's depth came from)
          host_depth_proportions_<seed>.csv  (each host's strain proportions)

Example:
  krsim vcf-host-depth --strain_vcf results/output_42_strain.vcf \\
      --ind_node_map_csv results/sim_ind_node_map_42.csv \\
      --samples SPT52544 PF1066-Cx RCN08168 --pf8_dir ~/pf8 \\
      --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--strain_vcf", type=str, required=True, metavar="FILE",
                          help="Strain VCF from `krsim vcf-strain` (output_<seed>_strain.vcf)")
    required.add_argument("--ind_node_map_csv", type=str, required=True, metavar="FILE",
                          help="Host-strain map from `krsim simulate` (sim_ind_node_map_<seed>.csv)")
    required.add_argument("--samples", type=str, nargs="+", required=True, metavar="ID",
                          help="Pf8 sample IDs, exactly one per host (1st sample -> 1st host)")
    required.add_argument("--pf8_dir", type=str, required=True, metavar="DIR",
                          help="Folder for Pf8 data (input/) and depth profiles (output/)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")


def run(args):
    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas)
    from krsim.output.vcf_host_depth import process_host_depth_vcf

    run = Run("vcf-host-depth", args.output_dir, args.random_seed)

    try:
        paths = process_host_depth_vcf(
            sim_id=args.random_seed,
            input_vcf=args.strain_vcf,
            ind_node_map_csv=args.ind_node_map_csv,
            samples=args.samples,
            pf8_dir=args.pf8_dir,
            output_dir=args.output_dir,
            random_seed=args.random_seed,
        )
    except ValueError as exc:
        raise KrsimError(str(exc)) from exc
    for path in paths:
        run.wrote(path)

    run.done()
