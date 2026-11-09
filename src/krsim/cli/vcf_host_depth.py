"""krsim vcf-host-depth — Stage 5 (alternative): per-host VCF with read depth from real depth profiles."""

from krsim.cli._common import KrsimError, Run

HELP = "Stage 5 (alternative): per-host VCF, read depth from a zip of depth files"

EPILOG = """\
Each host's read depth follows one real sample's depth profile (mean and SD of
depth in every 500 bp window), taken from --depth_zip. The zip needs at least one
file per host; if it has more, files are picked at random (set by --random_seed,
so the same seed always gives each host the same file), each used at most once.

The depth files must come from reads mapped to the 3D7 reference used by
MalariaGEN Pf8 (PlasmoDB-54), the genome krsim simulates. Make the zip from
mosdepth's per-base output with `krsim depth-zip`, or use an example zip
(100 sWGA or 100 gDNA Pf8 samples) from the krsim GitHub releases.

Outputs:
  always  output_<seed>_host_depth.vcf
          host_depth_profiles_<seed>.csv     (which depth file each host got)
          host_depth_proportions_<seed>.csv  (each host's strain proportions)

Example:
  krsim vcf-host-depth --strain_vcf results/output_42_strain.vcf \\
      --ind_node_map_csv results/sim_ind_node_map_42.csv \\
      --depth_zip pf8_sWGA_100.zip --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--strain_vcf", type=str, required=True, metavar="FILE",
                          help="Strain VCF from `krsim vcf-strain` (output_<seed>_strain.vcf)")
    required.add_argument("--ind_node_map_csv", type=str, required=True, metavar="FILE",
                          help="Host-strain map from `krsim simulate` (sim_ind_node_map_<seed>.csv)")
    required.add_argument("--depth_zip", type=str, required=True, metavar="FILE",
                          help="Zip of <sample>.window_stats.tsv depth files, at least one per host "
                               "(made by `krsim depth-zip`)")
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
            depth_zip=args.depth_zip,
            output_dir=args.output_dir,
            random_seed=args.random_seed,
        )
    except ValueError as exc:
        raise KrsimError(str(exc)) from exc
    for path in paths:
        run.wrote(path)

    run.done()
