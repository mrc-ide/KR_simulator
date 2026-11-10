"""krsim vcf-host-depth — Stage 5: per-host VCF with realistic read depth."""

import json
from importlib import resources

from krsim.cli._common import KrsimError, Run

HELP = "Stage 5: per-host VCF, realistic read depth from depth files or a model"


def _defaults_table():
    """The --depth_model defaults, read from the small JSON that ships with krsim (no numpy needed)."""
    try:
        with resources.files("krsim.data").joinpath("depth_model_defaults.json").open() as f:
            d = json.load(f)
    except (OSError, ValueError):
        return "  (defaults file missing: reinstall krsim)\n"

    def nums(v):
        return ", ".join(f"{x:.4g}" for x in v) if v else "none"

    lines = []
    for t in ("sWGA", "gDNA"):
        s = d[t]
        lines += [
            f"  {t}: {s['model']}",
            f"    ar           {nums(s['ar'])}",
            f"    ma           {nums(s['ma'])}",
            f"    mu {s['mu']:.4g}, scale {s['scale']:.4g}, nu {s['nu']:.4g}, offset {s['offset']:.4g}",
            f"    depth_mean   {s['depth_mean']:.1f}   depth_spread {s['depth_spread']:.3f}",
            f"    cap          {s['cap']:.1f}x the host's typical depth",
            f"    end_width    {s['end_width']:,.0f} bp at end_factor {s['end_factor']:.3g}",
        ]
    return "\n".join(lines) + "\n"


EPILOG = f"""\
Give one depth source:

--depth_zip: each host's depth follows one real sample's depth profile (mean
and SD of depth in every 500 bp window). The zip needs at least one file per
host; if it has more, files are picked at random (set by --random_seed, so the
same seed always gives each host the same file), each used at most once. Files
must come from reads mapped to the 3D7 reference used by MalariaGEN Pf8
(PlasmoDB-54). Make the zip with `krsim depth-zip`, or use an example zip
(100 sWGA or 100 gDNA Pf8 samples) from the krsim GitHub releases.

--depth_model sWGA|gDNA: depth is simulated from a model fitted to Pf8 (sWGA
AR(8), gDNA ARMA(2,2), both with Student t kicks). Each host's typical depth
is drawn once from a log-normal (--depth_mean, --depth_spread; spread 0 gives
every host the same typical depth); its own depth track along the genome
comes from the model, capped, and for sWGA lowered near chromosome ends. The
read depth at each site is Poisson around its window's depth. Any value below
can be changed in a JSON file given to --depth_params, e.g.
{{"nu": 3.0, "cap": 1000}}; values not listed keep their defaults.

Flat depth (a test control, not realistic): --depth_spread 0 with
--depth_params {{"scale": 0, "end_width": 0}} gives every host and every window
--depth_mean; each site's reads are still Poisson, so a site can get 0 reads.

Defaults (median over QC-pass Pf8 samples):
{_defaults_table()}
Outputs:
  always  output_<seed>_host_depth.vcf      (header lists the depth source and values used)
          host_depth_profiles_<seed>.csv     (each host's depth file, or its typical depth)
          host_depth_proportions_<seed>.csv  (each host's strain proportions)

Examples:
  krsim vcf-host-depth --strain_vcf results/output_42_strain.vcf \\
      --ind_node_map_csv results/sim_ind_node_map_42.csv \\
      --depth_zip pf8_sWGA_100.zip --random_seed 42 --output_dir results/

  krsim vcf-host-depth --strain_vcf results/output_42_strain.vcf \\
      --ind_node_map_csv results/sim_ind_node_map_42.csv \\
      --depth_model sWGA --depth_mean 30 --random_seed 42 --output_dir results/
"""


def add_arguments(parser):
    required = parser.add_argument_group("required arguments")

    required.add_argument("--strain_vcf", type=str, required=True, metavar="FILE",
                          help="Strain VCF from `krsim vcf-strain` (output_<seed>_strain.vcf)")
    required.add_argument("--ind_node_map_csv", type=str, required=True, metavar="FILE",
                          help="Host-strain map from `krsim simulate` (sim_ind_node_map_<seed>.csv)")
    required.add_argument("--random_seed", type=int, required=True, metavar="SEED",
                          help="Random seed")
    required.add_argument("--output_dir", type=str, required=True, metavar="DIR",
                          help="Output directory")

    source = parser.add_argument_group("depth source (give one)")
    source.add_argument("--depth_zip", type=str, metavar="FILE",
                        help="Zip of <sample>.window_stats.tsv depth files, at least one per host "
                             "(made by `krsim depth-zip`)")
    source.add_argument("--depth_model", type=str, choices=["sWGA", "gDNA"], metavar="TYPE",
                        help="Simulate depth from the Pf8 model for this sample type: sWGA or gDNA")

    model = parser.add_argument_group("--depth_model options")
    model.add_argument("--depth_mean", type=float, metavar="X",
                       help="Centre of the hosts' typical depths (default: Pf8, see below)")
    model.add_argument("--depth_spread", type=float, metavar="S",
                       help="SD of log typical depth between hosts; 0 = every host at --depth_mean "
                            "(default: Pf8)")
    model.add_argument("--depth_params", type=str, metavar="FILE",
                       help="JSON file changing any other model value, e.g. {\"nu\": 3.0}")


def run(args):
    if (args.depth_zip is None) == (args.depth_model is None):
        raise KrsimError("give one depth source: --depth_zip FILE or --depth_model sWGA|gDNA")
    model_only = [f"--{k}" for k in ("depth_mean", "depth_spread", "depth_params") if getattr(args, k) is not None]
    if args.depth_zip is not None and model_only:
        raise KrsimError(f"{', '.join(model_only)} can only be used with --depth_model, not --depth_zip")

    # imported here rather than at module level so `krsim --help` stays fast
    # (these pull in numpy/pandas/scipy)
    from krsim.output.vcf_host_depth import process_host_depth_vcf
    from krsim.output.depth_model import model_settings

    run = Run("vcf-host-depth", args.output_dir, args.random_seed)

    try:
        settings = None
        if args.depth_model is not None:
            settings = model_settings(args.depth_model, args.depth_mean, args.depth_spread, args.depth_params)
            run.detail(f"depth model settings: {settings}")
        paths = process_host_depth_vcf(
            sim_id=args.random_seed,
            input_vcf=args.strain_vcf,
            ind_node_map_csv=args.ind_node_map_csv,
            output_dir=args.output_dir,
            random_seed=args.random_seed,
            depth_zip=args.depth_zip,
            depth_model=settings,
        )
    except ValueError as exc:
        raise KrsimError(str(exc)) from exc
    for path in paths:
        run.wrote(path)

    run.done()
