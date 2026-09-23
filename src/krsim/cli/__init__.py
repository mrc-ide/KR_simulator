"""The `krsim` command: one entry point, one subcommand per pipeline step."""

import argparse

from krsim.cli import (
    simulate,
    mutate,
    coal,
    vcf_strain,
    vcf_host,
    relate_strain,
    relate_host,
    transmission_map,
    pairwise_coalescence,
)

# subcommand name -> module providing HELP, add_arguments(parser) and run(args)
COMMANDS = {
    "simulate": simulate,
    "mutate": mutate,
    "coal": coal,
    "vcf-strain": vcf_strain,
    "vcf-host": vcf_host,
    "relate-strain": relate_strain,
    "relate-host": relate_host,
    "transmission-map": transmission_map,
    "pairwise-coalescence": pairwise_coalescence,
}


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="krsim",
        description="KR Malaria Transmission Simulator",
    )
    subparsers = parser.add_subparsers(dest="command", required=True, metavar="<command>")

    for name, module in COMMANDS.items():
        sub = subparsers.add_parser(name, help=module.HELP, description=module.HELP)
        module.add_arguments(sub)
        sub.set_defaults(func=module.run)

    args = parser.parse_args(argv)
    args.func(args)
