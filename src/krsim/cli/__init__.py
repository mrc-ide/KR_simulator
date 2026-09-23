"""The `krsim` command: one entry point, one subcommand per pipeline step."""

import argparse

from rich_argparse import RichHelpFormatter

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
        formatter_class=RichHelpFormatter,
    )
    subparsers = parser.add_subparsers(dest="command", required=True, metavar="<command>")

    for name, module in COMMANDS.items():
        sub = subparsers.add_parser(name, help=module.HELP, description=module.HELP,
                                    formatter_class=RichHelpFormatter)
        module.add_arguments(sub)
        # list the "required arguments" section before argparse's built-in optional one
        sub._action_groups.sort(key=lambda group: group.title != "required arguments")
        sub.set_defaults(func=module.run)

    args = parser.parse_args(argv)
    args.func(args)
