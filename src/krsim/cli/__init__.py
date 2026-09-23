"""The `krsim` command: one entry point, one subcommand per pipeline step."""

import argparse

from rich_argparse import RawDescriptionRichHelpFormatter

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

# subcommand name -> module providing HELP, add_arguments(parser), run(args)
# and EPILOG (outputs + example, shown at the end of its --help)
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


class HelpFormatter(RawDescriptionRichHelpFormatter):
    """Raw* keeps the line layout of each command's EPILOG (outputs + example);
    the wider flag column stops most flags pushing their description onto a new line."""

    # colour the "Outputs:" / "Example:" / "Needs ..." headings at the end of each command's help
    highlights = [*RawDescriptionRichHelpFormatter.highlights,
                  r"(?m)^(?P<section>Outputs:|Example:|Needs)"]
    styles = {**RawDescriptionRichHelpFormatter.styles, "argparse.section": "bold green"}

    def __init__(self, prog, **kwargs):
        kwargs.setdefault("max_help_position", 30)
        super().__init__(prog, **kwargs)

    def add_argument(self, action):
        super().add_argument(action)
        # argparse sizes the column without the extra indent of the command list,
        # which pushes the longest command name's description onto its own line
        if isinstance(action, argparse._SubParsersAction):
            self._action_max_length += self._indent_increment


def _add_help(parser):
    parser.add_argument("-h", "--help", action="help", help="Show this help and exit")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="krsim",
        usage="krsim <command> [options]",
        description="KR Malaria Transmission Simulator",
        epilog="Run `krsim <command> --help` to see a command's options.",
        formatter_class=HelpFormatter,
        add_help=False,
    )
    _add_help(parser)
    subparsers = parser.add_subparsers(title="commands", dest="command", required=True, metavar="<command>")

    for name, module in COMMANDS.items():
        sub = subparsers.add_parser(name, help=module.HELP, description=module.HELP,
                                    usage=f"krsim {name} [options]",
                                    epilog=getattr(module, "EPILOG", None),
                                    formatter_class=HelpFormatter,
                                    add_help=False)
        _add_help(sub)
        module.add_arguments(sub)
        # list the "required arguments" section before argparse's built-in optional one
        sub._action_groups.sort(key=lambda group: group.title != "required arguments")
        sub.set_defaults(func=module.run)

    # list the commands before the -h/--help option
    parser._action_groups.sort(key=lambda group: group.title != "commands")

    args = parser.parse_args(argv)
    args.func(args)
