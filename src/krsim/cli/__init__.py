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
    strain_transmission_map,
    pairwise_strain_pedigree_map,
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
    "strain-transmission-map": strain_transmission_map,
    "pairwise-strain-pedigree-map": pairwise_strain_pedigree_map,
}


class HelpFormatter(RawDescriptionRichHelpFormatter):
    """Raw* keeps the line layout of each command's EPILOG (outputs + example);
    the wider flag column stops most flags pushing their description onto a new line."""

    # colour the "Outputs:" / "Example:" / "Needs ..." headings at the end of each command's help
    highlights = [*RawDescriptionRichHelpFormatter.highlights,
                  r"(?m)^(?P<section>Outputs:|Example:|Needs)"]
    styles = {**RawDescriptionRichHelpFormatter.styles, "argparse.section": "bold hot_pink"}

    def __init__(self, prog, **kwargs):
        kwargs.setdefault("max_help_position", 34)
        super().__init__(prog, **kwargs)

    def add_argument(self, action):
        super().add_argument(action)
        # argparse sizes the column without the extra indent of the command list,
        # which pushes the longest command name's description onto its own line
        if isinstance(action, argparse._SubParsersAction):
            self._action_max_length += self._indent_increment


class Parser(argparse.ArgumentParser):
    def format_help(self):
        # bold + underline the one-line description at the top; rich is imported here, not at
        # module level, so it only costs time when a help page is actually shown
        from rich.text import Text

        description = self.description
        if isinstance(description, str):
            self.description = Text(description, style="bold underline")
        try:
            # blank line above and below, so the page stands apart from the shell prompt
            # (error messages don't go through here and are unaffected)
            return "\n" + super().format_help() + "\n"
        finally:
            self.description = description


def _add_help(parser):
    parser.add_argument("-h", "--help", action="help", help="Show this help and exit")


def main(argv=None):
    parser = Parser(
        prog="krsim",
        usage="krsim <command> [options]",
        description="KR Malaria Transmission Simulator",
        epilog="Run `krsim <command> --help` to see a command's options.",
        formatter_class=HelpFormatter,
        add_help=False,
    )
    _add_help(parser)
    subparsers = parser.add_subparsers(title="commands", dest="command", required=True, metavar="<command>",
                                       parser_class=Parser)

    for name, module in COMMANDS.items():
        sub = subparsers.add_parser(name, prog=f"krsim {name}", help=module.HELP, description=module.HELP,
                                    usage=argparse.SUPPRESS,  # the Example at the end shows usage
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
