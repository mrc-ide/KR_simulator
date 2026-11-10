"""The `krsim` command: one entry point, one subcommand per pipeline step."""

import argparse
import sys

from rich_argparse import RawDescriptionRichHelpFormatter

from krsim.cli._common import KrsimError, Run

from krsim.cli import (
    simulate,
    mutate,
    coal,
    vcf_strain,
    vcf_host_depth,
    relate_strain,
    relate_host,
    strain_transmission_map,
    pairwise_strain_pedigree_map,
    depth_zip,
)

# subcommand name -> module providing HELP, add_arguments(parser), run(args)
# and EPILOG (outputs + example, shown at the end of its --help)
COMMANDS = {
    # main commands: the simulator itself, in pipeline order
    "simulate": simulate,
    "mutate": mutate,
    "vcf-strain": vcf_strain,
    "vcf-host-depth": vcf_host_depth,
    # helper commands: optional extras, not needed for the VCFs (listed in their own section, names in pink)
    "coal-mut": coal,
    "relate-strain": relate_strain,
    "relate-host": relate_host,
    "strain-transmission-map": strain_transmission_map,
    "pairwise-strain-pedigree-map": pairwise_strain_pedigree_map,
    "depth-zip": depth_zip,
}
HELPERS = {"coal-mut", "relate-strain", "relate-host", "strain-transmission-map", "pairwise-strain-pedigree-map",
           "depth-zip"}


class HelpFormatter(RawDescriptionRichHelpFormatter):
    """Raw* keeps the line layout of each command's EPILOG (outputs + example);
    the wider flag column stops most flags pushing their description onto a new line."""

    # colour the "Outputs:" / "Example:" / "Needs ..." headings at the end of each command's help
    highlights = [*RawDescriptionRichHelpFormatter.highlights,
                  r"(?m)^(?P<section>Outputs:|Example:|Needs)"]
    styles = {**RawDescriptionRichHelpFormatter.styles, "argparse.section": "bold hot_pink"}

    # section headings as written ("Main commands"), not Title Cased
    group_name_formatter = staticmethod(lambda title: title[:1].upper() + title[1:])

    def __init__(self, prog, **kwargs):
        kwargs.setdefault("max_help_position", 34)
        super().__init__(prog, **kwargs)

    def add_argument(self, action):
        super().add_argument(action)
        # argparse sizes the column without the extra indent of the command list,
        # which pushes the longest command name's description onto its own line
        if isinstance(action, argparse._SubParsersAction):
            self._action_max_length += self._indent_increment

    def _rich_format_action_invocation(self, action):
        text = super()._rich_format_action_invocation(action)
        # a helper command in the list (no --option strings): pink, not the usual blue
        if not action.option_strings and action.dest in HELPERS:
            text.stylize("bold hot_pink", 0, len(text))
        return text


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
    # -h/--help still works; it isn't listed, since anyone reading the help already used it
    parser.add_argument("-h", "--help", action="help", help=argparse.SUPPRESS)


def main(argv=None):
    parser = Parser(
        prog="krsim",
        usage="krsim <command> [options]",
        description="KR Malaria Transmission Simulator",
        epilog="Run `krsim <command> --help` to see a command's options.\n"
               "Every command records what it did in <output_dir>/KR_sim_<seed>.log.",
        formatter_class=HelpFormatter,
        add_help=False,
    )
    _add_help(parser)
    subparsers = parser.add_subparsers(dest="command", required=True, metavar="<command>",
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

    # two sections in `krsim --help`: argparse lists every subcommand under one heading, so their help lines are
    # moved into two plain groups and the subcommand block itself is hidden; parsing is unchanged
    main_group = parser.add_argument_group("main commands")
    helper_group = parser.add_argument_group("helper commands (optional, not needed for the VCFs)")
    for line in list(subparsers._choices_actions):
        (helper_group if line.dest in HELPERS else main_group)._group_actions.append(line)
    subparsers._choices_actions.clear()
    subparsers.help = argparse.SUPPRESS

    # list the commands before the -h/--help option
    order = {main_group.title: 0, helper_group.title: 1}
    parser._action_groups.sort(key=lambda group: order.get(group.title, 2))

    args = parser.parse_args(argv)
    try:
        args.func(args)
    except Exception as exc:
        _report_error(exc)
        sys.exit(1)


def _report_error(exc):
    """One line on screen; the traceback goes to the run's log file (if it got that far)."""
    if isinstance(exc, FileNotFoundError) and exc.filename:
        message = f"input file not found: {exc.filename}"
    elif isinstance(exc, KrsimError):
        message = str(exc)
    else:
        message = f"{type(exc).__name__}: {exc}"

    if Run.current is not None:
        Run.current.fail(message)
    else:
        print(f"\n  error: {message}\n", file=sys.stderr)
