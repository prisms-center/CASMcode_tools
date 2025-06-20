import argparse
import os
import sys

from casm.tools.shared import contexts

from .commands.search import make_search_parser
from .commands.write import make_write_parser


def make_parser():
    """Make a parser for the casm-map command line interface."""
    parser = argparse.ArgumentParser(
        description="CASM structure mapping CLI tool",
    )

    m = parser.add_subparsers(title="Select which method to use")
    make_search_parser(m)
    make_write_parser(m)

    return parser


def main(argv=None, working_dir=None):
    if argv is None:
        argv = sys.argv
    if working_dir is None:
        working_dir = os.getcwd()

    print("argv:", argv)

    parser = make_parser()

    # if "--desc" is in the arguments, print the description:
    if "--desc" in argv:

        if "search" in argv:
            from .commands.search import search_desc

            print(search_desc)
            return 0
        elif "write" in argv:
            from .commands.write import write_desc

            print(write_desc)
            return 0
        else:
            parser.print_help()
            return 1

    if len(argv) < 2:
        parser.print_help()
        return 1
    args = parser.parse_args(argv[1:])

    with contexts.working_dir(working_dir):
        code = args.func(args)

    return code
