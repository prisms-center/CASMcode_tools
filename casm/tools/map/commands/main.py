"""Implements ``casm-map ...``"""

import argparse
import os
import sys


def make_parser():
    """Make the ``casm-map ...`` argument parser

    Returns
    -------
    parser: argparse.ArgumentParser
        The argument parser for the `casm-map` program.

    """

    from .search import make_search_parser

    # from .search2 import make_search2_parser
    from .write import make_write_parser

    parser = argparse.ArgumentParser(
        description="CASM structure mapping CLI tool",
    )

    m = parser.add_subparsers(title="Select which method to use")
    make_search_parser(m)
    # make_search2_parser(m)
    make_write_parser(m)

    return parser


def main(argv=None, working_dir=None):
    """Implements ``casm-map ...``

    Parameters
    ----------
    argv : list of str, optional
        The command line arguments to parse. If None, uses `sys.argv`.
    working_dir : str, optional
        The working directory to use. If None, uses the current working directory.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    from casm.tools.shared import contexts

    if argv is None:
        argv = sys.argv
    if working_dir is None:
        working_dir = os.getcwd()

    print("argv:", argv)

    parser = make_parser()

    # if "--desc" is in the arguments, print the description:
    if "--desc" in argv:

        if "search" in argv:
            from .search import print_desc

            print_desc()
            return 0
        # elif "search2" in argv:
        #     from .search import print_desc
        #
        #     print_desc()
        #     return 0

        elif "write" in argv:
            from .write import print_desc

            print_desc()
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
