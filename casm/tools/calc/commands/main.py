"""Implements ``casm-calc ...``"""

import argparse
import os
import sys


def make_parser():
    """Make the ``casm-calc ...`` argument parser

    Returns
    -------
    parser: argparse.ArgumentParser
        The argument parser for the `casm-map` program.

    """

    from .status import make_status_subparser
    from .submit import make_submit_subparser
    from .vasp import make_vasp_subparser

    ### casm-calc ...
    parser = argparse.ArgumentParser(
        description="CASM calculation CLI tool",
    )
    c = parser.add_subparsers(title="Select which calculator to use")

    ### casm-calc ...
    make_vasp_subparser(c)
    make_status_subparser(c)
    make_submit_subparser(c)

    return parser


def main(argv=None, working_dir=None):
    """Implements ``casm-calc ...``

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

    parser = make_parser()

    # if "--desc" is in the arguments, print the description:
    if "--desc" in argv:

        if "vasp" in argv:
            from .vasp import print_desc

            print_desc(argv=argv)
            return 0
        elif "status" in argv:
            from .submit import print_desc

            print_desc(argv=argv)
            return 0

        elif "submit" in argv:

            from .submit import print_desc

            print_desc(argv=argv)

            return 0

        else:
            parser.print_help()
            return 1

    if len(argv) < 2:
        parser.print_help()
        return 1
    args = parser.parse_args(argv[1:])

    code = 0
    with contexts.working_dir(working_dir):
        code = args.func(args)

    return code
