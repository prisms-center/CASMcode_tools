import argparse
import os
import pathlib
import sys

from casm.tools.shared import contexts
from casm.tools.shared.io import read_structure, write_structure


def get_parent_format(args):
    if args.parent_format is not None:
        return args.parent_format
    if args.format is not None:
        return args.format
    return None


def get_child_format(args):
    if args.child_format is not None:
        return args.child_format
    if args.format is not None:
        return args.format
    return None


def make_parser():
    """Make a parser for the casm-map command line interface."""
    parser = argparse.ArgumentParser(
        description="CASM structure mapping CLI tool",
    )
    parser.add_argument("parent", type=pathlib.Path, help="Parent structure file")
    parser.add_argument("child", type=pathlib.Path, help="Child structure file")
    parser.add_argument(
        "--format",
        type=str,
        default=None,
        help=(
            "Set the format for reading the structure files. One of 'vasp', 'casm', "
            "or a format recognized by `ase.io.read`. "
            "The default is that no suffix or a '.vasp' suffix will read the file as a "
            "VASP POSCAR file, and a '.json' or '.casm' suffix will read the file as a "
            "CASM Structure JSON file. Otherwise, the file is read using ASE's "
            "`ase.io.read` method, if ASE is installed."
        ),
    )
    parser.add_argument(
        "--parent-format",
        type=str,
        default=None,
        help=(
            "Set the format for reading the parent structure file "
            "(overrides --format). "
        ),
    )
    parser.add_argument(
        "--child-format",
        type=str,
        default=None,
        help=(
            "Set the format for reading the child structure file "
            "(overrides --format). "
        ),
    )

    return parser


def main(argv=None, working_dir=None):
    if argv is None:
        argv = sys.argv
    if working_dir is None:
        working_dir = os.getcwd()

    parser = make_parser()
    if len(argv) < 2:
        parser.print_help()
        return 1
    args = parser.parse_args(argv[1:])

    with contexts.working_dir(working_dir):

        parent = read_structure(path=args.parent, format=get_parent_format(args))
        print("Parent:")
        print(parent)
        print()
        write_structure(
            casm_structure=parent,
            path=pathlib.Path("parent.xyz"),
            format=None,
        )

        child = read_structure(path=args.child, format=get_child_format(args))
        print("Child:")
        print(child)
        print()
        write_structure(
            casm_structure=parent,
            path=pathlib.Path("child.xyz"),
            format=None,
        )

    return 0
