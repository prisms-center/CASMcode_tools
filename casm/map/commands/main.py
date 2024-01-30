"""
This is the casm-map interactive mapping script.
"""
import argparse
from pathlib import Path

import casm.map.utils as utils

from casm.map.commands import apply, equiv, interp, search


def run_search(args):
    parent_path = Path(args.parent)
    child_path = Path(args.child)
    if not parent_path.is_file():
        print(f"missing parent at: {parent_path}")
        return
    elif not child_path.is_file():
        print(f"missing child at: {child_path}")
        return
    maps, parent, child = search.search(args)
    print(f"found {len(maps)} maps")
    utils.write_maps(maps, parent, child)


def run():
    """Run the interactive tool."""
    parser = argparse.ArgumentParser(
        description="Interface to libcasm.mapping utilities."
    )
    parser.add_argument(
        "--verbose", action="store_true", help="verbose output for debugging"
    )

    # choose method
    method = parser.add_subparsers(title="Select which method to use")
    apply_method = method.add_parser(
        "apply", help="apply mapping and construct resulting structure"
    )
    equiv_method = method.add_parser("equiv", help="find equivalent mappings")
    interp_method = method.add_parser(
        "interp", help="interpolate between two structures"
    )
    search_method = method.add_parser(
        "search", help="find a mapping between two structures"
    )

    # arguments for each method
    ## search
    search_method.set_defaults(func=run_search)
    search_method.add_argument("parent", help="path to the parent crystal structure")
    search_method.add_argument("child", help="path to the child crystal structure")
    search_method.add_argument(
        "--symmetrize",
        choices=["parent", "child", "both"],
        default=[],
        help="use Pymatgen to symmetrize the structure",
    )
    search_method.add_argument(
        "--symmetrize-if-necessary",
        choices=["parent", "child", "both"],
        default=[],
        help="use Pymatgen to symmetrize the structure if make_factor_group fails",
    )
    search_method.add_argument(
        "--primify",
        choices=["parent", "child", "both"],
        default=[],
        help="ensure that the structures are primitive",
    )
    search_method.add_argument(
        "--max-vol",
        type=int,
        default=None,
        help="maximum volume of the parent (after primifying if selected)",
    )
    search_method.add_argument(
        "--child-supercells",
        type=int,
        default=None,
        help="make supercells of the child up to this size",
    )
    search_method.add_argument(
        "--include-vacancies", action="store_true", help="allow vacancies"
    )
    search_method.add_argument(
        "--mask-occupants",
        action="store_true",
        help="treat all atoms as if they were the same species",
    )

    # parse
    args = parser.parse_args()
    return args


def main():
    args = run()
    args.func(args)


if __name__ == "__main__":
    main()
