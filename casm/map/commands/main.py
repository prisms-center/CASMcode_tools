"""
This is the casm-map interactive mapping script.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import libcasm.mapping.methods as mapmethods

import casm.map.utils as utils
from casm.map.commands import equiv, interp, search


def run_search(args):
    parent_path = Path(args.parent)
    child_path = Path(args.child)
    if not parent_path.is_file():
        print(f"missing parent at: {parent_path}")
        return
    elif not child_path.is_file():
        print(f"missing child at: {child_path}")
        return
    additional_data = []
    maps, parent, child = search.search(args)
    print(f"found {len(maps)} maps")
    if args.symmetry_adapted_strain:
        for m in maps:
            symmetry_adapted_strain = utils.strain_from_lattice_mapping(
                m.lattice_mapping()
            )
            additional_data.append(
                {"symmetry_adapted_strain": symmetry_adapted_strain.tolist()}
            )
    utils.write_maps(maps, parent, child, additional_data)


def run_interp(args):
    parent_path = Path(args.parent)
    child_path = Path(args.child)
    if not parent_path.is_file():
        print(f"missing parent at: {parent_path}")
        return
    elif not child_path.is_file():
        print(f"missing child at: {child_path}")
        return
    maps, parent, child = search.search(args)
    if Path("interp").exists():
        print("interp directory exists")
        return
    else:
        Path("interp").mkdir()
        steps = np.linspace(0, 1, args.steps)
        for i, m in enumerate(maps):
            p = Path("interp") / str(i)
            p.mkdir()
            structures = [
                mapmethods.make_mapped_structure(m.interpolated(s), child)
                for s in steps
            ]
            utils.write_structures(structures, p)


def parse_args(args):
    parser = argparse.ArgumentParser(
        description="Interface to libcasm.mapping utilities."
    )
    # parser.add_argument(
    #     "--verbose", action="store_true", help="verbose output for debugging"
    # )

    # choose method
    method = parser.add_subparsers(title="Select which method to use")
    deform_method = method.add_parser(
        "deform",
        help="apply mapping and construct resulting structure",
        description="Given a map and a structure, apply the deformation in the map.",
    )
    equiv_method = method.add_parser(
        "equiv",
        help="find equivalent mappings",
        description="Given a map between two structures, find all equivalent child structures.",
    )
    interp_method = method.add_parser(
        "interp",
        help="interpolate between two structures",
        description="Linearly interpolate along a mapping pathway between two structures.",
    )
    search_method = method.add_parser(
        "search",
        help="find a mapping between two structures",
        description="Given a parent and a child structure, find the best mapping from child to parent.",
    )

    # arguments for each method
    ## interpolate
    interp_method.set_defaults(func=run_interp)
    interp_method.add_argument(
        "parent",
        type=Path,
        help="path to the parent crystal structure (prim.json or POSCAR format)",
    )
    interp_method.add_argument(
        "child", type=Path, help="path to the child crystal structure (POSCAR format)"
    )
    interp_method.add_argument(
        "steps", type=int, help="number of interpolation steps (inclusive)"
    )
    interp_method.add_argument(
        "--xdatcar",
        action="store_true",
        help="write XDATCAR (append structures into one file)",
    )
    # interp_method.add_argument(
    #     "--symmetrize",
    #     choices=["parent", "child", "both"],
    #     default=[],
    #     help="use Pymatgen to symmetrize the structure",
    # )
    # interp_method.add_argument(
    #     "--symmetrize-if-necessary",
    #     choices=["parent", "child", "both"],
    #     default=[],
    #     help="use Pymatgen to symmetrize the structure if make_factor_group fails",
    # )
    interp_method.add_argument(
        "--primify",
        choices=["parent", "child", "both"],
        default=[],
        help="ensure that the structures are primitive",
    )
    interp_method.add_argument(
        "--max-vol",
        type=int,
        default=None,
        help="maximum volume of the parent (after primifying if selected)",
    )
    # interp_method.add_argument(
    #     "--child-supercells",
    #     type=int,
    #     default=None,
    #     help="make supercells of the child up to this size",
    # )
    # interp_method.add_argument(
    #     "--include-vacancies", action="store_true", help="allow vacancies"
    # )
    interp_method.add_argument(
        "--mask-occupants",
        action="store_true",
        help="treat all atoms as if they were the same species",
    )

    ## search
    search_method.set_defaults(func=run_search)
    search_method.add_argument(
        "parent",
        type=Path,
        help="path to the parent crystal structure (prim.json or POSCAR format)",
    )
    search_method.add_argument(
        "child", type=Path, help="path to the child crystal structure (POSCAR format)"
    )
    # search_method.add_argument(
    #     "--symmetrize",
    #     choices=["parent", "child", "both"],
    #     default=[],
    #     help="use Pymatgen to symmetrize the structure",
    # )
    # search_method.add_argument(
    #     "--symmetrize-if-necessary",
    #     choices=["parent", "child", "both"],
    #     default=[],
    #     help="use Pymatgen to symmetrize the structure if make_factor_group fails",
    # )
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
    # search_method.add_argument(
    #     "--child-supercells",
    #     type=int,
    #     default=None,
    #     help="make supercells of the child up to this size",
    # )
    # search_method.add_argument(
    #     "--include-vacancies", action="store_true", help="allow vacancies"
    # )
    search_method.add_argument(
        "--mask-occupants",
        action="store_true",
        help="treat all atoms as if they were the same species",
    )
    search_method.add_argument(
        "--symmetry-adapted-strain",
        action="store_true",
        help="write symmetry-adapted Hencky strain vector into mapping results",
    )

    # parse
    return parser.parse_args(args)


def main():
    args = parse_args(sys.argv[1:])
    if not hasattr(args, "func"):  # in case no input is passed to casm-map
        print("run 'casm-map --help' to print usage information")
        sys.exit()
    args.func(args)


if __name__ == "__main__":
    main()
