"""Implements ``casm-map search ...``"""

import argparse
import pathlib


# <-- max width = 80 characters                            --> #
################################################################
def print_desc():
    desc = """
# The `casm-map search` command:

## Method

The `casm-map search` command is intended for finding mappings 
between two crystal structures that have the same stoichiometry, 
and where no sites that allow vacancies or alloying exist in the
parent structure.

The `casm-map search` command reads a parent and child structure 
file, validates the structures have the same stoichiometry, 
searches for mappings, and writes the results to a 
`mappings.json` file in a specified results directory. 

Mappings are found using the search method described in Ref. 
[2], applied to the case that the parent structure specifies one
and only one atom type occupying each site. In summary, this does
the following:

1. Make superstructures of the child structure, for a specified 
   range of sizes. The default includes only the minimum size 
   superstructures which have the same number of atoms as a 
   superstructure of the parent structure.
2. For each child superstructure, make superstructures of the 
   parent structure which have the same number of atoms as the 
   child superstructure.
3. For each pair of parent and child superstructures, find the 
   `lattice_k_best` lowest lattice cost mappings. 
4. For each lattice mapping, use atoms from the minority type to 
   generate a minimal set of trial translations. For each trial 
   translation, find the best structure mapping and add it as a 
   "node" in the search queue.
5. While there are nodes in the search queue, find next-best 
   structure mappings, keeping the `k_best` mappings with the 
   lowest total mapping cost.
6. Deduplicate the results by comparing interpolated structures 
   on the path between the parent and mapped child.


## Results

Results are written to a JSON file in the specified results 
directory and a summary table is printed to the console.

Results are written to:

    <results_dir>/
    ├── mappings.json
    └── options_history.json


The `mappings.json` output file contains:

    "parent": libcasm.xtal.Structure
        The parent structure.
    "child": libcasm.xtal.Structure
        The child structure.
    "mappings": list[libcasm.mapping.info.ScoredStructureMapping]
        The list of scored structure mappings found between a 
        superstructure of the parent and a superstructure of the 
        child.
    "uuids": list[str]
        A list of UUIDs for the mappings, one per mapping.

The `options_history.json` output file is a JSON list with the 
history of options used for the search. When `casm-map search` 
is re-run with the `--merge` option new results are merged with 
existing results and the options used are appended to the list 
in this file.


## Mapping non-primitive structures

For a general mapping search, it is recommended to use primitive 
cells of the structures being mapped. The program will continue 
if the input structures are not primitive, but it will also 
print a notice and write files named `parent.primitive.json` and 
`child.primitive.json` containing the primitive structures. 
These can be used to re-run the search command with the 
primitive structures.


## Mapping relaxations from known starting structures

When mapping a relaxed structure (child) from a known starting 
structure (parent), the `--fix-parent` option can be used to 
skip searching over parent superstructures and lattice 
reorientations. The deformation gradient is still calculated and
atom mapping is still performed.


## Forcing and suppressing atom mappings

The `--forced-on` option can be used to constrain the atom
mapping search and force the mapping of specific atoms from
the child structure to specific sites in the parent structure.

The `--forced-off` option can be used to suppress the mapping
of specific atoms from the child structure to specific atoms in
the parent structure.


## Mapping relaxations with fixed atoms

By default, mean displacements are removed from the atom 
mappings. For relaxations where some atom positions are fixed, 
it may be desirable to suppress this with the
`--no-remove-mean-displacement` option. If this option is 
provided, then `--forced-on` must be specify at least one atom 
mapping. A trial translation resulting in an displacement of
zero is generated for each forced mapping.


## Parameters

Total mapping options:

--max-n-atoms: Optional[int]=None
    The maximum number of atoms in superstructures to include
    in the search. By default, the least common multiple of 
    the number of atoms in the child and parent structures is
    used.
--min-n-atoms: int=1
    The minimum number of atoms in superstructures to include
    in the search.
--min-total-cost: float=0.0
    Only mappings with a total cost greater than or equal to 
    this value are included in the final results.
--max-total-cost: float=0.3
    Only mappings with a total cost less than or equal to this 
    value are included in the final results. 
--k-best: int=100
    Keep the `k_best` mappings with lowest total cost that also 
    satisfy the min/max total cost criteria. Approximate ties 
    with the current `k_best`-ranked result are also kept.
--cost-tol: float=1e-5
    The cost tolerance for approximate ties. If the total cost 
    of a mapping is within this value of the current 
    `k_best`-ranked result, it is also kept in the results.
--lattice-cost-weight: float=0.5
    The fraction of the total cost that is due to the lattice 
    mapping cost. The remaining fraction is due to the atom 
    mapping cost.
--iso-cost: bool=False
    If given, use isotropic strain and displacement costs for
    the lattice and atom mapping costs, respectively. If not
    given, use symmetry-breaking strain and displacement costs.
    To set separately, use `--iso-strain-cost` or 
    `--iso-disp-cost`.
--no-remove-mean-displacement: bool=False
    If given, do not remove the mean displacement from the 
    structure mappings. By default, the mean displacement is
    removed from the structure mappings.
--fix-parent: bool=False
    If given, skip searching over parent superstructures and 
    lattice reorientations. The deformation gradient is still 
    calculated and atom mapping is still performed. The parent
    and child are required to have the same number of atoms.
    

Lattice mapping options:

--iso-strain-cost: bool=False
    If given, use the isotropic strain cost for the lattice 
    mapping cost.
--min-lattice-cost: float=0.0
    Only lattice mappings with a lattice cost greater than or 
    equal to this value are used to find structure mappings.
--max-lattice-cost: float=0.6
    Only lattice mappings with a lattice cost less than or equal
    to this value are used to find structure mappings.
--lattice-k-best: int=10
    Use the `lattice_k_best` lattice mappings with lowest 
    lattice cost (subject to the `min_lattice_cost` / 
    `max_lattice_cost` limits) for each parent/child 
    superstructure pair to find structure mappings.
--lattice-reorientation-range: int=1
    The absolute value of the maximum element in the lattice 
    mapping reorientation matrix, N. This determines how many 
    equivalent lattice vector reorientations are checked. 
    Increasing the value results in more checks. The value 1 is 
    generally expected to be sufficient because reduced cell 
    lattices are compared internally.

Atom mapping options:

--iso-disp-cost: bool=False
    If given, use the isotropic displacement cost method for 
    the atom mapping cost.
--forced-on: Optional[str]=None
    If given, a JSON representation of a list[tuple[int,int]]
    giving forced atom mapping. The first element of each tuple 
    is a parent structure site index and the second element is a
    child structure atom index. Indices start from 0. For 
    example, `--forced-on "[[0,0],[2,3]]"` specifies that the 
    first atom in the child structure (index 0) must be mapped 
    to the first site (index 0) in the parent structure, and the
    fourth atom (index 3) in the child structure must be mapped 
    to the third site (index 2) in the parent structure.
--forced-off: Optional[str]=None
    If given, a JSON representation of a list[tuple[int,int]] 
    giving suppressed atom mappings. The first element of each 
    tuple is a parent structure site index and the second 
    element is a child structure atom index. Indices start from 
    0. For example, `--forced-off "[[0,2],[0,3]]"` specifies 
    that the third atom (index 2) in the child structure must 
    not be mapped to the first site (index 0) in the parent 
    structure, and the fourth atom (index 3) in the child 
    structure must also not be mapped to the first site 
    (index 0) in the parent structure.
    
Deduplication options:

--dedup-interp-factors: Optional[str]=None
    If given, a JSON representation of a list[float] giving 
    interpolation factors to use for deduplication. A value of
    0.0 corresponds to the parent structure and a value of
    1.0 corresponds to the mapped child structure. If None, the 
    default, equivalent to ``"[0.5,1.0]"``, is used.

Input options:

--format: Optional[str]=None
    Specify the format for reading the structure files. If not 
    specified, the format is inferred from the file suffix. 
    Supported formats include 'vasp', 'casm', and any format 
    recognized by `ase.io.read` if ASE is installed.
--parent-format: Optional[str]=None
    Same as `--format`, but overrides it to specify the format 
    for reading the parent structure file.
--child-format: Optional[str]=None
    Same as `--format`, but overrides it to specify the format 
    for reading the child structure file.
--options: Optional[pathlib.Path]=None
    Path to a JSON file containing options to use instead of
    reading them from command line arguments. See (TODO) for
    the format of the options file.
    
Output options:

--results-dir: str="results"
    Directory where results are written. A new directory will be 
    created. If the directory already exists, the program exits 
    with an error unless --merge is used.
--merge: bool=False
    If given, read existing results and merge new results. If 
    not given, the program will exit with an error if the 
    results directory already exists.

Additional options:

--next: bool=False
    Shortcut which reads the last used options from the results
    directory (if results exist) and increases `--max-n-atoms`
    and `--min-n-atoms` to the next greatest common multiple of 
    the number of atoms in the child and parent structures while 
    keeping all other options the same. New results are merged 
    with existing results. If no results exist, only the least 
    common multiple of the number of atoms in the child and
    parent structures is searched.
    

## Citing

A suggested way to cite this program is as follows:

"Structure mappings were found by the method of Thomas et al. 
[1] using the `casm-map` program [2] provided by CASM [3]."

## References

[1] J. C. Thomas, A. R. Natarajan, and A. Van der Ven, Comparing 
    crystal structures with symmetry and geometry, npj 
    Computational Materials, 7 (2021), 164.
[2] B. Puchala, J. Thomas, and A. Van der Ven, "casm-map...".
[3] B. Puchala, J. C. Thomas, A. R. Natarajan, J. G. Goiri, 
    S. S. Behara, J. L. Kaufman, A. Van der Ven, CASM—A software
    package for first-principles based study of multicomponent 
    crystalline solids, Computational Materials Science 217 
    (2023) 111897.

"""
    print(desc)


_other_desc = """
--parent-superstructure: Optional[pathlib.Path]=None
    A file containing the parent superstructure. If given, skip 
    searching over superstructures and lattice reorientations. 
    The deformation gradient is still calculated and atom 
    mapping is still performed. The specified file must be a 
    valid structure file, but only the lattice is actually used.
--output-format: Optional[str]=None
    The format for writing structure files. If not specified, 
    the format is inferred from the child format. Supported 
    formats include 'vasp', 'casm', and any format recognized by
    `ase.io.write` if ASE is installed.
"""


def _get_parent_format(args):
    if args.parent_format is not None:
        return args.parent_format
    if args.format is not None:
        return args.format
    return None


def _get_child_format(args):
    if args.child_format is not None:
        return args.child_format
    if args.format is not None:
        return args.format
    return None


def run_search(args):
    """Implements ``casm-map search ...``

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """

    import math
    import sys

    import libcasm.configuration as casmconfig
    import libcasm.xtal as xtal
    from casm.tools.map import (
        StructureMappingSearch,
        StructureMappingSearchOptions,
    )
    from casm.tools.shared.json_io import read_optional, read_required
    from casm.tools.shared.structure_io import read_structure

    if args.desc:
        print_desc()
        return 0

    if args.prim:
        parent_prim = casmconfig.Prim.from_dict(data=read_required(args.prim))
        print("ParentPrim:")
        print(parent_prim)
        print()
    else:
        parent_prim = None

    parent = read_structure(path=args.parent, format=_get_parent_format(args))
    print("Parent:")
    print(parent)
    print()

    child = read_structure(path=args.child, format=_get_child_format(args))
    print("Child:")
    print(child)
    print()

    if args.options is not None:
        data = read_required(args.options)
        opt = StructureMappingSearchOptions.from_dict(data)
    else:

        lattice_mapping_cost_method = "symmetry_breaking_strain_cost"
        atom_mapping_cost_method = "symmetry_breaking_disp_cost"
        if args.iso_cost:
            lattice_mapping_cost_method = "isotropic_strain_cost"
            atom_mapping_cost_method = "isotropic_disp_cost"
        if args.iso_strain_cost:
            lattice_mapping_cost_method = "isotropic_strain_cost"
        if args.iso_disp_cost:
            atom_mapping_cost_method = "isotropic_disp_cost"

        opt = StructureMappingSearchOptions(
            max_n_atoms=args.max_n_atoms,
            min_n_atoms=args.min_n_atoms,
            child_transformation_matrix_to_super_list=None,
            parent_transformation_matrix_to_super_list=None,
            total_min_cost=args.min_total_cost,
            total_max_cost=args.max_total_cost,
            total_k_best=args.k_best,
            no_remove_mean_displacement=args.no_remove_mean_displacement,
            fix_parent=args.fix_parent,
            lattice_cost_weight=args.lattice_cost_weight,
            cost_tol=args.cost_tol,
            lattice_mapping_min_cost=args.min_lattice_cost,
            lattice_mapping_max_cost=args.max_lattice_cost,
            lattice_mapping_k_best=args.lattice_k_best,
            lattice_mapping_reorientation_range=args.lattice_reorientation_range,
            lattice_mapping_cost_method=lattice_mapping_cost_method,
            atom_mapping_cost_method=atom_mapping_cost_method,
            forced_on=args.forced_on,
            forced_off=args.forced_off,
            deduplication_interpolation_factors=args.dedup_interp_factors,
        )

    # --merge option
    merge = args.merge

    # --next option:
    #     run search with next greatest common multiple number of atoms and
    #     all other options the same; set merge=True
    if args.next:
        n_atoms_parent = len(parent.atom_type())
        n_atoms_child = len(child.atom_type())
        n_atoms_lcm = math.lcm(n_atoms_parent, n_atoms_child)

        results_dir = "results" if args.results_dir is None else args.results_dir

        data = read_optional(results_dir / "options_history.json", default=[])
        options_history = [
            StructureMappingSearchOptions.from_dict(data=x) for x in data
        ]
        last_max_n_atoms = None
        if len(options_history):
            if options_history[-1].max_n_atoms is None:
                last_max_n_atoms = n_atoms_lcm
            else:
                last_max_n_atoms = options_history[-1].max_n_atoms

        if last_max_n_atoms is None:
            next_max_n_atoms = n_atoms_lcm
        else:
            next_max_n_atoms = 0
            while next_max_n_atoms <= last_max_n_atoms:
                next_max_n_atoms += n_atoms_lcm

        opt.max_n_atoms = next_max_n_atoms
        opt.min_n_atoms = next_max_n_atoms
        merge = True

        print()
        print(
            f"""
--next: 
    Expanding search to next greatest common multiple number
    of atoms ({next_max_n_atoms} atoms) and merging results.
"""
        )
        sys.stdout.flush()

    print("Options:")
    print(xtal.pretty_json(opt.to_dict()))
    sys.stdout.flush()

    f = StructureMappingSearch(opt=opt)
    code = f(
        parent=parent,
        parent_prim=parent_prim,
        child=child,
        results_dir=args.results_dir,
        merge=merge,
    )

    return code


def _print(*x):
    import sys

    print(*x)
    sys.stdout.flush()


def _validate_forced_on(value):
    import json

    type_exception = argparse.ArgumentTypeError(
        "The --forced-on option must be a JSON list[tuple[int,int]]."
    )

    value_exception = argparse.ArgumentTypeError(
        "For the --forced-on option, "
        "all parent site indices (first index in each pair) must be unique, "
        "and all child atom indices (second index in each pair) must be unique."
    )

    try:
        forced_on = json.loads(value)
    except json.JSONDecodeError:
        raise type_exception

    # Ensure that forced_on is a list of tuples[int, int]
    if not isinstance(forced_on, list):
        raise type_exception

    for item in forced_on:
        if not isinstance(item, list) or len(item) != 2:
            raise type_exception
        if not isinstance(item[0], int) or not isinstance(item[1], int):
            raise type_exception

    # check that all first indices are unique:
    seen_indices = set()
    for x in forced_on:
        if x[0] in seen_indices:
            raise value_exception
        seen_indices.add(x[0])

    # check that all second indices are unique:
    seen_indices = set()
    for x in forced_on:
        if x[1] in seen_indices:
            raise value_exception
        seen_indices.add(x[1])

    return {x[0]: x[1] for x in forced_on}


def _validate_forced_off(value):
    import json

    exception = argparse.ArgumentTypeError(
        "The --forced-off option must be a JSON list[tuple[int,int]]."
    )

    try:
        forced_off = json.loads(value)
    except json.JSONDecodeError:
        raise exception

    # Ensure that forced_off is a list of tuples[int, int]
    if not isinstance(forced_off, list):
        raise exception

    for item in forced_off:
        if not isinstance(item, list) or len(item) != 2:
            raise exception
        if not isinstance(item[0], int) or not isinstance(item[1], int):
            raise exception

    return [tuple(item) for item in forced_off]


def _validate_dedup_interp_factors(value):

    import json

    exception = argparse.ArgumentTypeError(
        "The --dedup-interp-factors option must be a JSON list of floats."
    )

    try:
        dedup_interp_factors = json.loads(value)
    except json.JSONDecodeError:
        raise exception

    # Ensure that dedup_interp_factors is a list of floats
    if not isinstance(dedup_interp_factors, list):
        raise exception

    for factor in dedup_interp_factors:
        if not isinstance(factor, (int, float)):
            raise exception

    return dedup_interp_factors


def make_search_parser(m):
    ### casm-map search ...
    search = m.add_parser(
        "search",
        help="Search for structure mappings",
    )
    search.set_defaults(func=run_search)

    ### Positional arguments
    positional = search.add_argument_group("Positional arguments")
    positional.add_argument("parent", type=pathlib.Path, help="Parent structure file")
    positional.add_argument("child", type=pathlib.Path, help="Child structure file")

    ### Total mapping options
    total = search.add_argument_group("Total mapping options")
    total.add_argument(
        "--min-n-atoms",
        type=int,
        default=1,
        help="Minimum number of atoms in superstructures (default=1).",
    )
    total.add_argument(
        "--max-n-atoms",
        type=int,
        help=(
            "Maximum number of atoms in superstructures "
            "(default= lcm of parent/child )."
        ),
    )
    total.add_argument(
        "--min-total-cost",
        type=float,
        default=0.0,
        help="Minimum total cost (default=0.0).",
    )
    total.add_argument(
        "--max-total-cost",
        type=float,
        default=0.3,
        help="Maximum total cost (default=0.3).",
    )
    total.add_argument(
        "--k-best",
        type=int,
        default=100,
        help="Total number of mapping results.",
    )
    total.add_argument(
        "--cost-tol",
        type=float,
        default=1e-5,
        help=("Cost tolerance for approximate ties (default=1e-5)."),
    )
    total.add_argument(
        "--lattice-cost-weight",
        type=float,
        default=0.5,
        help=("Fraction of total cost due to lattice cost (default=0.5)."),
    )
    total.add_argument(
        "--iso-cost",
        action="store_true",
        default=False,
        help=(
            "Use isotropic strain and disp costs "
            "(default= use symmetry-breaking strain and disp costs)."
        ),
    )
    total.add_argument(
        "--no-remove-mean-displacement",
        action="store_true",
        default=False,
        help=(
            "Do not remove the mean displacement from the structure mappings. "
            "(default= remove mean displacement )."
        ),
    )
    # total.add_argument(
    #     "--parent-superstructure",
    #     metavar="PARENT_SUPERSTRUCTURE",
    #     type=pathlib.Path,
    #     default=None,
    #     help=(
    #         "Parent superstructure file. Fixes the parent superstructure using "
    #         "the lattice of the specified structure file "
    #         "(default= search over parent superstructures )."
    #     ),
    # )
    total.add_argument(
        "--fix-parent",
        action="store_true",
        default=False,
        help=(
            "Map to parent structure as provided; skip checking superstructures "
            "and lattice reorientations."
        ),
    )

    ### Lattice mapping options
    latmap = search.add_argument_group("Lattice mapping options")
    latmap.add_argument(
        "--iso-strain-cost",
        action="store_true",
        default=False,
        help=(
            "Use isotropic strain cost for the lattice mapping cost. "
            "(default= use symmetry-breaking strain cost )."
        ),
    )
    latmap.add_argument(
        "--min-lattice-cost",
        type=float,
        default=0.0,
        help="Minimum lattice cost (default=0.0). ",
    )
    latmap.add_argument(
        "--max-lattice-cost",
        type=float,
        default=1e20,
        help="Maximum lattice cost (default=1e20). ",
    )
    latmap.add_argument(
        "--lattice-k-best",
        type=int,
        default=10,
        help="Number of lattice mappings per parent/child superstructure (default=10).",
    )
    latmap.add_argument(
        "--lattice-reorientation-range",
        type=int,
        default=1,
        help="Max lattice reorientation matrix element (default=1).",
    )

    ### Atom mapping options
    atommap = search.add_argument_group("Atom mapping options")
    atommap.add_argument(
        "--iso-disp-cost",
        action="store_true",
        default=False,
        help=(
            "Use isotropic disp cost for the atom mapping cost. "
            "(default= use symmetry-breaking disp cost )."
        ),
    )
    atommap.add_argument(
        "--forced-on",
        type=_validate_forced_on,
        default=None,
        help=("Force specific atom mappings (JSON list[tuple[int,int]])."),
    )
    atommap.add_argument(
        "--forced-off",
        type=_validate_forced_off,
        default=None,
        help=("Suppress specific atom mappings (JSON list[tuple[int,int]])"),
    )

    ### Deduplication options
    dedup = search.add_argument_group("Deduplication options")
    dedup.add_argument(
        "--dedup-interp-factors",
        type=_validate_dedup_interp_factors,
        default=None,
        help="Interpolation factors for deduplication (JSON list[float]).",
    )

    ### Input options
    input = search.add_argument_group("Input options")
    input.add_argument(
        "--prim",
        type=str,
        help=(
            "Read a CASM Prim and map child structure onto allowed sites. It is "
            "required to also set --fix-parent)."
        ),
    )
    input.add_argument(
        "--format",
        type=str,
        default=None,
        help="Structure files format (default= inferred from file suffix ).",
    )
    input.add_argument(
        "--parent-format",
        type=str,
        default=None,
        help="Parent structure file format (overrides --format).",
    )
    input.add_argument(
        "--child-format",
        type=str,
        default=None,
        help="Child structure file format (overrides --format).",
    )
    input.add_argument(
        "--options",
        type=pathlib.Path,
        default=None,
        help=(
            "JSON file containing options to use instead of "
            "reading them from command line arguments."
        ),
    )

    ### Output options
    output = search.add_argument_group("Output options")
    output.add_argument(
        "--results-dir",
        type=pathlib.Path,
        default=pathlib.Path("results"),
        help="Directory where results are written (default=results).",
    )
    output.add_argument(
        "--merge",
        action="store_true",
        default=False,
        help="Merge new results into existing results.",
    )

    ### Additional options
    additional = search.add_argument_group("Additional options")
    additional.add_argument(
        "--next",
        action="store_true",
        help=(
            "Expand previous search to include next greatest common multiple number "
            "of atoms."
        ),
    )

    ### Other options:
    other = search.add_argument_group("Other options")
    other.add_argument(
        "--desc",
        action="store_true",
        help="Print an extended description of the method and parameters.",
    )
