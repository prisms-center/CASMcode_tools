import pathlib

from casm.tools.shared.io import read_structure, write_structure

# <-- max width = 80 characters                            --> #
################################################################
search_desc = """
Extended description of the `casm-map search` command:

Method
------

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
and only one atom type occupying the site. In summary, this does
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


Results
-------

Results are written to a JSON file in the specified results 
directory and a summary table is printed to the console.

Results are written to:

    <results_dir>/
    └── mappings.json


The `mappings.json` output file is containing JSON 
representations of the following:

    "parent": libcasm.xtal.Structure
        The parent structure.
    "child": libcasm.xtal.Structure
        The child structure.
    "options": casm.tools.map.search.StructureMappingSearchOptions
        The structure mapping search options used.
    "mappings": list[libcasm.mapping.info.ScoredStructureMapping]
        The list of scored structure mappings found between a 
        superstructure of the parent and a superstructure of the 
        child.


Notes
-----

1. Generally it is recommended to use primitive cells of the 
   structures being mapped. The program will continue if the 
   input structures are not primitive, but it will also print a 
   notice and write files named `parent.primitive.json` and 
   `child.primitive.json` containing the primitive structures. 
   These can be used to re-run the search command with the 
   primitive structures.


Parameters
----------

Total mapping options:

--child-max-supercell-size: Optional[int]=None
    The maximum supercell size of the child to search over, in 
    multiples of the input child structure. By default, the 
    maximum supercell size is set based on the least common 
    multiple of the number of atoms in the child and parent 
    structures.
--child-min-supercell-size: int=1
    The minimum supercell size of the child to search over, in 
    multiples of the input child structure.
--min-total-cost: float=0.0
    Only mappings with a total cost greater than or equal to 
    this value are included in the final results.
--max-total-cost: float=1e20
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

Lattice mapping options:

--lattice-cost-method: str="symmetry_breaking_strain_cost"
    Selects the method used to calculate the lattice mapping 
    cost. One of "isotropic_strain_cost" or 
    "symmetry_breaking_strain_cost".
--min-lattice-cost: float=0.0
    Only lattice mappings with a lattice cost greater than or 
    equal to this value are used to find structure mappings.
--max-lattice-cost: float=1e20
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

--atom-cost-method: str="symmetry_breaking_disp_cost"
    Selects the method used to calculate the atom mapping cost. 
    One of "isotropic_disp_cost" or 
    "symmetry_breaking_disp_cost".
    
Deduplication options:

--dedup-interp-factors: Optional[list[float]]=None
    Interpolation factors to use for deduplication. A value of
    0.0 corresponds to the parent structure and a value of
    1.0 corresponds to the mapped child structure. If None, the 
    default ``[0.5, 1.0]`` is used.

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
    
Output options:

--results-dir: str="results"
    Directory where results are written. A new directory will be 
    created. If the directory already exists, the program exits 
    with an error.
    

Citing
------

A suggested way to cite this program is as follows:

"Structure mappings were found with the `casm-map` program [1], 
using the method of Thomas et al. [2] implemented in CASM [3]."


References
----------

[1] B. Puchala, J. Thomas, and A. Van der Ven, "casm-map...".
[2] J. C. Thomas, A. R. Natarajan, and A. Van der Ven, Comparing 
    crystal structures with symmetry and geometry, npj 
    Computational Materials, 7 (2021), 164.
[3] B. Puchala, J. C. Thomas, A. R. Natarajan, J. G. Goiri, 
    S. S. Behara, J. L. Kaufman, A. Van der Ven, CASM—A software
    package for first-principles based study of multicomponent 
    crystalline solids, Computational Materials Science 217 
    (2023) 111897.

"""

_other_desc = """
--output-format: Optional[str]=None
    The format for writing structure files. If not specified, 
    the format is inferred from the child format. Supported 
    formats include 'vasp', 'casm', and any format recognized by
    `ase.io.write` if ASE is installed.
"""


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


def run_search(args):

    if args.desc:
        print(search_desc)
        return 0

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

    from casm.tools.map import (
        StructureMappingSearch,
        StructureMappingSearchOptions,
    )

    opt = StructureMappingSearchOptions(
        child_max_supercell_size=args.child_max_supercell_size,
        child_min_supercell_size=args.child_min_supercell_size,
        search_min_cost=args.min_total_cost,
        search_max_cost=args.max_total_cost,
        search_k_best=args.k_best,
        lattice_cost_weight=args.lattice_cost_weight,
        cost_tol=args.cost_tol,
        lattice_mapping_min_cost=args.min_lattice_cost,
        lattice_mapping_max_cost=args.max_lattice_cost,
        lattice_mapping_k_best=args.lattice_k_best,
        lattice_mapping_reorientation_range=args.lattice_reorientation_range,
        lattice_mapping_cost_method=args.lattice_cost_method,
        atom_mapping_cost_method=args.atom_cost_method,
        deduplication_interpolation_factors=args.dedup_interp_factors,
    )

    import sys

    import libcasm.xtal as xtal

    print(xtal.pretty_json(opt.to_dict()))
    sys.stdout.flush()

    f = StructureMappingSearch(opt=opt)
    code = f(parent=parent, child=child, results_dir=args.results_dir)

    return code


def make_search_parser(m):
    ### casm-calc vasp ...
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
        "--child-max-supercell-size",
        type=int,
        help=(
            "Maximum child supercell size "
            "(default= determine from lcm of parent/child )."
        ),
    )
    total.add_argument(
        "--child-min-supercell-size",
        type=int,
        default=1,
        help=("Minimum child supercell size (default=1)."),
    )
    total.add_argument(
        "--min-total-cost",
        type=int,
        default=0.0,
        help="Minimum total cost (default=0.0).",
    )
    total.add_argument(
        "--max-total-cost",
        type=int,
        default=1e20,
        help="Maximum total cost (default=1e20).",
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

    ### Lattice mapping options
    latmap = search.add_argument_group("Lattice mapping options")
    latmap.add_argument(
        "--lattice-cost-method",
        type=str,
        default="symmetry_breaking_strain_cost",
        choices=[
            "isotropic_strain_cost",
            "symmetry_breaking_strain_cost",
        ],
        help=(
            "Method to use for calculating the lattice mapping cost. "
            "(default=symmetry_breaking_strain_cost)."
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
        "--atom-cost-method",
        type=str,
        default="symmetry_breaking_disp_cost",
        choices=[
            "isotropic_disp_cost",
            "symmetry_breaking_disp_cost",
        ],
        help=(
            "Method to use for calculating the atom mapping cost. "
            "(default=symmetry_breaking_disp_cost)."
        ),
    )

    ### Deduplication options
    dedup = search.add_argument_group("Deduplication options")
    dedup.add_argument(
        "--dedup-interp-factors",
        type=float,
        nargs="*",
        default=None,
        help="Interpolation factors for deduplication (default= 0.5 1.0).",
    )

    ### Input options
    input = search.add_argument_group("Input options")
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

    ### Output options
    output = search.add_argument_group("Output options")
    output.add_argument(
        "--results-dir",
        type=pathlib.Path,
        default=pathlib.Path("results"),
        help="Directory where results are written (default=results).",
    )
    # output.add_argument(
    #     "--output-format",
    #     type=str,
    #     default=None,
    #     help="Format for writing structure files (default= child format).",
    # )

    ### Other options:
    other = search.add_argument_group("Other options")
    other.add_argument(
        "--desc",
        action="store_true",
        help="Print an extended description of the method and parameters.",
    )
