import math
import pathlib
import sys
from typing import Optional

import numpy as np
from tabulate import tabulate

import libcasm.configuration as casmconfig
import libcasm.mapping.info as mapinfo
import libcasm.mapping.mapsearch as mapsearch
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
from casm.tools.shared.json_io import (
    safe_dump,
)

from .methods import (
    make_child_transformation_matrix_to_super,
)


class StructureMappingSearchOptions:

    def __init__(
        self,
        child_max_supercell_size: Optional[int] = None,
        child_min_supercell_size: int = 1,
        search_min_cost: float = 0.0,
        search_max_cost: float = 1e20,
        search_k_best: int = 1,
        lattice_mapping_min_cost: Optional[float] = 0.0,
        lattice_mapping_max_cost: Optional[float] = 1e20,
        lattice_mapping_k_best: Optional[int] = 10,
        lattice_mapping_reorientation_range: Optional[int] = 1,
        lattice_mapping_cost_method: str = "symmetry_breaking_strain_cost",
        atom_mapping_cost_method: str = "symmetry_breaking_disp_cost",
        lattice_cost_weight: float = 0.5,
        cost_tol: Optional[float] = 1e-5,
        deduplication_interpolation_factors: Optional[list[float]] = None,
    ):
        """Options for `map_structures_v2`.

        Parameters
        ----------
        child_max_supercell_size : int
            The maximum supercell size of the child to search over. If None, the
            maximum supercell size is set based on the least common multiple of the
            number of atoms in the child and parent structures.
        child_min_supercell_size : int = 1
            The minimum supercell size of the child to search over. If None, the
            minimum supercell size is set to 1.
        child_scel_range : Optional[tuple[int, int]] = None
            Defines the range of supercell sizes for the child structure which are
            included in the search. If None, the range is set to only include the
            child superstructure which has the least common multiple of the number of
            atoms in the child and parent structures.
        search_min_cost : float = 0.0
            The minimum total cost mapping to include in search results.
        search_max_cost : float = 0.5
            The maximum total cost mapping to include in search results.
        search_k_best : int = 1
            Keep the `k_best` mappings with lowest total cost that also
            satisfy the min/max cost criteria. Approximate ties with the
            current `k_best`-ranked result are also kept.
        lattice_mapping_min_cost : float = 0.0
            Keep lattice mappings with cost >= min_cost. Used when
            `map_lattices_with_reorientation` is True.
        lattice_mapping_max_cost : float = 0.5
            Keep results with cost <= max_cost. Used when
            `map_lattices_with_reorientation` is True.
        lattice_mapping_k_best : int = 10
            If not None, then only keep the k-best results (i.e. k lattice mappings
            with minimum cost) satisfying the min_cost and max_cost constraints.
            If there are approximate ties, those will also be kept. Used when
            `map_lattices_with_reorientation` is True.
        lattice_mapping_reorientation_range : int = 1
            The absolute value of the maximum element in the lattice mapping
            reorientation matrix, :math:`N`. This determines how many equivalent
            lattice vector reorientations are checked. Increasing the value results in
            more checks. The value 1 is generally expected to be sufficient because
            reduced cell lattices are compared internally.
        lattice_mapping_cost_method : str = 'symmetry_breaking_strain_cost'
            Selects the method used to calculate lattice mapping costs. Used when
            `map_lattices_with_reorientation` is True. One of
            "isotropic_strain_cost" or "symmetry_breaking_strain_cost".
        atom_mapping_cost_method : str = 'symmetry_breaking_disp_cost'
            Selects the method used to calculate atom mapping costs. One of
            "isotropic_disp_cost" or "symmetry_breaking_disp_cost".
        lattice_cost_weight : float = 0.5
            The weight of the lattice cost in the total structure mapping cost.
        cost_tol : float = 1e-5
            Tolerance for checking if mapping costs are approximately equal.
        deduplication_interpolation_factors : Optional[list[float]] = None
            Interpolation factors to use for deduplication. If None, the default value
            ``[0.5, 1.0]`` is used.

        """
        self.child_min_supercell_size = child_min_supercell_size
        self.child_max_supercell_size = child_max_supercell_size
        self.search_min_cost = search_min_cost
        self.search_max_cost = search_max_cost
        self.search_k_best = search_k_best
        self.lattice_mapping_min_cost = lattice_mapping_min_cost
        self.lattice_mapping_max_cost = lattice_mapping_max_cost
        self.lattice_mapping_k_best = lattice_mapping_k_best
        self.lattice_mapping_reorientation_range = lattice_mapping_reorientation_range
        self.lattice_mapping_cost_method = lattice_mapping_cost_method
        self.atom_mapping_cost_method = atom_mapping_cost_method
        self.lattice_cost_weight = lattice_cost_weight
        self.cost_tol = cost_tol

        # Deduplication options
        if deduplication_interpolation_factors is None:
            deduplication_interpolation_factors = [0.5, 1.0]
        self.deduplication_interpolation_factors = deduplication_interpolation_factors

    def to_dict(self):
        return {
            "child_min_supercell_size": self.child_min_supercell_size,
            "child_max_supercell_size": self.child_max_supercell_size,
            "search_min_cost": self.search_min_cost,
            "search_max_cost": self.search_max_cost,
            "search_k_best": self.search_k_best,
            "lattice_mapping_min_cost": self.lattice_mapping_min_cost,
            "lattice_mapping_max_cost": self.lattice_mapping_max_cost,
            "lattice_mapping_k_best": self.lattice_mapping_k_best,
            "lattice_mapping_reorientation_range": self.lattice_mapping_reorientation_range,  # noqa: E501
            "lattice_mapping_cost_method": self.lattice_mapping_cost_method,
            "atom_mapping_cost_method": self.atom_mapping_cost_method,
            "lattice_cost_weight": self.lattice_cost_weight,
            "cost_tol": self.cost_tol,
            "deduplication_interpolation_factors": self.deduplication_interpolation_factors,  # noqa: E501
        }

    @staticmethod
    def from_dict(data):
        return StructureMappingSearchOptions(
            child_min_supercell_size=data["child_min_supercell_size"],
            child_max_supercell_size=data["child_max_supercell_size"],
            search_min_cost=data["search_min_cost"],
            search_max_cost=data["search_max_cost"],
            search_k_best=data["search_k_best"],
            lattice_mapping_min_cost=data["lattice_mapping_min_cost"],
            lattice_mapping_max_cost=data["lattice_mapping_max_cost"],
            lattice_mapping_k_best=data["lattice_mapping_k_best"],
            lattice_mapping_reorientation_range=data[
                "lattice_mapping_reorientation_range"
            ],
            lattice_mapping_cost_method=data["lattice_mapping_cost_method"],
            atom_mapping_cost_method=data["atom_mapping_cost_method"],
            lattice_cost_weight=["lattice_cost_weight"],
            cost_tol=data["cost_tol"],
            deduplication_interpolation_factors=data[
                "deduplication_interpolation_factors"
            ],
        )


class SearchResult:
    def __init__(self):
        self.parent_structure: Optional[xtal.Structure] = None
        """Optional[xtal.Structure]: The parent structure, with lattice 
        :math:`L_{1}`."""

        self.parent_prim: Optional[xtal.Prim] = None
        """Optional[xtal.Prim]: The `xtal.Prim` used to represent the parent structure
        in the search methods."""

        self.child_structure: Optional[xtal.Structure] = None
        """Optional[xtal.Structure]: The child structure, with lattice :math:`L_{2}`."""

        self.child_T: Optional[np.ndarray] = None
        """Optional[np.ndarray]: The transformation matrix :math:`T_{2}` used to create
        a superstructure of the child for input to the search methods."""

        self.child_superstructure: Optional[xtal.Structure] = None
        """Optional[xtal.Structure]: The child superstructure, with lattice
        :math:`L_{2} T_{2}`."""

        self.scored_structure_mapping: Optional[mapinfo.ScoredStructureMapping] = None
        """Optional[mapinfo.ScoredStructureMapping]: The scored structure mapping
        between the parent structure and the child superstructure."""

        self.dedup_chain_orbit: Optional[list[list[xtal.Structure]]] = None
        """Optional[list[list[xtal.Structure]]]: A list of lists of structures, where
        `chain_orbit[i][0]` is the parent and `chain_orbit[i][-1]` is the child in the
        :math:`i`-th chain of interpolated structures in the orbit of equivalent 
        chains, put into primitive, canonical form for use in deduplication."""


def results_dir_exists_error(results_dir: pathlib.Path) -> None:
    """Print an error message if the results directory already exists."""

    error = f"""
################################################################################
# Error: results directory already exists                                      #
#                                                                              #
# A directory already exists at the specified path.                            #
#                                                                              #

--results-dir={results_dir}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def primitive_parent_notice() -> None:
    """Write a notice to the console that the parent structure is not primitive."""

    notice = """
################################################################################
# Notice: parent is not primitive                                              #
# Writing primitive parent structure: parent.primitive.json                    #
#                                                                              #
# The parent structure is not primitive, and the search will continue with the #
# non-primitive parent structure. If you want to use the primitive parent,     #
# please use the file `parent.primitive.json` instead.                         #
################################################################################
"""
    print(notice)
    sys.stdout.flush()


def primitive_child_notice() -> None:
    """Write a notice to the console that the child structure is not primitive."""

    notice = """
################################################################################
# Notice: child is not primitive                                               #
# Writing primitive child structure: child.primitive.json                      #
#                                                                              #
# The child structure is not primitive, and the search will continue with the  #
# non-primitive child structure. If you want to use the primitive child,       #
# please use the file `child.primitive.json` instead.                           #
################################################################################
"""
    print(notice)
    sys.stdout.flush()


def invalid_child_min_supercell_size_error(child_min_supercell_size: int):
    """Print an error message for invalid child minimum supercell size."""

    error = f"""
################################################################################
# Error: Invalid child minimum supercell size                                  #
#                                                                              #
# The minimum supercell size for the child structure must be at least 1.       #
#                                                                              #

--child-min-supercell-size={child_min_supercell_size}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_child_max_supercell_size_error(
    child_min_supercell_size: int,
    child_max_supercell_size: int,
    computed_msg: str,
):
    """Print an error message for invalid child maximum supercell size."""

    error = f"""
################################################################################
# Error: Invalid child maximum supercell size                                  #
#                                                                              #
# The maximum supercell size for the child structure must be greater than or   #
# equal to the minimum.                                                        #

--child-min-supercell-size={child_min_supercell_size}
--child-max-supercell-size={child_max_supercell_size} {computed_msg}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def no_valid_parent_supercell_size_error(
    child_min_supercell_size: int,
    child_max_supercell_size: int,
):
    """Print an error message for no valid parent supercell sizes."""

    error = f"""
################################################################################
# Error: No valid parent supercell size                                        #
#                                                                              #
# In the range of child supercell sizes specified, there are no parent         #
# supercells with a matching number of atoms. To automatically determine the   #
# child supercell size from the lcm of the number of atoms in the parent and   #
# child structures, do not set the --child-max-supercell-size option.          #

--child-min-supercell-size={child_min_supercell_size}
--child-max-supercell-size={child_max_supercell_size} 

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_lattice_mapping_cost_method_error(method: str):
    """Print an error message for invalid lattice mapping cost method."""

    error = f"""
################################################################################
# Error: Invalid lattice mapping cost method                                   #
#                                                                              #
# The lattice mapping cost method must be one of:                              #
# - 'isotropic_strain_cost'                                                    #
# - 'symmetry_breaking_strain_cost'                                            #
#                                                                              #

--lattice-cost-method={method}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_atom_mapping_cost_method_error(method: str):
    """Print an error message for invalid atom mapping cost method."""

    error = f"""
################################################################################
# Error: Invalid atom mapping cost method                                      #
#                                                                              #
# The atom mapping cost method must be one of:                                 #
# - 'isotropic_disp_cost'                                                      #
# - 'symmetry_breaking_disp_cost'                                              #
#                                                                              #

--atom-cost-method={method}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_deduplication_interpolation_factors_error(dedup_factors):
    """Print an error message for invalid deduplication interpolation factors."""

    error = f"""
################################################################################
# Error: Invalid deduplication interpolation factors                           #
#                                                                              #
# The deduplication interpolation factors must be a list of floats.            #
#                                                                              #

--dedup-interp-factors={dedup_factors}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


class StructureMappingSearch:
    """Search for mappings between superstructures of parent and child structures.

    Find structure mappings of the type:

    - parent and child structures have matching atom types and fractions
    - no vacancies, no parent sites with >1 allowed atom type
    - enable mean displacement removal


    To do so, this class finds lattice mappings of the type:

    .. math::

        F * L_1 * T_1 * N = L_2 * T_2

    where:

    - :math:`T_2` is a shape=(3,3) integer transformation matrix that generates a
      superlattice of the child lattice :math:`L_2`
    - other variables are defined as in the class
      :class:`libcasm.mapping.info.LatticeMapping`


    Notes
    -----

    This mapping search is limited to the case where the parent and child structures:

    - have the same atom types
    - have the same atom fractions

    Constraints that can be applied to the search:

    - min / max child supercell size
    - min / max total cost of the mapping
    - min / max cost of the lattice mapping
    - lattice mapping reorientation range
    - k-best mappings to keep

    Other options include:

    - Choice of mapping cost methods:
      - Lattice mapping cost: "isotropic_strain_cost" or "symmetry_breaking_strain_cost"
      - Atom mapping cost: "isotropic_disp_cost" or "symmetry_breaking_disp_cost"
      - Lattice cost weight
    - Choice of interpolation factors used for deduplication

    """

    def __init__(
        self,
        opt: StructureMappingSearchOptions,
    ):
        self.opt: StructureMappingSearchOptions = opt
        """StructureMappingSearchOptions: Options for the search."""

    def _get_child_max_supercell_size(
        self, parent: xtal.Structure, child: xtal.Structure
    ):
        """Get the maximum supercell size of the child structure based on the parent
        structure.

        If `child_max_supercell_size` is not set, the maximum supercell size is set to
        the least common multiple of the number of atoms in the child and parent
        structures.
        """
        if self.opt.child_max_supercell_size is None:
            n_atoms_parent = len(parent.atom_type())
            n_atoms_child = len(child.atom_type())
            return math.lcm(n_atoms_parent, n_atoms_child) // n_atoms_child
        return self.opt.child_max_supercell_size

    def _get_parent_vol(
        self,
        parent: xtal.Structure,
        child: xtal.Structure,
    ):
        child_min_supercell_size = self.opt.child_min_supercell_size
        child_max_supercell_size = self._get_child_max_supercell_size(parent, child)
        child_n_atoms = len(child.atom_type())
        parent_n_atoms = len(parent.atom_type())

        parent_vol = {}
        for child_vol in range(child_min_supercell_size, child_max_supercell_size + 1):
            child_superstructure_n_atoms = child_n_atoms * child_vol
            _vol = child_superstructure_n_atoms / parent_n_atoms

            # if parent_vol is integer, then it is a valid supercell size:
            if _vol.is_integer():
                parent_vol[child_vol] = int(_vol)

        return parent_vol

    def _enable_symmetry_breaking_atom_cost(self):
        """Check if symmetry breaking atom cost is enabled based on the options."""
        return self.opt.atom_mapping_cost_method == "symmetry_breaking_disp_cost"

    def _atom_cost_f(self):
        """Get the atom cost function based on the options."""
        if self.opt.atom_mapping_cost_method == "isotropic_disp_cost":
            return mapsearch.IsotropicAtomCost()
        elif self.opt.atom_mapping_cost_method == "symmetry_breaking_disp_cost":
            return mapsearch.SymmetryBreakingAtomCost()
        else:
            raise ValueError(
                f"Unknown atom mapping cost method: {self.opt.atom_mapping_cost_method}"
            )

    def _total_cost_f(self):
        """Get the total cost function based on the options."""
        return mapsearch.WeightedTotalCost(
            lattice_cost_weight=self.opt.lattice_cost_weight
        )

    def validate(
        self,
        parent: xtal.Structure,
        child: xtal.Structure,
    ) -> None:
        """Raise if atom types or fractions differ between parent and child."""

        # Check atom types and stoichiometry
        parent_atom_types, parent_counts = np.unique(
            parent.atom_type(), return_counts=True
        )
        total_atoms = np.sum(parent_counts)
        parent_atom_frac = parent_counts / total_atoms

        child_atom_types, child_counts = np.unique(
            child.atom_type(), return_counts=True
        )
        total_atoms = np.sum(child_counts)
        child_atom_frac = child_counts / total_atoms

        if (parent_atom_types != child_atom_types).any():
            print("Error: Parent atom types differs from child atom types")
            print(f"- Parent atom types: {parent_atom_types}")
            print(f"- Child atom types: {child_atom_types}")
            print()
            print("Stopping")
            sys.exit(1)

        if not np.allclose(parent_atom_frac, child_atom_frac):
            print("Error: Parent and child structures have different atom fractions")
            print(f"- Atom types: {parent_atom_types}")
            print(f"- Parent atom fraction: {parent_atom_frac}")
            print(f"- Child atom fraction: {child_atom_frac}")
            print()
            print("Stopping")
            sys.exit(1)

        # Print notice if parent or child are not primitive, and write the
        # primitive structures
        primitive_parent = xtal.make_primitive_structure(parent)
        if len(primitive_parent.atom_type()) != len(parent.atom_type()):
            safe_dump(
                xtal.pretty_json(primitive_parent.to_dict()),
                path="parent.primitive.json",
                force=True,
                quiet=True,
            )
            primitive_parent_notice()

        primitive_child = xtal.make_primitive_structure(child)
        if len(primitive_child.atom_type()) != len(child.atom_type()):
            safe_dump(
                xtal.pretty_json(primitive_child.to_dict()),
                path="child.primitive.json",
                force=True,
                quiet=True,
            )
            primitive_child_notice()

        # Validate the child supercell size range
        if self.opt.child_min_supercell_size < 1:
            print("Error: --min-child-supercell-size must be at least 1")
            print(f"- --min-child-supercell-size={self.opt.child_min_supercell_size}")
            print("Stopping")
            sys.exit(1)

        _child_max_supercell_size = self._get_child_max_supercell_size(parent, child)
        if _child_max_supercell_size < self.opt.child_min_supercell_size:
            computed_msg = (
                "(computed from lcm of atom counts)"
                if self.opt.child_max_supercell_size is None
                else ""
            )
            invalid_child_max_supercell_size_error(
                child_min_supercell_size=self.opt.child_min_supercell_size,
                child_max_supercell_size=_child_max_supercell_size,
                computed_msg=computed_msg,
            )

        # Validate the parent supercell size
        parent_vol = self._get_parent_vol(parent, child)
        if len(parent_vol) == 0:
            no_valid_parent_supercell_size_error(
                child_min_supercell_size=self.opt.child_min_supercell_size,
                child_max_supercell_size=_child_max_supercell_size,
            )

        # Validate lattice mapping cost method
        if self.opt.lattice_mapping_cost_method not in [
            "isotropic_strain_cost",
            "symmetry_breaking_strain_cost",
        ]:
            invalid_lattice_mapping_cost_method_error(
                self.opt.lattice_mapping_cost_method
            )

        # Validate atom mapping cost method
        if self.opt.atom_mapping_cost_method not in [
            "isotropic_disp_cost",
            "symmetry_breaking_disp_cost",
        ]:
            invalid_atom_mapping_cost_method_error(self.opt.atom_mapping_cost_method)

        # Validate that deduplication_interpolation_factors is a list of floats:
        dedup_factors = self.opt.deduplication_interpolation_factors
        if not isinstance(dedup_factors, list) or not all(
            isinstance(factor, float) for factor in dedup_factors
        ):
            invalid_deduplication_interpolation_factors_error(dedup_factors)

    def __call__(
        self,
        parent: xtal.Structure,
        child: xtal.Structure,
        results_dir: pathlib.Path,
        force: bool = False,
    ):
        """Perform the structure mapping search.

        Parameters
        ----------
        parent : xtal.Structure
            The parent structure.
        child : xtal.Structure
            The child structure.
        results_dir : pathlib.Path
            The directory to write the results to. If the directory already exists,
            the program exits with an error.

        """
        if results_dir.exists():
            results_dir_exists_error(results_dir=results_dir)
            sys.exit(1)

        self.validate(parent, child)

        ## Parameters
        _child_max_supercell_size = self._get_child_max_supercell_size(parent, child)
        _child_min_supercell_size = self.opt.child_min_supercell_size
        _parent_vol = self._get_parent_vol(parent, child)
        _enable_symmetry_breaking_atom_cost = self._enable_symmetry_breaking_atom_cost()
        _search_min_cost = self.opt.search_min_cost
        _search_max_cost = self.opt.search_max_cost
        _search_k_best = self.opt.search_k_best
        _lattice_mapping_min_cost = self.opt.lattice_mapping_min_cost
        _lattice_mapping_max_cost = self.opt.lattice_mapping_max_cost
        _lattice_mapping_cost_method = self.opt.lattice_mapping_cost_method
        _lattice_mapping_k_best = self.opt.lattice_mapping_k_best
        _lattice_mapping_reorientation_range = (
            self.opt.lattice_mapping_reorientation_range
        )
        _atom_cost_f = self._atom_cost_f()
        _total_cost_f = self._total_cost_f()
        _cost_tol = self.opt.cost_tol

        ## Fixed parameters
        _infinity = 1e20
        _enable_remove_mean_displacement = True
        _atom_to_site_cost_f = mapsearch.make_atom_to_site_cost

        ## Create a parent structure search data object.
        prim = casmconfig.Prim(xtal.Prim.from_atom_coordinates(structure=parent))

        parent_search_data = mapsearch.PrimSearchData(
            prim=prim.xtal_prim,
            enable_symmetry_breaking_atom_cost=_enable_symmetry_breaking_atom_cost,
        )
        init_child_structure_data = mapsearch.StructureSearchData(
            lattice=child.lattice(),
            atom_coordinate_cart=child.atom_coordinate_cart(),
            atom_type=child.atom_type(),
            override_structure_factor_group=None,
        )

        # Create a MappingSearch object.
        # This will hold a queue of possible mappings,
        # sorted by cost, as we generate them.
        search = mapsearch.MappingSearch(
            min_cost=_search_min_cost,
            max_cost=_search_max_cost,
            k_best=_search_k_best,
            atom_cost_f=_atom_cost_f,
            total_cost_f=_total_cost_f,
            atom_to_site_cost_f=_atom_to_site_cost_f,
            enable_remove_mean_displacement=_enable_remove_mean_displacement,
            infinity=_infinity,
            cost_tol=_cost_tol,
        )

        parent_crystal_point_group = xtal.make_structure_crystal_point_group(parent)
        child_crystal_point_group = xtal.make_structure_crystal_point_group(child)
        child_superlattices = xtal.enumerate_superlattices(
            unit_lattice=child.lattice(),
            point_group=child_crystal_point_group,
            max_volume=_child_max_supercell_size,
            min_volume=_child_min_supercell_size,
        )

        for i_child_superlattice, child_superlattice in enumerate(child_superlattices):
            child_T = xtal.make_transformation_matrix_to_super(
                unit_lattice=child.lattice(),
                superlattice=child_superlattice,
            )
            child_vol = int(round(np.linalg.det(child_T)))

            # for each child, make lattice mapping solutions
            child_structure_data = mapsearch.make_superstructure_data(
                prim_structure_data=init_child_structure_data,
                transformation_matrix_to_super=child_T,
            )

            print(
                f"Mapping child supercell {i_child_superlattice + 1}"
                f"/{len(child_superlattices)} (supercell size={child_vol})... "
            )
            sys.stdout.flush()

            _vol = _parent_vol.get(child_vol, None)
            if _vol is None:
                print("- Skipping: (no matching parent supercell size)")
                continue
            print(f"- Mapping to parent supercell size={_vol}")
            sys.stdout.flush()

            parent_superlattices = xtal.enumerate_superlattices(
                unit_lattice=parent.lattice(),
                point_group=parent_crystal_point_group,
                max_volume=_vol,
                min_volume=_vol,
            )
            if len(parent_superlattices) == 0:
                raise Exception(
                    "Implementation error: no valid parent superlattices found "
                )
            elif len(parent_superlattices) == 1:
                print("- There is 1 parent superlattice at this size...")
            else:
                print(
                    f"- There are {len(parent_superlattices)} parent superlattices "
                    "at this size..."
                )
            sys.stdout.flush()

            for i_parent_superlattice, parent_superlattice in enumerate(
                parent_superlattices
            ):

                parent_T = xtal.make_transformation_matrix_to_super(
                    unit_lattice=parent.lattice(),
                    superlattice=parent_superlattice,
                )

                print(
                    f"  - Parent supercell {i_parent_superlattice + 1}"
                    f"/{len(parent_superlattices)}... "
                )
                sys.stdout.flush()

                lattice_mappings = mapmethods.map_lattices(
                    lattice1=parent.lattice(),
                    lattice2=child_structure_data.lattice(),
                    transformation_matrix_to_super=parent_T,
                    lattice1_point_group=parent_search_data.prim_crystal_point_group(),
                    lattice2_point_group=child_structure_data.structure_crystal_point_group(),
                    min_cost=_lattice_mapping_min_cost,
                    max_cost=_lattice_mapping_max_cost,
                    cost_method=_lattice_mapping_cost_method,
                    k_best=_lattice_mapping_k_best,
                    reorientation_range=_lattice_mapping_reorientation_range,
                    cost_tol=_cost_tol,
                )

                for scored_lattice_mapping in lattice_mappings:
                    lattice_mapping_data = mapsearch.LatticeMappingSearchData(
                        prim_data=parent_search_data,
                        structure_data=child_structure_data,
                        lattice_mapping=scored_lattice_mapping,
                    )

                    # for each lattice mapping, generate possible translations
                    trial_translations = mapsearch.make_trial_translations(
                        lattice_mapping_data=lattice_mapping_data,
                    )

                    # for each combination of lattice mapping and translation,
                    # make and insert a mapping solution (MappingNode)
                    for trial_translation in trial_translations:
                        search.make_and_insert_mapping_node(
                            lattice_cost=scored_lattice_mapping.lattice_cost(),
                            lattice_mapping_data=lattice_mapping_data,
                            trial_translation_cart=trial_translation,
                            forced_on={},
                            forced_off=[],
                        )

        print()
        print("# of initial mappings: ", search.size())
        print()
        sys.stdout.flush()

        print("Searching next best mappings...")
        sys.stdout.flush()

        while search.size():
            search.partition()

        print("DONE")
        print()
        sys.stdout.flush()

        search_results = search.results()

        if len(search_results) == 0:
            print("No valid mappings found.")
        print(f"# Results: {len(search_results)}\n")
        sys.stdout.flush()

        self.write_results(
            search_results=search_results,
            parent=parent,
            child=child,
            results_dir=results_dir,
        )

        self.tabulate_results(
            search_results=search_results,
            parent=parent,
            child=child,
        )

        # records = []
        # for i, scored_structure_mapping in enumerate(search_results):
        #     new_result = StructureMappingRecord(
        #         prim=parent,
        #         scored_structure_mapping=scored_structure_mapping,
        #         opt=opt,
        #     )
        #
        #     if opt.deduplicate:
        #         # Check for duplicates:
        #         found_duplicate = False
        #         duplicate_index = 0
        #         for existing in records:
        #             if new_result.is_duplicate_of(
        #                 existing,
        #                 f_chain=opt.deduplicate_using,
        #                 only_check_strain=opt.deduplicate_only_check_strain,
        #             ):
        #                 found_duplicate = True
        #                 break
        #             duplicate_index += 1
        #
        #         if found_duplicate:
        #             if new_result.volume >= records[duplicate_index].volume:
        #                 # print(
        #                 #     f"Result {i} is a duplicate of result "
        #                 #     f"{duplicate_index}... skipping"
        #                 # )
        #                 # sys.stdout.flush()
        #                 continue
        #             else:
        #                 # print(
        #                 #     f"Result {i} is a duplicate of result "
        #                 #     f"{duplicate_index}... replacing"
        #                 # )
        #                 # sys.stdout.flush()
        #                 records[duplicate_index] = new_result
        #         else:
        #             # print(f"Result {i} is a new result... adding")
        #             sys.stdout.flush()
        #             records.append(new_result)
        #     else:
        #         # Include all:
        #         records.append(new_result)

        return 0

    def write_results(
        self,
        search_results: list[mapinfo.ScoredStructureMapping],
        parent: xtal.Structure,
        child: xtal.Structure,
        results_dir: pathlib.Path,
    ) -> None:
        """Write the results of the search."""
        data = {
            "parent": parent.to_dict(),
            "child": child.to_dict(),
            "scored_structure_mappings": [smap.to_dict() for smap in search_results],
            "options": self.opt.to_dict(),
        }
        safe_dump(
            data,
            path=results_dir / "mappings.json",
            force=True,
            quiet=True,
        )

    def tabulate_results(
        self,
        search_results: list[mapinfo.ScoredStructureMapping],
        parent: xtal.Structure,
        child: xtal.Structure,
    ) -> str:
        """Tabulate the results of the search."""

        prec = 5
        headers = [
            "Index",
            "TotCost",
            "LatCost",
            "AtmCost",
            "Parent Vol.",
            "Child Vol.",
        ]

        data = []
        for i, scored_structure_mapping in enumerate(search_results):
            smap = scored_structure_mapping

            latmap = smap.lattice_mapping()
            T_parent = latmap.transformation_matrix_to_super()
            parent_volume = int(round(np.linalg.det(T_parent)))
            T_child = make_child_transformation_matrix_to_super(
                parent_lattice=parent.lattice(),
                child_lattice=child.lattice(),
                structure_mapping=scored_structure_mapping,
            )
            child_volume = int(round(np.linalg.det(T_child)))

            total_cost = f"{smap.total_cost():.{prec}f}"
            lattice_cost = f"{smap.lattice_cost():.{prec}f}"
            atom_cost = f"{smap.atom_cost():.{prec}f}"

            data.append(
                [
                    i,
                    total_cost,
                    lattice_cost,
                    atom_cost,
                    parent_volume,
                    child_volume,
                ]
            )

        print("Lattice cost method:", self.opt.lattice_mapping_cost_method)
        print("Atom cost method:", self.opt.atom_mapping_cost_method)
        print("Lattice cost weight:", self.opt.lattice_cost_weight)
        print(tabulate(data, headers=headers, tablefmt="grid"))
        print()
