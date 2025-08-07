import math
import pathlib
import sys
import uuid
from typing import Optional, Union

import numpy as np
from tabulate import tabulate

import libcasm.configuration as casmconfig
import libcasm.mapping.info as mapinfo
import libcasm.mapping.mapsearch as mapsearch
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
from casm.tools.shared.json_io import (
    read_optional,
    read_required,
    safe_dump,
)

from .methods import (
    chain_is_in_orbit,
    make_chain_orbit,
    make_child_supercell_info,
    make_child_transformation_matrix_to_super,
    make_parent_supercell_info,
    make_primitive_chain,
    make_primitive_chain_orbit,
    parent_supercell_size,
)


def ceildiv(a, b):
    return -(a // -b)


def floordiv(a, b):
    return a // b


def _get_max_n_atoms_for_parent_structure(
    max_n_atoms: Optional[int],
    parent_structure: xtal.Structure,
    child: xtal.Structure,
):
    """Get the maximum number of atoms to use when generating supercells of the child
    when the parent is provided as a structure.

    Parameters
    ----------
    max_n_atoms : Optional[int]
        If provided, use this value. Otherwise, use the least common multiple of the
        number of atoms in the child and parent structures.

    Returns
    -------
    max_n_atoms: int
        The maximum number of atoms to use when generating supercells of the child.
    """
    if max_n_atoms is not None:
        return max_n_atoms

    n_atoms_parent = len(parent_structure.atom_type())
    n_atoms_child = len(child.atom_type())
    return math.lcm(n_atoms_parent, n_atoms_child)


def _get_max_n_atoms_for_parent_prim(
    max_n_atoms: Optional[int],
    child: xtal.Structure,
):
    """Get the maximum number of atoms to use when generating supercells of the child
    when the parent is provided as a prim.

    Parameters
    ----------
    max_n_atoms : Optional[int]
        If provided, use the maximum of this value and the number of atoms in the child.
        Otherwise, use the number of atoms in the child structure.

    Returns
    -------
    max_n_atoms: int
        The maximum number of atoms to use when generating supercells of the child.
    """
    if max_n_atoms is not None:
        return max(max_n_atoms, len(child.atom_type()))
    else:
        return len(child.atom_type())


def _make_child_to_parent_vol(
    max_n_atoms: int,
    parent_structure: xtal.Structure,
    child: xtal.Structure,
):
    child_n_atoms = len(child.atom_type())
    parent_n_atoms = len(parent_structure.atom_type())

    child_to_parent_vol = {}
    child_vol = 1
    while child_vol * child_n_atoms <= max_n_atoms:
        child_superstructure_n_atoms = child_n_atoms * child_vol
        _vol = child_superstructure_n_atoms / parent_n_atoms

        # if parent_vol is integer, then it is a valid supercell size:
        if _vol.is_integer():
            child_to_parent_vol[child_vol] = int(_vol)

        child_vol += 1

    return child_to_parent_vol


def _make_T_pairs_for_parent_structure(
    parent_structure: xtal.Structure,
    child: xtal.Structure,
    parent_prim: casmconfig.Prim,
    min_n_atoms: int,
    max_n_atoms: int,
    child_T_list: Optional[list[np.ndarray]] = None,
    parent_T_list: Optional[list[np.ndarray]] = None,
):
    """Make a list of (T_child, T_parent) pairs for the search when a parent structure
    is given.

    Parameters
    ----------
    parent_structure : xtal.Structure
        The parent structure.
    child : xtal.Structure
        The child structure.
    parent_prim : casmconfig.Prim
        The primitive parent structure.
    min_n_atoms : int
        The minimum number of atoms in the superstructures that should be included
        in the search.
    max_n_atoms : Optional[int]
        The maximum number of atoms in the superstructures that should be included
        in the search.
    child_T_list : Optional[list[np.ndarray]] = None
        For the child superstructures, a list of transformation matrices
        :math:`T_{2}` to use. If None, the child superstructures are enumerated
        based on the `min_n_atoms` and `max_n_atoms` options.
    parent_T_list : Optional[list[np.ndarray]] = None
        For the parent superstructures, a list of transformation matrices
        :math:`T_{1}` to use. If None, the parent superstructures are enumerated
        based on the `min_n_atoms` and `max_n_atoms` options.

    Returns
    -------
    T_pairs: list[tuple[np.ndarray, np.ndarray]]
        List of (T_child, T_parent) pairs.

    """
    # Results, list of (T_child, T_parent) pairs
    T_pairs = []

    max_n_atoms = _get_max_n_atoms_for_parent_structure(
        max_n_atoms=max_n_atoms,
        parent_structure=parent_structure,
        child=child,
    )

    # Parameters
    child_crystal_point_group = xtal.make_structure_crystal_point_group(child)
    child_n_atoms = len(child.atom_type())
    child_to_parent_vol = _make_child_to_parent_vol(
        max_n_atoms=max_n_atoms,
        parent_structure=parent_structure,
        child=child,
    )

    # If child_T_list is not provided, enumerate the child supercells
    if child_T_list is None:
        child_T_list = []

        child_superlattices = xtal.enumerate_superlattices(
            unit_lattice=child.lattice(),
            point_group=child_crystal_point_group,
            max_volume=floordiv(max_n_atoms, child_n_atoms),
            min_volume=ceildiv(min_n_atoms, child_n_atoms),
        )
        for child_superlattice in child_superlattices:
            child_T_list.append(
                xtal.make_transformation_matrix_to_super(
                    unit_lattice=child.lattice(),
                    superlattice=child_superlattice,
                )
            )

    # For each child superstructure...
    for child_T in child_T_list:
        child_vol = int(round(np.linalg.det(child_T)))

        # If no valid parent volume, continue
        if child_vol not in child_to_parent_vol:
            continue
        parent_vol = child_to_parent_vol[child_vol]

        # Get the list of valid parent supercells
        restricted_parent_T_list = []

        # If parent_T_list is not provided, enumerate the parent supercells
        if parent_T_list is None:
            parent_superlattices = xtal.enumerate_superlattices(
                unit_lattice=parent_structure.lattice(),
                point_group=parent_prim.crystal_point_group.elements,
                max_volume=parent_vol,
                min_volume=parent_vol,
            )
            for parent_superlattice in parent_superlattices:
                restricted_parent_T_list.append(
                    xtal.make_transformation_matrix_to_super(
                        unit_lattice=parent_structure.lattice(),
                        superlattice=parent_superlattice,
                    )
                )

        # If parent_T_list is provided, filter the parent supercells
        else:
            for parent_T in parent_T_list:
                if int(round(np.linalg.det(parent_T))) == parent_vol:
                    restricted_parent_T_list.append(parent_T)

        # Add the (child_T, parent_T) pairs
        for parent_T in restricted_parent_T_list:
            T_pairs.append((child_T, parent_T))

    return T_pairs


class ParentVolumeSearchOptions:
    """Options for parent supercell volume to search over when the specific supercells
    have not been given.

    Options are:

    - "atoms-per-unitcell": Search a range of parent supercell sizes based on the
      range of atoms per parent unit cell in which a mapping is expected to be found.
    - "point-defect-count": Search a range of parent supercell sizes based on the
      number of vacancies and interstitials expected to be found.
    - "parent-volume-range": Explicitly give a range of parent supercell sizes to
      search.

    The default method is "atoms-per-unitcell" and the default range is
    `(n_expected_atoms_per_parent_unitcell, n_expected_atoms_per_parent_unitcell)`,
    where `n_expected_atoms_per_parent_unitcell` is by default the number of sites in
    the `parent_prim` which do not have a vacancy as the first occupant. This results
    in a single parent superstructure volume being searched, determined by the number
    of atoms in the child structure.

    """

    def __init__(
        self,
        method: str = "atoms-per-unitcell",
        expected_n_vacancy: Optional[int] = None,
        expected_n_interstitial: Optional[int] = None,
        parent_volume_range: Optional[tuple[int, int]] = None,
        atoms_per_unitcell_range: Optional[tuple[float, float]] = None,
        n_expected_atoms_per_parent_unitcell: Optional[int] = None,
    ):
        """

        .. rubric:: Constructor

        """

        self.method: str = method
        """str: Method used to determine the range of parent supercell sizes to search,
        if parent supercells are not explicitly provided.
        
        Options are:
        
        - "atoms-per-unitcell": Search a range of parent supercell sizes based on the
          number of atoms per parent unit cell expected to be found.
        - "point-defect-count": Search at a single parent supercell size determined 
          from the expected number of vacancies and interstitials.
        - "parent-volume-range": Explicitly give a range of parent supercell sizes.
        """

        self.expected_n_vacancy: Optional[int] = expected_n_vacancy
        """Optional[int]: The number of vacancies expected after mapping, if the 
        `method` is set to "point-defect-count"."""

        self.expected_n_interstitial: Optional[int] = expected_n_interstitial
        """Optional[int]: The number of interstitials expected after mapping, if the 
        `method` is set to "point-defect-count"."""

        self.parent_volume_range: Optional[tuple[int, int]] = parent_volume_range
        """Optional[tuple[int, int]]: The range of parent supercell sizes to search,
        if `method` is set to "parent-volume-range". 
        
        This is a tuple of ``(min, max)`` volume, as integer multiple of the parent 
        prim."""

        self.atoms_per_unitcell_range: Optional[tuple[float, float]] = (
            atoms_per_unitcell_range
        )
        """Optional[tuple[float, float]]: The range of atoms per parent unit cell to
        search, if `method` is set to "atoms-per-unitcell"."""

        self.n_expected_atoms_per_parent_unitcell: Optional[int] = (
            n_expected_atoms_per_parent_unitcell
        )
        """Optional[int]: The number of atoms expected in the parent unit cell if there
        are no point defects, used if `method` is set to "point-defect-count".

        By default, the number of sites in the `parent_prim` which do not have a 
        vacancy as the first occupant is used. If this is provided, it overrides the
        default value."""

    def to_dict(self) -> dict:
        """Convert the options to a dictionary."""
        return {
            "method": self.method,
            "expected_n_vacancy": self.expected_n_vacancy,
            "expected_n_interstitial": self.expected_n_interstitial,
            "parent_volume_range": self.parent_volume_range,
            "atoms_per_unitcell_range": self.atoms_per_unitcell_range,
            "n_expected_atoms_per_parent_unitcell": (
                self.n_expected_atoms_per_parent_unitcell
            ),
        }

    @staticmethod
    def from_dict(data: dict) -> "ParentVolumeSearchOptions":
        """Create a ParentVolumeSearchOptions object from a dictionary."""
        return ParentVolumeSearchOptions(
            method=data.get("method", "atoms-per-unitcell"),
            expected_n_vacancy=data.get("expected_n_vacancy"),
            expected_n_interstitial=data.get("expected_n_interstitial"),
            parent_volume_range=data.get("parent_volume_range"),
            atoms_per_unitcell_range=data.get("atoms_per_unitcell_range"),
            n_expected_atoms_per_parent_unitcell=data.get(
                "n_expected_atoms_per_parent_unitcell"
            ),
        )

    def parent_vol_range(
        self,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
        child_vol: int,
    ):
        """Return the parent volume range as a tuple of (min, max) integer volumes.

        Parameters
        ----------
        child: xtal.Structure
            The child structure
        parent_prim: casmconfig.Prim
            The parent prim
        child_vol: int
            The volume of the child supercell, as an integer multiple of the
            child structure.

        Returns
        -------
        min_parent_vol: int
            The minimum parent volume to search, as an integer multiple of the parent
            prim volume.

        max_parent_vol: int
            The maximum parent volume to search, as an integer multiple of the parent
            prim volume.

        """

        vacancy_names = ["Va", "va", "VA"]
        if self.n_expected_atoms_per_parent_unitcell is not None:
            _expected_per = self.n_expected_atoms_per_parent_unitcell
        else:
            _expected_per = 0
            occ_dof = parent_prim.xtal_prim.occ_dof()
            for site_dof in occ_dof:
                if len(site_dof):
                    if site_dof[0] not in vacancy_names:
                        _expected_per += 1

        n_child_atoms = len(child.atom_type()) * child_vol

        if self.method == "atoms-per-unitcell":
            atoms_per_unitcell_range = self.atoms_per_unitcell_range
            if atoms_per_unitcell_range is None:
                atoms_per_unitcell_range = (_expected_per, _expected_per)

            # validate that self.atoms_per_unitcell_range is a tuple[float, float]:
            if (
                not isinstance(atoms_per_unitcell_range, tuple)
                or len(atoms_per_unitcell_range) != 2
                or not all(
                    isinstance(x, (int, float)) for x in atoms_per_unitcell_range
                )
            ):
                raise ValueError(
                    "Error in ParentVolumeSearchOptions: "
                    "method is 'atoms-per-unitcell', "
                    "but atoms_per_unitcell_range is not a tuple of two numbers."
                )
            if (
                atoms_per_unitcell_range[0] <= 0.0
                or atoms_per_unitcell_range[1] <= 0.0
                or atoms_per_unitcell_range[0] > atoms_per_unitcell_range[1]
            ):
                raise ValueError(
                    "Error in ParentVolumeSearchOptions: "
                    "'atoms_per_unitcell_range' must be a tuple of two positive "
                    "numbers, where the first is less than or equal to the second."
                )
            _min_per, _max_per = self.atoms_per_unitcell_range
            min_parent_vol = int(math.floor(n_child_atoms / _max_per))
            max_parent_vol = int(math.ceil(n_child_atoms / _min_per))
            return (min_parent_vol, max_parent_vol)

        elif self.method == "point-defect-count":
            n_atoms = n_child_atoms
            n_atoms -= self.expected_n_interstitial
            n_atoms += self.expected_n_vacancy

            parent_vol = n_atoms / _expected_per

            # check if parent_vol is approximately integer:
            if not np.isclose(parent_vol, round(parent_vol), atol=1e-5):
                raise ValueError(
                    "Error in ParentVolumeSearchOptions: "
                    "method is 'point-defect-count', but the expected number of "
                    "vacancies and interstitials does not result in an integer "
                    "parent supercell volume."
                )

            parent_vol = int(round(parent_vol))

            return (parent_vol, parent_vol)

        elif self.method == "parent-volume-range":
            if self.parent_volume_range is None:
                raise ValueError(
                    "Error in ParentVolumeSearchOptions: "
                    "method is 'parent-volume-range', but no parent_volume_range "
                    "is provided."
                )
            # validate that self.parent_volume_range is a tuple[int, int]:
            if (
                not isinstance(self.parent_volume_range, tuple)
                or len(self.parent_volume_range) != 2
                or not all(isinstance(x, int) for x in self.parent_volume_range)
            ):
                raise ValueError(
                    "Error in ParentVolumeSearchOptions: "
                    "method is 'parent-volume-range', "
                    "but parent_volume_range is not a tuple of two integers."
                )
            if (
                self.parent_volume_range[0] < 1
                or self.parent_volume_range[1] < 1
                or self.parent_volume_range[0] > self.parent_volume_range[1]
            ):
                raise ValueError(
                    "Error in ParentVolumeSearchOptions: "
                    "'parent_volume_range' must be a tuple of two positive integers, "
                    "where the first is less than or equal to the second."
                )

            return self.parent_volume_range

        else:
            raise ValueError(
                "Error in ParentVolumeSearchOptions: "
                f"method={self.method} is not a valid option"
            )


def _make_T_pairs_for_parent_prim(
    child: xtal.Structure,
    parent_prim: casmconfig.Prim,
    child_atom_counts_of_parent_types: np.ndarray,
    min_atom_count_per_parent_unitcell: np.ndarray,
    max_atom_count_per_parent_unitcell: np.ndarray,
    min_n_atoms: int = 1,
    max_n_atoms: Optional[int] = None,
    child_T_list: Optional[list[np.ndarray]] = None,
    parent_T_list: Optional[list[np.ndarray]] = None,
    parent_vol_options: Optional[ParentVolumeSearchOptions] = None,
):
    """Make a list of (T_child, T_parent) pairs for the search when a parent structure
    is given.

    Method:

    1. Determine child supercells to search. If provided, use `child_T_list`. Otherwise,
       enumerate supercells based on the number of atoms in the child structure and the
       `min_n_atoms` and `max_n_atoms` options.
    2. For each child supercell, determine the parent supercells to try to map to. If
       provided, use `parent_T_list`. Otherwise, enumerate parent supercells based on
       `parent_vol_options`.
    3. Finally, filter out parent supercells that are impossible based on the
       `child_atom_counts_of_parent_types`, `min_atom_count_per_parent_unitcell`, and
       `max_atom_count_per_parent_unitcell` parameters.
       `

    Parameters
    ----------
    child : xtal.Structure
        The child structure.
    parent_prim : casmconfig.Prim
        The primitive parent structure.
    child_atom_counts_of_parent_types: np.ndarray
        The number of atoms of each parent type in the child structure.
    min_atom_count_per_parent_unitcell: np.ndarray
        The minimum number of each parent type per parent unit cell.
    max_atom_count_per_parent_unitcell: np.ndarray
        The maximum number of each parent type per parent unit cell.
    min_n_atoms : int
        The minimum number of atoms in the superstructures that should be included
        in the search.
    max_n_atoms : Optional[int]
        The maximum number of atoms in the superstructures that should be included
        in the search. If provided, use the maximum of this value and the number of
        atoms in the child. Otherwise, use the number of atoms in the child structure.
    child_T_list : Optional[list[np.ndarray]] = None
        For the child superstructures, a list of transformation matrices
        :math:`T_{2}` to use. If None, the child superstructures are enumerated
        based on the `min_n_atoms` and `max_n_atoms` options.
    parent_T_list : Optional[list[np.ndarray]] = None
        For the parent superstructures, a list of transformation matrices
        :math:`T_{1}` to use. If None, the parent superstructures are enumerated
        based on the `min_n_atoms` and `max_n_atoms` options.
    parent_vol_options : Optional[ParentVolumeSearchOptions] = None
        Options for the parent supercell volumes to search over when the specific
        supercells have not been given. If None, the default options are used,
        which is to search a single parent supercell size based on the number of
        atoms in the child structure.

    Returns
    -------
    T_pairs: list[tuple[np.ndarray, np.ndarray]]
        List of (T_child, T_parent) pairs.

    """
    # Results, list of (T_child, T_parent) pairs
    T_pairs = []

    max_n_atoms = _get_max_n_atoms_for_parent_prim(
        max_n_atoms=max_n_atoms,
        child=child,
    )

    # Parameters
    child_crystal_point_group = xtal.make_structure_crystal_point_group(child)
    child_n_atoms = len(child.atom_type())

    # If child_T_list is not provided, enumerate the child supercells
    if child_T_list is None:
        child_T_list = []

        child_superlattices = xtal.enumerate_superlattices(
            unit_lattice=child.lattice(),
            point_group=child_crystal_point_group,
            max_volume=floordiv(max_n_atoms, child_n_atoms),
            min_volume=ceildiv(min_n_atoms, child_n_atoms),
        )
        for child_superlattice in child_superlattices:
            child_T_list.append(
                xtal.make_transformation_matrix_to_super(
                    unit_lattice=child.lattice(),
                    superlattice=child_superlattice,
                )
            )

    # If parent_T_list is not provided, generate possible parent supercells for each
    # child supercell
    for child_T in child_T_list:
        child_vol = int(round(np.linalg.det(child_T)))
        superchild_atom_counts = child_atom_counts_of_parent_types * child_vol

        # Get the list of valid parent supercells

        # If parent_T_list is not provided,
        # then enumerate parent supercells based on `parent_vol_options`
        if parent_T_list is None:
            parent_T_list = []
            if parent_vol_options is None:
                parent_vol_options = ParentVolumeSearchOptions()
            parent_vol_range = parent_vol_options.parent_vol_range(
                child=child,
                parent_prim=parent_prim,
                child_vol=child_vol,
            )
            parent_superlattices = xtal.enumerate_superlattices(
                unit_lattice=parent_prim.xtal_prim.lattice(),
                point_group=parent_prim.crystal_point_group.elements,
                max_volume=parent_vol_range[0],
                min_volume=parent_vol_range[1],
            )
            for parent_superlattice in parent_superlattices:
                parent_T_list.append(
                    xtal.make_transformation_matrix_to_super(
                        unit_lattice=parent_prim.xtal_prim.lattice(),
                        superlattice=parent_superlattice,
                    )
                )

        # Filter parent_T_list based on the child atom counts and the
        # min/max atom counts per parent unit cell
        restricted_parent_T_list = []
        for parent_T in parent_T_list:
            parent_vol = int(round(np.linalg.det(parent_T)))

            # Check if the parent atom counts are within the min/max range
            if np.all(
                superchild_atom_counts
                >= min_atom_count_per_parent_unitcell * parent_vol
            ) and np.all(
                superchild_atom_counts
                <= max_atom_count_per_parent_unitcell * parent_vol
            ):
                restricted_parent_T_list.append(parent_T)

        # Add the (child_T, parent_T) pairs
        for parent_T in restricted_parent_T_list:
            T_pairs.append((child_T, parent_T))

    return T_pairs


class StructureMappingSearchOptions:
    """Options controlling the structure mapping search."""

    def __init__(
        self,
        max_n_atoms: Optional[int] = None,
        min_n_atoms: int = 1,
        child_transformation_matrix_to_super_list: Optional[list[np.ndarray]] = None,
        parent_transformation_matrix_to_super_list: Optional[list[np.ndarray]] = None,
        total_min_cost: float = 0.0,
        total_max_cost: float = 0.3,
        total_k_best: int = 1,
        no_remove_mean_displacement: bool = False,
        fix_parent: bool = False,
        lattice_mapping_min_cost: Optional[float] = 0.0,
        lattice_mapping_max_cost: Optional[float] = 1e20,
        lattice_mapping_k_best: Optional[int] = 10,
        lattice_mapping_reorientation_range: Optional[int] = 1,
        lattice_mapping_cost_method: str = "symmetry_breaking_strain_cost",
        atom_mapping_cost_method: str = "symmetry_breaking_disp_cost",
        forced_on: Optional[dict[int, int]] = None,
        forced_off: Optional[list[tuple[int, int]]] = None,
        lattice_cost_weight: float = 0.5,
        cost_tol: Optional[float] = 1e-5,
        deduplication_interpolation_factors: Optional[list[float]] = None,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        max_n_atoms : Optional[int] = None
            The maximum number of atoms in the superstructures that should be included
            in the search. If None, the least common multiple of the number of atoms
            in the child and parent structures.
        min_n_atoms : int = 1
            The minimum number of atoms in the superstructures that should be included
            in the search.
        child_transformation_matrix_to_super_list : Optional[list[np.ndarray]] = None
            If provided, overrides the `min_n_atoms` and `max_n_atoms` options to
            directly specify the transformation matrices to use for the child
            superstructures.
        parent_transformation_matrix_to_super_list : Optional[list[np.ndarray]] = None
            If provided, only use the specified transformation matrices to create
            parent superstructures. If None, the parent superstructures are
            enumerated.
        total_min_cost : float = 0.0
            The minimum total cost mapping to include in search results.
        total_max_cost : float = 0.3
            The maximum total cost mapping to include in search results.
        total_k_best : int = 1
            Keep the `k_best` mappings with lowest total cost that also
            satisfy the min/max cost criteria. Approximate ties with the
            current `k_best`-ranked result are also kept.
        no_remove_mean_displacement : bool = False
            If True, do not remove the mean displacement from the atom mapping.
        fix_parent : bool = False
            If True, map to the parent structure as provided and skip searching over
            parent superstructures and lattice reorientations. The deformation
            gradient is still calculated and atom mapping is still performed. Only
            allowed if the number of atoms in the parent structure is the same as
            the number of atoms in the child structure.
        lattice_mapping_min_cost : float = 0.0
            Keep lattice mappings with cost >= min_cost. Used when
            `map_lattices_with_reorientation` is True.
        lattice_mapping_max_cost : float = 1e20
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
        forced_on : Optional[dict[int, int]] = None
            A map of assignments `parent_atom_index: child_atom_index` that are forced
            on. Indices begin at 0. Requires that `fix_parent` is True.
        forced_off : Optional[list[tuple[int, int]]] = None
            A list of tuples of assignments `(parent_atom_index, child_atom_index) that
            are forced off. Indices begin at 0. Requires that `fix_parent` is True.
        lattice_cost_weight : float = 0.5
            The weight of the lattice cost in the total structure mapping cost.
        cost_tol : float = 1e-5
            Tolerance for checking if mapping costs are approximately equal.
        deduplication_interpolation_factors : Optional[list[float]] = None
            Interpolation factors to use for deduplication. If None, the default value
            ``[0.5, 1.0]`` is used.

        """
        self.min_n_atoms = min_n_atoms
        self.max_n_atoms = max_n_atoms
        self.child_transformation_matrix_to_super_list = (
            child_transformation_matrix_to_super_list
        )
        self.parent_transformation_matrix_to_super_list = (
            parent_transformation_matrix_to_super_list
        )

        self.total_min_cost = total_min_cost
        self.total_max_cost = total_max_cost
        self.total_k_best = total_k_best
        self.no_remove_mean_displacement = no_remove_mean_displacement
        self.fix_parent = fix_parent
        self.lattice_mapping_min_cost = lattice_mapping_min_cost
        self.lattice_mapping_max_cost = lattice_mapping_max_cost
        self.lattice_mapping_k_best = lattice_mapping_k_best
        self.lattice_mapping_reorientation_range = lattice_mapping_reorientation_range
        self.lattice_mapping_cost_method = lattice_mapping_cost_method
        self.atom_mapping_cost_method = atom_mapping_cost_method
        self.forced_on = forced_on
        self.forced_off = forced_off
        self.lattice_cost_weight = lattice_cost_weight
        self.cost_tol = cost_tol

        # Deduplication options
        if deduplication_interpolation_factors is None:
            deduplication_interpolation_factors = [0.5, 1.0]
        self.deduplication_interpolation_factors = deduplication_interpolation_factors

    def to_dict(self):
        return {
            "min_n_atoms": self.min_n_atoms,
            "max_n_atoms": self.max_n_atoms,
            "child_transformation_matrix_to_super_list": (
                [x.tolist() for x in self.child_transformation_matrix_to_super_list]
                if self.child_transformation_matrix_to_super_list is not None
                else None
            ),
            "parent_transformation_matrix_to_super_list": (
                [x.tolist() for x in self.parent_transformation_matrix_to_super_list]
                if self.parent_transformation_matrix_to_super_list is not None
                else None
            ),
            "total_min_cost": self.total_min_cost,
            "total_max_cost": self.total_max_cost,
            "total_k_best": self.total_k_best,
            "no_remove_mean_displacement": self.no_remove_mean_displacement,
            "fix_parent": self.fix_parent,
            "lattice_mapping_min_cost": self.lattice_mapping_min_cost,
            "lattice_mapping_max_cost": self.lattice_mapping_max_cost,
            "lattice_mapping_k_best": self.lattice_mapping_k_best,
            "lattice_mapping_reorientation_range": self.lattice_mapping_reorientation_range,  # noqa: E501
            "lattice_mapping_cost_method": self.lattice_mapping_cost_method,
            "atom_mapping_cost_method": self.atom_mapping_cost_method,
            "forced_on": (
                [[key, value] for key, value in self.forced_on.items()]
                if self.forced_on is not None
                else None
            ),
            "forced_off": self.forced_off,
            "lattice_cost_weight": self.lattice_cost_weight,
            "cost_tol": self.cost_tol,
            "deduplication_interpolation_factors": self.deduplication_interpolation_factors,  # noqa: E501
        }

    @staticmethod
    def from_dict(data: dict):
        return StructureMappingSearchOptions(
            max_n_atoms=data["max_n_atoms"],
            min_n_atoms=data["min_n_atoms"],
            child_transformation_matrix_to_super_list=(
                [np.array(x) for x in data["child_transformation_matrix_to_super_list"]]
                if data["child_transformation_matrix_to_super_list"] is not None
                else None
            ),
            parent_transformation_matrix_to_super_list=(
                [
                    np.array(x)
                    for x in data["parent_transformation_matrix_to_super_list"]
                ]
                if data["parent_transformation_matrix_to_super_list"] is not None
                else None
            ),
            total_min_cost=data["total_min_cost"],
            total_max_cost=data["total_max_cost"],
            total_k_best=data["total_k_best"],
            no_remove_mean_displacement=data["no_remove_mean_displacement"],
            fix_parent=data["fix_parent"],
            lattice_mapping_min_cost=data["lattice_mapping_min_cost"],
            lattice_mapping_max_cost=data["lattice_mapping_max_cost"],
            lattice_mapping_k_best=data["lattice_mapping_k_best"],
            lattice_mapping_reorientation_range=data[
                "lattice_mapping_reorientation_range"
            ],
            lattice_mapping_cost_method=data["lattice_mapping_cost_method"],
            atom_mapping_cost_method=data["atom_mapping_cost_method"],
            forced_on=(
                {x[0]: x[1] for x in data["forced_on"]}
                if data["forced_on"] is not None
                else None
            ),
            forced_off=(
                [tuple(x) for x in data["forced_off"]]
                if data["forced_off"] is not None
                else None
            ),
            lattice_cost_weight=data["lattice_cost_weight"],
            cost_tol=data["cost_tol"],
            deduplication_interpolation_factors=data[
                "deduplication_interpolation_factors"
            ],
        )


class MappingSearchData:
    def __init__(
        self,
        parent: Union[xtal.Structure, casmconfig.Prim],
        child: xtal.Structure,
        options: Optional[StructureMappingSearchOptions] = None,
        mappings: list[mapinfo.ScoredStructureMapping] = [],
        uuids: list[str] = [],
        options_history: list[StructureMappingSearchOptions] = [],
    ):

        if isinstance(parent, xtal.Structure):
            parent_structure = parent
            xtal_prim = xtal.Prim.from_atom_coordinates(structure=parent)
            parent_prim = casmconfig.Prim(xtal_prim)
        elif isinstance(parent, casmconfig.Prim):
            parent_structure = None
            parent_prim = parent
        else:
            raise TypeError(
                "Error in MappingSearchData: `parent` must be either a "
                "libcasm.xtal.Structure or a libcasm.configuration.Prim."
            )

        self.parent_structure: Optional[xtal.Structure] = parent_structure
        """Optional[xtal.Structure]: The parent structure, with lattice 
        :math:`L_{1}`, if mapping to a particular structure."""

        self.parent_prim: casmconfig.Prim = parent_prim
        """casmconfig.Prim: The :class:`~libcasm.configuration.Prim` for the parent, 
        which determines allowed occupants on each basis site."""

        self.child: xtal.Structure = child
        """Optional[xtal.Structure]: The child structure, with lattice :math:`L_{2}`."""

        self.options: Optional[StructureMappingSearchOptions] = options
        """Optional[StructureMappingSearchOptions]: Options for the current search."""

        self.mappings: list[mapinfo.ScoredStructureMapping] = mappings
        """list[mapinfo.ScoredStructureMapping]: The list of scored structure mappings
        between the parent structure and the child superstructure."""

        self.uuids: list[str] = uuids
        """list[str]: A list of UUIDs for the mappings."""

        self.options_history: list[StructureMappingSearchOptions] = options_history
        """list[StructureMappingSearchOptions]: A history of options used by previous
        searches."""

    @property
    def parent_atom_types(self):
        """list[str]: The list of atom types in the parent structure, sorted."""
        if self.parent_structure is not None:
            parent_atom_types = set(self.parent_structure.atom_type())
        else:
            occ_dof = self.parent_prim.xtal_prim.occ_dof()
            parent_atom_types = {name for site_dof in occ_dof for name in site_dof}
        return sorted(list(parent_atom_types))

    @property
    def parent_atom_frac(self):
        """Optional[np.ndarray]: The fraction of each atom type in the parent, in
        order corresponding to `parent_atom_types`.

        If `parent_structure` is None, or has no atoms, the value is None."""
        if self.parent_structure is None:
            return None
        if len(self.parent_structure.atom_type()) == 0:
            return None
        _atom_types = self.parent_atom_types
        _atom_count = [0] * len(_atom_types)
        for atom_type in self.parent_structure.atom_type():
            _atom_count[_atom_types.index(atom_type)] += 1
        _atom_counts = np.array(_atom_count)
        total = np.sum(_atom_counts)
        return _atom_counts / total

    @property
    def min_atom_count_per_parent_unitcell(self):
        """np.array: The minimum number of atoms per parent unit cell of each type, in
        order corresponding to `parent_atom_types`."""
        _atom_types = self.parent_atom_types
        _atom_count = [0] * len(_atom_types)
        if self.parent_structure is not None:
            for atom_type in self.parent_structure.atom_type():
                _atom_count[_atom_types.index(atom_type)] += 1
        else:
            occ_dof = self.parent_prim.xtal_prim.occ_dof()
            for site_dof in occ_dof:
                if len(site_dof) == 1:
                    _atom_count[_atom_types.index(site_dof[0])] += 1
        return np.array(_atom_count)

    @property
    def max_atom_count_per_parent_unitcell(self):
        """np.array: The maximum number of atoms per parent unit cell of each type, in
        order corresponding to `parent_atom_types`."""
        _atom_types = self.parent_atom_types
        _atom_count = [0] * len(_atom_types)
        if self.parent_structure is not None:
            for atom_type in self.parent_structure.atom_type():
                _atom_count[_atom_types.index(atom_type)] += 1
        else:
            occ_dof = self.parent_prim.xtal_prim.occ_dof()
            for site_dof in occ_dof:
                for name in site_dof:
                    _atom_count[_atom_types.index(name)] += 1
        return np.array(_atom_count)

    @property
    def child_atom_count(self):
        """np.array: The number of atoms in the child of each type, in
        order corresponding to `child_atom_types`."""
        _atom_types = self.child_atom_types
        _atom_count = [0] * len(_atom_types)
        for atom_type in self.child_structure.atom_type():
            index = _atom_types.index(atom_type)
            if index >= 0:
                _atom_count[index] += 1
        return np.array(_atom_count)

    @property
    def child_atom_count_of_parent_types(self):
        """np.array: The number of atoms in the child of each parent type, in
        order corresponding to `parent_atom_types`."""
        _atom_types = self.parent_atom_types
        _atom_count = [0] * len(_atom_types)
        for atom_type in self.child_structure.atom_type():
            index = _atom_types.index(atom_type)
            if index >= 0:
                _atom_count[index] += 1
        return np.array(_atom_count)

    @property
    def child_atom_types(self):
        """list[str]: The list of atom types in the child structure, sorted."""
        child_atom_types = set(self.child.atom_type())
        return sorted(list(child_atom_types))

    @property
    def child_atom_frac(self):
        """Optional[np.ndarray]: The fraction of each atom type in the child, in
        order corresponding to `child_atom_types`.

        If `child` has no atoms, return None."""
        if len(self.child.atom_type()) == 0:
            return None
        _atom_types = self.child_atom_types
        _atom_count = [0] * len(_atom_types)
        for atom_type in self.child.atom_type():
            _atom_count[_atom_types.index(atom_type)] += 1
        _atom_counts = np.array(_atom_count)
        total = np.sum(_atom_counts)
        return _atom_counts / total

    def validate_atom_types(self):
        """Validate parent and child atom types are consistent

        If `parent_structure` is not None, check that the atom types in the parent
        structure and child structure are the same.

        If `parent_structure` is None, check that the atom types in the child structure
        are a subset of the atom types in the parent prim.

        If the check fails, print an error message and exit.

        """
        if self.parent_structure is not None:
            if set(self.parent_atom_types) != set(self.child_atom_types):
                atom_types_mismatch_error(self.parent_atom_types, self.child_atom_types)
        else:
            if not set(self.child_atom_types).issubset(set(self.parent_atom_types)):
                atom_types_mismatch_error(self.parent_atom_types, self.child_atom_types)

    def validate_atom_frac(self):
        """Validate the parent and child atom fractions are consistent.

        If `parent_structure` is not None, check that the atom fractions in the parent
        structure and child structure are the same.
        """
        if self.parent_structure is not None:
            if not np.allclose(self.parent_atom_frac, self.child_atom_frac, atol=1e-5):
                atom_fraction_mismatch_error(
                    self.parent_atom_frac, self.child_atom_frac
                )

    def validate_forced_on(self):
        """Validate `forced_on` values, if the parent_structure is given

        If `parent_structure` is not None, check that the `--forced-on` option maps
        parent and child atoms of the same type.

        If `parent_structure` is None, this method currently does nothing.

        The `forced_on` values are also validated before atom mapping.

        """
        if self.parent_structure is None:
            return
        _allowed = [list([x]) for x in self.parent_structure.atom_type()]
        _child_types = self.child.atom_type()
        for parent_site_index, child_atom_index in self.opt.forced_on.items():
            child_type = _child_types[child_atom_index]
            if child_type not in _allowed[parent_site_index]:
                invalid_forced_on_values_error(
                    parent_site_index=parent_site_index,
                    child_atom_index=child_atom_index,
                    child_type=child_type,
                    allowed_types=_allowed[parent_site_index],
                )

    def validate_fix_parent(self):
        """Validate `fix_parent` option, if the parent_structure is given

        If `parent_structure` is not None, check that the `--fix-parent` option is used
        only when the number of atoms in the parent structure is the same as the number
        of atoms in the child structure.

        If `parent_structure` is None, this method currently does nothing.

        """
        if self.parent_structure is None:
            return
        if self.opt.fix_parent:
            child_n_atoms = len(self.child.atom_type())
            parent_n_atoms = len(self.parent_structure.atom_type())
            if child_n_atoms != parent_n_atoms:
                invalid_fix_parent_error()

    def notify_if_non_primitive(self):
        """Print a notice if the parent or child is not primitive, and write
        the primitive form to a file.

        If the parent is not primitive, write its primitive form to
        `parent.primitive.json`. If parent_structure is not None, it is checked.
        Otherwise, the parent_prim is checked.

        If the child is not primitive, write its primitive form to
        `child.primitive.json`.
        """
        if self.parent_structure is not None:
            parent = self.parent_structure
            primitive_parent = xtal.make_primitive_structure(parent)
            if len(primitive_parent.atom_type()) != len(parent.atom_type()):
                safe_dump(
                    xtal.pretty_json(primitive_parent.to_dict()),
                    path="parent.primitive.json",
                    force=True,
                    quiet=True,
                )
                primitive_parent_notice()
        else:
            parent = self.parent_structure
            primitive_parent = xtal.make_primitive_prim(parent)
            if len(primitive_parent.occ_dof()) != len(parent.occ_dof()):
                safe_dump(
                    xtal.pretty_json(primitive_parent.to_dict()),
                    path="parent.primitive.json",
                    force=True,
                    quiet=True,
                )
                primitive_parent_notice()

        primitive_child = xtal.make_primitive_structure(self.child)
        if len(primitive_child.atom_type()) != len(self.child.atom_type()):
            safe_dump(
                xtal.pretty_json(primitive_child.to_dict()),
                path="child.primitive.json",
                force=True,
                quiet=True,
            )
            primitive_child_notice()

    def validate_n_atoms(self):
        """Validate the min_n_atoms and max_n_atoms options.

        If `child_transformation_matrix_to_super_list` is None, this does nothing
        because the user has requested which child supercells to try mapping.

        Otherwise, it checks that the `min_n_atoms` >= 1 and that `max_n_atoms` is
        greater than or equal to `min_n_atoms`. If the user does not specify
        `max_n_atoms` explicitly, it is computed using the least common multiple of the
        number of atoms in the parent and child structures.

        """
        if self.options.child_transformation_matrix_to_super_list is None:

            min_n_atoms = self.options.min_n_atoms
            max_n_atoms = _get_max_n_atoms_for_parent_structure(
                max_n_atoms=self.options.max_n_atoms,
                parent_structure=self.parent_structure,
                child=self.child,
            )

            # Validate the min/max number of atoms
            if min_n_atoms < 1:
                invalid_min_n_atoms_error(min_n_atoms=min_n_atoms)

            if self.parent_structure is not None:

                if max_n_atoms < min_n_atoms:
                    computed_msg = (
                        "(computed from lcm of atom counts)"
                        if self.options.max_n_atoms is None
                        else ""
                    )
                    invalid_max_n_atoms_error(
                        min_n_atoms=min_n_atoms,
                        max_n_atoms=max_n_atoms,
                        computed_msg=computed_msg,
                    )

    def move_options_to_history(self):
        """Move the current options to the options history."""
        if self.options is not None:
            self.options_history.append(self.options)
            self.options = None

    def to_dict(self):
        """Convert the search data to a Python dictionary.

        Notes
        -----

        This does not move the current options to the options history. Use
        :func:`move_options_to_history` before calling this method if you want to
        include the current options in the history.

        Returns
        -------
        data: dict
            A Python dict representation of the search data.

        """
        return {
            "parent_structure": (
                self.parent_structure.to_dict() if self.parent_structure else None
            ),
            "parent_prim": self.parent_prim.to_dict(),
            "child": self.child.to_dict(),
            "options": self.options.to_dict() if self.options else None,
            "mappings": [mapping.to_dict() for mapping in self.mappings],
            "uuids": self.uuids,
            "options_history": [opt.to_dict() for opt in self.options_history],
        }

    @staticmethod
    def from_dict(
        self,
        data: dict,
    ):
        """Create a MappingSearchData object from a Python dictionary.

        Parameters
        ----------
        data: dict
            A Python dict representation of the search data.

        Returns
        -------
        search_data: MappingSearchData
            The MappingSearchData object created from the dictionary.

        """
        parent_structure = (
            xtal.Structure.from_dict(data["parent_structure"])
            if data["parent_structure"] is not None
            else None
        )
        parent_prim = casmconfig.Prim.from_dict(data["parent_prim"])
        child = xtal.Structure.from_dict(data["child"])
        options = (
            StructureMappingSearchOptions.from_dict(data["options"])
            if data["options"] is not None
            else None
        )
        mappings = [
            mapinfo.ScoredStructureMapping.from_dict(data=x, prim=parent_prim.xtal_prim)
            for x in data["mappings"]
        ]
        uuids = data["uuids"]
        options_history = [
            StructureMappingSearchOptions.from_dict(data=x)
            for x in data["options_history"]
        ]

        return MappingSearchData(
            parent=parent_structure or parent_prim,
            child=child,
            options=options,
            mappings=mappings,
            uuids=uuids,
            options_history=options_history,
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


def update_options_to_next_n_atoms(
    self,
    options: StructureMappingSearchOptions,
    results_dir: pathlib.Path,
):
    n_atoms_parent = len(self.parent.atom_type())
    n_atoms_child = len(self.child.atom_type())
    n_atoms_lcm = math.lcm(n_atoms_parent, n_atoms_child)

    last_max_n_atoms = None
    if len(self.options_history):
        if self.options_history[-1].max_n_atoms is None:
            last_max_n_atoms = n_atoms_lcm
        else:
            last_max_n_atoms = self.options_history[-1].max_n_atoms

    if last_max_n_atoms is None:
        next_max_n_atoms = n_atoms_lcm
    else:
        next_max_n_atoms = 0
        while next_max_n_atoms <= last_max_n_atoms:
            next_max_n_atoms += n_atoms_lcm

    options.max_n_atoms = next_max_n_atoms
    options.min_n_atoms = next_max_n_atoms


def results_dir_exists_error(results_dir: pathlib.Path) -> None:
    """Print an error message if the results directory already exists."""

    error = f"""
################################################################################
# Error: Results directory already exists                                      #
#                                                                              #
# A directory already exists at the specified path.                            #
#                                                                              #
# To merge new results, use --merge. Otherwise, delete the existing directory  #
# or specify a new one.                                                        #

--results-dir={results_dir}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def atom_types_mismatch_error(parent_atom_types, child_atom_types) -> None:
    """Print an error message if the parent and child atom types do not match."""
    error = f"""
    ################################################################################
    # Error: Child cannot map to parent due to atom types mismatch                 #
    #                                                                              #
    
    - Parent atom types: {parent_atom_types}
    - Child atom types: {child_atom_types}
    
    # Stopping...                                                                  #
    ################################################################################
    """
    print(error)
    sys.exit(1)


def atom_fraction_mismatch_error(parent_atom_frac, child_atom_frac) -> None:
    """Print an error message if the parent and child atom types do not match."""
    error = f"""
    ################################################################################
    # Error: Parent and child structures have different atom fractions             #
    #                                                                              #

    - Parent atom fraction: {parent_atom_frac}
    - Child atom fraction: {child_atom_frac}

    # Stopping...                                                                  #
    ################################################################################
    """
    print(error)
    sys.exit(1)


def invalid_forced_on_values_error(
    parent_site_index: int,
    child_atom_index: int,
    child_type: str,
    allowed_types: list[str],
) -> None:
    """Print an error message if the `--forced-on` option is used with invalid
    values."""

    error = f"""
################################################################################
# Error: Invalid --forced-on values                                            #
#                                                                              #
# The `--forced-on` option requires that the child atom type is allowed on the #
# parent site.                                                                 #

child_atom_index={child_atom_index}
child_type={child_type}
parent_site_index={parent_site_index}
allowed_types={allowed_types}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_fix_parent_error() -> None:
    """Print an error message if the `--fix-parent` option is used with a parent and
    child structure that have different numbers of atoms."""

    error = """
################################################################################
# Error: --fix-parent requires parent and child w/ same number of atoms.       #
#                                                                              #
# The `--fix-parent` option is used to map to the parent structure as          #
# provided, without searching over parent superstructures and lattice          #
# reorientations. It is only allowed if the number of atoms in the parent      #
# structure is the same as the number of atoms in the child structure.         #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_parent_error(results_dir: pathlib.Path) -> None:
    """When merging, print an error message if the parent has changed."""

    error = """
################################################################################
# Error: parent structure has changed                                          #
#                                                                              #
# When using the --merge option, the parent structure must remain the same.    #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_child_error(results_dir: pathlib.Path) -> None:
    """When merging, print an error message if the child has changed."""

    error = """
################################################################################
# Error: child structure has changed                                           #
#                                                                              #
# When using the --merge option, the child structure must remain the same.     #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_lattice_mapping_cost_method_error() -> None:
    """When merging, print an error message if the lattice mapping cost method has
    changed."""

    error = """
################################################################################
# Error: lattice mapping cost method has changed                               #
#                                                                              #
# When using the --merge option, the lattice mapping cost method must remain   #
# the same.                                                                    #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_atom_mapping_cost_method_error() -> None:
    """When merging, print an error message if the atom mapping cost method has
    changed."""

    error = """
################################################################################
# Error: atom mapping cost method has changed                                  #
#                                                                              #
# When using the --merge option, the atom mapping cost method must remain      #
# the same.                                                                    #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def different_lattice_cost_weight_error() -> None:
    """When merging, print an error message if the lattice cost weight has changed."""

    error = """
################################################################################
# Error: lattice cost weight has changed                                       #
#                                                                              #
# When using the --merge option, the lattice cost weight must remain the same. #
#                                                                              #
# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def primitive_parent_notice() -> None:
    """Write a notice to the console that the parent is not primitive."""

    notice = """
################################################################################
# Notice: parent is not primitive                                              #
# Writing primitive parent: parent.primitive.json                              #
#                                                                              #
# The parent is not primitive, and the search will continue with the           #
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


def invalid_min_n_atoms_error(min_n_atoms: int):
    """Print an error message for invalid min_n_atoms."""

    error = f"""
################################################################################
# Error: Invalid min_n_atoms                                                   #
#                                                                              #
# The value of min_n_atoms must be at least 1.                                 #
#                                                                              #

--min-n-atoms={min_n_atoms}

# Stopping...                                                                  #
################################################################################
"""
    print(error)
    sys.exit(1)


def invalid_max_n_atoms_error(
    min_n_atoms: int,
    max_n_atoms: int,
    computed_msg: str,
):
    """Print an error message for invalid max_n_atoms."""

    error = f"""
################################################################################
# Error: Invalid max_n_atoms                                                   #
#                                                                              #
# The value of max_n_atoms must be greater than or equalt to min_n_atoms.      #
# equal to the minimum.                                                        #

--min-n-atoms={min_n_atoms}
--max-n-atoms={max_n_atoms} {computed_msg}

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


    To do so, this method finds lattice mappings of the type:

    .. math::

        F L_1 T_1 N = L_2 T_2

    where:

    - :math:`T_2` is a shape=(3,3) integer transformation matrix that generates a
      superlattice of the child lattice :math:`L_2`
    - other variables are defined as in the class
      :class:`libcasm.mapping.info.LatticeMapping`, using :math:`T_1` for :math:`T`.


    Notes
    -----

    This mapping search is limited to the case where the parent and child structures:

    - have the same atom types
    - have the same atom fractions

    Constraints that can be applied to the search:

    - min / max number of atoms
    - parent / child supercells used
    - min / max total cost of the mapping
    - min / max cost of the lattice mapping
    - lattice mapping reorientation range
    - k-best mappings to keep

    Other options include:

    - Choice of mapping cost methods:

      - Lattice mapping cost: "isotropic_strain_cost" or "symmetry_breaking_strain_cost"
      - Atom mapping cost: "isotropic_disp_cost" or "symmetry_breaking_disp_cost"
      - Lattice cost weight: The fraction of the total cost that is due to the lattice
        mapping cost. The remainder is due to the atom mapping cost.

    - Choice of interpolation factors used for deduplication

    """

    def __init__(
        self,
        opt: StructureMappingSearchOptions,
    ):
        self.opt: StructureMappingSearchOptions = opt
        """StructureMappingSearchOptions: Options for the search."""

    def _get_max_n_atoms(self, parent: xtal.Structure, child: xtal.Structure):
        """Get the maximum supercell size of the child structure based on the parent
        structure.

        If `child_max_supercell_size` is not set, the maximum supercell size is set to
        the least common multiple of the number of atoms in the child and parent
        structures.
        """
        if self.opt.max_n_atoms is not None:
            return self.opt.max_n_atoms

        n_atoms_parent = len(parent.atom_type())
        n_atoms_child = len(child.atom_type())
        return math.lcm(n_atoms_parent, n_atoms_child)

    def _get_child_to_parent_vol(
        self,
        parent: xtal.Structure,
        child: xtal.Structure,
    ):
        max_n_atoms = self._get_max_n_atoms(parent, child)
        child_n_atoms = len(child.atom_type())
        parent_n_atoms = len(parent.atom_type())

        child_to_parent_vol = {}
        child_vol = 1
        while child_vol * child_n_atoms <= max_n_atoms:
            child_superstructure_n_atoms = child_n_atoms * child_vol
            _vol = child_superstructure_n_atoms / parent_n_atoms

            # if parent_vol is integer, then it is a valid supercell size:
            if _vol.is_integer():
                child_to_parent_vol[child_vol] = int(_vol)

            child_vol += 1

        return child_to_parent_vol

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

        if self.opt.forced_on is not None:
            _allowed = [list([x]) for x in parent.atom_type()]
            _child_types = child.atom_type()
            for parent_site_index, child_atom_index in self.opt.forced_on.items():
                child_type = _child_types[child_atom_index]
                if child_type not in _allowed[parent_site_index]:
                    invalid_forced_on_values_error(
                        parent_site_index=parent_site_index,
                        child_atom_index=child_atom_index,
                        child_type=child_type,
                        allowed_types=_allowed[parent_site_index],
                    )

        if self.opt.fix_parent:
            child_n_atoms = len(child.atom_type())
            parent_n_atoms = len(parent.atom_type())
            if child_n_atoms != parent_n_atoms:
                invalid_fix_parent_error()

        else:
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

        if self.opt.child_transformation_matrix_to_super_list is None:
            # Validate the min/max number of atoms
            if self.opt.min_n_atoms < 1:
                invalid_min_n_atoms_error(min_n_atoms=self.opt.min_n_atoms)

            _max_n_atoms = self._get_max_n_atoms(parent, child)
            if _max_n_atoms < self.opt.min_n_atoms:
                computed_msg = (
                    "(computed from lcm of atom counts)"
                    if self.opt.max_n_atoms is None
                    else ""
                )
                invalid_max_n_atoms_error(
                    min_n_atoms=self.opt.min_n_atoms,
                    max_n_atoms=_max_n_atoms,
                    computed_msg=computed_msg,
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

    def _make_T_pairs(
        self,
        parent: xtal.Structure,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
        min_n_atoms: int,
        max_n_atoms: int,
        child_T_list: Optional[list[np.ndarray]] = None,
        parent_T_list: Optional[list[np.ndarray]] = None,
    ):
        """Make a list of (T_child, T_parent) pairs for the search.

        Parameters
        ----------
        parent : xtal.Structure
            The parent structure.
        child : xtal.Structure
            The child structure.
        parent_prim : casmconfig.Prim
            The primitive parent structure.
        min_n_atoms : int
            The minimum number of atoms in the superstructures that should be included
            in the search.
        max_n_atoms : int
            The maximum number of atoms in the superstructures that should be included
            in the search.
        child_T_list : Optional[list[np.ndarray]] = None
            For the child superstructures, a list of transformation matrices
            :math:`T_{2}` to use. If None, the child superstructures are enumerated
            based on the `min_n_atoms` and `max_n_atoms` options.
        parent_T_list : Optional[list[np.ndarray]] = None
            For the parent superstructures, a list of transformation matrices
            :math:`T_{1}` to use. If None, the parent superstructures are enumerated
            based on the `min_n_atoms` and `max_n_atoms` options.

        Returns
        -------
        T_pairs: list[tuple[np.ndarray, np.ndarray]]
            List of (T_child, T_parent) pairs.

        """
        # Results, list of (T_child, T_parent) pairs
        T_pairs = []

        # Parameters
        child_crystal_point_group = xtal.make_structure_crystal_point_group(child)
        child_n_atoms = len(child.atom_type())
        child_to_parent_vol = self._get_child_to_parent_vol(
            parent=parent,
            child=child,
        )

        # If child_T_list is not provided, enumerate the child supercells
        if child_T_list is None:
            child_T_list = []

            child_superlattices = xtal.enumerate_superlattices(
                unit_lattice=child.lattice(),
                point_group=child_crystal_point_group,
                max_volume=floordiv(max_n_atoms, child_n_atoms),
                min_volume=ceildiv(min_n_atoms, child_n_atoms),
            )
            for child_superlattice in child_superlattices:
                child_T_list.append(
                    xtal.make_transformation_matrix_to_super(
                        unit_lattice=child.lattice(),
                        superlattice=child_superlattice,
                    )
                )

        # For each child superstructure...
        for child_T in child_T_list:
            child_vol = int(round(np.linalg.det(child_T)))

            # If no valid parent volume, continue
            if child_vol not in child_to_parent_vol:
                continue
            parent_vol = child_to_parent_vol[child_vol]

            # Get the list of valid parent supercells
            restricted_parent_T_list = []

            # If parent_T_list is not provided, enumerate the parent supercells
            if parent_T_list is None:
                parent_superlattices = xtal.enumerate_superlattices(
                    unit_lattice=parent.lattice(),
                    point_group=parent_prim.crystal_point_group.elements,
                    max_volume=parent_vol,
                    min_volume=parent_vol,
                )
                for parent_superlattice in parent_superlattices:
                    restricted_parent_T_list.append(
                        xtal.make_transformation_matrix_to_super(
                            unit_lattice=parent.lattice(),
                            superlattice=parent_superlattice,
                        )
                    )

            # If parent_T_list is provided, filter the parent supercells
            else:
                for parent_T in parent_T_list:
                    if int(round(np.linalg.det(parent_T))) == parent_vol:
                        restricted_parent_T_list.append(parent_T)

            # Add the (child_T, parent_T) pairs
            for parent_T in restricted_parent_T_list:
                T_pairs.append((child_T, parent_T))

        return T_pairs

    def __call__(
        self,
        parent: xtal.Structure,
        parent_prim: Optional[casmconfig.Prim],
        child: xtal.Structure,
        results_dir: pathlib.Path,
        merge: bool = False,
    ):
        """Perform the structure mapping search.

        Parameters
        ----------
        parent : xtal.Structure
            The parent structure.
        parent_prim : Optional[casmconfig.Prim]
            The parent primitive structure. If both `parent` and `parent_prim` are
            provided, the `parent` structure is used to create the fix the supercell
            being mapped to with the "fix_parent" option.
        child : xtal.Structure
            The child structure.
        results_dir : pathlib.Path
            The directory to write the results to. If the directory already exists,
            the program exits with an error.
        merge: bool = False
            If True, merge the results with existing results in the directory. If False,
            exit with error if the directory already exists.
        """
        alloy = False
        if parent_prim is None:
            parent_prim = casmconfig.Prim(
                xtal.Prim.from_atom_coordinates(structure=parent)
            )
        else:
            alloy = True
            if not self.opt.fix_parent:
                raise NotImplementedError(
                    "Mapping to a prim is only supported with the --fix-parent option."
                )

        self._supercell_set = casmconfig.SupercellSet(prim=parent_prim)

        search_results = []
        uuids = []
        chain_orbits = []

        if results_dir.exists():
            if merge is False:
                results_dir_exists_error(results_dir=results_dir)
                sys.exit(1)
            else:
                data = read_required(results_dir / "mappings.json")

                # Validate same parent and child structures:
                if alloy is False:
                    _last_parent = xtal.Structure.from_dict(data.get("parent"))
                    if not parent.is_equivalent_to(_last_parent):
                        different_parent_error()
                    _last_child = xtal.Structure.from_dict(data.get("child"))
                    if not child.is_equivalent_to(_last_child):
                        different_child_error()
                else:
                    # TODO validation
                    pass

                search_results = [
                    mapinfo.ScoredStructureMapping.from_dict(
                        data=x, prim=parent_prim.xtal_prim
                    )
                    for x in data["mappings"]
                ]
                uuids = data.get("uuids", [])

                options_data = read_required(results_dir / "options_history.json")
                last_options = StructureMappingSearchOptions.from_dict(options_data[-1])

                if (
                    self.opt.lattice_mapping_cost_method
                    != last_options.lattice_mapping_cost_method
                ):
                    different_lattice_mapping_cost_method_error()
                if (
                    self.opt.atom_mapping_cost_method
                    != last_options.atom_mapping_cost_method
                ):
                    different_atom_mapping_cost_method_error()
                if not math.isclose(
                    self.opt.lattice_cost_weight,
                    last_options.lattice_cost_weight,
                    abs_tol=1e-5,
                ):
                    different_lattice_cost_weight_error()

        if alloy is False:
            self.validate(parent, child)
        else:
            # TODO
            pass

        ## Parameters
        if alloy is False:
            _max_n_atoms = self._get_max_n_atoms(parent, child)
        else:
            _max_n_atoms = _get_max_n_atoms_for_parent_prim(
                max_n_atoms=self.opt.max_n_atoms,
                child=child,
            )
        _min_n_atoms = self.opt.min_n_atoms
        _child_T_list = self.opt.child_transformation_matrix_to_super_list
        _parent_T_list = self.opt.parent_transformation_matrix_to_super_list
        _enable_symmetry_breaking_atom_cost = self._enable_symmetry_breaking_atom_cost()
        _total_min_cost = self.opt.total_min_cost
        _total_max_cost = self.opt.total_max_cost
        _total_k_best = self.opt.total_k_best
        _enable_remove_mean_displacement = not self.opt.no_remove_mean_displacement
        _lattice_cost_weight = self.opt.lattice_cost_weight
        _lattice_mapping_min_cost = self.opt.lattice_mapping_min_cost
        _lattice_mapping_max_cost = self.opt.lattice_mapping_max_cost
        _lattice_mapping_cost_method = self.opt.lattice_mapping_cost_method
        _lattice_mapping_k_best = self.opt.lattice_mapping_k_best
        _lattice_mapping_reorientation_range = (
            self.opt.lattice_mapping_reorientation_range
        )
        _atom_cost_f = self._atom_cost_f()
        _forced_on = self.opt.forced_on if self.opt.forced_on is not None else {}
        _forced_off = self.opt.forced_off if self.opt.forced_off is not None else []
        _total_cost_f = self._total_cost_f()
        _cost_tol = self.opt.cost_tol

        ## Fixed parameters
        _infinity = 1e20
        _atom_to_site_cost_future_f = mapsearch.make_atom_to_site_cost_future

        ## Create a parent structure search data object.
        parent_search_data = mapsearch.PrimSearchData(
            prim=parent_prim.xtal_prim,
            enable_symmetry_breaking_atom_cost=_enable_symmetry_breaking_atom_cost,
        )
        init_child_structure_data = mapsearch.StructureSearchData(
            lattice=child.lattice(),
            atom_coordinate_cart=child.atom_coordinate_cart(),
            atom_type=child.atom_type(),
            override_structure_factor_group=None,
        )

        if self.opt.fix_parent:
            if alloy is False:
                I_matrix = np.eye(3, dtype="int")
                T_pairs = [(I_matrix, I_matrix)]
            else:
                # If fixing the parent, we only need one pair of transformation matrices
                # (identity for both child and parent).
                I_matrix = np.eye(3, dtype="int")
                T_parent = xtal.make_transformation_matrix_to_super(
                    superlattice=parent.lattice(),
                    unit_lattice=parent_prim.xtal_prim.lattice(),
                )
                T_pairs = [(I_matrix, T_parent)]
        else:
            if alloy is False:
                # Get a list of (T_child, T_parent) pairs
                T_pairs = self._make_T_pairs(
                    parent=parent,
                    child=child,
                    parent_prim=parent_prim,
                    min_n_atoms=_min_n_atoms,
                    max_n_atoms=_max_n_atoms,
                    child_T_list=_child_T_list,
                    parent_T_list=_parent_T_list,
                )
            else:
                raise NotImplementedError(
                    "Alloy structures are only supported with the --fix-parent option."
                )

        total = len(T_pairs)
        print(f"Beginning search over {total} parent / child superstructure pairs...")
        print()
        print("Search results so far:")
        print()
        print()
        sys.stdout.flush()

        last_child_T = None
        child_structure_data = None

        n_atoms = 0
        if len(T_pairs):
            # Get the number of atoms in the child structure
            child_vol = int(round(np.linalg.det(T_pairs[0][0])))
            child_n_atoms = len(child.atom_type())
            n_atoms = child_n_atoms * child_vol
        min_total_cost = 0.0
        max_total_cost = 0.0

        for i_pair, _pair in enumerate(T_pairs):

            child_T, parent_T = _pair
            if last_child_T is None or not np.allclose(child_T, last_child_T):
                child_structure_data = mapsearch.make_superstructure_data(
                    prim_structure_data=init_child_structure_data,
                    transformation_matrix_to_super=child_T,
                )

            child_vol = int(round(np.linalg.det(child_T)))
            child_n_atoms = len(child.atom_type())
            n_atoms = child_n_atoms * child_vol

            # Create a MappingSearch object.
            # This will hold a queue of possible mappings,
            # sorted by cost, as we generate them.
            search = mapsearch.MappingSearch(
                min_cost=_total_min_cost,
                max_cost=_total_max_cost,
                k_best=_total_k_best,
                atom_cost_f=_atom_cost_f,
                total_cost_f=_total_cost_f,
                atom_to_site_cost_future_f=_atom_to_site_cost_future_f,
                enable_remove_mean_displacement=_enable_remove_mean_displacement,
                infinity=_infinity,
                cost_tol=_cost_tol,
            )

            if self.opt.fix_parent:
                lattice_mapping = mapmethods.map_lattices_without_reorientation(
                    lattice1=parent_search_data.prim_lattice(),
                    lattice2=child_structure_data.lattice(),
                    transformation_matrix_to_super=parent_T,
                )
                F = lattice_mapping.deformation_gradient()
                if _lattice_mapping_cost_method == "isotropic_strain_cost":
                    lattice_cost = mapinfo.isotropic_strain_cost(
                        deformation_gradient=F,
                    )
                elif _lattice_mapping_cost_method == "symmetry_breaking_strain_cost":
                    lattice_cost = mapinfo.symmetry_breaking_strain_cost(
                        deformation_gradient=F,
                        lattice1_point_group=parent_search_data.prim_crystal_point_group(),
                    )
                else:
                    raise ValueError(
                        f"Unknown lattice mapping cost method: "
                        f"{_lattice_mapping_cost_method}"
                    )
                lattice_mappings = [
                    mapinfo.ScoredLatticeMapping(
                        lattice_cost=lattice_cost,
                        lattice_mapping=lattice_mapping,
                    )
                ]

            else:
                # Might be able to tighten lattice max cost limit:
                _curr_search_max = _total_max_cost
                if len(search_results) >= _total_k_best:
                    _curr_search_max = search_results[-1].total_cost()

                _curr_lattice_max = min(
                    _lattice_mapping_max_cost,
                    _curr_search_max / _lattice_cost_weight,
                )

                lattice_mappings = mapmethods.map_lattices(
                    lattice1=parent.lattice(),
                    lattice2=child_structure_data.lattice(),
                    transformation_matrix_to_super=parent_T,
                    lattice1_point_group=parent_search_data.prim_crystal_point_group(),
                    lattice2_point_group=child_structure_data.structure_crystal_point_group(),
                    min_cost=_lattice_mapping_min_cost,
                    max_cost=_curr_lattice_max,
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

                # Check if 'forced_on' values are valid.
                if len(_forced_on) > 0:
                    _allowed = lattice_mapping_data.supercell_allowed_atom_types()
                    _child_types = child_structure_data.atom_type()
                    for parent_site_index, child_atom_index in _forced_on.items():
                        child_type = _child_types[child_atom_index]
                        if child_type not in _allowed[parent_site_index]:
                            raise ValueError(
                                f"Invalid --forced-on values: "
                                f"child atom {child_atom_index} (type={child_type}) "
                                f"is not allowed to map to "
                                f"parent site {parent_site_index} "
                                f"(allowed types: {_allowed[parent_site_index]})"
                            )

                # for each lattice mapping, generate possible translations
                if not _enable_remove_mean_displacement:
                    # If mean displacement removal is disabled, then we need info
                    # on which parent/atom mappings to force on. (We could also allow
                    # generating every combination here.)
                    if len(_forced_on) == 0:
                        raise ValueError(
                            "If --no-remove-mean-displacement is set, "
                            "the --forced-on option must be set."
                        )
                    # If forced_on is set, also use parent/child pairs to generate
                    # trial translations
                    trial_translations = []
                    parent_cart = parent_search_data.prim_site_coordinate_cart()
                    child_cart = (
                        lattice_mapping_data.atom_coordinate_cart_in_supercell()
                    )
                    for parent_index, child_index in _forced_on.items():
                        trial_translations.append(
                            parent_cart[:, parent_index] - child_cart[:, child_index]
                        )
                else:
                    # Make a minimal set of trial translations
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
                        forced_on=_forced_on,
                        forced_off=_forced_off,
                    )

            while search.size():
                search.partition()

            search_results, uuids, chain_orbits = self.add_new_results(
                new_results=search.results().data(),
                existing_results=search_results,
                uuids=uuids,
                chain_orbits=chain_orbits,
                parent=parent,
                child=child,
                parent_prim=parent_prim,
                k_best=_total_k_best,
                cost_tol=_cost_tol,
            )

            self.write_results(
                search_results=search_results,
                uuids=uuids,
                parent=parent,
                child=child,
                parent_prim=parent_prim,
                results_dir=results_dir,
            )

            if len(search_results) > 0:
                min_total_cost = search_results[0].total_cost()
                max_total_cost = search_results[-1].total_cost()

            # Delete the last line
            sys.stdout.write("\033[F")  # Move cursor up one line
            sys.stdout.write("\033[K")  # Clear the line
            sys.stdout.flush()

            print(
                (
                    f"Pair: {i_pair + 1} / {total} (#atoms: {n_atoms}), "
                    f"MinTotalCost: {min_total_cost:.5f}, "
                    f"MaxTotalCost: {max_total_cost:.5f}, "
                    f"#mappings: {len(search_results)}"
                )
            )
            sys.stdout.flush()
            # pbar.update(1)

        print("DONE")
        print()
        sys.stdout.flush()
        print(f"# Results: {len(search_results)}\n")
        sys.stdout.flush()

        self.tabulate_results(
            search_results=search_results,
            uuids=uuids,
            parent=parent,
            child=child,
            parent_prim=parent_prim,
        )

        # Write the options history
        self.write_options_history(results_dir=results_dir)

        return 0

    def add_new_results(
        self,
        new_results: list[mapinfo.ScoredStructureMapping],
        existing_results: list[mapinfo.ScoredStructureMapping],
        uuids: list[str],
        chain_orbits: list[list[xtal.Structure]],
        parent: xtal.Structure,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
        k_best: int,
        cost_tol: float,
    ) -> tuple[
        list[mapinfo.ScoredStructureMapping],
        list[str],
        list[list[xtal.Structure]],
    ]:
        """Add new results to the existing search results, deduplicating them.

        Parameters
        ----------
        new_results : list[libcasm.mapping.info.ScoredStructureMapping]
            The new results to add to the existing search results.
        existing_results : list[libcasm.mapping.info.ScoredStructureMapping]
            The existing search results to which the new results will be added.
        uuids : list[str]
            The UUIDs of the existing search results.
        chain_orbits : list[list[xtal.Structure]]
            The chain orbits of the existing search results.
        parent : xtal.Structure
            The parent structure.
        child : xtal.Structure
            The child structure.
        parent_prim : casmconfig.Prim
            The parent structure, as a Prim.
        k_best : int
            The number of best results to keep after deduplication. Any approximate ties
            will also be kept.
        cost_tol : float
            The tolerance for comparing costs.

        Returns
        -------
        search_results : list[libcasm.mapping.info.ScoredStructureMapping]
            The updated list of search results after deduplication.
        uuids : list[str]
            The updated list of UUIDs corresponding to the search results.
        chain_orbits : list[list[xtal.Structure]]
            The updated list of chain orbits corresponding to the search results.

        """
        search_results = existing_results

        # Deduplicate the new results
        f_chain = self.opt.deduplication_interpolation_factors

        def make_chain(structure_mapping):
            return make_primitive_chain(
                parent_lattice=parent_prim.xtal_prim.lattice(),
                child=child,
                structure_mapping=structure_mapping,
                f_chain=f_chain,
            )

        def make_orbit(chain_prototype):
            return make_chain_orbit(
                chain_prototype=chain_prototype,
                parent_prim=parent_prim,
            )

        while len(chain_orbits) < len(search_results):
            smap = search_results[len(chain_orbits)]
            chain_orbits.append(make_orbit(make_chain(smap)))
            uuids.append(str(uuid.uuid4()))

        if len(new_results) == 0:
            return search_results, uuids, chain_orbits

        for i, smap_new in enumerate(new_results):
            primitive_chain = make_chain(smap_new)

            # Check for duplicates:
            found_duplicate = False
            i_duplicate = 0
            for smap_existing, chain_orbit_existing in zip(
                search_results, chain_orbits
            ):
                if chain_is_in_orbit(primitive_chain, chain_orbit_existing):
                    found_duplicate = True
                    break
                i_duplicate += 1

            if found_duplicate:
                smap_existing = search_results[i_duplicate]
                scel_size_new = parent_supercell_size(smap_new)
                scel_size_existing = parent_supercell_size(smap_existing)

                prefer_new = False
                if scel_size_new < scel_size_existing:
                    prefer_new = True

                # prefer smaller volume mappings
                if prefer_new:
                    search_results[i_duplicate] = smap_new
                    uuids[i_duplicate] = str(uuid.uuid4())
                    chain_orbits[i_duplicate] = make_orbit(primitive_chain)
                else:
                    continue
            else:
                search_results.append(smap_new)
                uuids.append(str(uuid.uuid4()))
                chain_orbits.append(make_orbit(primitive_chain))

        # Sort the search results and chain orbits, by total cost
        sys.stdout.flush()
        isorted = [
            x[0]
            for x in sorted(enumerate(search_results), key=lambda x: x[1].total_cost())
        ]
        search_results = [search_results[i] for i in isorted]
        uuids = [uuids[i] for i in isorted]
        chain_orbits = [chain_orbits[i] for i in isorted]

        # Keep only the k-best results
        if len(search_results) > k_best:
            next_index = k_best
            while next_index < len(search_results):
                next_cost = search_results[next_index].total_cost()
                if math.isclose(
                    search_results[k_best - 1].total_cost(), next_cost, abs_tol=cost_tol
                ):
                    next_index += 1
                else:
                    break

            search_results = search_results[:(next_index)]
            uuids = uuids[:(next_index)]
            chain_orbits = chain_orbits[:(next_index)]

        return search_results, uuids, chain_orbits

    def write_results(
        self,
        search_results: list[mapinfo.ScoredStructureMapping],
        uuids: list[str],
        parent: xtal.Structure,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
        results_dir: pathlib.Path,
    ) -> None:
        """Write the results of the search."""
        data = {
            "parent": parent.to_dict(),
            "child": child.to_dict(),
            "parent_prim": parent_prim.to_dict(),
            "mappings": [smap.to_dict() for smap in search_results],
            "uuids": [x for x in uuids],
        }
        safe_dump(
            data,
            path=results_dir / "mappings.json",
            force=True,
            quiet=True,
        )

    def write_options_history(
        self,
        results_dir: pathlib.Path,
    ) -> None:
        options = read_optional(results_dir / "options_history.json", default=[])
        options.append(self.opt.to_dict())
        safe_dump(
            options,
            path=results_dir / "options_history.json",
            force=True,
            quiet=True,
        )

    def tabulate_results(
        self,
        search_results: list[mapinfo.ScoredStructureMapping],
        uuids: list[str],
        parent: xtal.Structure,
        child: xtal.Structure,
        parent_prim: casmconfig.Prim,
    ) -> str:
        """Tabulate the results of the search."""

        prec = 5
        headers = [
            "Index",
            "TotCost",
            "LatCost",
            "AtmCost",
            "Parent Vol., Grp., #Ops",
            "Child Vol., Grp., #Ops",
            "Mult.",
            "UUID",
        ]
        f_chain = self.opt.deduplication_interpolation_factors
        child_prim = casmconfig.Prim(xtal.Prim.from_atom_coordinates(structure=child))

        data = []
        for i, scored_structure_mapping in enumerate(search_results):
            smap = scored_structure_mapping

            latmap = smap.lattice_mapping()
            T_parent = latmap.transformation_matrix_to_super()
            parent_volume = abs(int(round(np.linalg.det(T_parent))))
            T_child = make_child_transformation_matrix_to_super(
                parent_lattice=parent_prim.xtal_prim.lattice(),
                child_lattice=child.lattice(),
                structure_mapping=scored_structure_mapping,
            )
            child_volume = abs(int(round(np.linalg.det(T_child))))

            total_cost = f"{smap.total_cost():.{prec}f}"
            lattice_cost = f"{smap.lattice_cost():.{prec}f}"
            atom_cost = f"{smap.atom_cost():.{prec}f}"

            chain_orbit = make_primitive_chain_orbit(
                parent_prim=parent_prim,
                child=child,
                structure_mapping=smap,
                f_chain=f_chain,
            )
            mult = len(chain_orbit)

            parent_info = make_parent_supercell_info(
                structure_mapping=smap,
                parent_prim=parent_prim,
            )
            parent_grp = parent_info["spacegroup_type"]["international_short"]
            fg_size = parent_info["factor_group_size"]

            child_info = make_child_supercell_info(
                T_child=T_child,
                child_prim=child_prim,
            )
            child_grp = child_info["spacegroup_type"]["international_short"]
            child_fg_size = child_info["factor_group_size"]

            data.append(
                [
                    i,
                    total_cost,
                    lattice_cost,
                    atom_cost,
                    str(parent_volume) + ", " + parent_grp + ", " + str(fg_size),
                    str(child_volume) + ", " + child_grp + ", " + str(child_fg_size),
                    mult,
                    uuids[i],
                ]
            )

        print("Lattice cost method:", self.opt.lattice_mapping_cost_method)
        print("Atom cost method:", self.opt.atom_mapping_cost_method)
        print("Lattice cost weight:", self.opt.lattice_cost_weight)
        print(tabulate(data, headers=headers, tablefmt="grid"))
        print()
