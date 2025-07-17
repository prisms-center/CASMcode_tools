"""Methods used to implement ``casm-map``"""

import os
import sys
import typing

import numpy as np

import libcasm.configuration as casmconfig
import libcasm.configuration.io as config_io
import libcasm.mapping.info as mapinfo
import libcasm.mapping.methods as mapmethods
import libcasm.sym_info as sym_info
import libcasm.xtal as xtal


def _suppress_output(func, *args, **kwargs):
    with open(os.devnull, "w") as devnull:
        # Save the original file descriptors for stdout and stderr
        old_stdout_fd = os.dup(sys.stdout.fileno())
        old_stderr_fd = os.dup(sys.stderr.fileno())
        try:
            # Redirect stdout and stderr to /dev/null
            os.dup2(devnull.fileno(), sys.stdout.fileno())
            os.dup2(devnull.fileno(), sys.stderr.fileno())
            # Call the function
            return func(*args, **kwargs)
        finally:
            # Restore the original file descriptors
            os.dup2(old_stdout_fd, sys.stdout.fileno())
            os.dup2(old_stderr_fd, sys.stderr.fileno())
            os.close(old_stdout_fd)
            os.close(old_stderr_fd)


def _get_symgroup_classification(
    obj: typing.Any,
    symgroup: sym_info.SymGroup,
):
    """Call symgroup_to_dict_with_group_classification with suppressed spglib output."""
    return _suppress_output(
        config_io.symgroup_to_dict_with_group_classification,
        obj,
        symgroup,
    )


def make_child_transformation_matrix_to_super(
    parent_lattice: xtal.Lattice,
    child_lattice: xtal.Lattice,
    structure_mapping: mapinfo.StructureMapping,
) -> xtal.Structure:
    """Create the transformation matrix to a child superstructure from the parent and
    structure mapping.

    Parameters
    ----------
    parent_lattice : libcasm.xtal.Lattice
        The parent lattice.
    child_lattice : libcasm.xtal.Lattice
        The child lattice.
    structure_mapping : libcasm.mapping.info.StructureMapping
        The structure mapping between the parent and child structures.

    Returns
    -------
    child_transformation_matrix_to_super : np.ndarray
        The transformation matrix from the child structure to the superstructure
        that was mapped.
    """
    # F * L_1 * T_1 * N = L_2 * T_2

    smap = structure_mapping
    lmap = smap.lattice_mapping()

    F = lmap.deformation_gradient()
    L1 = parent_lattice.column_vector_matrix()
    L2 = child_lattice.column_vector_matrix()
    T1 = lmap.transformation_matrix_to_super()
    N = lmap.reorientation()

    return np.linalg.solve(L2, F @ L1 @ T1 @ N)


def make_child_superstructure(
    parent_lattice: xtal.Lattice,
    child: xtal.Structure,
    structure_mapping: mapinfo.StructureMapping,
) -> xtal.Structure:
    """Create a child superstructure from the parent and structure mapping.

    Parameters
    ----------
    parent_lattice : libcasm.xtal.Lattice
        The parent lattice.
    child : libcasm.xtal.Structure
        The child structure.
    structure_mapping : libcasm.mapping.info.StructureMapping
        The structure mapping between the parent and child structures.

    Returns
    -------
    libcasm.xtal.Structure
        The child superstructure.
    """
    T2 = make_child_transformation_matrix_to_super(
        parent_lattice=parent_lattice,
        child_lattice=child.lattice(),
        structure_mapping=structure_mapping,
    )

    T2_int = np.round(T2).astype(int)

    return xtal.make_superstructure(
        transformation_matrix_to_super=T2_int,
        structure=child,
    )


def make_mapped_child_superstructure(
    parent_lattice: xtal.Lattice,
    child: xtal.Structure,
    structure_mapping: mapinfo.StructureMapping,
    f: float,
) -> xtal.Structure:
    """Create a mapped child superstructure from the parent and structure mapping.

    Parameters
    ----------
    parent_lattice : xtal.Lattice
        The parent lattice.
    child : xtal.Structure
        The child structure.
    structure_mapping : libcasm.mapping.info.StructureMapping
        The structure mapping between the parent and child structures.
    f : float
        The interpolation factor, where a value of 0.0 corresponds to the ideal parent
        superstructure and a value of 1.0 corresponds to the mapped child
        superstructure.

    Returns
    -------
    xtal.Structure
        The mapped child superstructure.
    """
    s = mapmethods.make_mapped_structure(
        structure_mapping=structure_mapping.interpolated(f),
        unmapped_structure=make_child_superstructure(
            parent_lattice=parent_lattice,
            child=child,
            structure_mapping=structure_mapping,
        ),
    )
    return xtal.make_structure_within(s)


def _make_primitive_structure(init_structure: xtal.Structure):
    xtal_prim = xtal.make_primitive_prim(
        xtal.Prim.from_atom_coordinates(structure=init_structure)
    )
    prim_structure = xtal.Structure(
        lattice=xtal_prim.lattice(),
        atom_coordinate_frac=xtal_prim.coordinate_frac(),
        atom_type=[occ[0] for occ in xtal_prim.occ_dof()],
    )
    return prim_structure


def is_equivalent_chain(
    chain_A: list[xtal.Structure],
    chain_B: list[xtal.Structure],
) -> bool:
    """Check if two chains of structures are equivalent.

    Parameters
    ----------
    chain_A : list[libcasm.xtal.Structure]
        The first chain of structures. Expects all structures to be in
        canonical form.

    chain_B : list[libcasm.xtal.Structure]
        The second chain of structures. Expects all structures to be in
        canonical form.

    Returns
    -------
    True if the chains are equivalent, False otherwise.
    """
    if len(chain_A) != len(chain_B):
        return False
    for s_A, s_B in zip(chain_A, chain_B):
        if not s_A.is_equivalent_to(s_B):
            return False
    return True


def chain_is_in_orbit(
    chain: list[xtal.Structure],
    orbit: list[list[xtal.Structure]],
) -> bool:
    """Check if a chain of structures is in a given orbit.

    Parameters
    ----------
    chain : list[libcasm.xtal.Structure]
        The chain of structures. Expects all structures to be in
        canonical form.

    orbit : list[list[libcasm.xtal.Structure]]
        The orbit of chains of structures. Expects all structures to be in
        canonical form.

    Returns
    -------
    True if the chain is in the orbit, False otherwise.
    """
    for existing in orbit:
        if is_equivalent_chain(chain, existing):
            return True
    return False


def make_primitive_chain(
    parent_lattice: xtal.Lattice,
    child: xtal.Structure,
    structure_mapping: mapinfo.StructureMapping,
    f_chain: list[float],
) -> list[xtal.Structure]:
    """Make a chain of primitive structures interpolating between the parent and child.

    Parameters
    ----------
    parent_lattice : xtal.Lattice
        The parent lattice.

    child : xtal.Structure
        The child structure.

    structure_mapping : libcasm.mapping.info.StructureMapping
        The structure mapping between the parent and child structures.

    f_chain : list[float]
        The interpolation factors for the structures in the chain, where a value of 0.0
        corresponds to the ideal parent structure and a value of 1.0 corresponds
        to the mapped child structure.

    Returns
    -------
    primitive_chain : list[xtal.Structure]
        A list of primitive structures interpolating between the parent and child.

    """
    unmapped_structure = make_child_superstructure(
        parent_lattice=parent_lattice,
        child=child,
        structure_mapping=structure_mapping,
    )

    def _make_interpolated_primitive_structure(f):
        """Make an interpolated, primitive, prototype structure."""
        s = mapmethods.make_mapped_structure(
            structure_mapping=structure_mapping.interpolated(f),
            unmapped_structure=unmapped_structure,
        )
        return _make_primitive_structure(xtal.make_structure_within(s))

    primitive_chain = []
    for f in f_chain:
        primitive_chain.append(_make_interpolated_primitive_structure(f))
    return primitive_chain


def make_chain_orbit(
    chain_prototype: list[xtal.Structure],
    parent_prim: casmconfig.Prim,
) -> list[list[xtal.Structure]]:
    """Make an orbit of chains of structures from a prototype chain.

    Parameters
    ----------
    chain_prototype : list[libcasm.xtal.Structure]
        A prototype chain of structures.
    parent_prim : casmconfig.Prim
        Parent structure, as a Prim. The factor group of this Prim will be used
        to generate the orbit of equivalent chains.

    Returns
    -------
    chain_orbit : list[list[libcasm.xtal.Structure]]
        A list of chains, where each element is a distinct chain that is equivalent by a
        symmetry operation in the factor group of the parent Prim.
    """
    chain_orbit = []
    for op in parent_prim.factor_group.elements:
        chain = []
        for structure in chain_prototype:
            chain.append(xtal.make_canonical_structure(op * structure))
        if not chain_is_in_orbit(chain, chain_orbit):
            chain_orbit.append(chain)

    return chain_orbit


def make_primitive_chain_orbit(
    parent_prim: casmconfig.Prim,
    child: xtal.Structure,
    structure_mapping: mapinfo.StructureMapping,
    f_chain: list[float],
) -> list[list[xtal.Structure]]:
    """Make an orbit of chains of primitive structures mapping between the parent and
    child.

    Parameters
    ----------
    parent_prim : casmconfig.Prim
        Parent structure, as a Prim. The factor group of this Prim will be used
        to generate the orbit of equivalent chains.
    child : xtal.Structure
        The child structure.
    structure_mapping : libcasm.mapping.info.StructureMapping
        The structure mapping between the parent and child structures.
    f_chain : list[float]
        The interpolation factors for the structures in the chain, where a value of 0.0
        corresponds to the ideal parent structure and a value of 1.0 corresponds
        to the mapped child structure.

    Returns
    -------
    primitive_chain_orbit : list[list[xtal.Structure]]
        A list of chains, where each element is a distinct chain that is equivalent by a
        symmetry operation in the factor group of the parent Prim. The structures in
        each chain are made primitive.

    """
    return make_chain_orbit(
        chain_prototype=make_primitive_chain(
            parent_lattice=parent_prim.xtal_prim.lattice(),
            child=child,
            structure_mapping=structure_mapping,
            f_chain=f_chain,
        ),
        parent_prim=parent_prim,
    )


def parent_supercell_factor_group_size(
    structure_mapping: mapinfo.StructureMapping,
    supercell_set: casmconfig.SupercellSet,
) -> int:
    """Get the size of the factor group of the parent supercell.

    Parameters
    ----------
    structure_mapping : libcasm.mapping.info.StructureMapping
        The structure mapping between the parent and child structures.
    supercell_set : casmconfig.SupercellSet
        A :class:`~libcasm.configuration.SupercellSet` object that allows re-using
        supercell information.

    Returns
    -------
    factor_group_size : int
        The size of the factor group of the supercell corresponding to the parent
        superstructure.

    """
    T = structure_mapping.lattice_mapping().transformation_matrix_to_super()
    T_int = np.round(T).astype(int)
    record = supercell_set.add_by_transformation_matrix_to_super(
        transformation_matrix_to_super=T_int,
    )
    return len(record.supercell.factor_group.elements)


def make_supercell_info(
    supercell: casmconfig.Supercell,
) -> dict:
    """Make a dictionary with information about the supercell.

    Parameters
    ----------
    supercell : casmconfig.Supercell
        A supercell.

    Returns
    -------
    data: dict
        A dictionary containing information about the supercell, including:

        - volume: The volume of the supercell relative to the prim.
        - factor_group_size: The size of the factor group of the supercell.
        - spacegroup_type: The space group type of the supercell, determined by spglib
          from the symmetry determined by CASM.

    """
    info = {}
    T_int = supercell.transformation_matrix_to_super
    info["volume"] = round(int(np.linalg.det(T_int)))
    default_config = casmconfig.Configuration(supercell=supercell)
    symgroup = supercell.factor_group
    data = _get_symgroup_classification(
        obj=default_config,
        symgroup=symgroup,
    )
    info["factor_group_size"] = len(symgroup.elements)

    # This handles getting the supercell crystal symmetry
    classification = data["group_classification"]
    if classification.get("spacegroup_type_from_casm_symmetry") is None:
        info["spacegroup_type"] = data["group_classification"]["spacegroup_type"]
    else:
        info["spacegroup_type"] = classification["spacegroup_type_from_casm_symmetry"]
    return info


def make_parent_supercell_info(
    structure_mapping: mapinfo.StructureMapping,
    parent_prim: casmconfig.Prim,
) -> dict:
    """Make a dictionary with information about the parent supercell.

    Parameters
    ----------
    structure_mapping : libcasm.mapping.info.StructureMapping
        The structure mapping between the parent and child structures.
    parent_prim : casmconfig.Prim
        The parent structure, as a Prim.

    Returns
    -------
    data: dict
        A dictionary containing information about the supercell, including:

        - volume: The volume of the supercell relative to `child_prim`.
        - factor_group_size: The size of the factor group of the supercell.
        - spacegroup_type: The space group type of the supercell, determined by spglib
          from the symmetry determined by CASM.

    """
    T = structure_mapping.lattice_mapping().transformation_matrix_to_super()
    T_int = np.round(T).astype(int)
    supercell = casmconfig.Supercell(
        prim=parent_prim,
        transformation_matrix_to_super=T_int,
    )
    return make_supercell_info(supercell=supercell)


def make_child_supercell_info(
    T_child: np.ndarray,
    child_prim: casmconfig.Prim,
) -> dict:
    """Make a dictionary with information about the child supercell.

    Parameters
    ----------
    T_child : np.ndarray[np.int[3,3]]
        The transformation matrix to the child superstructure, :math:`T_2`.
    child_prim : casmconfig.Prim
        The child structure, as a Prim.

    Returns
    -------
    data: dict
        A dictionary containing information about the supercell, including:

        - volume: The volume of the supercell relative to `child_prim`.
        - factor_group_size: The size of the factor group of the supercell.
        - spacegroup_type: The space group type of the supercell, determined by spglib
          from the symmetry determined by CASM.

    """
    T_int = np.round(T_child).astype(int)
    supercell = casmconfig.Supercell(
        prim=child_prim,
        transformation_matrix_to_super=T_int,
    )
    return make_supercell_info(supercell=supercell)


def make_chain_info(
    chain: list[xtal.Structure],
    parent_prim: casmconfig.Prim,
):
    """Make a dictionary with information about the structures in a chain.

    Parameters
    ----------
    chain: list[libcasm.xtal.Structure]
        A chain of structures.
    parent_prim : casmconfig.Prim
        The parent structure, as a Prim.

    Returns
    -------
    chain_info: list[dict]
        A ``list[dict]``, with one dictionary for each structure in the chain,
        including:

        - volume: The volume of the structure relative to `parent_prim`.
        - factor_group_size: The size of the factor group of the structure.
        - spacegroup_type: The space group type of the structure, determined by spglib
          from the symmetry determined by CASM.

    """
    chain_info = []
    vol_parent = abs(parent_prim.xtal_prim.lattice().volume())
    for structure in chain:
        info = {}

        # Get the volume w.r.t. the prim of the interpolated structure
        vol_structure = abs(structure.lattice().volume())
        info["volume"] = vol_structure / vol_parent

        # Get the space group of the interpolated structure
        tmp_prim = casmconfig.Prim(
            xtal_prim=xtal.Prim.from_atom_coordinates(structure=structure)
        )
        data = _get_symgroup_classification(
            obj=tmp_prim,
            symgroup=tmp_prim.factor_group,
        )
        info["factor_group_size"] = len(tmp_prim.factor_group.elements)

        classification = data["group_classification"]
        if classification.get("spacegroup_type_from_casm_symmetry") is None:
            info["spacegroup_type"] = data["group_classification"]["spacegroup_type"]
        else:
            info["spacegroup_type"] = classification[
                "spacegroup_type_from_casm_symmetry"
            ]

        chain_info.append(info)
    return chain_info


def parent_supercell_size(
    structure_mapping: mapinfo.StructureMapping,
) -> int:
    """Return the size of the parent superstructure the child superstructure was mapped
    to, relative to the size of the parent structure."""
    T = structure_mapping.lattice_mapping().transformation_matrix_to_super()
    return abs(int(round(np.linalg.det(T))))


def child_supercell_size(
    parent_lattice: xtal.Lattice,
    child: xtal.Structure,
    structure_mapping: mapinfo.StructureMapping,
) -> int:
    """Return the size of the child superstructure that was mapped to the parent
    superstructure, relative to the size of the child structure."""
    T_child = make_child_transformation_matrix_to_super(
        parent_lattice=parent_lattice,
        child_lattice=child.lattice(),
        structure_mapping=structure_mapping,
    )
    return abs(int(round(np.linalg.det(T_child))))
