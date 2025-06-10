import numpy as np

import libcasm.mapping.info as mapinfo
import libcasm.xtal as xtal


def make_child_transformation_matrix_to_super(
    parent_lattice: xtal.Lattice,
    child_lattice: xtal.Lattice,
    structure_mapping: mapinfo.StructureMapping,
) -> xtal.Structure:
    """Create a child superstructure from the parent and structure mapping.

    Parameters
    ----------
    parent_lattice : xtal.Lattice
        The parent lattice.
    child_lattice : xtal.Lattice
        The child lattice.
    structure_mapping : mapinfo.StructureMapping
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
    parent_lattice : xtal.Lattice
        The parent lattice.
    child : xtal.Structure
        The child structure.
    structure_mapping : mapinfo.StructureMapping
        The structure mapping between the parent and child structures.

    Returns
    -------
    xtal.Structure
        The child superstructure.
    """
    T2 = make_child_transformation_matrix_to_super(
        parent_lattice=parent_lattice,
        child_lattice=child.lattice(),
        structure_mapping=structure_mapping,
    )

    return xtal.make_superstructure(
        transformation_matrix_to_super=T2,
        structure=child,
    )
