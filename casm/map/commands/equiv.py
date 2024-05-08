"""Find equivalent structures.
Mostly from libcasm.configuration tests/configuration/test_structure_conversions.py
"""
import libcasm.configuration as casmconfig
import libcasm.mapping.info as mapinfo
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal
import numpy as np


def make_equivalent_structures(
    structure: xtal.Structure, point_group: list[xtal.SymOp]
):
    """Given a structure and a point group, construct all equivalent structures.
    If this function is used to generate equivalent mapped structures, then the structure should be the mapped child structure from `mapmethods.make_mapped_structure` and the point group should be the point group of the parent. This function will only return structures with positive determinants, and will raise an error if a structure with negative determinant is unique under the point group symmetry operations.

    Parameters
    ----------
    structure : xtal.Structure
        The structure used to generate equivalents.
    point_group : list[xtal.SymOp]
        The point group symmetry operations to apply to the structure. To generate equivalent mapped child structures, use the point group of the parent.

    Returns
    -------
    equivalent_structures : list[xtal.Structure]
    """
    equivalent_structures = []
    negative_determinant_structures = []
    for S in point_group:
        _test = S * structure
        found = False
        if np.linalg.det(_test.lattice().column_vector_matrix()) > 0:
            for e in equivalent_structures:
                if e.is_equivalent_to(_test):
                    found = True
                    break
        else:
            # store negative determinant structures separately
            negative_determinant_structures.append(_test)
            continue
        if not found:
            equivalent_structures.append(_test)

    # check negative determinant structures
    for neg_struc in negative_determinant_structures:
        found = False
        for equiv_struc in equivalent_structures:
            if equiv_struc.is_equivalent_to(neg_struc):
                found = True
                break
        if not found:
            raise ValueError(
                f"There is a equivalent structure with a negative determinant: {i.lattice().column_vector_matrix()}"
            )
    return equivalent_structures


def make_equivalent_paths(
    path: list[xtal.Structure],
    point_group: list[xtal.SymOp],
):
    """Given a list of structures along a deformation path and a point group, construct all equivalent deformation paths.
    This function shares its logic with `make_equivalent_structures`. The point group should be the point group of the parent. This function will only return structures with positive determinants, and will raise an error if a structure with negative determinant is unique under the point group symmetry operations.

    Parameters
    ----------
    path : list[xtal.Structure]
        A deformation pathway. The structure used to generate and compare equivalents is the last item in the list.
    point_group : list[xtal.SymOp]
        The point group symmetry operations to apply to the path. To generate equivalent deformation pathways from a parent structure to a child structure, use the point group of the parent.

    Returns
    -------
    equivalent_paths : list[list[xtal.Structure]]
        The outer list contains each path and the inner list contains the structures in the path.
    """
    equivalent_paths = []
    negative_determinant_paths = []
    for S in point_group:
        _test = [S * structure for structure in path]
        found = False
        if np.linalg.det(_test[-1].lattice().column_vector_matrix()) > 0:
            for equiv_path in equivalent_paths:
                if all(
                    list(
                        structure.is_equivalent_to(_test[i])
                        for i, structure in enumerate(equiv_path)
                    )
                ):
                    found = True
                    break
        else:
            # store paths with negative determinant structures separately
            negative_determinant_paths.append(_test)
            continue
        if not found:
            equivalent_paths.append(_test)

    # check negative determinant paths
    for neg_path in negative_determinant_paths:
        found = False
        for equiv_path in equivalent_paths:
            if np.all(
                structure.is_equivalent_to(neg_path[i])
                for i, structure in enumerate(equiv_path)
            ):
                found = True
                break
        if not found:
            raise ValueError(
                "There is an equivalent path with negative determinant structures"
            )
    return equivalent_paths


def check_equiv(s1: xtal.Structure, s2: xtal.Structure, point_group: list[xtal.SymOp]):
    """Check if two structures are equivalent, including symmetry operations.
    This function checks if there are any point group operations which map s1 onto s2. If comparing two mapped child structures, the point group should be the point group of the parent crystal. This will not check for translational equivalence.

    Parameters
    ----------
    s1 : xtal.Structure
        First structure to compare.
    s2 : xtal.Structure
        Second structure to compare.
    point_group : list[xtal.SymOp]
        List of symmetry operations to apply to s2. If comparing mapped child structures, use the point group of the parent crystal.

    Returns
    -------
    equivalent : bool
        True if the structures are related by a point group symmetry operation.
    """
    equivalent = False
    for S in point_group:
        if s1.is_equivalent_to(S * s2):
            equivalent = True
            break
    return equivalent


# def make_unique_mapped_structures(
#     unmapped_structure: xtal.Structure,
#     structure_mappings: mapinfo.StructureMappingResults,
#     prim_factor_group: list[xtal.SymOp],
# ) -> list[xtal.Structure]:
#     """From libcasm.configuration tests/configuration/test_structure_conversions.py"""
#     ### Find unique mapped structures
#     unique_mapped_structures = []
#     for i, smap in enumerate(structure_mappings):
#         # print(f"~~~ {i} ~~~")
#         # print(xtal.pretty_json(smap.to_dict()))
#         # print()

#         # Construct mapped_strucutre from structure mapping
#         mapped_structure = mapmethods.make_mapped_structure(
#             structure_mapping=smap,
#             unmapped_structure=unmapped_structure,
#         )

#         # Append symmetrically unique mapped structures to list
#         if len(unique_mapped_structures) == 0:
#             unique_mapped_structures.append(mapped_structure)
#             continue
#         found = False
#         for S in prim_factor_group:
#             transformed_mapped_structure = S * mapped_structure
#             if mapped_structure.is_equivalent_to(transformed_mapped_structure):
#                 found = True
#                 break
#         if not found:
#             unique_mapped_structures.append(mapped_structure)
#     return unique_mapped_structures


# def construct_equivalent_mapped_structures(prototype_structure, prim_factor_group):
#     """From libcasm.configuration tests/configuration/test_structure_conversions.py"""
#     ### Construct equivalent mapped structures
#     equivalent_mapped_structures = []
#     for S in prim_factor_group:
#         test = S * prototype_structure
#         found = False
#         for existing in equivalent_mapped_structures:
#             if test.is_equivalent_to(existing):
#                 found = True
#                 break
#         if not found:
#             equivalent_mapped_structures.append(test)
#     return equivalent_mapped_structures


# def construct_equivalent_configurations(prim, prototype_structure):
#     """From libcasm.configuration tests/configuration/test_structure_conversions.py"""
#     ### Construct equivalent configuration in a Prim with strain and displacement DoF
#     # BCC_strain_disp_prim = casmconfig.Prim(
#     #     xtal_prims.BCC(
#     #         r=1.0,
#     #         occ_dof=["A", "B", "Va"],
#     #         local_dof=[xtal.DoFSetBasis("disp")],
#     #         global_dof=[xtal.DoFSetBasis("Hstrain")],
#     #     )
#     # )
#     prim_with_dof = xtal.Prim(
#         prim.lattice(),
#         prim.coordinate_frac(),
#         prim.occ_dof(),
#         [[xtal.DoFSetBasis("disp")]],
#         [xtal.DoFSetBasis("Hstrain")],
#         prim.occupants(),
#         "title",
#     )
#     prim = casmconfig.Prim(prim_with_dof)
#     print(prim.to_dict())
#     mapped_configuration = casmconfig.make_canonical_configuration(
#         configuration=casmconfig.Configuration.from_structure(
#             prim=prim, structure=prototype_structure
#         ),
#         in_canonical_supercell=True,
#     )

#     ### Make the fully commensurate superdupercell
#     superduperlattice = xtal.make_superduperlattice(
#         lattices=[mapped_configuration.supercell.superlattice],
#         mode="fully_commensurate",
#         point_group=prim.crystal_point_group.elements,
#     )
#     T = xtal.make_transformation_matrix_to_super(
#         superlattice=superduperlattice,
#         unit_lattice=prim.xtal_prim.lattice(),
#     )
#     superdupercell = casmconfig.make_canonical_supercell(
#         casmconfig.Supercell(
#             prim=prim,
#             transformation_matrix_to_super=T,
#         ),
#     )

#     ### Put the mapped configuration into the fully commensurate superdupercell
#     prototype_configuration = casmconfig.copy_configuration(
#         motif=mapped_configuration,
#         supercell=superdupercell,
#     )

#     ### Generate equivalent configurations
#     equivalent_configurations = casmconfig.make_equivalent_configurations(
#         configuration=prototype_configuration,
#     )
#     return equivalent_configurations


def main():
    print("equiv")


if __name__ == "__main__":
    main()
