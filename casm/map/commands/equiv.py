"""Find equivalent structures.
Mostly from libcasm.configuration tests/configuration/test_structure_conversions.py
"""
import libcasm.configuration as casmconfig
import libcasm.mapping.info as mapinfo
import libcasm.mapping.methods as mapmethods
import libcasm.xtal as xtal


def make_unique_mapped_structures(
    unmapped_structure: xtal.Structure,
    structure_mappings: mapinfo.StructureMappingResults,
    prim_factor_group: list[xtal.SymOp],
) -> list[xtal.Structure]:
    """From libcasm.configuration tests/configuration/test_structure_conversions.py"""
    ### Find unique mapped structures
    unique_mapped_structures = []
    for i, smap in enumerate(structure_mappings):
        # print(f"~~~ {i} ~~~")
        # print(xtal.pretty_json(smap.to_dict()))
        # print()

        # Construct mapped_strucutre from structure mapping
        mapped_structure = mapmethods.make_mapped_structure(
            structure_mapping=smap,
            unmapped_structure=unmapped_structure,
        )

        # Append symmetrically unique mapped structures to list
        if len(unique_mapped_structures) == 0:
            unique_mapped_structures.append(mapped_structure)
            continue
        found = False
        for S in prim_factor_group:
            transformed_mapped_structure = S * mapped_structure
            if mapped_structure.is_equivalent_to(transformed_mapped_structure):
                found = True
                break
        if not found:
            unique_mapped_structures.append(mapped_structure)
    return unique_mapped_structures


def construct_equivalent_mapped_structures(prototype_structure, prim_factor_group):
    """From libcasm.configuration tests/configuration/test_structure_conversions.py"""
    ### Construct equivalent mapped structures
    equivalent_mapped_structures = []
    for S in prim_factor_group:
        test = S * prototype_structure
        found = False
        for existing in equivalent_mapped_structures:
            if test.is_equivalent_to(existing):
                found = True
                break
        if not found:
            equivalent_mapped_structures.append(test)
    return equivalent_mapped_structures


def construct_equivalent_configurations(prim, prototype_structure):
    """From libcasm.configuration tests/configuration/test_structure_conversions.py"""
    ### Construct equivalent configuration in a Prim with strain and displacement DoF
    # BCC_strain_disp_prim = casmconfig.Prim(
    #     xtal_prims.BCC(
    #         r=1.0,
    #         occ_dof=["A", "B", "Va"],
    #         local_dof=[xtal.DoFSetBasis("disp")],
    #         global_dof=[xtal.DoFSetBasis("Hstrain")],
    #     )
    # )
    prim_with_dof = xtal.Prim(
        prim.lattice(),
        prim.coordinate_frac(),
        prim.occ_dof(),
        [[xtal.DoFSetBasis("disp")]],
        [xtal.DoFSetBasis("Hstrain")],
        prim.occupants(),
        "title",
    )
    prim = casmconfig.Prim(prim_with_dof)
    print(prim.to_dict())
    mapped_configuration = casmconfig.make_canonical_configuration(
        configuration=casmconfig.Configuration.from_structure(
            prim=prim, structure=prototype_structure
        ),
        in_canonical_supercell=True,
    )

    ### Make the fully commensurate superdupercell
    superduperlattice = xtal.make_superduperlattice(
        lattices=[mapped_configuration.supercell.superlattice],
        mode="fully_commensurate",
        point_group=prim.crystal_point_group.elements,
    )
    T = xtal.make_transformation_matrix_to_super(
        superlattice=superduperlattice,
        unit_lattice=prim.xtal_prim.lattice(),
    )
    superdupercell = casmconfig.make_canonical_supercell(
        casmconfig.Supercell(
            prim=prim,
            transformation_matrix_to_super=T,
        ),
    )

    ### Put the mapped configuration into the fully commensurate superdupercell
    prototype_configuration = casmconfig.copy_configuration(
        motif=mapped_configuration,
        supercell=superdupercell,
    )

    ### Generate equivalent configurations
    equivalent_configurations = casmconfig.make_equivalent_configurations(
        configuration=prototype_configuration,
    )
    return equivalent_configurations


def main():
    print("equiv")


if __name__ == "__main__":
    main()
