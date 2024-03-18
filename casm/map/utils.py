"""Helper functions for the mapping tool."""
import json

import libcasm.xtal as xtal


def read_prim(prim_path, primify=False, symmetrize=False):
    if symmetrize:
        raise NotImplementedError("symmetrize prim not implemented")
    if prim_path.suffix == ".json":
        with open(prim_path, "r") as f:
            prim = xtal.Prim.from_dict(json.load(f))
    else:
        prim = xtal.Prim.from_poscar(prim_path.as_posix())
    if primify:
        prim = xtal.make_primitive(prim)
    return prim


def read_structure(structure_path, primify=False, symmetrize=False):
    if symmetrize:
        raise NotImplementedError("symmetrize structure not implemented")
    if structure_path.suffix == ".json":
        raise NotImplementedError("prim.json format not supported structure")
    structure = xtal.Structure.from_poscar(structure_path.as_posix())
    if primify:
        structure = xtal.make_primitive(structure)
    return structure


def write_maps(maps, parent, child, additional_data=[]):
    """Writes mapping data to json files.

    For each map, the parent and child are also written. Additional
    data can be added to each map by passing a list of dictionaries,
    where the index corresponds to the map index. An empty dictionary
    will be skipped.
    """
    for i, m in enumerate(maps):
        print(f"writing map_{i}")
        data = {
            "map": m.to_dict(),
            "parent": parent.to_dict(),
            "child": child.to_dict(),
        }
        # add any additional data to the output file
        if len(additional_data) > 0:
            if len(additional_data[i]) > 0:
                for k, v in additional_data[i].items():
                    data[k] = v
        with open(f"map_{i}.json", "w") as f:
            json.dump(data, f)


def write_structures(structures, xdatcar=False):
    if xdatcar is False:
        for i, s in enumerate(structures):
            print(f"writing structure_{i}")
            with open(f"structure_{i}.vasp", "w") as f:
                f.write(s.to_poscar_str())
    else:
        print(
            "writing XDATCAR"
        )  # needs some work, not sure if actual xdatcar format+appends to existing
        for i, s in enumerate(structures):
            with open(f"XDATCAR", "a") as f:
                f.write(s.to_poscar_str())


def strain_from_lattice_mapping(lattice_mapping, strain_metric="Hstrain"):
    """Given a lattice mapping, return the symmetry adapted strain values."""
    basis = xtal.make_symmetry_adapted_strain_basis()
    converter = xtal.StrainConverter(strain_metric, basis)
    strain_vector = converter.from_F(lattice_mapping.deformation_gradient())
    return strain_vector
