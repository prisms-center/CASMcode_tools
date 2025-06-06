"""Helper functions for the mapping tool."""

import json
from pathlib import Path

import numpy as np

import libcasm.mapping.methods as mapmethods  # remove
import libcasm.xtal as xtal


def read_prim(
    prim_path: Path, primify: bool = False, symmetrize: bool = False
) -> xtal.Prim:
    """Read a prim file into an xtal.Prim.

    If a prim.json file is used, all occupational degrees of freedom
    are preserved. If a POSCAR is used, occupants on each site will
    be fixed to the species in the POSCAR.

    Parameters
    ----------
    prim_path : Path
        The path to a prim.json or POSCAR.
    primify : bool, default=False
        If true, call xtal.make_primitive on the prim.
    symmetrize : bool, default=False
        If true, symmetrize the prim.

    Returns
    --------
    xtal.Prim

    """
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


def read_structure(
    structure_path: Path, primify: bool = False, symmetrize: bool = False
) -> xtal.Structure:
    """Read a POSCAR file into an xtal.Structure.

    Parameters
    ----------
    structure_path : Path
        The path to a POSCAR.
    primify : bool, default=False
        If true, call xtal.make_primitive on the structure.
    symmetrize : bool, default=False
        If true, symmetrize the structure.

    Returns
    -------
    xtal.Structure

    """

    if symmetrize:
        raise NotImplementedError("symmetrize structure not implemented")
    if structure_path.suffix == ".json":
        raise NotImplementedError("json format not supported for structures")
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
            f.write(xtal.pretty_json(data))


def pretty_print_maps(maps, parent, child, additional_data=[]):
    """Writes mapping data to text files.

    For each map, the parent and child are also written. Additional
    data can be added to each map by passing a list of dictionaries,
    where the index corresponds to the map index. An empty dictionary
    will be skipped.
    """
    for i, m in enumerate(maps):
        print(f"pretty printing map_{i}")
        with open(f"map_{i}.out", "w") as f:
            f.write("mapping scores\n")
            f.write("--------------\n")
            f.write(
                f"atom_cost: {round(m.atom_cost(), 3)}, "
                f"lattice_cost: {round(m.lattice_cost(), 3)} "
                f"total_cost: {round(m.total_cost(), 3)}\n"
            )
            f.write("\n")
            f.write("deformation gradient\n")
            f.write("--------------------\n")
            f.write(f"{np.round(m.lattice_mapping().deformation_gradient(), 5)}")
            f.write("\n")
            f.write("displacements\n")
            f.write("-------------\n")
            f.write(f"{np.round(m.atom_mapping().displacement(), 5)}")
            f.write("\n")
            f.write("parent\n")
            f.write("------\n")
            f.write(parent.to_json())
            f.write("\n")
            f.write("child\n")
            f.write("-----\n")
            f.write(child.to_poscar_str())
            f.write("\n")
            f.write("mapped child\n")
            f.write("-----\n")
            mapped_child = mapmethods.make_mapped_structure(m, child)
            f.write(mapped_child.to_poscar_str())
            f.write("\n")
            # just hardcode strain here for the time being
            if len(additional_data) != 0:
                f.write("strain\n")
                f.write("------\n")
                f.write(f"metric: {list(additional_data[i].keys())[0]}\n")
                f.write(f"{np.round(list(additional_data[i].values())[0], 3)}\n")


def write_structures(structures, path=".", xdatcar=False):
    path = Path(path)
    if xdatcar is False:
        for i, s in enumerate(structures):
            print(f"writing structure_{i}")
            with open(path / f"structure_{i}.vasp", "w") as f:
                f.write(s.to_poscar_str())
    else:
        print(
            "writing XDATCAR"
        )  # needs some work, not sure if actual xdatcar format+appends to existing
        for i, s in enumerate(structures):
            with open(path / "XDATCAR", "a") as f:
                f.write(s.to_poscar_str())


def strain_from_lattice_mapping(lattice_mapping, strain_metric="Hstrain"):
    """Given a lattice mapping, return the symmetry adapted strain values."""
    basis = xtal.make_symmetry_adapted_strain_basis()
    converter = xtal.StrainConverter(strain_metric, basis)
    strain_vector = converter.from_F(lattice_mapping.deformation_gradient())
    return strain_vector
