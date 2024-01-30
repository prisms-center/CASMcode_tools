import json

import libcasm.xtal as xtal
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


# helper functions
## it is a bit confusing about which structure is actually being used with primify/symmetrize
def read_prim(prim_path, primify=False, symmetrize=False):
    p = Structure.from_file(prim_path)  # does not work for CASM prim
    if symmetrize:
        sg = SpacegroupAnalyzer(p)
        p = sg.get_refined_structure()
    prim = xtal.Prim.from_poscar_str(p.to(fmt="poscar"))
    if primify:
        prim = xtal.make_primitive(
            prim
        )  # removes occ_dof (add example to libcasm.xtal) (done)
    return prim


## add make_primitive for structures to xtal package? (already exists?)
## no title in Structure for POSCAR output?
def read_structure(structure_path, primify=False, symmetrize=False):
    prim = read_prim(structure_path, primify, symmetrize)
    structure = xtal.Structure.from_dict(
        {
            "atom_coords": prim.coordinate_frac().T.tolist(),
            "atom_type": [i[0] for i in prim.occ_dof()],
            "coordinate_mode": "Fractional",
            "lattice_vectors": prim.lattice().column_vector_matrix().tolist(),
        }
    )
    return structure


def write_maps(maps, parent, child):
    for i, m in enumerate(maps):
        print(f"writing map_{i}")
        with open(f"map_{i}.json", "w") as f:
            json.dump(
                {
                    "map": m.to_dict(),
                    "parent": parent.to_dict(),
                    "child": child.to_dict(),
                },
                f,
            )


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
