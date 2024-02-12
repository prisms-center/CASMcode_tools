import json

import libcasm.xtal as xtal


# helper functions
## it is a bit confusing about which structure is actually being used with primify/symmetrize
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


## add make_primitive for structures to xtal package? (already exists?)
## no title in Structure for POSCAR output?
def read_structure(structure_path, primify=False, symmetrize=False):
    if symmetrize:
        raise NotImplementedError("symmetrize structure not implemented")
    if primify:
        raise NotImplementedError("primify structure not implemented")
    if structure_path.suffix == ".json":
        raise NotImplementedError("prim.json format not supported structure")
    structure = xtal.Structure.from_poscar(structure_path.as_posix())
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
