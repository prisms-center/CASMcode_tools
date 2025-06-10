import shutil
import subprocess

from casm.tools.map.main import main
from casm.tools.shared.json_io import read_required
from libcasm.mapping.info import ScoredStructureMapping
from libcasm.xtal import Prim, Structure


def test_1(get_structure, tmp_path):

    main(
        argv=[
            "casm-map",
            get_structure("BCC_Li.vasp"),
            get_structure("FCC_Li.vasp"),
        ],
        working_dir=tmp_path,
    )

    files = [file.name for file in tmp_path.iterdir()]
    assert "map_0.json" in files

    data = read_required(tmp_path / "map_0.json")

    assert "child" in data
    assert "parent" in data
    assert "map" in data

    prim = Prim.from_dict(data["parent"])
    assert isinstance(prim, Prim)

    child = Structure.from_dict(data["child"])
    assert isinstance(child, Structure)

    mapping = ScoredStructureMapping.from_dict(data=data["map"], prim=prim)
    assert isinstance(mapping, ScoredStructureMapping)


def test_example_map_1(examples_dir, tmp_path):

    # copy examples_dir / "map_1" contents to tmp_path:
    shutil.copytree(examples_dir / "map_1", tmp_path)

    # Run `casm-map`
    result = subprocess.run(
        ["casm-map", "BCC_Li.vasp", "FCC_Li.vasp"],
        check=True,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )

    print("---stdout---")
    print(result.stdout)
    print("------------")
    print("---stderr---")
    print(result.stderr)
    print("------------")

    files = [file.name for file in tmp_path.iterdir()]
    assert "map_0.json" in files

    data = read_required(tmp_path / "map_0.json")

    assert "child" in data
    assert "parent" in data
    assert "map" in data

    prim = Prim.from_dict(data["parent"])
    assert isinstance(prim, Prim)

    child = Structure.from_dict(data["child"])
    assert isinstance(child, Structure)

    mapping = ScoredStructureMapping.from_dict(data=data["map"], prim=prim)
    assert isinstance(mapping, ScoredStructureMapping)


# def test_mask_occupants_1(get_structure, tmp_path):
#     main(
#         argv=[
#             "casm-map",
#             "search",
#             get_structure("HCP_Ag.vasp"),
#             get_structure("BCC_Li_conv.vasp"),
#         ],
#         working_dir=tmp_path,
#     )
#
#     files = [file.name for file in tmp_path.iterdir()]
#     assert "map_0.json" not in files
#
#
# def test_mask_occupants_2(get_structure, tmp_path):
#     main(
#         argv=[
#             "casm-map",
#             "search",
#             get_structure("HCP_Ag.vasp"),
#             get_structure("BCC_Li_conv.vasp"),
#             "--mask-occupants",
#         ],
#         working_dir=tmp_path,
#     )
#
#     files = [file.name for file in tmp_path.iterdir()]
#     assert "map_0.json" in files
