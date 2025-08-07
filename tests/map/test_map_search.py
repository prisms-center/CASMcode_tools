import shutil
import subprocess

from casm.tools.map._StructureMappingSearch import (
    StructureMappingSearchOptions,
)
from casm.tools.shared.json_io import read_required
from libcasm.mapping.info import ScoredStructureMapping
from libcasm.xtal import Prim, Structure


def test_example_map_1_bcc_fcc(examples_dir, tmp_path):
    # copy examples_dir / "map_1" contents to tmp_path:
    shutil.copytree(examples_dir / "map_1", tmp_path, dirs_exist_ok=True)

    results_dir = tmp_path / "results"

    # Run `casm-map`
    result = subprocess.run(
        [
            "casm-map",
            "search",
            "--max-n-atoms=4",
            "--max-total-cost=0.2",
            "BCC_Li.vasp",
            "FCC_Li.vasp",
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0

    print("---stdout---")
    print(result.stdout)
    print("------------")
    print("---stderr---")
    print(result.stderr)
    print("------------")

    files = [file.name for file in results_dir.iterdir()]
    assert "mappings.json" in files

    data = read_required(results_dir / "mappings.json")

    assert "child" in data
    assert "parent" in data
    assert "mappings" in data
    assert "uuids" in data

    parent = Structure.from_dict(data["parent"])
    assert isinstance(parent, Structure)

    parent_xtal_prim = Prim.from_atom_coordinates(structure=parent)

    child = Structure.from_dict(data["child"])
    assert isinstance(child, Structure)

    mappings = [
        ScoredStructureMapping.from_dict(data=data, prim=parent_xtal_prim)
        for data in data["mappings"]
    ]
    assert len(mappings) == 20
    for mapping in mappings:
        assert isinstance(mapping, ScoredStructureMapping)

    # Check options file
    data = read_required(results_dir / "options_history.json")
    assert isinstance(data, list)
    assert len(data) == 1
    options = [StructureMappingSearchOptions.from_dict(x) for x in data]
    for opt in options:
        assert isinstance(opt, StructureMappingSearchOptions)


def test_example_map_1_bcc_hcp(examples_dir, tmp_path):
    # copy examples_dir / "map_1" contents to tmp_path:
    shutil.copytree(examples_dir / "map_1", tmp_path, dirs_exist_ok=True)

    results_dir = tmp_path / "results"

    # Run `casm-map`
    result = subprocess.run(
        [
            "casm-map",
            "search",
            "--max-n-atoms=4",
            "BCC_Li.vasp",
            "HCP_Li.vasp",
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0

    print("---stdout---")
    print(result.stdout)
    print("------------")
    print("---stderr---")
    print(result.stderr)
    print("------------")

    files = [file.name for file in results_dir.iterdir()]
    assert "mappings.json" in files

    data = read_required(results_dir / "mappings.json")

    assert "child" in data
    assert "parent" in data
    assert "mappings" in data
    assert "uuids" in data

    parent = Structure.from_dict(data["parent"])
    assert isinstance(parent, Structure)

    parent_xtal_prim = Prim.from_atom_coordinates(structure=parent)

    child = Structure.from_dict(data["child"])
    assert isinstance(child, Structure)

    mappings = [
        ScoredStructureMapping.from_dict(data=data, prim=parent_xtal_prim)
        for data in data["mappings"]
    ]
    assert len(mappings) == 52
    for mapping in mappings:
        assert isinstance(mapping, ScoredStructureMapping)

    # Check options file
    data = read_required(results_dir / "options_history.json")
    assert isinstance(data, list)
    assert len(data) == 1
    options = [StructureMappingSearchOptions.from_dict(x) for x in data]
    for opt in options:
        assert isinstance(opt, StructureMappingSearchOptions)
