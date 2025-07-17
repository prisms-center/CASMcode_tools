import os
import subprocess

import numpy as np
import pytest

import casm.tools.shared.structure_io as structure_io
import libcasm.xtal as xtal


def print_run_output(result):
    """Prints the stdout and stderr of a subprocess run."""
    print("---stdout---")
    print(result.stdout)
    print("------------")
    print("---stderr---")
    print(result.stderr)
    print("------------")


@pytest.mark.requires_ase
def test_calc_vasp_report_1(examples_dir, tmp_path):
    vasp_pp_path = examples_dir / "dummy_vasp_potentials/"
    os.environ["VASP_PP_PATH"] = str(vasp_pp_path.resolve())

    calc_dir = examples_dir / "calc_vasp_report_1" / "run.1"

    # Run `casm-map`
    result = subprocess.run(
        [
            "casm-calc",
            "vasp",
            "report",
            calc_dir.as_posix(),
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0

    files = [file.name for file in tmp_path.iterdir()]
    print(files)
    expected_files = ["structure_with_properties.json"]
    for expected_file in expected_files:
        assert expected_file in files

    structure = structure_io.read_structure(tmp_path / "structure_with_properties.json")
    assert isinstance(structure, xtal.Structure)
    assert structure.atom_type() == ["Cd", "Cd", "Mg", "Mg"]

    global_properties = structure.global_properties()
    assert "energy" in global_properties
    assert isinstance(global_properties["energy"], np.ndarray)
    assert global_properties["energy"].shape == (1, 1)

    atom_properties = structure.atom_properties()
    assert "force" in atom_properties
    assert isinstance(atom_properties["force"], np.ndarray)
    assert atom_properties["force"].shape == (3, 4)


@pytest.mark.requires_ase
def test_calc_vasp_report_1_traj(examples_dir, tmp_path):
    vasp_pp_path = examples_dir / "dummy_vasp_potentials/"
    os.environ["VASP_PP_PATH"] = str(vasp_pp_path.resolve())

    calc_dir = examples_dir / "calc_vasp_report_1" / "run.1"

    # Run `casm-map`
    result = subprocess.run(
        [
            "casm-calc",
            "vasp",
            "report",
            calc_dir.as_posix(),
            "--traj",
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0
    # print_run_output(result)

    files = [file.name for file in tmp_path.iterdir()]
    print(files)
    expected_files = ["structure_with_properties.traj.json"]
    for expected_file in expected_files:
        assert expected_file in files

    structure_traj = structure_io.read_structure_traj(
        tmp_path / "structure_with_properties.traj.json"
    )
    assert isinstance(structure_traj, list)
    for structure in structure_traj:
        assert isinstance(structure, xtal.Structure)
        assert structure.atom_type() == ["Cd", "Cd", "Mg", "Mg"]

        global_properties = structure.global_properties()
        assert "energy" in global_properties
        assert isinstance(global_properties["energy"], np.ndarray)
        assert global_properties["energy"].shape == (1, 1)

        atom_properties = structure.atom_properties()
        assert "force" in atom_properties
        assert isinstance(atom_properties["force"], np.ndarray)
        assert atom_properties["force"].shape == (3, 4)
