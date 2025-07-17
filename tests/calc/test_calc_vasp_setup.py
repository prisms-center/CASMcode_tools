import os
import shutil
import subprocess

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


def _run_calc_vasp_setup_1(structure_file, examples_dir, tmp_path):
    # copy examples_dir / "map_1" contents to tmp_path:
    shutil.copytree(examples_dir / "calc_vasp_setup_1", tmp_path, dirs_exist_ok=True)

    vasp_pp_path = examples_dir / "dummy_vasp_potentials/"
    os.environ["VASP_PP_PATH"] = str(vasp_pp_path.resolve())

    calc_dir = tmp_path / "run.1"
    settings_dir = examples_dir / "settings" / "calctype.MgCd_pbe"

    # Run `casm-map`
    result = subprocess.run(
        [
            "casm-calc",
            "vasp",
            "setup",
            structure_file.as_posix(),
            calc_dir.as_posix(),
            settings_dir.as_posix(),
        ],
        check=True,
        capture_output=True,
        text=True,
        cwd=tmp_path,
    )
    assert result.returncode == 0

    assert calc_dir.exists()

    files = [file.name for file in calc_dir.iterdir()]
    expected_files = ["INCAR", "POSCAR", "POTCAR", "ase-sort.dat", "KPOINTS"]
    for expected_file in expected_files:
        assert expected_file in files

    structure = structure_io.read_structure(calc_dir / "POSCAR")
    assert isinstance(structure, xtal.Structure)
    assert structure.atom_type() == ["Cd", "Cd", "Mg", "Mg"]


@pytest.mark.requires_ase
def test_calc_vasp_setup_1_from_casm_structure(examples_dir, tmp_path):
    structure_file = tmp_path / "structure.json"
    _run_calc_vasp_setup_1(structure_file, examples_dir, tmp_path)


@pytest.mark.requires_ase
def test_calc_vasp_setup_1_from_POSCAR(examples_dir, tmp_path):
    structure_file = tmp_path / "POSCAR"
    _run_calc_vasp_setup_1(structure_file, examples_dir, tmp_path)


@pytest.mark.requires_ase
def test_calc_vasp_setup_1_from_xyz(examples_dir, tmp_path):
    structure_file = tmp_path / "structure.xyz"
    _run_calc_vasp_setup_1(structure_file, examples_dir, tmp_path)
