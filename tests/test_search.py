import argparse
from pathlib import Path

import pytest
import numpy as np

from casm.map.commands.main import parse_args
import casm.map.commands.search as search
import casm.map.utils as utils


@pytest.fixture
def hcp_ag_primitive_poscar():
    return Path(__file__).parent / "structures/HCP_Ag.vasp"


@pytest.fixture
def bcc_li_conventional_poscar():
    return Path(__file__).parent / "structures/BCC_Li_conv.vasp"


@pytest.fixture
def bcc_li_primitive_poscar():
    return Path(__file__).parent / "structures/BCC_Li.vasp"


@pytest.fixture
def fcc_li_primitive_poscar():
    return Path(__file__).parent / "structures/FCC_Li.vasp"


@pytest.fixture
def hcp_li_primitive_poscar():
    return Path(__file__).parent / "structures/HCP_Li.vasp"


def test_mask_occupants(hcp_ag_primitive_poscar, bcc_li_conventional_poscar):
    args = parse_args(
        [
            "search",
            hcp_ag_primitive_poscar.as_posix(),
            bcc_li_conventional_poscar.as_posix(),
        ]
    )
    maps_without_masking, _, _ = search.search(args)
    assert len(maps_without_masking) == 0

    args = parse_args(
        [
            "search",
            hcp_ag_primitive_poscar.as_posix(),
            bcc_li_conventional_poscar.as_posix(),
            "--mask-occupants",
        ]
    )
    maps_with_masking, _, _ = search.search(args)
    assert len(maps_with_masking) >= 1


def test_high_symmetry_strain(
    bcc_li_primitive_poscar, fcc_li_primitive_poscar, hcp_li_primitive_poscar
):
    args = parse_args(
        [
            "search",
            bcc_li_primitive_poscar.as_posix(),
            fcc_li_primitive_poscar.as_posix(),
            "--symmetry-adapted-strain",
        ]
    )
    maps, _, _ = search.search(args)
    for m in maps:
        strain = utils.strain_from_lattice_mapping(m.lattice_mapping())
        assert np.allclose(strain, [-0.00058, 0, 0.28298, 0, 0, 0], atol=1e-5)

    # HCP test failing
    args = parse_args(
        [
            "search",
            bcc_li_primitive_poscar.as_posix(),
            hcp_li_primitive_poscar.as_posix(),
            "--symmetry-adapted-strain",
        ]
    )
    maps, _, _ = search.search(args)
    for m in maps:
        strain = utils.strain_from_lattice_mapping(m.lattice_mapping())
        assert np.allclose(strain, [-0.00290, 0, -0.13974, 0, 0, 0.04468], atol=1e-5)
