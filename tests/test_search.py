import argparse
from pathlib import Path

import pytest

from casm.map.commands.main import parse_args
import casm.map.commands.search as search

@pytest.fixture
def hcp_ag_primitive_poscar():
    return Path(__file__).parent/'structures/HCP_Ag.vasp'

@pytest.fixture
def bcc_li_conventional_poscar():
    return Path(__file__).parent/'structures/BCC_Li_conv.vasp'

def test_mask_occupants(hcp_ag_primitive_poscar, bcc_li_conventional_poscar):
    args = parse_args([
        "search",
        hcp_ag_primitive_poscar.as_posix(),
        bcc_li_conventional_poscar.as_posix()
    ])
    maps_without_masking, _, _ = search.search(args)
    assert len(maps_without_masking) == 0

    args = parse_args([
        "search",
        hcp_ag_primitive_poscar.as_posix(),
        bcc_li_conventional_poscar.as_posix(),
        "--mask-occupants"])
    maps_with_masking, _, _ = search.search(args)
    assert len(maps_with_masking) >= 1
