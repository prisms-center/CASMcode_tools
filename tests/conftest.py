import pathlib

import pytest


@pytest.fixture
def data_dir():
    tests_dir = pathlib.Path(__file__).parent
    return tests_dir / "data"


@pytest.fixture
def get_structure(data_dir):
    def _get_structure(name):
        return (data_dir / "structures" / name).as_posix()

    return _get_structure


@pytest.fixture
def examples_dir(data_dir):
    return data_dir / "examples"


@pytest.fixture
def settings_dir(examples_dir):
    return dir / "settings"


@pytest.fixture
def dummy_vasp_potentials_dir(examples_dir):
    return examples_dir / "dummy_vasp_potentials"
