import pytest


def structure_file(shared_datadir, name):
    return shared_datadir / "structures" / name


@pytest.fixture
def get_structure(shared_datadir):
    def _get_structure(name):
        return structure_file(shared_datadir, name).as_posix()

    return _get_structure
