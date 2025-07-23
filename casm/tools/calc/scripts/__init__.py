"""Distribution of calculation scripts"""

import pathlib


def get_script_path(
    name: str,
) -> pathlib.Path:
    """Get the path to a calculation script.

    Parameters
    ----------
    name: str
        The name of the script to get the path for. For example,
        'vasp/relaxandstatic.sh'.

    Returns
    -------
    script_path: pathlib.Path
        The path to the script.
    """
    script_path = pathlib.Path(__file__).parent / pathlib.Path(name)
    if not script_path.exists():
        raise FileNotFoundError(f"Script {name} does not exist at {script_path}.")
    return script_path
