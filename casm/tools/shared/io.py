import pathlib
from typing import Optional

import libcasm.xtal as xtal
from casm.tools.shared.json_io import (
    read_required,
    safe_dump,
)
from casm.tools.shared.text_io import (
    safe_write,
)


def read_structure(
    path: pathlib.Path,
    format: Optional[str] = None,
) -> xtal.Structure:
    """Read a structure from a file.

    .. attention::

        This method does not read magnetic moments.

    Notes
    -----

    This method reads a structure from a file. For CASM and VASP files, ASE does not
    need to be installed, but for other formats, ASE is required.


    Parameters
    ----------
    path : pathlib.Path
        The path to the structure file. If the file has a suffix, it will be used to
        determine how to read the file. If the file has no suffix, or the suffix is
        ".vasp", it is read as a VASP POSCAR file, using CASM. If the suffix is
        ".json" or ".casm", it is read as a CASM Structure JSON file, using CASM.
        Otherwise, if `ASE <https://wiki.fysik.dtu.dk/ase/>`_ is installed, then the
        `ase.io.read` method will be used to read the structure.
    format : Optional[str]=None
        If not None, ignore the path suffix and read the structure using the specified
        format. If the format is "vasp" or "casm", a VASP POSCAR or CASM Structure is
        read, using CASM. For any other value, use `ase.io.read` to read the structure,
        using the specified format.

    Returns
    -------
    structure: libcasm.xtal.Structure
        A CASM Structure read from the file.

    """
    if not path.exists():
        raise FileNotFoundError(f"Structure file '{path}' does not exist.")

    def _read_structure_using_ase(path, format):
        from casm.tools.shared.ase_utils import read_structure_using_ase

        return read_structure_using_ase(path=path, format=format)

    if format is None:
        if path.suffix in ("", ".vasp"):
            # Read as a VASP POSCAR file
            structure = xtal.Structure.from_poscar_str(path.read_text())
        elif path.suffix in (".json", ".casm"):
            structure = xtal.Structure.from_dict(read_required(path))
        else:
            structure = _read_structure_using_ase(path, format)
    else:
        if format in ("vasp",):
            # Read as a VASP POSCAR file
            structure = xtal.Structure.from_poscar_str(path.read_text())
        elif format in ("casm",):
            structure = xtal.Structure.from_dict(read_required(path))
        else:
            structure = _read_structure_using_ase(path, format)

    return structure


def write_structure(
    path: pathlib.Path,
    casm_structure: xtal.Structure,
    format: Optional[str] = None,
    force: bool = False,
    quiet: bool = False,
) -> None:
    """Write a structure to a file.

    Parameters
    ----------
    path : pathlib.Path
        The path to the structure file. If the file has no suffix, or the suffix is
        ".vasp", it will be written as a VASP POSCAR file, using CASM. If the suffix
        is ".json" or ".casm", it will be written as a CASM Structure JSON file, using
        CASM. Otherwise, if `ASE <https://wiki.fysik.dtu.dk/ase/>`_ is installed, then
        the `ase.io.write` method will be used to write the structure.
    structure : libcasm.xtal.Structure
        The CASM Structure to write.
    format : Optional[str]=None
        If not None, ignore the path suffix and write the structure with specified
        format. If the format is "vasp" or "casm", a VASP POSCAR or CASM Structure is
        written, using CASM. For any other value, use `ase.io.write`
        to write the structure with the specified format.
    force : bool=False
        By default, if the file already exists, it will not be overwritten. If `force`
        is True, the file will be overwritten.
    quiet : bool=False
        By default, messages about writing the file will be printed. If `quiet` is
        True, no messages will be printed.
    """
    if not path.parent.exists():
        raise FileNotFoundError(f"Parent directory '{path.parent}' does not exist.")

    def _write_structure_using_ase(path, casm_structure, format):
        from casm.tools.shared.ase_utils import write_structure_using_ase

        return write_structure_using_ase(
            path=path,
            casm_structure=casm_structure,
            format=format,
        )

    if format is None:
        if path.suffix in ("", ".vasp"):
            text = casm_structure.to_poscar_str()
            safe_write(text, path=path, force=force, quiet=quiet)
        elif path.suffix in (".json", ".casm"):
            data = casm_structure.to_dict()
            safe_dump(data, path=path, force=force, quiet=quiet)
        else:
            _write_structure_using_ase(path, casm_structure, format)
    else:
        if format in ("vasp",):
            text = casm_structure.to_poscar_str()
            safe_write(text, path=path, force=force, quiet=quiet)
        elif format in ("casm",):
            data = casm_structure.to_dict()
            safe_dump(data, path=path, force=force, quiet=quiet)
        else:
            _write_structure_using_ase(path, casm_structure, format)
