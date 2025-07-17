"""Structure input/output, optionally with `ASE <https://wiki.fysik.dtu.dk/ase/>`_"""

import pathlib
import sys
import typing

import casm.tools.shared.json_io as json_io
import casm.tools.shared.text_io as text_io
import libcasm.xtal as xtal


def read_structure(
    path: pathlib.Path,
    format: typing.Optional[str] = None,
) -> xtal.Structure:
    """Read a structure from a file.

    .. attention::

        This method does not read magnetic moments.

    Notes
    -----

    This method reads a structure from a file. For CASM and VASP files, ASE does not
    need to be installed, but for other formats, ASE is required. If ASE is not
    installed, or ASE does not recognize the file format, an error message will be
    printed and the program will exit.


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
        try:
            import ase
        except ImportError:
            print(
                f"""Cannot read {path}

CASM does not support the file format and ASE is not installed.
Install ASE for additional conversions with:

    pip install ase"""
            )
            sys.exit(1)

        from casm.tools.shared.ase_utils import read_structure_using_ase

        try:
            return read_structure_using_ase(path=path, format=format)
        except ase.io.formats.UnknownFileTypeError:
            print(
                f"""Cannot read {path}
                
Neither CASM nor ASE recognize the file format."""
            )
            sys.exit(1)

    if format is None:
        if path.suffix in ("", ".vasp"):
            # Read as a VASP POSCAR file
            structure = xtal.Structure.from_poscar_str(path.read_text())
        elif path.suffix in (".json", ".casm"):
            structure = xtal.Structure.from_dict(json_io.read_required(path))
        else:
            structure = _read_structure_using_ase(path, format)
    else:
        if format in ("vasp",):
            # Read as a VASP POSCAR file
            structure = xtal.Structure.from_poscar_str(path.read_text())
        elif format in ("casm",):
            structure = xtal.Structure.from_dict(json_io.read_required(path))
        else:
            structure = _read_structure_using_ase(path, format)

    return structure


def read_structure_traj(
    path: pathlib.Path,
    format: typing.Optional[str] = None,
) -> list[xtal.Structure]:
    """Read a structure trajectory from a file.

    .. attention::

        This method does not read magnetic moments.

    Notes
    -----

    This method reads a structure from a file. For CASM and VASP files, ASE does not
    need to be installed, but for other formats, ASE is required. If ASE is not
    installed, or ASE does not recognize the file format, an error message will be
    printed and the program will exit.


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
    structure_traj: list[libcasm.xtal.Structure]
        A CASM Structure trajectory read from the file.

    """
    if not path.exists():
        raise FileNotFoundError(f"Structure file '{path}' does not exist.")

    def _read_structure_traj_using_ase(path, format):
        try:
            import ase
        except ImportError:
            print(
                f"""Cannot read {path}

CASM does not support the file format and ASE is not installed.
Install ASE for additional conversions with:

    pip install ase"""
            )
            sys.exit(1)

        from casm.tools.shared.ase_utils import read_structure_traj_using_ase

        try:
            return read_structure_traj_using_ase(path=path, format=format)
        except ase.io.formats.UnknownFileTypeError:
            print(
                f"""Cannot read {path}

Neither CASM nor ASE recognize the file format."""
            )
            sys.exit(1)

    if format is None:
        if path.suffix in ("", ".vasp"):
            # Read as a VASP POSCAR file
            structure_traj = [xtal.Structure.from_poscar_str(path.read_text())]
        elif path.suffix in (".json", ".casm"):
            structure_traj = [
                xtal.Structure.from_dict(x) for x in json_io.read_required(path)
            ]
        else:
            structure_traj = _read_structure_traj_using_ase(path, format)
    else:
        if format in ("vasp",):
            # Read as a VASP POSCAR file
            structure_traj = [xtal.Structure.from_poscar_str(path.read_text())]
        elif format in ("casm",):
            structure_traj = [
                xtal.Structure.from_dict(x) for x in json_io.read_required(path)
            ]
        else:
            structure_traj = _read_structure_traj_using_ase(path, format)

    return structure_traj


def write_structure(
    path: pathlib.Path,
    casm_structure: xtal.Structure,
    format: typing.Optional[str] = None,
    force: bool = False,
    quiet: bool = False,
) -> None:
    """Write a structure to a file.

    Notes
    -----

    This method writes a structure to a file. For CASM and VASP files, ASE does not
    need to be installed, but for other formats, ASE is required. If ASE is not
    installed, or ASE does not recognize the file format, an error message will be
    printed and the program will exit.


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
        try:
            import ase
        except ImportError:
            print(
                f"""Cannot write {path}

CASM does not support the file format and ASE is not installed. 
Install ASE for additional conversions with:

    pip install ase"""
            )
            sys.exit(1)

        from casm.tools.shared.ase_utils import write_structure_using_ase

        try:

            return write_structure_using_ase(
                path=path,
                casm_structure=casm_structure,
                format=format,
            )
        except ase.io.formats.UnknownFileTypeError:
            print(
                f"""Cannot write {path}

Neither CASM nor ASE support the file format."""
            )
            sys.exit(1)

    if format is None:
        if path.suffix in ("", ".vasp"):
            text = casm_structure.to_poscar_str()
            text_io.safe_write(text, path=path, force=force, quiet=quiet)
        elif path.suffix in (".json", ".casm"):
            data = casm_structure.to_dict()
            json_io.safe_dump(data, path=path, force=force, quiet=quiet)
        else:
            _write_structure_using_ase(path, casm_structure, format)
    else:
        if format in ("vasp",):
            text = casm_structure.to_poscar_str()
            text_io.safe_write(text, path=path, force=force, quiet=quiet)
        elif format in ("casm",):
            data = casm_structure.to_dict()
            json_io.safe_dump(data, path=path, force=force, quiet=quiet)
        else:
            _write_structure_using_ase(path, casm_structure, format)


def write_structure_traj(
    path: pathlib.Path,
    casm_structure_traj: list[xtal.Structure],
    format: typing.Optional[str] = None,
    force: bool = False,
    quiet: bool = False,
) -> None:
    """Write a structure to a file.

    Notes
    -----

    This method writes a structure to a file. For CASM and VASP files, ASE does not
    need to be installed, but for other formats, ASE is required. If ASE is not
    installed, or ASE does not recognize the file format, an error message will be
    printed and the program will exit.


    Parameters
    ----------
    path : pathlib.Path
        The path to the structure file. If the file has no suffix, or the suffix is
        ".vasp", it will be written as a VASP POSCAR file, using CASM. If the suffix
        is ".json" or ".casm", it will be written as a CASM Structure JSON file, using
        CASM. Otherwise, if `ASE <https://wiki.fysik.dtu.dk/ase/>`_ is installed, then
        the `ase.io.write` method will be used to write the structure.
    casm_structure_traj : list[libcasm.xtal.Structure]
        The CASM Structure trajectory to write.
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

    def _write_structure_traj_using_ase(path, casm_structure_traj, format):
        try:
            import ase
        except ImportError:
            print(
                f"""Cannot write {path}

CASM does not support the file format and ASE is not installed. 
Install ASE for additional conversions with:

    pip install ase"""
            )
            sys.exit(1)

        from casm.tools.shared.ase_utils import write_structure_traj_using_ase

        try:

            return write_structure_traj_using_ase(
                path=path,
                casm_structure_traj=casm_structure_traj,
                format=format,
            )
        except ase.io.formats.UnknownFileTypeError:
            print(
                f"""Cannot write {path}

Neither CASM nor ASE support the file format."""
            )
            sys.exit(1)

    if format is None:
        if path.suffix in ("", ".vasp"):
            if len(casm_structure_traj) != 1:
                raise ValueError(
                    "Error in write_structure_traj: "
                    "VASP POSCAR files can only contain one structure."
                )
            text = casm_structure_traj[0].to_poscar_str()
            text_io.safe_write(text, path=path, force=force, quiet=quiet)
        elif path.suffix in (".json", ".casm"):
            data = [x.to_dict() for x in casm_structure_traj]
            json_io.safe_dump(data, path=path, force=force, quiet=quiet)
        else:
            _write_structure_traj_using_ase(path, casm_structure_traj, format)
    else:
        if format in ("vasp",):
            if len(casm_structure_traj) != 1:
                raise ValueError(
                    "Error in write_structure_traj: "
                    "VASP POSCAR files can only contain one structure."
                )
            text = casm_structure_traj[0].to_poscar_str()
            text_io.safe_write(text, path=path, force=force, quiet=quiet)
        elif format in ("casm",):
            data = [x.to_dict() for x in casm_structure_traj]
            json_io.safe_dump(data, path=path, force=force, quiet=quiet)
        else:
            _write_structure_traj_using_ase(path, casm_structure_traj, format)
