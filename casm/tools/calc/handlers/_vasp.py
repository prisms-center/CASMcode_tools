"""Vasp calculation report handlers"""

import gzip
import json
import os
import pathlib
import sys
import tarfile
import typing
from abc import ABC, abstractmethod

from casm.tools.shared.ase_utils import AseVaspTool
from casm.tools.shared.json_io import read_required, safe_dump


class ReportHandlerBase(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def is_calcdir(self, path: typing.Union[pathlib.Path, tarfile.TarInfo]) -> bool:
        """Check if a path is a calculation directory.

        Parameters
        ----------
        path: Union[pathlib.Path, tarfile.TarInfo]
            The path or archive member to check.

        Returns
        -------
        value: bool
            True if the path is a calculation directory, False otherwise.
        """
        pass

    @abstractmethod
    def is_complete(self, calcdir: pathlib.Path) -> bool:
        """Check if the calculation is complete

        Parameters
        ----------
        calcdir: Union[pathlib.Path, tarfile.TarInfo]
            The calculation directory, as a path or tarfile member.

        Returns
        -------
        value: bool
            True if the calculation is complete, False otherwise.

        """
        pass

    @abstractmethod
    def status(self, calcdir: pathlib.Path):
        """Get the status of a calculation directory

        Parameters
        ----------
        calcdir: pathlib.Path
            The calculation directory.

        Returns
        -------
        status: str
            The status of the calculation.

        """
        pass

    @abstractmethod
    def report(self, calcdir: pathlib.Path) -> dict:
        """Report calculation results as a Python dict

        Parameters
        ----------
        calcdir: pathlib.Path
            The calculation directory.

        Returns
        -------
        data: dict
            The calculated structure with properties, as a Python dict

        """
        pass

    @abstractmethod
    def output_path(
        self,
        target: pathlib.Path,
        suffix: str,
    ) -> pathlib.Path:
        """Make the path for an output file

        Parameters
        ----------
        target: pathlib.Path
            The target directory or file for which the output is constructed.
        suffix: str
            A suffix to append to the target, e.g. ".complete.json" or ".results.json",
            to create the output path.

        Returns
        -------
        path: pathlib.Path
            An output path.
        """
        pass


class ArchiveReportHandler(ReportHandlerBase):
    """A base class for handling reporting results from archive files
    (either .tar.gz or .tgz)

    This class is not meant to be used directly, but rather as a base class for
    specific report handlers that implement the methods for specific layouts.
    """

    def __init__(self):
        """

        .. rubric:: Constructor

        """
        self.archive = None
        """Optional[tarfile.TarFile]: If this handler is used to report results from a
        tar archive, this will be set to the opened tarfile. If this handler is used
        to report results from a directory, this will be None."""

    def _getmember(self, path: pathlib.Path):
        """Get a tarfile member by path or return None if it does not exist"""
        if self.archive is None:
            raise ValueError(
                "Error in CasmV1VaspReportHandler._getmember: " "archive is not set."
            )
        try:
            return self.archive.getmember(str(path))
        except KeyError:
            return None

    def _member_exists(self, path: pathlib.Path):
        """Check by path if a tarfile member exists"""
        return self._getmember(path) is not None

    def _member_isdir(self, path: pathlib.Path):
        """Check by path if a tarfile member exists and is a directory"""
        member = self._getmember(path)
        return member is not None and member.isdir()

    def _member_isfile(self, path: pathlib.Path):
        """Check by path if a tarfile member exists and is a file"""
        member = self._getmember(path)
        return member is not None and member.isfile()

    def _read_member(self, path: pathlib.Path):
        """Read a tarfile member by path and return its content as a string"""
        member = self._getmember(path)
        if member is None:
            raise FileNotFoundError(f"Member {path} not found in archive.")
        with self.archive.extractfile(member) as f:
            return f.read().decode("utf-8")

    def _load_json_member(self, path: pathlib.Path):
        """Load a tarfile JSON file member by path"""
        member = self._getmember(path)
        if member is None:
            raise FileNotFoundError(f"Member {path} not found in archive.")
        with self.archive.extractfile(member) as f:
            return json.load(f)


class VaspArchiveReportHandler(ArchiveReportHandler):
    """A base class for handling reporting results from VASP calculations stored in
    a directory or an archive file (either .tar.gz or .tgz)

    This class is not meant to be used directly, but rather as a base class for
    specific report handlers that implement the methods for specific layouts.
    """

    def __init__(
        self,
        update: bool,
        tool: typing.Optional[typing.Any],
        structure_relpath: pathlib.Path,
        config_relpath: pathlib.Path,
    ):
        """

        .. rubric:: Constructor

        Notes
        -----
        This class is a base class for handling reporting results from VASP archive
        files. It is not meant to be used directly, but rather as a base class for
        specific VASP archive handlers that implement the methods for specific VASP
        archive formats.

        Parameters
        ----------

        update: bool
            If True, the report will be run even if the
            "structure_with_properties.json" file already exists in the calculation
            directory. If False, the report will only be run if that file does not
            exist.
        tool : Any
            A tool to use for reporting results, expected to have a `report` method
            with signature ``tool.report(calc_dir: pathlib.Path)``. By default,
            an :class:`~casm.tools.shared.ase_utils.AseVaspTool` will be created.
        structure_relpath: pathlib.Path
            The relative path to the structure file from the calculation directory.
        config_relpath: pathlib.Path
            The relative path to the configuration file from the calculation directory.

        """
        super().__init__()

        if tool is None:
            tool = AseVaspTool()

        self.update = update
        """bool: If True, the report will be run even if the 
        "structure_with_properties.json" file already exists in the calculation
        directory. If False, the report will only be run if that file does not
        exist."""

        self.tool = tool
        """Any: A tool to use for reporting results, expected to have a `report` method
        with signature ``tool.report(calc_dir: pathlib.Path)``."""

        self.structure_relpath = structure_relpath
        """str: The relative path to the structure file from the calculation 
        directory."""

        self.config_relpath = config_relpath
        """str: The relative path to the configuration file from the calculation
        directory."""

        self.sentinel_name = "status.json"
        """str: The name of the sentinel file that indicates the parent
        directory is a calculation directory."""

    def is_complete(self, calcdir: pathlib.Path):
        """Check if the calculation is complete

        Notes
        -----
        A calculation directory is considered complete if it contains a "status.json"
        file with attribute ``"status": True``. Completed calculations are expected
        to have an "OUTCAR" or "OUTCAR.gz" file in a "run.final" directory.

        Parameters
        ----------
        calcdir: Union[pathlib.Path, tarfile.TarInfo]
            The calculation directory, as a path or tarfile member.

        Returns
        -------
        value: bool
            True if the calculation is complete, False otherwise.

        """
        if isinstance(calcdir, tarfile.TarInfo):
            path = pathlib.Path(calcdir.name) / self.sentinel_name
            data = self._load_json_member(path)

        else:
            data = read_required(calcdir / self.sentinel_name)

        return isinstance(data, dict) and data.get("status") == "complete"

    def status(self, calcdir: pathlib.Path):
        """Get the status of a calculation directory

        Notes
        -----
        This method reads the "status.json" file in the calculation directory and
        returns its `"status"` attribute value.

        Parameters
        ----------
        calcdir: pathlib.Path
            The calculation directory.

        Returns
        -------
        status: str
            The value of the `"status"` attribute in the "status.json" file. If the
            file does not exist a FileNotFoundError exception is raised. If the file
            cannot be read, or does not contain a `"status"` attribute, a ValueError
            exception is raised.

        """
        if isinstance(calcdir, tarfile.TarInfo):
            path = pathlib.Path(calcdir.name) / self.sentinel_name
            data = self._load_json_member(path)
        else:
            data = read_required(calcdir / self.sentinel_name)
        status = data.get("status")
        if status is None:
            raise ValueError(
                f"Invalid status for {calcdir}: 'status' attribute not found in "
                f"{self.sentinel_name}"
            )
        elif not isinstance(status, str):
            raise ValueError(
                f"Invalid status value type for {calcdir}: expected str, "
                f"got {type(status)}"
            )
        return status

    def _report_archive(self, calcdir: pathlib.Path):
        """Report the structure with properties from an archived calculation directory

        Parameters
        ----------
        calcdir: pathlib.Path
            The calculation directory.

        Returns
        -------
        data: dict
            The calculation data, as a Python dict, containing:

            - "structure_with_properties": ``dict``, The final structure with
              properties.
            - "structure": ``Optional[dict]``, The initial structure, if present.
            - "config": ``Optional[dict]``, The initial configuration, if present.


        """
        data = dict()

        calcdir_path = pathlib.Path(calcdir.name)

        structure_path = pathlib.Path(
            os.path.normpath(calcdir_path / self.structure_relpath)
        )
        print("Checking structure path:", structure_path)
        sys.stdout.flush()
        member = self._getmember(structure_path)
        if member:
            print("Found structure file")
            sys.stdout.flush()
            with self.archive.extractfile(member) as f:
                data["structure"] = json.load(f)
        print()

        config_path = pathlib.Path(os.path.normpath(calcdir_path / self.config_relpath))
        print("Checking config path:", config_path)
        sys.stdout.flush()
        member = self._getmember(config_path)
        if member:
            print("Found config file")
            sys.stdout.flush()
            with self.archive.extractfile(member) as f:
                data["config"] = json.load(f)

        # If not force updating, check if the
        # structure_with_properties.json file already exists
        if not self.update:
            structure_with_properties_path = (
                calcdir_path / "structure_with_properties.json"
            )
            member = self._getmember(structure_with_properties_path)
            if member:
                print("Found existing structure_with_properties.json in", calcdir_path)
                print()
                sys.stdout.flush()
                with self.archive.extractfile(member) as f:
                    data["structure_with_properties"] = json.load(f)
                return data

        # If force updating, or the file does not exist, find the OUTCAR and
        # report the structure
        try:
            run_final = calcdir_path / "run.final"
            outcar_file = run_final / "OUTCAR"
            outcar_member = self._getmember(outcar_file)
            if outcar_member:
                extract_dir = pathlib.Path("extract")
                extract_dir.mkdir(exist_ok=True, parents=True)

                with self.archive.extractfile(outcar_member) as f_in:
                    with open(extract_dir / "OUTCAR", "w") as f_out:
                        f_out.write(f_in.read().decode("utf-8"))
                casm_structure = self.tool.report(calc_dir=extract_dir)

                # remove extract_dir after use:
                for item in extract_dir.iterdir():
                    item.unlink()
                extract_dir.rmdir()
            else:
                outcar_gz_file = run_final / "OUTCAR.gz"
                outcar_gz_member = self._getmember(outcar_gz_file)
                if outcar_gz_member:
                    # This work pretty well:
                    extract_dir = pathlib.Path("extract")
                    extract_dir.mkdir(exist_ok=True, parents=True)

                    with self.archive.extractfile(outcar_gz_member) as g:
                        with gzip.open(g, "r") as f_in:
                            with open(extract_dir / "OUTCAR", "w") as f_out:
                                f_out.write(f_in.read().decode("utf-8"))
                    casm_structure = self.tool.report(calc_dir=extract_dir)

                    # remove extract_dir after use:
                    for item in extract_dir.iterdir():
                        item.unlink()
                    extract_dir.rmdir()
                else:
                    raise FileNotFoundError(
                        f"Neither OUTCAR nor OUTCAR.gz found in {run_final}"
                    )

            data["structure_with_properties"] = casm_structure.to_dict()
        except Exception as e:
            print(f"{calcdir}: failed to report")
            print()
            raise e

        return data

    def _report_directory(self, calcdir: pathlib.Path):
        """Report the structure with properties from a calculation directory

        Parameters
        ----------
        calcdir: pathlib.Path
            The calculation directory.

        Returns
        -------
        data: dict
            The calculation data, as a Python dict, containing:

            - "structure_with_properties": ``dict``, The final structure with
              properties.
            - "structure": ``Optional[dict]``, The initial structure, if present.
            - "config": ``Optional[dict]``, The initial configuration, if present.

        """
        data = dict()

        structure_path = calcdir / self.structure_relpath
        if structure_path.exists():
            data["structure"] = read_required(structure_path)

        config_path = calcdir / self.config_relpath
        if config_path.exists():
            data["config"] = read_required(config_path)

        structure_with_properties_path = calcdir / "structure_with_properties.json"

        if not self.update:
            if structure_with_properties_path.exists():
                data["structure_with_properties"] = read_required(
                    structure_with_properties_path
                )
                return data

        try:
            run_final = calcdir / "run.final"
            casm_structure = self.tool.report(calc_dir=run_final)
            casm_structure_data = casm_structure.to_dict()
            data["structure_with_properties"] = casm_structure_data
            safe_dump(
                data=casm_structure_data,
                path=structure_with_properties_path,
                force=True,
                quiet=True,
            )
        except Exception as e:
            print(f"{calcdir}: failed to report")
            print()
            sys.stdout.flush()
            raise e

        return data

    def report(self, calcdir: pathlib.Path):
        """Write the structure_with_properties.json file, and return its value as a
        Python dict

        Notes
        -----
        This method also writes `calcdir / structure_with_properties.json`

        Parameters
        ----------
        calcdir: pathlib.Path
            The calculation directory.

        Returns
        -------
        data: dict
            The calculation data, as a Python dict, containing:

            - "structure_with_properties": ``dict``, The final structure with
              properties.
            - "structure": ``Optional[dict]``, The initial structure, if present.
            - "config": ``Optional[dict]``, The initial configuration, if present.

        """
        if isinstance(calcdir, tarfile.TarInfo):
            return self._report_archive(calcdir=calcdir)
        else:
            return self._report_directory(calcdir=calcdir)


class CasmV1VaspReportHandler(VaspArchiveReportHandler):
    """Handler for VASP calculations stored in the CASM v1 directory structure

    This works for calculations with the following directory structure:

    .. code-block:: text

        dir/
        ├── <configuration_name>/
        │   ├── config.json
        │   ├── structure.json
        │   ├── calctype.<calc_id>/
        │   │   ├── ...
        │   │   ├── run.final/OUTCAR
        │   │   ├── run.final/OUTCAR.gz
        │   │   └── structure_with_properties.json
        ...

    When the handler is used is run on `dir`, it assumes that:

    1. All subdirectories named "calctype.<calc_id>" are calculation directories.
    2. A calculation directory with a "run.final" directory containing either
       an "OUTCAR" or "OUTCAR.gz" file is complete. Otherwise, the calculation is
       incomplete.
    3. The "config.json" and "structure.json" files are optional, but may be present
       if the calculation was setup from a CASM configuration or structure.

    If :func:`CasmV1VaspReportHandler.report` is run on a completed calculation
    directory, it will parse the results from the "OUTCAR" or "OUTCAR.gz" file, and
    store the resulting "structure_with_properties.json" file in the calculation
    directory. It will also collect the "config.json" and "structure.json" values if
    they are present.

    """

    def __init__(
        self,
        calc_id: str,
        update: bool = False,
        tool: typing.Optional[typing.Any] = None,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        calc_id : str
            The calculation ID. This is used to identify calculation directories, which
            are expected to be named "calctype." + calc_id. For example, if calc_id is
            "gga", then the calculation directories are expected to be named
            "calctype.gga".
        update: bool = False
            If True, the report will be run even if the
            "structure_with_properties.json" file already exists in the calculation
            directory. If False, the report will only be run if that file does not
            exist.
        tool : Optional[Any] = None
            A tool to use for reporting results, expected to have a `report` method
            with signature ``tool.report(calc_dir: pathlib.Path)``. By default,
            a :class:`~casm.tools.shared.ase_utils.AseVaspTool` will be created.

        """
        super().__init__(
            update=update,
            tool=tool,
            structure_relpath=pathlib.Path("..") / "structure.json",
            config_relpath=pathlib.Path("..") / "config.json",
        )

        self.calc_id = calc_id
        """str: The calculation ID. This is used to identify calculation directories,
        which are expected to be named "calctype." + calc_id. For example, if calc_id is
        "parameter_set_2", then the calculation directories are expected to be named
        "calctype.parameter_set_2"."""

    def is_calcdir(
        self,
        path: typing.Union[pathlib.Path, tarfile.TarInfo],
    ):
        """Check if a path is a calculation directory

        Notes
        -----
        If `path` is a directory and its name is ``"calctype." + self.calc_id``, then
        it is a calculation directory.

        Parameters
        ----------
        path: Union[pathlib.Path, tarfile.TarInfo]
            The path or archive member to check.

        Returns
        -------
        value: bool
            True if the path is a calculation directory, False otherwise.

        """
        if isinstance(path, tarfile.TarInfo):
            # If path is a tarfile member, check its name
            return (
                path.isdir()
                and "calctype." + self.calc_id in path.name
                and self._member_exists(pathlib.Path(path.name) / self.sentinel_name)
            )
        else:
            return (
                path.is_dir()
                and "calctype." + self.calc_id in str(path)
                and (path / self.sentinel_name).exists()
            )

    def output_path(
        self,
        target: pathlib.Path,
        suffix: str,
    ):
        """Make the path for an output file

        Parameters
        ----------
        target: pathlib.Path
            The target directory or file for which the output is constructed.
        suffix: str
            A suffix to append to the target, e.g. ".complete.json" or ".results.json",
            to create the output path.

        Returns
        -------
        path: pathlib.Path
            An output path. If `target` is ``"path/to/training_data"` and `suffix` is
            ``".results.json"``, then the output path will be
            ``"path/to/training_data.<calc_id>.results.json"``.


        """
        return target.parent / (target.name + "." + self.calc_id + suffix)


class CasmVaspReportHandler(VaspArchiveReportHandler):
    """Handler for reporting VASP calculations

    This works for calculations which can be identified by a "status.json" file
    in the calculation directory:

    .. code-block:: text

        dir/
        ├── <configuration_name>/
        │   ├── config.json
        │   ├── structure.json
        │   ├── status.json
        │   ├── ...
        │   ├── run.final/OUTCAR
        │   ├── run.final/OUTCAR.gz
        │   └── structure_with_properties.json
        ...

    When the handler is used is run on `dir`, it assumes that:

    1. Any subdirectory containing a file named "status.json" is a calculation
       directory.
    2. A calculation directory with a "status.json" file containing
       ``"status": "complete"`` in a top-level JSON attribute is complete. Otherwise,
       the calculation is incomplete. Completed calculations are expected to have
       an "OUTCAR" or "OUTCAR.gz" file in a "run.final" directory.

    If :func:`CasmVaspReportHandler.report` is run on a completed calculation
    directory, it will parse the results from the "OUTCAR" or "OUTCAR.gz" file, and
    store the resulting "structure_with_properties.json" file in the calculation
    directory.

    """

    def __init__(
        self,
        update: bool = False,
        tool: typing.Optional[typing.Any] = None,
    ):
        """

        .. rubric:: Constructor

        Parameters
        ----------
        update: bool = False
            If True, the report will be run even if the
            "structure_with_properties.json" file already exists in the calculation
            directory. If False, the report will only be run if that file does not
            exist.
        tool : Optional[Any] = None
            A tool to use for reporting results, expected to have a `report` method
            with signature ``tool.report(calc_dir: pathlib.Path)``. By default,
            an `AseVaspTool` will be created.
        """
        super().__init__(
            update=update,
            tool=tool,
            structure_relpath=pathlib.Path("structure.json"),
            config_relpath=pathlib.Path("config.json"),
        )

    def is_calcdir(
        self,
        path: typing.Union[pathlib.Path, tarfile.TarInfo],
    ):
        """Check if a path is a calculation directory

        Notes
        -----
        If `path` is a directory and its name is ``"calctype." + self.calc_id``, then
        it is a calculation directory.

        Parameters
        ----------
        path: Union[pathlib.Path, tarfile.TarInfo]
            The path or archive member to check.

        Returns
        -------
        value: bool
            True if the path is a calculation directory, False otherwise.

        """
        if isinstance(path, tarfile.TarInfo):
            return self._member_exists(pathlib.Path(path.name) / self.sentinel_name)
        else:
            return (path / self.sentinel_name).exists()

    def output_path(
        self,
        target: pathlib.Path,
        suffix: str,
    ):
        """Make the path for an output file

        Parameters
        ----------
        target: pathlib.Path
            The target directory or file for which the output is constructed.
        suffix: str
            A suffix to append to the target, e.g. ".complete.json" or ".results.json",
            to create the output path.

        Returns
        -------
        path: pathlib.Path
            An output path. If `target` is ``"path/to/training_data"` and `suffix` is
            ``".results.json"``, then the output path will be
            ``"path/to/training_data.results.json"``.


        """
        return target.parent / (target.name + suffix)
