"""Methods for casm-calc"""

import pathlib
import sys
import tarfile
import typing

from casm.tools.shared.json_io import read_optional, safe_dump


def delete_last_line():
    # Delete the last line
    sys.stdout.write("\033[F")  # Move cursor up one line
    sys.stdout.write("\033[K")  # Clear the line
    sys.stdout.flush()


def import_directory(
    dir: pathlib.Path,
    handler: typing.Any,
):
    """Walk a directory and report results in CASM structure format for all calculation
    subdirectories.

    Notes
    -----

    This method writes three files:

    - ``<output_base>.complete.json``: list[pathlib.Path], A list of relative paths to
      completed calculation directories.
    - ``<output_base>.incomplete.json``: list[pathlib.Path], A list of relative paths
      to incomplete calculation directories.
    - ``<output_base>.results.json`` : dict[str, dict], A dictionary where keys are
      relative paths to the calculation directories, and values are dictionaries
      containing:

      - "structure_with_properties": ``dict``, The final structure with properties.
      - "structure": ``dict``, The initial structure, if present.
      - "config": ``dict``, The initial configuration, if present.


    Parameters
    ----------
    dir : pathlib.Path
        The directory to walk. This should contain subdirectories that may be
        calculation directories. Calculation directories should not contain other
        calculation directories.
    handler : typing.Any
        An instance of a handler class that implements the methods
        `is_calcdir`, `is_complete`, `report_is_needed`, and `run_report`. This
        handler will be used to determine if a directory is a calculation directory,
        if it is complete, if a report is needed, and to run the report, creating
        the "structure_with_properties.json" file in the calculation directory. See
        :class:`~casm.tools.calc.handlers.CasmV1VaspReportHandler` for an example
        implementation.

    """
    if not dir.exists():
        raise FileNotFoundError(
            f"Error in import_directory: " f"Directory {dir} does not exist."
        )
    if not dir.is_dir():
        raise NotADirectoryError(
            f"Error in import_directory: " f"{dir} is not a directory."
        )

    start = dir

    # Walk `dir` recursively, and if a directory is a calculation directory for
    # which a report is needed, run the report in that directory.
    def _walk_directory(
        handler: typing.Any,
        dir: pathlib.Path,
        complete: typing.List[pathlib.Path],
        incomplete: typing.List[pathlib.Path],
        results: typing.Dict[str, dict],
    ):
        for path in dir.iterdir():
            if not path.is_dir():
                continue
            if path.name == "settings":
                # Skip settings directories
                continue
            if handler.is_calcdir(path=path):
                if not handler.is_complete(calcdir=path):
                    incomplete.append(path)
                else:
                    data = handler.report(calcdir=path)
                    complete.append(path)
                    results[str(path.relative_to(start))] = data

                # If a directory is a calculation directory,
                # no need to walk into its subdirectories.
                continue
            else:
                # Recursively walk subdirectories
                _walk_directory(
                    dir=path,
                    handler=handler,
                    complete=complete,
                    incomplete=incomplete,
                    results=results,
                )

        delete_last_line()
        print(f"#Complete: {len(complete)}, #Incomplete: {len(incomplete)}")
        sys.stdout.flush()

    complete = []
    incomplete = []
    results = dict()
    print(f"Looking for calculations in {dir} ...")
    print(f"#Complete: {len(complete)}, #Incomplete: {len(incomplete)}")
    sys.stdout.flush()
    _walk_directory(
        dir=dir,
        handler=handler,
        complete=complete,
        incomplete=incomplete,
        results=results,
    )

    safe_dump(
        data=[str(path.relative_to(start)) for path in complete],
        path=handler.output_path(target=dir, suffix=".complete.json"),
        force=True,
    )
    safe_dump(
        data=[str(path.relative_to(start)) for path in incomplete],
        path=handler.output_path(target=dir, suffix=".incomplete.json"),
        force=True,
    )
    safe_dump(
        data=results,
        path=handler.output_path(target=dir, suffix=".results.json"),
        force=True,
    )


# This is a work in progress...
# TODO:
# - Fix output as relative paths
# - Currently, this testing shows this slows down as it progresses.
def _import_archive(
    archive_path: pathlib.Path,
    handler: typing.Any,
):
    """Walk a tar, gzipped archive and report results in CASM structure format for all
    calculation subdirectories.

    Notes
    -----
    If there are existing results in `<archiv_path>.results.json`, these will be
    read and used to skip re-reporting calculation directories, unless updates are
    forced with `handler.update` is True.

    This method writes three files:

    - ``<output_base>.complete.json``: list[pathlib.Path], A list of completed
      calculation directories.
    - ``<output_base>.incomplete.json``: list[pathlib.Path], A list of incomplete
      calculation directories.
    - ``<output_base>.results.json``: dict[str, dict], A dictionary where keys are
      relative paths to the calculation directories, and values are dictionaries
      containing:

      - "structure_with_properties": ``dict``, The final structure with properties.
      - "structure": ``dict``, The initial structure, if present.
      - "config": ``dict``, The initial configuration, if present.


    Parameters
    ----------
    archive_path : pathlib.Path
        The archive file to walk. This should contain subdirectories that may be
        calculation directories. Calculation directories should not contain other
        calculation directories.
    handler : typing.Any
        An instance of a handler class that implements the methods
        `is_calcdir`, `is_complete`, and `report`. This handler will be used to
        determine if a directory is a calculation directory, if it is complete,
        and will run the report creating the "structure_with_properties.json" data.
        See :class:`~casm.tools.calc.handlers.CasmV1VaspReportHandler` for an example
        implementation.

    Returns
    -------
    complete : list[pathlib.Path]
        A list of completed calculation directories.
    incomplete : list[pathlib.Path]
        A list of incomplete calculation directories.
    results : dict[str, dict]
        A dictionary where keys are relative paths to the reported
        "structure_with_properties.json" files, and values are the file values,
        Python dict representations of the CASM structures with properties.

    """
    complete = []
    incomplete = []
    results = dict()

    existing_results = read_optional(
        pathlib.Path(archive_path).with_suffix(".results.json"), default=[]
    )

    print(f"Opening: {archive_path}")
    sys.stdout.flush()
    archive = tarfile.open(archive_path, "r:gz")
    handler.archive = archive
    first = True
    for member in archive.getmembers():
        if first:
            print("Looking for calculations...")
            print(f"#Complete: {len(complete)}, #Incomplete: {len(incomplete)}")
            sys.stdout.flush()
            first = False
        member_path = pathlib.Path(member.name)
        if not member.isdir():
            continue
        if "settings" in member_path.parts:
            # Skip directories containing "settings" in their path
            continue
        if not handler.is_calcdir(path=member):
            continue
        if not handler.is_complete(calcdir=member):
            incomplete.append(member_path)
            continue
        else:
            if not handler.update and str(member_path) in existing_results:
                print("found existing result for", member_path)
                print()
                complete.append(member_path)
                results[str(member_path)] = existing_results[str(member_path)]
            else:
                data = handler.report(calcdir=member)
                complete.append(member_path)
                results[str(member_path)] = data

        delete_last_line()
        print(f"#Complete: {len(complete)}, #Incomplete: {len(incomplete)}")
        sys.stdout.flush()
    # Close the archive after use:
    archive.close()

    safe_dump(
        data=[str(x) for x in complete],
        path=handler.output_path(target=archive_path, suffix=".complete.json"),
        force=True,
    )
    safe_dump(
        data=[str(x) for x in incomplete],
        path=handler.output_path(target=archive_path, suffix=".incomplete.json"),
        force=True,
    )
    safe_dump(
        data=results,
        path=handler.output_path(target=archive_path, suffix=".results.json"),
        force=True,
    )
