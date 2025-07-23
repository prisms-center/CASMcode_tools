"""Implements ``casm-calc vasp ...``"""

import importlib.util
import os
import pathlib
import sys
import typing


def _get_format(args):
    if args.format is not None:
        return args.format
    return None


def _validate_ase():
    if importlib.util.find_spec("ase") is None:
        print(
            """ASE is not installed. To use `casm-calc vasp` install ASE with:

    pip install ase"""
        )
        sys.exit(1)


def _validate_vasp_pp_path():
    if "VASP_PP_PATH" not in os.environ:
        print(
            """Please set the environment variable VASP_PP_PATH to the directory
containing the VASP POTCAR files. It should contain the directories
potpaw_PBE, potpaw, and potpaw_GGA, as necessary."""
        )
        sys.exit(1)


def vasp_setup(args):
    """Implements ``casm-calc vasp setup ...``

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line. Uses:

        - `args.settings`: pathlib.Path, The directory containing the VASP settings
          (expected to contain a `calc.json` file with settings to be passed to the
          ase.calculators.vasp.Vasp constructor, and template INCAR and KPOINTS files).
        - `args.structure`: pathlib.Path, The structure file to use for the calculation.
          For CASM and VASP files, ASE does not need to be installed, but for other
          formats, ASE is required.
        - `args.format`: str, optional, The format for reading the structure file. If
          not specified, the suffix of the file is used to determine the format. If the
          file has no suffix or a '.vasp' suffix, it is read as a VASP POSCAR file. If
          it has a '.json' or '.casm' suffix, it is read as a CASM Structure JSON file.
          Otherwise, if ASE is installed, it is read using `ase.io.read`.
        - `args.calcdir`: pathlib.Path, The directory to write the calculation input
          files to. This directory will be created if it does not exist.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    _validate_ase()
    _validate_vasp_pp_path()

    from casm.tools.shared.ase_utils import AseVaspTool
    from casm.tools.shared.structure_io import read_structure

    tool = AseVaspTool(calctype_settings_dir=args.settings)
    tool.setup(
        casm_structure=read_structure(path=args.structure, format=_get_format(args)),
        calc_dir=args.calcdir,
    )
    return 0


def vasp_report(args):
    """Implements ``casm-calc vasp report ...``

    Notes
    -----

    On success, this writes a file named `structure_with_properties.json` containing
    a CASM structure for the final state of the calculation results. If the option
    `--traj` is used, this instead writes a file named
    `structure_with_properties.traj.json` containing a list of CASM Structures, one for
    each step in the calculation.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line. Uses:

        - `args.calcdir`: pathlib.Path, The directory containing the VASP calculation
          results.
        - `args.traj`: If set, saves the trajectory of structures instead of just the
          final structure.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    _validate_ase()

    from casm.tools.shared.ase_utils import AseVaspTool
    from casm.tools.shared.json_io import safe_dump

    tool = AseVaspTool()
    if args.traj:
        traj = tool.report(calc_dir=args.calcdir, index=":")
        safe_dump(
            [x.to_dict() for x in traj],
            path=pathlib.Path("structure_with_properties.traj.json"),
            force=True,
        )
    else:
        casm_structure = tool.report(calc_dir=args.calcdir)
        safe_dump(
            casm_structure.to_dict(),
            path=pathlib.Path("structure_with_properties.json"),
            force=True,
        )
    return 0


def _vasp_report_all(
    target: pathlib.Path,
    handler: typing.Any,
):
    """Walk a directory or archive file and report results in CASM structure format for
    all calculation subdirectories.

    Notes
    -----
    If there are existing results in `<target>.results.json`, these will be read and
    used to skip re-reporting calculation directories, unless updates are forced with
    `handler.update` is True.

    This method writes three files:

    - `<target>.complete.json`: list[pathlib.Path], A list of completed calculation
      directories.
    - `<target>.incomplete.json`: list[pathlib.Path], A list of incomplete calculation
      directories.
    - `<target>.results.json`: dict[str, dict], A dictionary where keys are relative
      paths to the calculation directories, and values are dictionaries containing:

      - "structure_with_properties": ``dict``, The final structure with properties.
      - "structure": ``dict``, The initial structure, if present.
      - "config": ``dict``, The initial configuration, if present.

    Parameters
    ----------
    target : pathlib.Path
        The directory or archive file to walk. This should contain subdirectories that
        may be calculation directories. Calculation directories should not contain other
        calculation directories.
    handler : typing.Any
        An instance of a handler class that implements the methods `is_calcdir`,
        `is_complete`, and `report`. This handler will be used to determine if a
        directory is a calculation directory, if it is complete, and will run the
        report creating the "structure_with_properties.json" data. See
        :class:`~casm.tools.calc.handlers.CasmV1VaspReportHandler` for an example
        implementation.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    from casm.tools.calc.methods import report, report_archive

    if not target.exists():
        print(f"{target}: does not exist")
        return 1

    try:
        if str(target).endswith((".tar.gz", ".tgz")):
            report_archive(archive_path=target, handler=handler)
            return 0
        else:
            report(dir=target, handler=handler)
            return 0
    except Exception as e:
        print(f"{target}: Error processing archive")
        print(e)
        return 1


def vasp_report_all_v1(args):
    """Implements ``casm-calc vasp report-all-v1 ...``

    Notes
    -----
    This method assumes the CASM v1 organization of VASP calculations. See
    :class:`~casm.tools.calc.handlers.CasmV1VaspReportHandler` for details.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line. Uses:

        - `args.target`: pathlib.Path, The directory or archive file to search for VASP
          calculations.
        - `args.calc_id`: str, optional, The calculation type ID (i.e. gga in
          'calctype.gga'), used to identify calculation directories.
        - `args.update`: If set, forces re-parsing VASP output.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    from casm.tools.calc.handlers import CasmV1VaspReportHandler

    target = args.target
    calc_id = args.calc_id
    update = args.update

    handler = CasmV1VaspReportHandler(calc_id=calc_id, update=update)
    return _vasp_report_all(
        target=target,
        handler=handler,
    )


def vasp_report_all(args):
    """Implements ``casm-calc vasp report-all ...``

    Notes
    -----
    This method assumes the VASP calculations have sentinel file or directory
    which can be used to identify them. See
    :class:`~casm.tools.calc.handlers.CasmVaspReportHandler` for details.

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line. Uses:

        - `args.target`: pathlib.Path, The directory or archive file to search for VASP
          calculations.
        - `args.update`: If set, forces re-parsing VASP output.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    from casm.tools.calc.handlers import CasmVaspReportHandler

    target = args.target
    sentinel_name = args.sentinel
    update = args.update

    handler = CasmVaspReportHandler(sentinel_name=sentinel_name, update=update)
    return _vasp_report_all(
        target=target,
        handler=handler,
    )


################################################################################
vasp_setup_desc = """
# The `casm-calc vasp setup` command:

Setup a VASP calculation directory.
"""

vasp_report_desc = """
# The `casm-calc vasp report` command:

Read VASP calculation results and write them to a file in the current directory.
 
By default, this writes a file named `structure_with_properties.json` containing 
a CASM structure for the final state of the calculation results. 

If the option --traj is used, this writes a file named 
`structure_with_properties.traj.json` containing a list of CASM Structures, one 
for each step in the calculation.
"""

vasp_report_all_v1_desc = """
# The `casm-calc vasp report-all-v1` command:

This works for calculations with the following directory structure:

    .. code-block:: text

        <target>/
        ├── <configuration_name>/
        │   ├── config.json
        │   ├── structure.json
        │   ├── calctype.<calc_id>/
        │   │   ├── ...
        │   │   ├── run.final/OUTCAR
        │   │   ├── run.final/OUTCAR.gz
        │   │   └── structure_with_properties.json
        ...

When the command is run on `target`, it assumes that:

1. All subdirectories named "calctype.<calc_id>" are calculation directories.
2. A calculation directory with a "run.final" directory containing either
   an "OUTCAR" or "OUTCAR.gz" file is complete. Otherwise, the calculation is
   incomplete.
3. The "config.json" and "structure.json" files are optional, but may be present
   if the calculation was setup from a CASM configuration or structure.

For each completed calculation directory, it will parse the results from the 
"OUTCAR" or "OUTCAR.gz" file, and store the resulting 
"structure_with_properties.json" file in the calculation directory (this is not
done if the target is an archive). It will also collect the "config.json" and 
"structure.json" file contents if they are present.

This method writes three files:

- ``<target>.<calc_id>.complete.json``: list[pathlib.Path], A list of completed 
  calculation directories.
- ``<target>.<calc_id>.incomplete.json``: list[pathlib.Path], A list of incomplete 
  calculation directories.
- ``<target>.<calc_id>.results.json`` : dict[str, dict], A dictionary where keys 
  are relative paths to the calculation directories, and values are dictionaries 
  containing:

  - "structure_with_properties": ``dict``, The final structure with properties.
  - "structure": ``dict``, The initial structure, if present.
  - "config": ``dict``, The initial configuration, if present.

"""

vasp_report_all_desc = """
# The `casm-calc vasp report-all` command:

This works for calculations with the following directory structure:

.. code-block:: text

    <target>/
    ├── <configuration_name>/
    │   ├── config.json
    │   ├── structure.json
    │   ├── ...
    │   ├── run.final/OUTCAR
    │   ├── run.final/OUTCAR.gz
    │   └── structure_with_properties.json
    ...

When the commaned is run on `target`, it assumes that:

1. All subdirectories containing a file or directory named `sentinel_name`
   is a calculation directory. By default, the sentinel name is "run.final".
2. A calculation directory with a "run.final" directory containing either
   an "OUTCAR" or "OUTCAR.gz" file is complete. Otherwise, the calculation is
   incomplete.

For each completed calculation directory,  it will parse the results from the 
"OUTCAR" or "OUTCAR.gz" file, and store the resulting 
"structure_with_properties.json" file in the calculation directory (this is 
not done if the target is an archive). It will also collect the "config.json" 
and "structure.json" file contents if they are present.

This method writes three files:

- ``<target>.complete.json``: list[pathlib.Path], A list of completed 
  calculation directories.
- ``<target>.incomplete.json``: list[pathlib.Path], A list of incomplete 
  calculation directories.
- ``<target>.results.json`` : dict[str, dict], A dictionary where keys are 
  relative paths to the calculation directories, and values are dictionaries 
  containing:

  - "structure_with_properties": ``dict``, The final structure with properties.
  - "structure": ``dict``, The initial structure, if present.
  - "config": ``dict``, The initial configuration, if present.

"""


def print_desc(argv=None):
    if "setup" in argv:
        print(vasp_setup_desc)
    elif "report" in argv:
        print(vasp_report_desc)
    elif "report-all-v1" in argv:
        print(vasp_report_all_v1_desc)
    elif "report-all" in argv:
        print(vasp_report_all_desc)
    else:
        print("No extended description available.")


def make_vasp_subparser(c):
    """Constructs the ``casm-calc vasp ...`` argument parser, and attaches the methods
    for running the subcommands.

    Parameters
    ----------
    c: argparse._SubParsersAction
        The output from ``parser.add_subparsers`` to which ``casm-calc vasp`` arguments
        are added.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    vasp = c.add_parser(
        "vasp",
        help="VASP calculations",
        description="Setup, run, and collect VASP calculations",
    )

    ### casm-calc vasp setup ....
    vasp_m = vasp.add_subparsers(title="Select which method to use")
    m = vasp_m.add_parser(
        "setup",
        help="Setup calculation input",
        description="Setup input files for a VASP calculation.",
    )
    m.set_defaults(func=vasp_setup)

    m.add_argument("structure", type=pathlib.Path, help="Structure file")
    m.add_argument(
        "calcdir",
        type=pathlib.Path,
        help=("Directory to write the calculation input files to. "),
    )
    m.add_argument(
        "settings",
        type=pathlib.Path,
        help=(
            "Settings directory. Expected to contain a 'calc.json' file with settings "
            "to be passed to the ase.calculators.vasp.Vasp constructor, and template "
            "INCAR and KPOINTS files."
        ),
    )
    m.add_argument(
        "--format",
        type=str,
        default=None,
        help=(
            "Set the format for reading the structure file. One of 'vasp', 'casm', "
            "or a format recognized by `ase.io.read`. "
            "The default is that no suffix or a '.vasp' suffix will read the file as a "
            "VASP POSCAR file, and a '.json' or '.casm' suffix will read the file as a "
            "CASM Structure JSON file. Otherwise, the file is read using ASE's "
            "`ase.io.read` method, if ASE is installed."
        ),
    )

    ### casm-calc vasp report ....
    m = vasp_m.add_parser(
        "report",
        help="Report VASP calculation output",
        description=(
            "Report VASP calculation results. Results are read from an OUTCAR or "
            "OUTCAR.gz file present in the calculation directory."
        ),
    )
    m.set_defaults(func=vasp_report)

    m.add_argument(
        "calcdir",
        type=pathlib.Path,
        help="Calculation directory",
    )
    m.add_argument(
        "--traj",
        action="store_true",
        help=("Save the trajectory of structures instead of the final structure. "),
    )

    ### casm-calc vasp report-all ....
    m = vasp_m.add_parser(
        "report-all",
        help="Report calculations from a directory or archive",
        description=(
            "Report all VASP calculation results from a directory or archive file "
            "(*.tar.gz or *.tgz), identified by a sentinel file or directory (use "
            "--desc for details)."
        ),
    )
    m.set_defaults(func=vasp_report_all)

    m.add_argument(
        "target",
        type=pathlib.Path,
        help=(
            "Directory or archive file (*.tar.gz or *.tgz) to search for "
            "VASP calculations"
        ),
    )
    m.add_argument(
        "--sentinel",
        type=str,
        default="run.final",
        help=(
            "Name of file or directory used to identify calculation directories "
            "(default is 'run.final')"
        ),
    )
    m.add_argument(
        "--update",
        action="store_true",
        help=("Force re-parsing VASP output. "),
    )
    other = m.add_argument_group("Other options")
    other.add_argument(
        "--desc",
        action="store_true",
        help="Print an extended description of the method and parameters.",
    )

    ### casm-calc vasp report-all-v1 ....
    m = vasp_m.add_parser(
        "report-all-v1",
        help="Report calculations from a directory or archive (CASM v1 layout)",
        description=(
            "Report all VASP calculation results from a directory or archive file "
            "(*.tar.gz or *.tgz) with CASM v1 layout (use --desc for details)."
        ),
    )
    m.set_defaults(func=vasp_report_all_v1)

    m.add_argument(
        "target",
        type=pathlib.Path,
        help=(
            "Directory or archive file (*.tar.gz or *.tgz) to search for "
            "VASP calculations"
        ),
    )
    m.add_argument(
        "--calc-id",
        type=str,
        help="Calculation type ID (i.e. gga in 'calctype.gga'), used to identify "
        "calculation directories",
    )
    m.add_argument(
        "--update",
        action="store_true",
        help=("Force re-parsing VASP output. "),
    )
    other = m.add_argument_group("Other options")
    other.add_argument(
        "--desc",
        action="store_true",
        help="Print an extended description of the method and parameters.",
    )
