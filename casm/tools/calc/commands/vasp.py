"""Implements ``casm-calc vasp ...``"""

import importlib.util
import os
import pathlib
import sys


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


################################################################################
vasp_report_desc = """
Read VASP calculation results and write them to a file in the current directory. 
By default, this writes a file named `structure_with_properties.json` containing 
a CASM structure for the final state of the calculation results. If the option 
--traj is used, this writes a file named `structure_with_properties.traj.json` 
containing a list of CASM Structures, one for each step in the calculation.
"""


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
        description=vasp_report_desc,
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
