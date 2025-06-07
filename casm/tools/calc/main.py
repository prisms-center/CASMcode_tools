import argparse
import os
import pathlib
import sys


def get_format(args):
    if args.format is not None:
        return args.format
    return None


def vasp_setup(args):
    from casm.tools.shared.ase_utils import AseVaspTool
    from casm.tools.shared.io import read_structure

    tool = AseVaspTool(calctype_settings_dir=args.settings)
    tool.setup(
        casm_structure=read_structure(path=args.structure, format=get_format(args)),
        calc_dir=args.calcdir,
    )
    return 0


def vasp_report(args):
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


def make_parser():
    """Make a parser for the casm-map command line interface."""

    ### casm-calc ...
    parser = argparse.ArgumentParser(
        description="CASM calculation CLI tool",
    )
    c = parser.add_subparsers(title="Select which calculator to use")

    ### casm-calc vasp ...
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
        description="Setup input files for a calculation.",
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
            "for the calculation, and template INCAR and KPOINTS files."
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

    return parser


def main(argv=None, working_dir=None):

    from casm.tools.shared import contexts

    if argv is None:
        argv = sys.argv
    if working_dir is None:
        working_dir = os.getcwd()

    parser = make_parser()
    if len(argv) < 2:
        parser.print_help()
        return 1
    args = parser.parse_args(argv[1:])

    code = 0
    with contexts.working_dir(working_dir):
        code = args.func(args)

    return code
