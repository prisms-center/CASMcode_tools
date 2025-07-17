"""Implements ``casm-map write ...``"""

import pathlib


# <-- max width = 80 characters                            --> #
################################################################
def print_desc():
    desc = """
# The `casm-map write` command:

TODO
"""
    print(desc)


def run_write(args):
    """Implements ``casm-map write ...``

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line.

    """
    print("Running write command with arguments:", args)
    # Implement the logic for the write command here
    return 0


def make_write_parser(m):
    """Adds the ``casm-map write ...`` argument sub-parser.

    Parameters
    ----------
    m : argparse._SubParsersAction
        The output from ``parser.add_subparsers`` to which ``casm-map write`` arguments
        are added.

    """
    write = m.add_parser(
        "write",
        help=(
            "Write selected structure mappings, mapped structures, and interpolated "
            "paths."
        ),
    )
    write.set_defaults(func=run_write)

    ### Select mapping to write
    select = write.add_argument_group("Mapping selection options")
    select.add_argument(
        "-b",
        "--best",
        action="store_true",
        help=("Select the lowest-cost mapping. "),
    )
    select.add_argument(
        "--uuid",
        type=pathlib.Path,
        nargs="+",
        help="Select one or more mappings by UUID",
    )
    select.add_argument(
        "--index",
        type=pathlib.Path,
        nargs="+",
        help="Select one or more mappings by index",
    )

    ### Options
    output = write.add_argument_group("Output options")
    output.add_argument(
        "--orbit",
        action="store_true",
        help="Write orbit of equivalent mapped structures.",
    )
    output.add_argument(
        "--path",
        metavar="IMAGES",
        type=int,
        help=(
            "Write interpolated path between parent and mapped child structures "
            "with specified number of images."
        ),
    )
    output.add_argument(
        "--path-orbit",
        metavar="ORBIT_IMAGES",
        type=int,
        help="Write orbit of equivalent interpolated paths. ",
    )

    ### Location options:
    location = write.add_argument_group("Location")
    location.add_argument(
        "--results-dir",
        type=pathlib.Path,
        default=pathlib.Path("results"),
        help="Directory containing mapping search results (default=results).",
    )
    location.add_argument(
        "--prefix",
        type=pathlib.Path,
        help=(
            "Specify a custom location where output files are written "
            "(default= current working directory )."
        ),
    )
    location.add_argument(
        "-s",
        "--save",
        action="store_true",
        help=(
            "Write the mapping and output files in a subdirectory of the results "
            "directory, named by the UUID. (i.e. --save is equivalent to "
            "--prefix=<RESULTS_DIR>/<UUID>). "
        ),
    )

    ### Output format options:
    format = write.add_argument_group("Format options")
    format.add_argument(
        "--one-file-per-orbit",
        action="store_true",
        help=(
            "Write orbit of equivalent mapped structures in one file "
            "(default= write individual files)."
        ),
    )
    format.add_argument(
        "--one-file-per-path",
        action="store_true",
        help="Write path in one file (default= write individual files ). ",
    )
    format.add_argument(
        "--casm",
        action="store_true",
        help=("Write mapped structures as CASM structure JSON files. "),
    )
    format.add_argument(
        "--config",
        action="store_true",
        help=("Write mapped structures as CASM configuration JSON files. "),
    )
    format.add_argument(
        "--vasp",
        action="store_true",
        help=("Write mapped structures as VASP POSCAR files."),
    )
    format.add_argument(
        "--vasp-neb",
        action="store_true",
        help=(
            "Write path images as VASP POSCAR in sequentially "
            "numbered directories for NEB calculation input "
            "(default= files in same directory). "
        ),
    )
    format.add_argument(
        "--format",
        type=str,
        default="casm",
        help=("Output structure file format (default=vasp)."),
    )

    ### Other options:
    other = write.add_argument_group("Other options")
    other.add_argument(
        "--desc",
        action="store_true",
        help="Print an extended description of the method and parameters.",
    )
