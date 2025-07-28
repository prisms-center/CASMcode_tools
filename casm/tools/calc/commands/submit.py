"""Implements ``casm-calc submit ...``"""


def run_submit(args):
    """Implements ``casm-calc submit ...``

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line. Uses:

        - `args.enum`: str, Enumeration ID
        - `args.clex`: str, Cluster expansion ID to specify the calctype
        - `args.selection`: str, optional, Name of a selection in the enumeration with
          selected configurations to submit calculations for. Select all configurations
          if not specified.
        - `args.dry_run`: If set, prints configurations that would be submitted and
          their current calc status, but does not submit.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    from .status import (
        _validate_casm_project,
        list_options,
        make_config_selection,
    )

    _validate_casm_project()

    from casm.project import Project, project_path

    # --- Load project ---
    path = project_path()
    if path is None:
        print("No CASM project found.")
        return 1
    project = Project(path=path)

    # --- list options ---
    if args.list:
        list_options(
            project=project,
            args=args,
        )
        return 0
    elif args.enum is None:
        print("No enumeration ID provided. Use --list to see available options.")
        return 1

    # --- Variables ---
    if args.dry_run:
        submit_jobs = False
    else:
        submit_jobs = True

    # --- Select configurations ---
    config_selection = make_config_selection(
        project=project,
        args=args,
    )
    if isinstance(config_selection, int):
        # An error occurred in make_config_selection
        return config_selection

    # --- Submit jobs or print status ---
    for record in config_selection:
        if not record.is_selected:
            continue

        if not submit_jobs:
            # --- Print current status ---
            print(f"(dry-run) {record.name}: {record.calc_status}")
        else:
            # --- Submit job ---
            completed_processes = record.run_subprocess(
                args=["sbatch", "submit.sh"],
                capture_output=True,
                text=True,
            )

            if completed_processes[0].returncode != 0:
                print(
                    f"Error submitting job for {record.name}: "
                    f"{completed_processes[0].stderr}"
                )
                return 1

            stdout = completed_processes[0].stdout

            # --- Extract job ID from output ---
            # Use string operations only
            if "Submitted batch job" in stdout:
                jobid = stdout.strip().split()[-1]
                print("Submitted job ID:", jobid)
            else:
                print(
                    "Warning: Could not parse jobid\n"
                    f"Unexpected output format: \n"
                    f"{stdout}"
                )
                jobid = None

            # match = re.search(r"Submitted batch job (\d+)", result.stdout)

    return 0


################################################################################


def print_desc(argv=None):
    print("No extended description available.")


def make_submit_subparser(c):
    """Constructs the ``casm-calc submit ...`` argument parser, and attaches the methods
    for running the subcommands.

    Parameters
    ----------
    c: argparse._SubParsersAction
        The output from ``parser.add_subparsers`` to which ``casm-calc submit``
        arguments are added.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    submit = c.add_parser(
        "submit",
        help="Submit CASM project calculations",
        description="Submit CASM project calculations",
    )

    ### casm-calc submit ....
    submit.set_defaults(func=run_submit)

    submit.add_argument(
        "enum",
        type=str,
        help=("Enumeration ID"),
        nargs="?",
    )
    submit.add_argument(
        "-c",
        "--selection",
        type=str,
        help=(
            "Name of a selection in the enumeration with selected configurations to "
            "submit calculations for. Select all configurations if not specified. "
        ),
    )
    submit.add_argument(
        "--calctype",
        type=str,
        help=(
            "A calctype ID, used to indicate which calctype to submit"
            "calculations for. Uses --clex if not specified. "
        ),
    )
    submit.add_argument(
        "--clex",
        type=str,
        help=(
            "A cluster expansion key, used to indicate which calctype to submit "
            "calculations for. Uses the default clex if not specified. "
        ),
    )
    submit.add_argument(
        "--dry-run",
        action="store_true",
        help=(
            "If given, print configurations that would be submitted and their current "
            "calc status, but do not submit."
        ),
    )
    submit.add_argument(
        "-l",
        "--list",
        action="store_true",
        help=("If given, list enumeration, calctype, clex, and selection options."),
    )
