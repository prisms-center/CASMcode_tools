"""Implements ``casm-calc status ...``"""

import importlib.util
import sys
import typing

if typing.TYPE_CHECKING:
    import argparse

    from casm.project import ConfigSelection, Project


def _validate_casm_project():
    if importlib.util.find_spec("casm.project") is None:
        print(
            """
casm.project is not installed. To use `casm-calc status` install casm-project with:

    pip install casm-project"""
        )
        sys.exit(1)


def _print_available(
    items: list[str],
    item_type_desc: str,
    default_item: typing.Optional[str] = None,
):
    """Prints available items of a given type."""
    print(f"Available {item_type_desc}:")
    for item in items:
        suffix = ""
        if default_item is not None and item == default_item:
            suffix = " (default)"
        print(f"- {item}{suffix}")
    print()


def make_config_selection(
    project: "Project",
    args: "argparse.Namespace",
) -> "ConfigSelection":
    """Constructs a configuration selection based on the provided arguments.

    Parameters
    ----------
    project: Project
        The CASM project to work with.
    args: argparse.Namespace
        The parsed arguments from the command line. Uses:
        - `args.enum`: str, Enumeration ID
        - `args.calctype`: str, A calctype ID to specify the calculation type.
        - `args.clex`: str, Cluster expansion ID to specify the calctype. This is the
          fallback, used if `args.calctype` is not specified. If this is also not
          specified, the default cluster expansion is used.
        - `args.selection`: str, optional, Name of a selection in the enumeration
            to select configurations from.


    Returns
    -------
    config_selection: Union[casm.project.enum.ConfigSelection, int]
        A ConfigSelection, or a return code indicating failure (non-zero).

    """

    # --- Variables ---
    enum_id = args.enum
    clex_id = args.clex
    calc_id = args.calctype
    selection_name = args.selection

    # --- Construct enum ---
    if enum_id not in project.enum.all():
        print(f"Enumeration ID '{enum_id}' not found in project.")
        _print_available(
            items=project.enum.all(),
            item_type_desc="enumerations",
        )
        return 1
    enum = project.enum.get(enum_id)

    # --- Construct config_selection ---
    if calc_id is not None:
        if calc_id not in project.calc.all():
            print(f"Calculation type '{calc_id}' not found in project.")
            _print_available(
                items=project.calc.all(),
                item_type_desc="calculation types",
                default_item=project.settings.default_clex.calctype,
            )
            return 1
        clex = project.settings.default_clex
        clex.calctype = calc_id

    elif clex_id is not None:
        if clex_id not in project.settings.cluster_expansions:
            print(f"Cluster expansion '{clex_id}' not found in project.")
            _print_available(
                items=project.settings.cluster_expansions.keys(),
                item_type_desc="cluster expansions",
                default_item=project.settings.default_clex_name,
            )
            return 1
        clex = project.settings.get_clex(name=clex_id)
    else:
        clex = project.settings.default_clex

    all_selections = enum.all_config_selections()
    if selection_name is None:
        i = 0
        selection_name = f"all-{i}"
        while selection_name in all_selections:
            i += 1
            selection_name = f"all-{i}"
    elif selection_name not in all_selections:
        print(f"Selection '{selection_name}' not found in enumeration '{enum_id}'.")
        _print_available(
            items=all_selections,
            item_type_desc="selections",
        )
        return 1

    # --- Select configurations ---
    config_selection = enum.config_selection(
        name=selection_name,
        clex=clex,
    )
    print(config_selection)
    print()
    return config_selection


def list_options(
    project: "Project",
    args: "argparse.Namespace",
):
    """List available options for enumeration, calctype, clex, and selection."""
    print("Available enumerations:")
    for enum_id in project.enum.all():
        print(f"- {enum_id}")

    default_clex_name = project.settings.default_clex_name
    default_clex = project.settings.default_clex

    print("\nAvailable cluster expansions (--clex):")
    for clex_name, clex_desc in project.settings.cluster_expansions.items():
        is_default = ""
        if clex_name == default_clex_name:
            is_default = " (default)"
        print(f"- {clex_name} (calctype={clex_desc.calctype}){is_default}")

    print("\nAvailable calculation types (--calctype):")
    for calc_id in project.calc.all():
        is_default = ""
        if default_clex is not None:
            if calc_id == default_clex.calctype:
                is_default = " (default)"
        print(f"- {calc_id}{is_default}")

    if args.enum is not None and args.enum in project.enum.all():
        enum = project.enum.get(args.enum)
        print(f"\nAvailable selections (--selection) for --enum={args.enum}:")
        for selection_name in enum.all_config_selections():
            print(f"- {selection_name}")


def run_status(args):
    """Implements ``casm-calc status ...``

    Parameters
    ----------
    args : argparse.Namespace
        The parsed arguments from the command line. Uses:

        - `args.enum`: str, Enumeration ID
        - `args.clex`: str, Cluster expansion ID to specify the calctype
        - `args.selection`: str, optional, Name of a selection in the enumeration
          to check the status of. Select all configurations if not specified.
        - `args.tabulate`: If set, prints a table of job status counts.
        - `args.none`, `args.setup`, `args.started`, `args.stopped`, `args.complete`,
          `args.other`: If set, prints configurations with the corresponding status. If
          multiple are set, prints configurations with any of the specified statuses.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    from tabulate import tabulate

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

    config_selection = make_config_selection(
        project=project,
        args=args,
    )
    if isinstance(config_selection, int):
        # An error occurred in make_config_selection
        return config_selection

    status_count = dict()
    for record in config_selection:
        if not record.is_selected:
            continue
    details = []
    standard_status = ["none", "setup", "started", "stopped", "complete"]

    def add_details(record):
        details.append(
            [record.name, record.calc_status, record.calc_jobid, record.calc_runtime]
        )

    for record in config_selection:
        if args.all is None and not record.is_selected:
            continue

        status = record.calc_status
        if status in status_count:
            status_count[status] += 1
        else:
            status_count[status] = 1

        if (
            args.details
            or (args.none and status == "none")
            or (args.setup and status == "setup")
            or (args.started and status == "started")
            or (args.stopped and status == "stopped")
            or (args.complete and status == "complete")
            or (args.other and status not in standard_status)
        ):
            add_details(record)

    if args.tabulate:
        table = []
        for status, count in status_count.items():
            table.append([status, count])
        print(tabulate(table, headers=["Status", "Count"]))
        print()

    if len(details):
        headers = ["Name", "Status", "Job ID", "Runtime"]
        print(tabulate(details, headers=headers))
        print()

    return 0


################################################################################


def print_desc(argv=None):
    print("No extended description available.")


def make_status_subparser(c):
    """Constructs the ``casm-calc status ...`` argument parser, and attaches the methods
    for running the subcommands.

    Parameters
    ----------
    c: argparse._SubParsersAction
        The output from ``parser.add_subparsers`` to which ``casm-calc status``
        arguments are added.

    Returns
    -------
    code: int
        A return code indicating success (0) or failure (non-zero).

    """
    status = c.add_parser(
        "status",
        help="Chck the status of CASM project calculations",
        description="Check the status of CASM project calculations",
    )

    ### casm-calc status ....
    status.set_defaults(func=run_status)

    status.add_argument(
        "enum",
        type=str,
        help=("Enumeration ID"),
        nargs="?",
    )
    status.add_argument(
        "-c",
        "--selection",
        type=str,
        help=(
            "Name of a selection in the enumeration with selected configurations to "
            "check the calculation status. Select all configurations if not specified."
        ),
    )
    status.add_argument(
        "-a",
        "--all",
        action="store_true",
        help=(
            "If given, check the status of all configurations in the enumeration, "
            "regardless of selection. If not given, only check selected configurations."
        ),
    )
    status.add_argument(
        "--calctype",
        type=str,
        help=(
            "A calctype ID, used to indicate which calctype to check "
            "the calculation status. Uses --clex if not specified. "
        ),
    )
    status.add_argument(
        "--clex",
        type=str,
        help=(
            "A cluster expansion key, used to indicate which calctype to check "
            "the calculation status. Uses the default clex if not specified. "
        ),
    )
    status.add_argument(
        "-t",
        "--tabulate",
        action="store_true",
        help=("Print table of job status counts."),
    )
    status.add_argument(
        "-d",
        "--details",
        action="store_true",
        help=("Print configurations with any status."),
    )
    status.add_argument(
        "--none",
        action="store_true",
        help=('Print configurations with status="none".'),
    )
    status.add_argument(
        "--setup",
        action="store_true",
        help=('Print configurations with status="setup".'),
    )
    status.add_argument(
        "--started",
        action="store_true",
        help=('Print configurations with status="started".'),
    )
    status.add_argument(
        "--stopped",
        action="store_true",
        help=('Print configurations with status="stopped".'),
    )
    status.add_argument(
        "--complete",
        action="store_true",
        help=('Print configurations with status="complete".'),
    )
    status.add_argument(
        "--other",
        action="store_true",
        help=("Print jobs with any other status."),
    )
    status.add_argument(
        "-l",
        "--list",
        action="store_true",
        help=("If given, list enumeration, calctype, clex, and selection options."),
    )
