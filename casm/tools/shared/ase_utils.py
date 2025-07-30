"""Interface with `ASE <https://wiki.fysik.dtu.dk/ase/>`_"""

import pathlib
import typing

import ase
import ase.calculators.vasp
import ase.io
import numpy as np

import casm.tools.shared.json_io as json_io
import libcasm.configuration as casmconfig
import libcasm.xtal as xtal


def make_ase_atoms(casm_structure: xtal.Structure) -> ase.Atoms:
    """Given a CASM Structure, convert it to an ASE Atoms

    .. attention::

        This method only works for non-magnetic atomic structures. If the structure
        contains molecular information, an error will be raised.

    Notes
    -----

    This method converts a CASM Structure object to an ASE Atoms object. It includes:

    - the lattice vectors
    - the atomic positions
    - the atomic types

    Parameters
    ----------
    casm_structure : libcasm.xtal.Structure

    Returns
    -------
    ase.Atoms

    """
    if len(casm_structure.mol_type()):
        raise ValueError(
            "Error: only non-magnetic atomic structures may be converted using "
            "to_ase_atoms"
        )

    symbols = casm_structure.atom_type()
    positions = casm_structure.atom_coordinate_cart().transpose()
    cell = casm_structure.lattice().column_vector_matrix().transpose()

    return ase.Atoms(
        symbols=symbols,
        positions=positions,
        cell=cell,
        pbc=True,
    )


def make_casm_structure(ase_atoms: ase.Atoms) -> xtal.Structure:
    """Given an ASE Atoms, convert it to a CASM Structure

    .. attention::

        This method only works for non-magnetic atomic structures.

    Notes
    -----

    This method converts an ASE Atoms object to a CASM Structure object. It includes:

    - the lattice vectors
    - atomic positions
    - atomic types, using the ASE chemical symbols.
    - Optional properties, if available from an ASE calculator:

      - Forces are added as atom properties named `"force"`.
      - Potential energy is added as the global property name `"energy"`.

    Parameters
    ----------
    ase_atoms : ase.Atoms
        A :class:`ase.Atoms` object

    Returns
    -------
    casm_structure: libcasm.xtal.Structure
        A :class:`~libcasm.xtal.Structure` object

    """

    lattice = xtal.Lattice(
        column_vector_matrix=ase_atoms.get_cell().transpose(),
    )
    atom_coordinate_frac = ase_atoms.get_scaled_positions().transpose()
    atom_type = ase_atoms.get_chemical_symbols()

    atom_properties = {}
    global_properties = {}
    if ase_atoms._calc is not None:
        try:
            forces = ase_atoms.get_forces()
            atom_properties["force"] = forces.transpose()
        except Exception:
            pass

        try:
            energy = ase_atoms.get_potential_energy()
            global_properties["energy"] = np.array([[energy]])
        except Exception:
            pass

    return xtal.Structure(
        lattice=lattice,
        atom_coordinate_frac=atom_coordinate_frac,
        atom_type=atom_type,
        atom_properties=atom_properties,
        global_properties=global_properties,
    )


def write_structure_using_ase(
    casm_structure: xtal.Structure,
    path: pathlib.Path,
    format: typing.Optional[str] = None,
    make_ase_atoms_f: typing.Optional[
        typing.Callable[[xtal.Structure], ase.Atoms]
    ] = None,
) -> None:
    """Write a structure using ASE's write function.

    .. attention::

        This method does not write magnetic moments.

    Parameters
    ----------
    casm_structure : libcasm.xtal.Structure
        The CASM Structure to write.
    path : pathlib.Path
        The path to the file where the structure will be written.
    format : Optional[str]=None
        The format to use for writing the file. If None, ASE will try to infer the
        format from the file extension.
    make_ase_atoms_f : Optional[Callable[[libcasm.xtal.Structure], ase.Atoms]] = None
        A function to convert the CASM structure to an ASE Atoms object. If None,
        the default function, :func:`make_ase_atoms` is used, which works for
        non-magnetic atomic structures.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. "
            "Please install ASE to write this structure format."
        )
    if make_ase_atoms_f is None:
        make_ase_atoms_f = make_ase_atoms

    ase_atoms = make_ase_atoms_f(casm_structure)
    ase.io.write(path.as_posix(), ase_atoms, format=format)


def write_structure_traj_using_ase(
    casm_structure_traj: list[xtal.Structure],
    path: pathlib.Path,
    format: typing.Optional[str] = None,
    make_ase_atoms_f: typing.Optional[
        typing.Callable[[xtal.Structure], ase.Atoms]
    ] = None,
) -> None:
    """Write a structure using ASE's write function.

    .. attention::

        This method does not write magnetic moments.

    Parameters
    ----------
    casm_structure_traj : list[libcasm.xtal.Structure]
        The CASM Structure trajectory to write.
    path : pathlib.Path
        The path to the file where the structure will be written.
    format : Optional[str]=None
        The format to use for writing the file. If None, ASE will try to infer the
        format from the file extension.
    make_ase_atoms_f : Optional[Callable[[libcasm.xtal.Structure], ase.Atoms]] = None
        A function to convert the CASM structure to an ASE Atoms object. If None,
        the default function, :func:`make_ase_atoms` is used, which works for
        non-magnetic atomic structures.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. "
            "Please install ASE to write this structure format."
        )
    if make_ase_atoms_f is None:
        make_ase_atoms_f = make_ase_atoms

    ase_atoms = [make_ase_atoms_f(x) for x in casm_structure_traj]
    ase.io.write(path.as_posix(), ase_atoms, format=format)


def read_structure_using_ase(
    path: pathlib.Path,
    format: typing.Optional[str] = None,
    make_casm_structure_f: typing.Optional[
        typing.Callable[[ase.Atoms], xtal.Structure]
    ] = None,
) -> xtal.Structure:
    """Read a structure using ASE's read function.

    .. attention::

        This method does not read magnetic moments.

    Parameters
    ----------
    path : pathlib.Path
        The path to the structure file.
    format : Optional[str]=None
        The format to use for reading the file. If None, ASE will try to infer the
        format from the file extension.
    make_casm_structure_f : typing.Callable[[ase.Atoms], libcasm.xtal.Structure]
        A function to convert an ASE Atoms object to a CASM structure. If None, the
        default function, :func:`make_casm_structure`, is used, which works for
        non-magnetic atomic structures.

    Returns
    -------
    casm_structure: libcasm.xtal.Structure
        A CASM Structure read from the file.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. Please install ASE to read this structure format."
        )
    if make_casm_structure_f is None:
        make_casm_structure_f = make_casm_structure

    return make_casm_structure_f(ase.io.read(path.as_posix(), format=format))


def read_structure_traj_using_ase(
    path: pathlib.Path,
    format: typing.Optional[str] = None,
    make_casm_structure_f: typing.Optional[
        typing.Callable[[ase.Atoms], xtal.Structure]
    ] = None,
) -> list[xtal.Structure]:
    """Read a structure trajectory using ASE's read function.

    .. attention::

        This method does not read magnetic moments.

    Parameters
    ----------
    path : pathlib.Path
        The path to the structure file.
    format : Optional[str]=None
        The format to use for reading the file. If None, ASE will try to infer the
        format from the file extension.
    make_casm_structure_f : typing.Callable[[ase.Atoms], libcasm.xtal.Structure]
        A function to convert an ASE Atoms object to a CASM structure. If None, the
        default function, :func:`make_casm_structure`, is used, which works for
        non-magnetic atomic structures.

    Returns
    -------
    casm_structure: libcasm.xtal.Structure
        A CASM Structure read from the file.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. Please install ASE to read this structure format."
        )
    if make_casm_structure_f is None:
        make_casm_structure_f = make_casm_structure

    return [
        make_casm_structure_f(x)
        for x in ase.io.read(path.as_posix(), format=format, index=":")
    ]


class AseVaspTool:
    def __init__(
        self,
        calctype_settings_dir: typing.Optional[pathlib.Path] = None,
        make_ase_atoms_f: typing.Optional[
            typing.Callable[[xtal.Structure], ase.Atoms]
        ] = None,
        make_casm_structure_f: typing.Optional[
            typing.Callable[[ase.Atoms], xtal.Structure]
        ] = None,
    ):
        """Setup, run, and collect VASP calculations using ASE.

        For details on the parameters and environment configuration, see the
        `ASE Vasp calculator documentation <https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html>`_.

        .. attention::

            ASE assumes that POTCAR files exist in one of `potpaw_PBE`, `potpaw`, or
            `potpaw_GGA`, located at the path specified by the environment
            variable VASP_PP_PATH.

        Parameters
        ----------
        calctype_settings_dir: Optional[pathlib.Path] = None
            Path to the directory containing settings, including a `calc`.json file for
            VASP calculator constructor arguments, and INCAR and KPOINTS
            template files. All files are optional, but if they exist, they will be
            used to set up the VASP calculator.
        make_ase_atoms_f: Optional[Callable[[libcasm.xtal.Structure], ase.Atoms]] = None
            A function to convert the CASM structure to an ASE Atoms object. The
            default function, :func:`make_ase_atoms` works for non-magnetic atomic
            structures.
        make_casm_structure_f: \
        Optional[Callable[[ase.Atoms], libcasm.xtal.Structure]] = None
            A function to convert an ASE Atoms object to a CASM structure. The
            default function, :func:`make_casm_structure` works for non-magnetic atomic
            structures.
        """
        if make_ase_atoms_f is None:
            make_ase_atoms_f = make_ase_atoms
        if make_casm_structure_f is None:
            make_casm_structure_f = make_casm_structure

        ### Read settings from the calctype_settings_dir if provided

        self.calctype_settings_dir = calctype_settings_dir
        """Optional[pathlib.Path]: Path to the directory containing settings, 
        including a `calc.json` file for 
        `ASE VASP calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#module-ase.calculators.vasp>`_
        constructor arguments, and template INCAR and KPOINTS files. All files are 
        optional, but if they exist, they will be used to set up the VASP 
        calculations."""

        incar_path = None
        kpoints_path = None
        settings = {}
        if self.calctype_settings_dir is not None:
            settings = json_io.read_required(
                path=self.calctype_settings_dir / "calc.json"
            )

            incar_path = calctype_settings_dir / "INCAR"
            if not incar_path.exists():
                incar_path = None

            kpoints_path = calctype_settings_dir / "KPOINTS"
            if not kpoints_path.exists():
                kpoints_path = None

        self.settings = settings
        """Optional[dict]: Settings read from the calc.json file in the 
        calctype_settings_dir, which will be passed to the 
        :class:`ase.calculators.vasp.Vasp` calculator constructor."""

        self.incar_path = incar_path
        """Optional[pathlib.Path]: Path to the template INCAR file, if it exists."""

        self.kpoints_path = kpoints_path
        """Optional[pathlib.Path]: Path to the template KPOINTS file, if it exists."""

        ### Functions to convert between CASM Structure and ASE Atoms

        self._make_ase_atoms_f = make_ase_atoms_f
        """Callable[[libcasm.xtal.Structure], ase.Atoms]: Function to convert CASM 
        Structure to ASE Atoms."""

        self._make_casm_structure_f = make_casm_structure_f
        """Callable[[ase.Atoms], libcasm.xtal.Structure]: Function to convert an 
        ASE Atoms to CASM Structure."""

    def make_calculator(
        self,
        ase_atoms: ase.Atoms,
        calc_dir: pathlib.Path,
    ) -> ase.calculators.vasp.Vasp:
        """Construct an ASE VASP calculator.

        Parameters
        ----------
        ase_atoms: ase.Atoms
            The ASE Atoms object to use for the calculation.
        calc_dir: pathlib.Path
            The directory in which to store the calculation files. This directory will
            be created if it does not exist.

        Returns
        -------
        vasp_calculator: ase.calculators.vasp.Vasp
            The ASE VASP calculator object, configured with the provided settings and
            paths to INCAR and KPOINTS files if they exist.

        """
        calc_dir.mkdir(parents=True, exist_ok=True)

        vasp_calculator = ase.calculators.vasp.Vasp(
            atoms=ase_atoms,
            directory=calc_dir,
            **self.settings,
        )

        if self.incar_path is not None:
            vasp_calculator.read_incar(self.incar_path)

        if self.kpoints_path is not None:
            vasp_calculator.read_kpoints(self.kpoints_path)

        return vasp_calculator

    def setup(
        self,
        casm_structure: xtal.Structure,
        calc_dir: pathlib.Path,
        config: typing.Optional[casmconfig.Configuration] = None,
    ) -> ase.calculators.vasp.Vasp:
        """Setup a VASP calculation for a given structure.

        Parameters
        ----------
        casm_structure: libcasm.xtal.Structure
            The structure to calculate. The structure is written to the calculation
            directory as `structure.json`.
        calc_dir: pathlib.Path
            The directory in which to store the calculation files.
        config: Optional[libcasm.configuration.Configuration] = None
            If provided, the configuration object associated with the structure is
            printed to the calculation directory as `config.json`.

        Returns
        -------
        vasp_calculator: ase.calculators.vasp.Vasp
            The ASE VASP calculator object.
        """
        ase_atoms = self._make_ase_atoms_f(casm_structure)
        vasp_calculator = self.make_calculator(ase_atoms=ase_atoms, calc_dir=calc_dir)

        # Write INCAR, KPOINTS, POTCAR, POSCAR
        vasp_calculator.write_input(atoms=ase_atoms)

        return vasp_calculator

    def report(
        self,
        calc_dir: typing.Optional[pathlib.Path] = None,
        index: typing.Any = None,
    ) -> typing.Union[xtal.Structure, list[xtal.Structure]]:
        """Report the results of a VASP calculation.

        Parameters
        ----------
        calc_dir: pathlib.Path
            The directory containing the VASP calculation files.
        index: int, slice or str
            Indicates the structures to return. By default, only the last structure
            is returned. Use `":"` to return all structures.

        Returns
        -------
        results: Union[libcasm.xtal.Structure, list[libcasm.xtal.Structure]]
            A CASM Structure or a list of CASM Structures, as specified by `index`.
        """

        outcar_file = calc_dir / "OUTCAR"
        if outcar_file.exists():
            value = ase.io.read(outcar_file, format="vasp-out", index=index)
        else:
            outcar_gz_file = calc_dir / "OUTCAR.gz"
            if outcar_gz_file.exists():
                import gzip

                with gzip.open(outcar_gz_file, "rt") as f:
                    value = ase.io.read(f, format="vasp-out", index=index)
            else:
                raise FileNotFoundError(
                    f"Error in AseVaspTool.report: "
                    f"Neither OUTCAR nor OUTCAR.gz found in {calc_dir.as_posix()}"
                )

        if isinstance(value, ase.Atoms):
            results = self._make_casm_structure_f(value)
        elif isinstance(value, list):
            results = [self._make_casm_structure_f(x) for x in value]
        else:
            raise Exception(f"Unrecognized type {type(value)} from ase.io.read")

        return results
