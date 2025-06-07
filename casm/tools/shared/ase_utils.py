import os
import pathlib
from typing import Any, Callable, Optional, Union

import ase
import ase.calculators.vasp
import numpy as np

import libcasm.xtal as xtal
from casm.tools.shared.json_io import read_required


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
    casm_structure : xtal.Structure

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
    format: Optional[str] = None,
    make_ase_atoms_f: Callable[[xtal.Structure], ase.Atoms] = make_ase_atoms,
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
    make_ase_atoms_f : Callable[[xtal.Structure], ase.Atoms]
        A function to convert the CASM structure to an ASE Atoms object. The
        default function, :func:`make_ase_atoms` works for non-magnetic atomic
        structures.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. "
            "Please install ASE to try to write this structure format."
        )

    ase_atoms = make_ase_atoms_f(casm_structure)
    ase.io.write(path.as_posix(), ase_atoms, format=format)


def read_structure_using_ase(
    path: pathlib.Path,
    format: Optional[str] = None,
    make_casm_structure_f: Callable[[ase.Atoms], xtal.Structure] = make_casm_structure,
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
    make_casm_structure_f : Callable[[ase.Atoms], xtal.Structure]
        A function to convert an ASE Atoms object to a CASM structure. The
        default function, :func:`make_casm_structure` works for non-magnetic atomic
        structures.

    Returns
    -------
    casm_structure: libcasm.xtal.Structure
        A CASM Structure read from the file.

    """
    try:
        import ase.io
    except ImportError:
        raise ImportError(
            "ASE is not installed. "
            "Please install ASE to try to read this structure format."
        )

    return make_casm_structure_f(ase.io.read(path.as_posix(), format=format))


class AseVaspTool:
    def __init__(
        self,
        calctype_settings_dir: Optional[pathlib.Path] = None,
        make_ase_atoms_f: Callable[[xtal.Structure], ase.Atoms] = make_ase_atoms,
        make_casm_structure_f: Callable[
            [ase.Atoms], xtal.Structure
        ] = make_casm_structure,
    ):
        """Setup, run, and collect VASP calculations using ASE.

        For details on the parameters, see the `ase documentation for the vasp
        calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html>`.

        .. attention::

            ase assumes that POTCAR files exist in one of `potpaw_PBE`, `potpaw`, or
            `potpaw_GGA`, located at the path specified by the environment
            variable VASP_PP_PATH.

        Parameters
        ----------
        calctype_settings_dir: Optional[pathlib.Path] = None
            Path to the directory containing settings, including a `calc`.json file for
            VASP calculator constructor arguments, and INCAR and KPOINTS
            template files. All files are optional, but if they exist, they will be
            used to set up the VASP calculator.
        setups: Optional[dict] = None
            Dictionary of pseudopotential setups for each element.
        xc: Optional[str] = None
            Exchange-correlation type, one of: "pbe", "lda", or "pw91"
        make_ase_atoms_f: Callable[[xtal.Structure], ase.Atoms] = make_ase_atoms
            A function to convert the CASM structure to an ASE Atoms object. The
            default function, :func:`make_ase_atoms` works for non-magnetic atomic
            structures.
        make_casm_structure_f: \
        Callable[[ase.Atoms], xtal.Structure] = make_casm_structure
            A function to convert an ASE Atoms object to a CASM structure. The
            default function, :func:`make_casm_structure` works for non-magnetic atomic
            structures.
        """
        if "VASP_PP_PATH" not in os.environ:
            raise ValueError(
                "Please set the environment variable VASP_PP_PATH to the directory "
                "containing the VASP POTCAR files. It should contain the directories "
                "potpaw_PBE, potpaw, and potpaw_GGA."
            )

        ### Read settings from the calctype_settings_dir if provided

        self.calctype_settings_dir = calctype_settings_dir
        """Optional[pathlib.Path]: Path to the directory containing settings, 
        including a `calc`.json file for 
        `ASE VASP calculator <https://wiki.fysik.dtu.dk/ase/ase/calculators/vasp.html#module-ase.calculators.vasp>`_
        constructor arguments, and template INCAR and KPOINTS files. All files are 
        optional, but if they exist, they will be used to set up the VASP 
        calculations."""

        incar_path = None
        kpoints_path = None
        settings = {}
        if self.calctype_settings_dir is not None:
            settings = read_required(path=self.calctype_settings_dir / "calc.json")

            incar_path = calctype_settings_dir / "INCAR"
            if not incar_path.exists():
                incar_path = None

            kpoints_path = calctype_settings_dir / "KPOINTS"
            if not kpoints_path.exists():
                kpoints_path = None

        self.settings = settings
        """Optional[dict]: Settings read from the calc.json file in the 
        calctype_settings_dir, which will be passed to the Vasp calculator
        constructor."""

        self.incar_path = incar_path
        """Optional[pathlib.Path]: Path to the template INCAR file, if it exists."""

        self.kpoints_path = kpoints_path
        """Optional[pathlib.Path]: Path to the template KPOINTS file, if it exists."""

        ### Functions to convert between CASM Structure and ASE Atoms

        self._make_ase_atoms_f = make_ase_atoms_f
        """Callable[[xtal.Structure], ase.Atoms]: Function to convert CASM Structure to 
        ASE Atoms."""

        self._make_casm_structure_f = make_casm_structure_f
        """Callable[[ase.Atoms], xtal.Structure]: Function to convert ASE Atoms to
        CASM Structure."""

    def make_calculator(
        self,
        ase_atoms: ase.Atoms,
        calc_dir: pathlib.Path,
    ) -> ase.calculators.vasp.Vasp:
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
    ) -> ase.calculators.vasp.Vasp:
        """Setup a VASP calculation for a given structure.

        Parameters
        ----------
        casm_structure: libcasm.xtal.Structure
            The structure to calculate.
        calc_dir: pathlib.Path
            The directory in which to store the calculation files.
        make_ase_atoms_f: Callable[[xtal.Structure], ase.Atoms] = make_ase_atoms
            A function to convert the CASM structure to an ASE Atoms object. The
            default function, :func:`make_ase_atoms` works for non-magnetic atomic
            structures.

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
        calc_dir: pathlib.Path,
        index: Any = None,
    ) -> Union[xtal.Structure, list[xtal.Structure]]:
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
        results: Union[xtal.Structure, list[xtal.Structure]]
            A CASM Structure or a list of CASM Structures, as specified by `index`.
        """
        value = ase.io.read(calc_dir / "OUTCAR", format="vasp-out", index=index)

        if isinstance(value, ase.Atoms):
            results = self._make_casm_structure_f(value)
        elif isinstance(value, list):
            results = [self._make_casm_structure_f(x) for x in value]
        else:
            raise Exception(f"Unrecognized type {type(value)} from ase.io.read")

        return results
