<img alt="Shows the CASM logo" src="https://raw.githubusercontent.com/prisms-center/CASMcode_global/main/python/doc/_static/logo.svg" width="600" />

#### casm-map

The casm-map package provides pure Python tools for structure mapping and import. This includes:

- A library and command line tool that enables:
  - searching for structure mappings
  - creating mapped structures
  - creating interpolated structures
  - generating structures with equivalent relative orientations
  - importing structures into CASM projects.

This package makes extensive use of lower-level methods which are implemented in [`libcasm-xtal`](https://github.com/prisms-center/CASMcode_crystallography), [`libcasm-mapping`](https://github.com/prisms-center/CASMcode_mapping), and [`libcasm-configuration`](https://github.com/prisms-center/CASMcode_configuration). 

Methods for searching for low-cost lattice, atom, and structure mappings, taking into account symmetry are based on the approach described in the paper `Thomas, Natarajan, and Van der Ven, npj Computational Materials, 7 (2021), 164 <https://doi.org/10.1038/s41524-021-00627-0>`_.


#### Install

    pip install casm-map


#### Usage

See the [casm docs](https://prisms-center.github.io/CASMcode_pydocs/casm/overview/latest/).


#### Instal Zsh Completions

To install Zsh completions, cp the completions function to your custom completions directory:

    cp completions/zsh/_casm_map ~/.oh-my-zsh/custom/completions/_casm-map

Then reload your shell or run:

    compinit


#### About CASM

The casm-map package is part of the [CASM](https://prisms-center.github.io/CASMcode_docs/) open source software package, which is designed to perform first-principles statistical mechanical studies of multi-component crystalline solids.

CASM is developed by the Van der Ven group, originally at the University of Michigan and currently at the University of California Santa Barbara.

For more information, see the [CASM homepage](https://prisms-center.github.io/CASMcode_docs/).


#### License

GNU Lesser General Public License (LGPL). Please see the file LICENSE for details.

