.. image:: _static/logo_outline.svg
  :alt: CASM logo
  :width: 600
  :class: only-light

.. image:: _static/logo_dark_outline.svg
  :alt: CASM logo
  :width: 600
  :class: only-dark


casm-map
========

The casm-map package provides pure Python tools for structure mapping and import. This includes:

- A library and command line tool that enables:
  - searching for structure mappings
  - creating mapped structures
  - creating interpolated structures
  - generating structures with equivalent relative orientations
  - importing structures into CASM projects.

This package makes extensive use of lower-level methods which are implemented in `libcasm-xtal <https://github.com/prisms-center/CASMcode_crystallography>`_, `libcasm-mapping <https://github.com/prisms-center/CASMcode_mapping>`_, and `libcasm-configuration <https://github.com/prisms-center/CASMcode_configuration>`_.

Methods for searching for low-cost lattice, atom, and structure mappings, taking into account symmetry are based on the approach described in :cite:t:`THOMAS2021a`.


About CASM
==========

The casm-map package is part of the CASM_ open source software package, which is designed to perform first-principles statistical mechanical studies of multi-component crystalline solids.

CASM is developed by the Van der Ven group, originally at the University of Michigan and currently at the University of California Santa Barbara.

For more information, see the `CASM homepage <CASM_>`_.


License
=======

GNU Lesser General Public License (LGPL). Please see the LICENSE file available on GitHub_.


Documentation
=============

.. toctree::
    :maxdepth: 2

    Installation <installation>
    Usage <usage>
    Reference <reference/casm/index>
    Bibliography <bibliography>

casm-map is available on GitHub_.

.. _CASM: https://prisms-center.github.io/CASMcode_docs/
.. _GitHub: https://github.com/prisms-center/CASMcode_crystallography
