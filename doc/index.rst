.. image:: _static/logo_outline.svg
  :alt: CASM logo
  :width: 600
  :class: only-light

.. image:: _static/logo_dark_outline.svg
  :alt: CASM logo
  :width: 600
  :class: only-dark

casm-tools
==========

The casm-tools package provides CLI tools based on capabilities implemented in CASM.
This includes:

- casm.tools.calc: Setup, run, and report results of structure calculations
- casm.tools.convert: Convert structures and configurations between CASM and other
  formats
- casm.tools.map: Structure mapping and import


casm-map
--------

The casm-map command line program performs structure mapping and import, including:

- Searching for structure mappings
- Creating mapped structures
- Creating interpolated structures
- Generating structures with equivalent relative orientations
- Importing structures into CASM projects as CASM configurations with calculated
  properties.

This program makes extensive use of lower-level methods which are implemented in
`libcasm-xtal <https://github.com/prisms-center/CASMcode_crystallography>`_,
`libcasm-mapping <https://github.com/prisms-center/CASMcode_mapping>`_, and
`libcasm-configuration <https://github.com/prisms-center/CASMcode_configuration>`_.

Methods for searching for low-cost lattice, atom, and structure mappings, taking into
account symmetry are based on the approach described in :cite:t:`THOMAS2021a`.

A suggested way to cite this program is as follows:

.. code-block:: text

    "Structure mappings were found with the `casm-map` program [1],
    using the method of Thomas et al. [2] implemented in CASM [3]."

    1. B. Puchala, J. Thomas, and A. Van der Ven, "casm-map...".
    2. J. C. Thomas, A. R. Natarajan, and A. Van der Ven, Comparing
        crystal structures with symmetry and geometry, npj
        Computational Materials, 7 (2021), 164.
    2. B. Puchala, J. C. Thomas, A. R. Natarajan, J. G. Goiri,
        S. S. Behara, J. L. Kaufman, A. Van der Ven, CASM—A software
        package for first-principles based study of multicomponent
        crystalline solids, Computational Materials Science 217
        (2023) 111897.


About CASM
==========

The casm-tools package is part of the CASM_ open source software suite, which is
designed to perform first-principles statistical mechanical studies of multi-component
crystalline solids.

CASM is developed by the Van der Ven group, originally at the University of Michigan
and currently at the University of California Santa Barbara.

For more information, see the `CASM homepage <CASM_>`_.


License
=======

GNU Lesser General Public License (LGPL). Please see the LICENSE file available on
GitHub_.


Documentation
=============

.. toctree::
    :maxdepth: 2

    Installation <installation>
    Usage <usage>
    Reference <reference/casm/index>
    Bibliography <bibliography>

casm-tools is available on GitHub_.

.. _CASM: https://prisms-center.github.io/CASMcode_docs/
.. _GitHub: https://github.com/prisms-center/CASMcode_tools
