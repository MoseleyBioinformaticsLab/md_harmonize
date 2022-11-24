User Guide
==========

Description
~~~~~~~~~~~

MDH was created to achieve creation and harmonization of compounds and reactions in
the public metabolic databases. It is a command line tool that provides


Installation
~~~~~~~~~~~~

The MDH package runs under Python 3.7+. Use pip_ to install.

Install on Linux, Mac OS X
--------------------------

.. code:: bash

   python3 -m pip install MDH


Get the source code
~~~~~~~~~~~~~~~~~~~

Code is available on GitHub: https://github.com/MoseleyBioinformaticsLab/MDH.git

You can clone the public repository:

.. code:: bash

   $ https://github.com/MoseleyBioinformaticsLab/MDH.git


Dependencies
~~~~~~~~~~~~

MDH requires the following Python libraries:

    * docopt_ for creating the command-line interface.
    * jsonpickle_ for saving Python objects in a JSON serializable form and outputting to a file.
    * numpy_ and cython_ for speeding optimization.
    * ctfile_ for parsing compound molfile representation.
    * indigo_ for detecting aromatic atoms in the compound.
    * pebble_ for multiprocessing of cythonized calculation.


Data
~~~~

The raw data from KEGG and MetaCyc databases can be accessed from this URL.


Basic usage
~~~~~~~~~~~

MDH provides functions to achieve compound and reaction harmonization across public metabolic databases. Details about
the usages are in the :doc:`tutorial`.

.. code-block:: console

 Usage:
    mdh -h | --help
    mdh --version
    mdh download <database_names> <working_directory>
    mdh standardize <database_names> <working_directory>
    mdh aromatize <database_names> <working_directory> <save_file> [--aromatic_manager=<aromatic_manager_file>]
    mdh initialize_compound <database_names> <working_directory> <aromatic_manager_file> [--parse_kegg_atom]
    mdh initialize_reaction <database_names> <working_directory>
    mdh harmonize <database_names> <working_directory>

 Options:
    -h, --help          Show this screen.
    --version           Show version.


.. _GitHub: https://github.com/MoseleyBioinformaticsLab/MDH
.. _jsonpickle: https://github.com/jsonpickle/jsonpickle
.. _pip: https://pip.pypa.io/
.. _docopt: https://pypi.org/project/docopt/
.. _cython: https://github.com/cython/cython
.. _numpy: https://github.com/numpy/numpy
.. _ctfile: https://github.com/MoseleyBioinformaticsLab/ctfile
.. _indigo: https://github.com/epam/Indigo
.. _pebble: https://pypi.org/project/Pebble/