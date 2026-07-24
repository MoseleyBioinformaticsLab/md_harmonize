User Guide
==========

Description
~~~~~~~~~~~

md_harmonize was created to harmonize compounds and reactions in the public metabolic databases. It is a command line tool
directly supporting harmonization of KEGG, MetaCyc, and HMDB, but has an API for expanding to other metabolic databases.


Installation
~~~~~~~~~~~~

The md_harmonize package runs under Python 3.7+. Use pip_ to install.

Install on Linux, Mac OS X
--------------------------

.. code:: bash

   python3 -m pip install md-harmonize


Get the source code
~~~~~~~~~~~~~~~~~~~

Code is available on GitHub: https://github.com/MoseleyBioinformaticsLab/md_harmonize.git

You can clone the public repository:

.. code:: bash

   $ https://github.com/MoseleyBioinformaticsLab/md_harmonize.git


Dependencies
~~~~~~~~~~~~

md_harmonize requires the following Python libraries:

    * docopt_ for creating the command-line interface.
    * jsonpickle_ for saving Python objects in a JSON serializable form and outputting to a file.
    * numpy_ and cython_ for speeding optimization.
    * ctfile_ for parsing compound molfile representation.
    * indigo_ for detecting aromatic atoms in the compound.
    * pebble_ for multiprocessing of cythonized calculation.


Data
~~~~

The raw data from KEGG and MetaCyc databases can be accessed from this URL: https://doi.org/10.5281/zenodo.7384576


Basic usage
~~~~~~~~~~~

md_harmonize provides functions to achieve compound and reaction harmonization across public metabolic databases. Details about
the usages are in the :doc:`tutorial`.

.. code-block:: console

 Usage:
    md_harmonize -h | --help
    md_harmonize --version
    md_harmonize download <database_names> <working_directory>
    md_harmonize standardize <database_names> <working_directory>
    md_harmonize aromatize <database_names> <working_directory> <save_file> [--pickle] [--aromatic_manager=<aromatic_manager_file>]
    md_harmonize initialize_compound <database_names> <working_directory> <aromatic_manager_file> [--parse_kegg_atom] [--pickle] [--split=k]
    md_harmonize initialize_reaction <database_names> <working_directory> [--pickle]
    md_harmonize harmonize_compound <database_names> <working_directory> [--pickle]
    md_harmonize harmonize_reaction <database_names> <working_directory> [--pickle]

 Options:
    -h, --help           Show this screen.
    --version            Show version.
    --aromatic_manager=<aromatic_manager_file>   A pre-constructed aromatic manager is provided.
    --pickle             Use pickle to save the results, otherwise use jsonpickle.
    --parse_kegg_atom    To parse KEGG atom mapping between compounds.
    --split=k              Split compounds to speed up construction.


.. _GitHub: https://github.com/MoseleyBioinformaticsLab/MDH
.. _jsonpickle: https://pypi.org/project/jsonpickle/
.. _pip: https://pip.pypa.io/
.. _docopt: https://pypi.org/project/docopt/
.. _cython: https://pypi.org/project/Cython/
.. _numpy: https://pypi.org/project/numpy/
.. _ctfile: https://pypi.org/project/ctfile/
.. _indigo: https://pypi.org/project/epam.indigo/
.. _pebble: https://pypi.org/project/Pebble/
