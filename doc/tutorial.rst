The MDH Tutorial
================

The `MDH` is intended to harmonize compounds and reactions across public metabolic databases. It provides functionalities to:
    * Update a local copy of the KEGG metabolic databases.
    * Curate compound molfile representations.
    * Construct compound/reaction entities.
    * To detect aromatic substructures in a compound.
    * Parse KEGG RClass representations to generate atom mappings between compounds.
    * Harmonize compounds/reactions across databases.

In this document, each usage will be explained in details.

The MDH API tutorial
~~~~~~~~~~~~~~~~~~~~

Data preparation
----------------
The raw KEGG, MetaCyc, and HMDB is available on zenodo https://doi.org/10.5281/zenodo.7384576. Please download the data first and keep the hierarchy of directory. The following data is included in the directory:
    * KEGG compound molfiles, compound kcf's, rclasses, and reactions.
    * MetaCyc compound molfiles, atom-mapping and reactions.
    * HMDB compound molfiles.

Using MDH to download KEGG databases
-------------------------------------

Updates the KEGG metabolic databases, including compound kcf/molfile, reaction, and rclass.

.. code-block:: console

    mdh download KEGG <working_directory>

Using MDH to curate compound molfile representation
----------------------------------------------------

Curates molfile representations, eg: add H.

.. code-block:: console

    mdh standardize <database_names> <working_directory>

Note: Multiple database names can be provided with comma separation, eg: KEGG,MetaCyc.

Using MDH to construct compound/reaction entity
-----------------------------------------------

Constructs compound/reaction entities.

Compound construction:

.. code-block:: console

    mdh initialize_compound <database_names> <working_directory> <aromatic_manager_file> [--parse_kegg_atom]


Options
-------

--parse-kegg-atom:
The parse-kegg-atom option is used for parsing KEGG atom mapping between compounds based on KEGG RClass definitions.


Reaction construction:

.. code-block:: console

    mdh initialize_reaction <database_names> <working_directory>

Using MDH to extract aromatic substructures from compounds
----------------------------------------------------------

Extracts aromatic substructures from compounds using either Indigo or BASS to construct
aromatic manager that can be used to detect aromatic substructures in the compound.

.. code-block:: console

    mdh aromatize <database_names> <working_directory> <save_file> [--aromatic_manager=<aromatic_manager_file>]

Options
-------
--aromatic_manager option to indicate a pre-constructed aromatic manager is provided and newly discovered aromatic substructures will be added to aromatic manager.

Using MDH to harmonize compounds/reactions across databases
-----------------------------------------------------------

Harmonizes compounds/reactions across metabolic databases.

.. code-block:: console

    mdh harmonize_compound <database_names> <working_directory>
    mdh harmonize_reaction <database_names> <working_directory>



