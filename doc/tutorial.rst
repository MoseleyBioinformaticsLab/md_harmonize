The MDH Tutorial
================

The `MDH` is intended to harmonize compound and reaction across public metabolic databases. It provides functionalities to:
    * Update KEGG database.
    * Curate compound molfile representation.
    * Construct compound/reaction entity.
    * To detect aromatic substructures in a compound.
    * Parse KEGG RClass representation to generate atom mappings between compounds.
    * Harmonize compounds/reactions across databases.

In this document, each usage will be explained in details.

The MDH API tutorial
~~~~~~~~~~~~~~~~~~~~

Data preparation
----------------
The raw KEGG, MetaCyc and HMD is available on zenodo https://doi.org/10.5281/zenodo.7384576. Please download the data first and keep the hierarchy of directory. The following data is included in the directory:
    * KEGG compound molfile, compound kcf, rclass, and reaction.
    * MetaCyc compound molfile, atom-mapping and reactions.
    * HMD compound molfile.

Using MDH to download KEGG databases
-------------------------------------

The `MDH` provides function to update the KEGG databases, including compound kcf/molfile, reaction, and rclass.

.. code-block:: console

    python3 -m mdh download KEGG <working_directory>

Using MDH to curate compound molfile representation
----------------------------------------------------

The `MDH` provides function to curate molfile representations, eg: add H.

.. code-block:: console

    python3 -m mdh standardize <database_names> <working_directory>

Note: Multiple database names can be provided with "," separation, eg: KEGG,MetaCyc.

Using MDH to construct compound/reaction entity
-----------------------------------------------

The `MDH` provides function to construct compound/reaction entities.

Compound construction:
.. code-block:: console

    python3 -m mdh initialize_compound <database_names> <working_directory> <aromatic_manager_file> [--parse_kegg_atom]


Options
-------

--parse-kegg-atom:
The parse-kegg-atom option is used for parsing KEGG atom mapping between compounds based on KEGG RClass definitions.

Reaction construction:
.. code-block:: console

    python3 -m mdh initialize_reaction <database_names> <working_directory>

Using MDH to extract aromatic substructures from compounds
----------------------------------------------------------

The `MDH` provides function to extract aromatic substructures from compounds using either indigo/BASS to construct
aromatic manager that can be used to detect aromatic substructures in the compound.

.. code-block:: console

    python3 -m mdh aromatize <database_names> <working_directory> <save_file> [--aromatic_manager=<aromatic_manager_file>]

Options
-------
--aromatic_manager option to indicate a pre-constructed aromatic manager is provided and newly discovered aromatic substructures will be added to aromatic manager.

Using MDH to harmonize compounds/reactions across databases
-----------------------------------------------------------

The `MDH` provides function to harmonize compounds/reactions across metabolic databases.

.. code-block:: console

    python3 -m mdh harmonize_compound <database_names> <working_directory>
    python3 -m mdh harmonize_reaction <database_names> <working_directory>



