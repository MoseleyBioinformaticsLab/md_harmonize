#!/usr/bin/python3
"""
MDH command-line interface

Usage:
    mdh -h | --help
    mdh --version
    mdh download <database_names> <working_directory>
    mdh standardize <database_names> <working_directory>
    mdh aromatize <database_names> <working_directory> <save_file> [--aromatic_manager=<aromatic_manager_file>]
    mdh initialize <database_names> <working_directory> <aromatic_manager_file>
    mdh harmonize <database_names> <working_directory>

Options:
    -h, --help          Show this screen.
    --version           Show version.

"""
import glob
import os
from . import KEGG_database_scraper
from . import KEGG_parser
from . import MetaCyc_parser
from . import aromatics
from . import tools
from . import compound
from . import openbabel_utils
from . import harmonization

scraper_dict = {"KEGG": KEGG_database_scraper}
parser_dict = {"KEGG": KEGG_parser, "MetaCyc": MetaCyc_parser}

def cli(args):

    if args['download']:
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        for database_name in database_names:
            path = os.path.join(working_directory + "/sources", database_name)
            os.makedirs(path, exist_ok=True)
            try:
                scraper_dict[database_name].download(path)
            except KeyError:
                raise KeyError("Scraper for {0} database does not exits".format(database_name))

    elif args['standardize']:
        # here to search the database directory, find the molfile subdirectory and create standardized_molfile subdirectory
        # and store the standardized molfile in it.
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        for database_name in database_names:
            from_path = working_directory + "sources/{0}/molfile".format(database_name)
            if not os.path.exists(from_path):
                raise OSError("The directory {0} does not exist.".format(from_path))
            to_path = working_directory + "/standardized/{0}/molfile".format(database_name)
            os.makedirs(to_path, exist_ok=True)
            molfiles = glob.glob(from_path + "/*")
            for molfile in molfiles:
                openbabel_utils.standardize_molfile(molfile, to_path)

    elif args['aromatize']:
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        save_file = args['<save_file>']
        if args['--aromatic_manager']:
            aromatic_manager = tools.open_jsonpickle(args['--aromatic_manager'])
        else:
            aromatic_manager = aromatics.AromaticManager()
        for database_name in database_names:
            from_path = working_directory + "standardized/{0}/molfile".format(database_name)
            if not os.path.exists(from_path):
                raise OSError("The directory {0} does not exist.".format(from_path))
            molfiles = glob.glob(from_path + "/*")
            for molfile in molfiles:
                aromatic_manager.indigo_aromatize(molfile)
        if "KEGG" in database_names:
            kcf_path = working_directory + "/sources/KEGG/kcf"
            if os.path.exists(kcf_path):
                kcf_files = glob.glob(kcf_path + "/*")
                parser = parser_dict['KEGG']
                for kcf_file in kcf_files:
                    kcf_cpd = parser.create_compound_kcf(kcf_file)
                    aromatic_manager.kegg_aromatize(kcf_cpd)
        tools.save_to_jsonpickle(aromatic_manager, save_file)

    elif args['initialize']:
        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        aromatic_manager = tools.open_jsonpickle(args['<aromatic_manager_file>'])
        to_directory = working_directory + "/initialized"
        for database_name in database_names:
            from_path =  working_directory + "standardized/{0}/molfile".format(database_name)
            if not os.path.exists(from_path):
                raise OSError("The directory {0} does not exist.".format(from_path))
            parser = parser_dict[database_name]
            # construct compounds
            compound_list = {}
            molfiles = glob.glob(from_path + "/*")
            for molfile in molfiles:
                cpd = compound.Compound.create(molfile)
                compound_list[cpd.compound_name] = cpd
            # add the kat to kegg compound
            if database_name == "KEGG":
                kcf_files = glob.glob(working_directory + "/sources/KEGG/kcf/*")
                for kcf_file in kcf_files:
                    kcf_compound = parser.create_compound_kcf(kcf_file)
                    parser.add_kat(compound_list["cpd:" + kcf_compound.cpd.compound_name], kcf_compound)
            # update compound
            for cpd_name in compound_list:
                cpd = compound_list[cpd_name]
                aromatic_manager.detect_aromatic_substructures(cpd)
                cpd.define_bond_stereochemistry()

            reaction_list = []
            # construct reactions
            if database_name == "KEGG":
                reaction_directory = working_directory + "/sources/KEGG/reaction/"
                rclass_directory = working_directory + "/sources/KEGG/rclass/"
                if not os.path.exists(reaction_directory):
                    raise OSError("The directory {0} does not exist.".format(reaction_directory))
                if not os.path.exists(rclass_directory):
                    raise OSError("The directory {0} does not exist.".format(rclass_directory))
                atom_mappings = parser.create_atom_mappings(rclass_directory, compound_list)
                reaction_list = parser.create_reactions(reaction_directory, compound_list, atom_mappings)

            elif database_name == "MetaCyc":
                reaction_file = working_directory + "/sources/MetaCyc/reaction.dat"
                atom_mapping_file = working_directory + "/sources/MetaCyc/atom-mapping.dat"
                if not os.path.exists(reaction_file):
                    raise OSError("The file {0} does not exist.".format(reaction_file))
                if not os.path.exists(atom_mapping_file):
                    raise OSError("The file {0} does not exist.".format(atom_mapping_file))
                reaction_list = parser.create_reactions(reaction_file, atom_mapping_file, compound_list)

            save_directory = to_directory + "/{0}".format(database_name)
            os.makedirs(save_directory, exist_ok=True)
            tools.save_to_jsonpickle(compound_list, save_directory + "/compounds.json")
            tools.save_to_jsonpickle(reaction_list, save_directory + "/reactions.json")

    elif args['harmonize']:
        # create harmonized compound manager first, then reaction harmonization manager.

        database_names = args['<database_names>'].split(",")
        working_directory = args['<working_directory>']
        compound_list = []
        reaction_list = []
        for database_name in database_names:
            from_directory = working_directory + "/initialized/{0}/".format(database_name)
            if not os.path.exists(from_directory + "compounds.json"):
                raise OSError("The file {0} does not exist.".format(from_directory + "compounds.json"))
            if not os.path.exists(from_directory + "reactions.json"):
                raise OSError("The file {0} does not exist.".format(from_directory + "reactions.json"))
            compound_list.append(tools.open_jsonpickle(from_directory + "compounds.json"))
            reaction_list.append(tools.open_jsonpickle(from_directory + "reactions.json"))

        compound_harmonization_manager = harmonization.harmonize_compound_list(compound_list)
        # while we do reaction harmonization, we need to pay attention to compound harmonization without same structural representations.
        # this includes: R group, linear-circular-transformation, resonance.
        reaction_harmonization_manager = harmonization.harmonize_reaction_list(reaction_list, compound_harmonization_manager)
        tools.save_to_jsonpickle(reaction_harmonization_manager, working_directory + "/harmonized/{0}.json".format("_".join(database_names)))




























