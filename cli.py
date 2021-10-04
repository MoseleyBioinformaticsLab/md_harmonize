#!/usr/bin/python3
"""
MDH command-line interface

Usage:
    mdh -h | --help
    mdh --version
    mdh download <database_name> <working_directory>
    mdh standardize <working_directory>
    mdh aromatize <working_directory> <save_file> [--aromatic_manager=<aromatic_manager_file>]
    mdh initialize <working_directory>
    mdh harmonize <working_directory>

Options:
    -h, --help          Show this screen.
    --version           Show version.

"""
import glob
import ctfile
import os
from . import KEGG_database_scraper
from . import KEGG_parser
from . import MetaCyc_parser
from . import aromatics
from . import tools
from . import compound
from pybel import *


scraper_dict = {"KEGG": KEGG_database_scraper}
parser_dict = {"KEGG": KEGG_parser, "MetaCyc": MetaCyc_parser}

def cli(args):

    if args['download']:
        database_name = args['<database_name>']
        working_directory = args['<working_directory>'] + "/sources"
        path = os.path.join(working_directory, database_name)
        os.makedirs(path, exist_ok=True)
        scraper_dict[database_name].download(path)

    elif args['standardize']:
        # here to search the database directory, find the molfile subdirectory and create standardized_molfile subdirectory
        # and store the standardized molfile in it.
        working_directory = args['<working_directory>']
        subdirectories = [directory[0] for directory in os.walk(working_directory + "/sources")]
        for subdirectory in subdirectories:
            path, base = os.path.split(subdirectory)
            if base == "molfile":
                database_name = path.split("/")[-1]
                standardized_path = os.path.join(working_directory + "/standardized", database_name)
                os.makedirs(standardized_path, exist_ok=True)
                files = glob.glob(subdirectory + "*")
                for file in files:
                    # need to standardize with ctfile.
                    with open(file, "r") as infile:
                        ct_object = ctfile.load(infile)


    elif args['aromatize']:
        working_directory = args['<working_directory>']
        subdirectories = [directory[0] for directory in os.walk(working_directory + "/standardized")]
        kcf_path = working_directory + "/sources/KEGG/kcf"
        save_file = args['<save_file>']
        if args['--aromatic_manager']:
            aromatic_manager = tools.open_jsonpickle(args['--aromatic_manager'])
        else:
            aromatic_manager = aromatics.AromaticManager()
        for subdirectory in subdirectories:
            molfiles = glob.glob(subdirectory + "/*")
            for molfile in molfiles:
                aromatic_manager.indigo_aromatize(molfile)
        if os.path.exists(kcf_path):
            kcf_files = glob.glob(kcf_path + "/*")
            parser = parser_dict['KEGG']
            for kcf_file in kcf_files:
                kcf_cpd = parser.create_compound_kcf(kcf_file)
                aromatic_manager.kegg_aromatize(kcf_cpd)
        tools.save_to_jsonpickle(aromatic_manager, save_file)

    elif args['initialize']:
        working_directory = args['<working_directory>']
        subdirectories = [directory[0] for directory in os.walk(working_directory + "/standardized")]
        aromatic_manager = tools.open_jsonpickle(working_directory + "/aromatic_manager.json")
        save_directory = working_directory + "/initialized"
        os.makedirs(save_directory)
        for subdirectory in subdirectories:
            path, database_name = os.path.split(subdirectory)
            parser = parser_dict[database_name]

            # construct compounds
            compounds = {}
            molfiles = glob.glob(subdirectory + "/*")
            for molfile in molfiles:
                cpd = compound.Compound.create(molfile)
                compounds[cpd.compound_name] = cpd
            # add the kat to kegg compound
            if database_name == "KEGG":
                kcf_files = glob.glob(working_directory + "/sources/KEGG/kcf/*")
                for kcf_file in kcf_files:
                    kcf_compound = parser.create_compound_kcf(kcf_file)
                    parser.add_kat(compounds[kcf_compound.cpd.compound_name], kcf_compound)
            # update compound
            for cpd_name in compounds:
                cpd = compounds[cpd_name]
                aromatic_manager.detect_aromatic_substructures(cpd)
                cpd.define_bond_stereochemistry()

            reactions = []
            # construct reactions
            if database_name == "KEGG":
                reaction_directory = working_directory + "/sources/KEGG/reaction/"
                rclass_directory = path + "rclass"
                atom_mappings = parser.create_atom_mappings(rclass_directory, compounds)
                reactions = parser.create_reactions(reaction_directory, compounds, atom_mappings)

            elif database_name == "MetaCyc":
                reaction_file = working_directory + "/sources/MetaCyc/reaction.dat"
                atom_mapping_file = working_directory + "/sources/MetaCyc/atom-mapping.dat"
                reactions = parser.create_reactions(reaction_file, atom_mapping_file, compounds)

            tools.save_to_jsonpickle(compounds, save_directory + "/" + database_name + "compounds.json")
            tools.save_to_jsonpickle(reactions, save_directory + "/" + database_name + "reactions.json")

    elif args['harmonize']:
        # create harmonized compound manager first, then reaction harmonization manager.



























