#!/usr/bin/python3
"""
MDH command-line interface

Usage:
    mdh -h | --help
    mdh --version
    mdh download <database_name> <save_directory>
    mdh standardize <file_directory> <save_directory>
    mdh aromatize <molfile> <save_directory> [--aromatic_manager=<aromatic_manager_file>] [--kcf=<kcf_directory>]
    mdh initialize <database_name> <molfile> <reaction> <atom_mapping> <save_directory> <aromatic_manager_file> [--kcf=<kcf_directory>]
    mdh harmonize <

Options:
    -h, --help          Show this screen.
    --version           Show version.

"""
import glob
import ctfile
from . import KEGG_database_scraper
from . import KEGG_parser
from . import MetaCyc_parser
from . import aromatics
from . import tools
from . import compound


scraper_dict = {"KEGG": KEGG_database_scraper}
parser_dict = {"KEGG": KEGG_parser, "MetaCyc": MetaCyc_parser}

def cli(args):

    if args['download']:
        database_name = args['<database_name>']
        save_directory = args['<save_directory']
        scraper_dict[database_name].download(save_directory)

    elif args['standardize']:
        file_directory = args['<file_directory>']
        save_directory = args['<save_directory>']
        files = glob.glob(file_directory + "*")
        for file in files:
            # need to standardize with ctfile.
            with open(file, "r") as infile:
                ct_object = ctfile.load(infile)


    elif args['aromatize']:
        molfiles = glob.glob(args['<molfile>'] + "*")
        # to save the updated aromatic manager
        save_directory = args['<save_directory>']
        if args['--aromatic_manager']:
            aromatic_manager = tools.open_jsonpickle(args['--aromatic_manager'])
        else:
            aromatic_manager = aromatics.AromaticManager()
        for molfile in molfiles:
            aromatic_manager.indigo_aromatize(molfile)
        if args['--kcf']:
            kcf_files = glob.glob(args['--kcf'] + "*")
            parser = parser_dict['KEGG']
            for kcf_file in kcf_files:
                kcf_cpd = parser.create_compound_kcf(kcf_file)
                aromatic_manager.kegg_aromatize(kcf_cpd)
        tools.save_to_jsonpickle(aromatic_manager, save_directory)

    elif args['initialize']:
        database_name = args['<database_name>']
        parser = parser_dict[database_name]
        molfiles = glob.glob(args['<molfile>'] + "*")
        aromatic_manager = tools.open_jsonpickle(args['--aromatic_manager'])
        compounds = {}
        for molfile in molfiles:
            cpd = compound.Compound.create(molfile)
            compounds[cpd.compound_name] = cpd
        # add the kat to kegg compound
        if args['--kcf']:
            kcf_files = glob.glob(args['--kcf'] + "*")
            for kcf_file in kcf_files:
                kcf_compound = parser.create_compound_kcf(kcf_file)
                parser.add_kat(compounds[kcf_compound.cpd.compound_name], kcf_compound)
        # update compound
        for cpd_name in compounds:
            cpd = compounds[cpd_name]
            aromatic_manager.detect_aromatic_substructures(cpd)
            cpd.define_bond_stereochemistry()

    elif args['harmonize']:

        # create harmonized compound manager first, then reaction harmonization manager.

























