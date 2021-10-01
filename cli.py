#!/usr/bin/python3
"""
MDH command-line interface

Usage:
    mdh -h | --help
    mdh --version
    mdh download <database_name> <save_directory>
    mdh standardize <file_directory> <save_directory>
    mdh aromatize <molfile> <save_directory> [--aromatic_manager=<aromatic_manager_file>] [--kcf=<kcf_directory>]
    mdh initialize <database_name> <molfile> <reaction> <save_directory> <aromatic_manager_file> [--kcf=<kcf_directory>]
    mdh harmonize

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

        molfile_directory = args['<molfile>']
        molfiles = glob.glob(molfile_directory + "*")
        # to save the updated aromatic manager
        save_directory = args['<save_directory>']
        if args['--aromatic_manager']:
            aromatic_manager = tools.open_jsonpickle(args['--aromatic_manager'])
        else:
            aromatic_manager = aromatics.AromaticManager()
        for molfile in molfiles:

        if args['--kcf']:




    elif args['initialize']:
        database_name = args['<database_name>']












