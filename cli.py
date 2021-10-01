#!/usr/bin/python3
"""
MDH command-line interface

Usage:
    mdh -h | --help
    mdh --version
    mdh download <database_name> <save_directory>
    mdh standardize <file_directory> <save_directory>
    mdh aromatize <database_name> <molfile> <save_directory> [--kcf=<kcf_directory>] [--aromatic_manager=<aromatic_manager_file>]
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


def cli(args):

    if args['download']:
        scraper_dict = {"KEGG": KEGG_database_scraper}
        database_name = args['<database_name>']
        save_directory = args['<save_directory']
        scraper_dict[database_name].download(save_directory)

    elif args['standardize']:
        file_directory = args['<file_directory>']
        save_directory = args['<save_directory>']
        files = glob.glob(file_directory + "*")
        for file in files:

    elif args['aromatize']:
        database_name = args['<database_name>']
        if database_name == "KEGG":

        else:

    elif args['initialize']:
        database_name = args['<database_name>']












