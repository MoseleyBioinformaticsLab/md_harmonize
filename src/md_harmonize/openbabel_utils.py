#!/usr/bin/python3


import ctfile
from pathlib import Path
import os


def standardize_molfile(molfile: str, to_path: str) -> None:
    """
    To standardize molfile using openbabel.

    :param molfile: the filename to the original molfile.
    :param to_path: the path to store the standardized molfile.
    :return: None.
    """
    compound_name = Path(molfile).stem
    tofile = "{0}/{1}.mol".format(to_path, compound_name)
    try:
        with open(molfile, "r") as infile:
            ct_object = ctfile.load(infile)
            ct_object.write(open(tofile, "w"), "ctfile")
    except UnicodeDecodeError:
        # Some MetaCyc molfiles require this type of encoding.
        with open(molfile, "r", encoding='cp1252') as infile:
            ct_object = ctfile.load(infile)
            ct_object.write(open(tofile, "w"), "ctfile")
    os.system("obabel {0} -O {0} -h -b".format(tofile))
