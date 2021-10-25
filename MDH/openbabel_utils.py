#!/usr/bin/python3


import ctfile
from pathlib import Path
import os

def standardize_molfile(molfile, to_path):

    with open(molfile, "r") as infile:
        ct_object = ctfile.load(infile)
        compound_name = Path(molfile).stem
        tofile = "{0}/{1}.mol".format(to_path, compound_name)
        ct_object.write(open(tofile, "w"), "ctfile")
        os.system("obabel {0} -O {0} -h".format(tofile))


