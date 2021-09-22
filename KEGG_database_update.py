#!/usr/bin/python3

"""
update.py is used to download KEGG compound and reaction data from web.
"""


import requests
import os
import urllib
import tools

KEGG_compound_list_URL = "http://rest.kegg.jp/list/compound"
KEGG_reaction_list_URL = "http://rest.kegg.jp/list/reaction"
KEGG_rclass_list_URL = "http://rest.kegg.jp/list/rclass"
this_directory = os.path.abspath(os.path.dirname(__file__))
default_directory = '{}/data/KEGG/'.format(this_directory)


def entry_list(dir):

    file = urllib.request.urlopen(dir)
    the_list = []
    for line in file:
        the_list.append(line.decode('utf-8').split()[0])
    return the_list
	
def update_entity(entries, sub_directory, dir=default_directory, suffix=""):
    """

    """
    for entry in entries:
        data = requests.get("http://rest.kegg.jp/get/{0}/{1}".format(entry, suffix))
        tools.save_to_text(data.text, dir+sub_directory+entry)
       

# def multiprocess(entry_list, function):

#     with multiprocessing.Pool() as pool:
#         results = pool.map(function, (entry for entry in entry_list))


if __name__ == '__main__':

    update_entity(entry_list(KEGG_reaction_list_URL), "reaction/")
    update_entity(entry_list(KEGG_rclass_list_URL), "rclass/")
    update_entity(entry_list(KEGG_compound_list_URL), "compound_molfile/", suffix="mol")
    update_entity(entry_list(KEGG_compound_list_URL), "compound_kcf/", suffix="kcf")
    # compound_list()
    # multiprocess(compound_list(), update_compound)
    # multiprocess(reaction_list(), update_reaction)


