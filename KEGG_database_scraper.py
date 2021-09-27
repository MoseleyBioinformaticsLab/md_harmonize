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


def entry_list(target_url):
    """
    To get the list of entity name to download.

    :param target_url: the url to fetch.
    :return:
    """

    file = urllib.request.urlopen(target_url)
    the_list = []
    for line in file:
        the_list.append(line.decode('utf-8').split()[0])
    return the_list
	
def update_entity(entries, sub_directory, directory=default_directory, suffix=""):
    """
    To download the KEGG entity (compound, reaction, or rclass) and save it into a file.

    :param entries: the list of entry name
    :param sub_directory: the subdirectory to save the downloaded file.
    :param directory: the main directory to save the downloaded file
    :param suffix: the suffix needed for download, like the mol for compound molfile and kcf for compound kcf file.
    :return:
    """
    for entry in entries:
        data = requests.get("http://rest.kegg.jp/get/{0}/{1}".format(entry, suffix))
        tools.save_to_text(data.text, directory + sub_directory + entry)
       

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


