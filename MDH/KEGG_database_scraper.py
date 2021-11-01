#!/usr/bin/python3

"""
KEGG_database_scraper.py is used to download KEGG data (including compound, reaction, kcf, and rclass) from web.

This URL can change.
"""

import requests
import urllib
import os

from . import tools

KEGG_compound_list_URL = "http://rest.kegg.jp/list/compound"
KEGG_reaction_list_URL = "http://rest.kegg.jp/list/reaction"
KEGG_rclass_list_URL = "http://rest.kegg.jp/list/rclass"

def entry_list(target_url):
    """
    To get the list of entity name to download.

    :param target_url: the url to fetch.
    :type target_url: :py:class:`str`.
    :return: the list of entry names.
    :rtype: :py:obj:`list`.
    """
    file = urllib.request.urlopen(target_url)
    the_list = []
    for line in file:
        the_list.append(line.decode('utf-8').split()[0])
    return the_list
	
def update_entity(entries, sub_directory, directory, suffix=""):
    """
    To download the KEGG entity (compound, reaction, or rclass) and save it into a file.

    :param entries: the list of entry name to download.
    :type entries: :py:obj:`list`.
    :param sub_directory: the subdirectory to save the downloaded file.
    :type sub_directory: :py:obj:`str`.
    :param directory: the main directory to save the downloaded file.
    :type directory: :py:obj:`str`.
    :param suffix: the suffix needed for download, like the mol for compound molfile and kcf for compound kcf file.
    :type suffix: :py:obj:`str`.
    :return: None.
    :rtype: :py:obj:`None`.
    """
    path = os.path.join(directory, sub_directory)
    if not os.path.exists(path):
        os.mkdir(path)
    for entry in entries:
        data = requests.get("http://rest.kegg.jp/get/{0}/{1}".format(entry, suffix))
        tools.save_to_text(data.text, path + "/" + entry)

def download(directory):
    """
    To down load all the KEGG required files.

    :param directory: the directory to store the data.
    :type directory: :py:obj:`str`.
    :return: None.
    :rtype: :py:obj:`None`.
    """
    update_entity(entry_list(KEGG_reaction_list_URL), "reaction", directory)
    update_entity(entry_list(KEGG_rclass_list_URL), "rclass", directory)
    update_entity(entry_list(KEGG_compound_list_URL), "molfile", directory, suffix="mol")
    update_entity(entry_list(KEGG_compound_list_URL), "kcf", directory, suffix="kcf")




