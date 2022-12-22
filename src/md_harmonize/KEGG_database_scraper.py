#!/usr/bin/python3

"""
md_harmonize.KEGG_database_scraper
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This module provides functions to download KEGG data (including compound, reaction, kcf, and rclass) from the KEGG (REST) API.

The URLs can change.
"""

import requests
import urllib
import os

from . import tools

KEGG_compound_list_URL = "http://rest.kegg.jp/list/compound"
KEGG_reaction_list_URL = "http://rest.kegg.jp/list/reaction"
KEGG_rclass_list_URL = "http://rest.kegg.jp/list/rclass"


def entry_list(target_url: str) -> list:
    """
    To get the list of entity name to download.

    :param target_url: the url to fetch.
    :return: the list of entry names.
"""
    file = urllib.request.urlopen(target_url)
    the_list = []
    for line in file:
        the_list.append(line.decode('utf-8').split()[0])
    return the_list


def update_entity(entries: list, sub_directory: str, directory: str, suffix: str = "") -> None:
    """
    To download the KEGG entity (compound, reaction, or rclass) and save it into a file.

    :param entries: the list of entry names to download.
    :param sub_directory: the subdirectory to save the downloaded file.
    :param directory: the main directory to save the downloaded file.
    :param suffix: the suffix needed for download, like the mol for compound molfile and kcf for compound kcf file.
    :return: None.
"""
    path = os.path.join(directory, sub_directory)
    if not os.path.exists(path):
        os.mkdir(path)
    for entry in entries:
        data = requests.get("http://rest.kegg.jp/get/{0}/{1}".format(entry, suffix))
        text = data.text
        tools.save_to_text(text, path + "/" + entry)
        if sub_directory == "molfile":
            curate_molfile(path + "/" + entry)


def curate_molfile(file_path: str) -> None:
    """
    To curate the molfile representation.

    :param file_path: the path to the molfile.
    :return: None.
"""
    new_text = ""
    with open(file_path, "r") as infile:
        for line in infile:
            new_text += line
            if line.startswith("M  END"):
                break
    tools.save_to_text(new_text, file_path)
    return 


def download(directory: str) -> None:
    """
    To download all the KEGG required files.

    :param directory: the directory to store the data.
    :return: None.
"""
    update_entity(entry_list(KEGG_reaction_list_URL), "reaction", directory)
    update_entity(entry_list(KEGG_rclass_list_URL), "rclass", directory)
    update_entity(entry_list(KEGG_compound_list_URL), "molfile", directory, suffix="mol")
    update_entity(entry_list(KEGG_compound_list_URL), "kcf", directory, suffix="kcf")
    



