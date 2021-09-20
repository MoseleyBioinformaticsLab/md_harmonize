#!/usr/bin/python3

"""
update.py is used to download KEGG compound and reaction data from web.
"""

import multiprocessing
import requests
import os
import urllib

KEGG_compound_list_URL = "http://rest.kegg.jp/list/compound"
KEGG_reaction_list_URL = "http://rest.kegg.jp/list/reaction"
this_directory = os.path.abspath(os.path.dirname(__file__))
default_directory = '{}/data/KEGG/'.format(this_directory)

def save_to_file(content, dir, filename):

	with open(dir+filename, "w") as outfile:
		outfile.write(content)

def compound_list(dir=KEGG_compound_list_URL):
	
	file = urllib.request.urlopen(dir)
	the_list = []
	for line in file:
		the_list.append(line.decode("utf-8").split()[0])
	return the_list
	
def update_compound(compound, dir=default_directory):
	"""

	"""
	compound_molfile = requests.get("http://rest.kegg.jp/get/{0}/mol".format(compound))
	compound_kcf = requests.get("http://rest.kegg.jp/get/{0}/kcf".format(compound))
	save_to_file(compound_molfile.text, dir+"compound_molfile/", "{0}.mol".format(compound))
	save_to_file(compound_kcf.text, dir+"compound_kcf/", "{0}.kcf".format(compound))


def reaction_list(dir=KEGG_reaction_list_URL):

	file = urllib.request.urlopen(dir)
	the_list = []
	for line in file:
		the_list.append(line.decode("utf-8").split()[0])
	return the_list

def update_reaction(reaction, dir=default_directory):

	reaction_content = requests.get("http://rest.kegg.jp/get/{0}".format(reaction))
	save_to_file(reaction_content.text, dir+"reaction/", reaction)

def multiprocess(entry_list, function):

	with multiprocessing.Pool() as pool:
		results = pool.map(function, (entry for entry in entry_list))


if __name__ == '__main__':
	
	# compound_list()
	multiprocess(compound_list(), update_compound)
	multiprocess(reaction_list(), update_reaction)


