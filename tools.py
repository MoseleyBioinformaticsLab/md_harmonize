#!/usr/bin/python3


import json
import jsonpickle


def save_to_text(data, filename):
    
    with open(filename, 'w') as outfile:
        outfile.write(data)

def save_to_jsonpickle(data, filename):
	
	with open(filename, 'w') as outfile:
        outfile.write(jsonpickle.encode(data, keys=True))

def save_to_json(data, filename):

	with open(filename, 'w') as outfile:
        json.dump(data, outfile, indent=4)

def open_text(filename, encoding='utf-8'):

	with open(filename, 'r', encoding=encoding) as infile:
		data = infile.read()
	return data

def open_jsonpickle(filename):

	with open(filename, 'r') as infile:
        data = jsonpickle.decode(infile.read(), keys=True)
    return data

def open_json(filename):

	with open(filename, 'r') as infile:
        data = json.load(infile)
    return data




