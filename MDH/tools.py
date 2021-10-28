#!/usr/bin/python3


import json
import jsonpickle

jsonpickle.set_encoder_options('json', sort_keys=True, indent=4)

def save_to_text(data, filename):
    """
    To save the data in a text file.

    :param data:
    :param filename:
    :return:
    """
    with open(filename, 'w') as outfile:
        outfile.write(data)

def save_to_jsonpickle(data, filename):
    """
    To save the data via jsonpickle.
    :param data:
    :param filename:
    :return:
    """
    with open(filename, 'w') as outfile:
        outfile.write(jsonpickle.encode(data, keys=True))

def save_to_json(data, filename):
    """
    To save the data into json.
    :param data:
    :param filename:
    :return:
    """
    with open(filename, 'w') as outfile:
        json.dump(data, outfile, indent=4)

def open_text(filename, encoding='utf-8'):
    """
    To open the text file.
    :param filename:
    :param encoding:
    :return:
    """
    with open(filename, 'r', encoding=encoding) as infile:
        data = infile.read()
    return data

def open_jsonpickle(filename):
    """
    To open the data via jsonpickle.
    :param filename:
    :return:
    """
    with open(filename, 'r') as infile:
        data = jsonpickle.decode(infile.read(), keys=True)
    return data

def open_json(filename):
    """
    To open the data via json.
    :param filename:
    :return:
    """
    with open(filename, 'r') as infile:
        data = json.load(infile)
    return data


if __name__ == '__main__':

    data = open_text("/mlab/data/hji236/projects/metabolic_database_harmonization/MetaCyc_parser.py").split("\n")
    for line in data:
        print(line)

