#!/usr/bin/python3

"""
tools.py provides functions to open file or save data to file.

"""

import json
import jsonpickle

jsonpickle.set_encoder_options('json', sort_keys=True, indent=4)

def save_to_text(data, filename):
    """
    To save the data in a text file.

    :param data: data to be saved.
    :type data: :py:class:`str`.
    :param filename: the file to save the data.
    :type filename: :py:class:`str`.
    :return: None.
    :rtype: :py:obj:`None`.
    """
    with open(filename, 'w') as outfile:
        outfile.write(data)

def save_to_jsonpickle(data, filename):
    """
    To save the data via jsonpickle.

    :param data: data to be saved.
    :type data: :py:class:`obj`.
    :param filename: the file to save the data.
    :type filename: :py:class:`str`.
    :return: None.
    :rtype: :py:obj:`None`.
    """
    with open(filename, 'w') as outfile:
        outfile.write(jsonpickle.encode(data, keys=True))

def save_to_json(data, filename):
    """
    To save the data into json.

    :param data: data to be saved.
    :type data: :py:class:`obj`.
    :param filename: the file to save the data.
    :type filename: :py:class:`str`.
    :return: None.
    :rtype: :py:obj:`None`.
    """
    with open(filename, 'w') as outfile:
        json.dump(data, outfile, indent=4)

def open_text(filename, encoding='utf-8'):
    """
    To load text file.

    :param filename: the file to be loaded.
    :type filename: :py:class:`str`.
    :param encoding: The name of the encoding used to decode the stream’s bytes into strings.
    :type encoding: :py:class:`str`.
    :return: the decoded data from the file.
    :rtype: :py:class:`str`.
    """
    with open(filename, 'r', encoding=encoding) as infile:
        data = infile.read()
    return data

def open_jsonpickle(filename):
    """
    To load data via jsonpickle.

    :param filename: the file to be loaded.
    :type filename: :py:class:`str`.
    :return: the decoded data from the file.
    :rtype: :py:class:`obj`.
    """
    with open(filename, 'r') as infile:
        data = jsonpickle.decode(infile.read(), keys=True)
    return data

def open_json(filename):
    """
    To load data via json.

    :param filename: the file to be loaded.
    :type filename: :py:class:`str`.
    :return: the decoded data from the file.
    :rtype: :py:class:`obj`.
    """
    with open(filename, 'r') as infile:
        data = json.load(infile)
    return data

