
#!/usr/bin/python3

"""
tools.py provides functions to open file or save data to file.

"""

import json
import jsonpickle
import multiprocessing
import multiprocessing.pool

jsonpickle.set_encoder_options('json', sort_keys=True, indent=4)


def timeout(seconds=10):
    def decorator(func):
        def wrapper(*args, **kwargs):
            pool = multiprocessing.pool.ThreadPool(processes=1)
            results = pool.apply_async(func, args, kwargs)
            pool.close()
            try:
                return results.get(seconds)
            except multiprocessing.TimeoutError:
                raise OSError('Timeout expired after: %s' % seconds)
            finally:
                pool.terminate()

        return wrapper

def save_to_text(data: str, filename: str) -> None:
    """
    To save the data in a text file.

    :param data: data to be saved.
    :param filename: the file to save the data.
    :return: None.
    """
    with open(filename, 'w') as outfile:
        outfile.write(data)


def save_to_jsonpickle(data, filename: str) -> None:
    """
    To save the data via jsonpickle.

    :param data: data to be saved.
    :param filename: the file to save the data.
    :return: None.
    """
    with open(filename, 'w') as outfile:
        outfile.write(jsonpickle.encode(data, keys=True))


def save_to_json(data, filename: str) -> None:
    """
    To save the data into json.

    :param data: data to be saved.
    :param filename: the file to save the data.
    :return: None.
    """
    with open(filename, 'w') as outfile:
        json.dump(data, outfile, indent=4)


def open_text(filename: str, encoding: str = 'utf-8') -> str:
    """
    To load text file.

    :param filename: the file to be loaded.
    :param encoding: The name of the encoding used to decode the stream’s bytes into strings.
    :return: the decoded data from the file.
    """
    with open(filename, 'r', encoding=encoding) as infile:
        data = infile.read()
    return data


def open_jsonpickle(filename: str):
    """
    To load data via jsonpickle.

    :param filename: the file to be loaded.
    :return: the decoded data from the file.
    """
    with open(filename, 'r') as infile:
        data = jsonpickle.decode(infile.read(), keys=True)
    return data


def open_json(filename: str):
    """
    To load data via json.

    :param filename: the file to be loaded.
    :return: the decoded data from the file.
    """
    with open(filename, 'r') as infile:
        data = json.load(infile)
    return data

