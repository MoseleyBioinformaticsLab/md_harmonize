
#!/usr/bin/python3

"""
tools.py provides functions to open file or save data to file.

"""

import json
import jsonpickle
import errno
import os
import signal
from functools import wraps

jsonpickle.set_encoder_options('json', sort_keys=True, indent=4)


class TimeoutError(Exception):
    pass


def timeout(seconds=10, error_message=os.strerror(errno.ETIME)):
    def decorator(func):
        def _handle_timeout(signum, frame):
            raise TimeoutError(error_message)

        def wrapper(*args, **kwargs):
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.setitimer(signal.ITIMER_REAL,seconds) #used timer instead of alarm
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result
        return wraps(func)(wrapper)
    return decorator


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

