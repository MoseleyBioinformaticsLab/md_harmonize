#!/usr/bin/python3

"""
tools.py provides functions to open file or save data to file.

"""

import json
import jsonpickle
import signal
import multiprocessing
import typing as t
import pickle

jsonpickle.set_encoder_options('json', sort_keys=True, indent=4)


def timeout(func:t.Callable[...,t.Any], args:t.Optional[tuple] = None, seconds:int = 1) -> None:
    """Executes a given function but only for a certain specified time.
    :param func: function to call
    :param args: function arguments to use
    :param seconds: length of timeout
    :param default: default value if function does not execute properly.
    """

    def function_with_queue(queue:multiprocessing.Queue, *args):
        queue.put(func(*(args or ())))

    queue = multiprocessing.Queue()
    process = multiprocessing.Process(target=function_with_queue, args=(queue, *(args or ())))
    process.start()
    process.join(seconds)

    if process.is_alive():
        process.terminate()
        raise TimeoutError("{func_name} timeout taking more than {seconds} seconds.".format(func_name=func.__name__,
                                                                                            seconds=seconds))
    else:
        return queue.get()


class timeout_context:
    """Implements a timeout context manager to work with a with statement."""

    def __init__(self, seconds=1, error_message='Timeout'):
        """
        :param seconds:
        :type seconds: :py:class:`int`
        :param error_message:
        :type error_message: :py:class:`str`
        """
        self.seconds = seconds
        self.error_message = error_message

    def handle_timeout(self, signum, frame):
        raise TimeoutError(self.error_message)

    def __enter__(self):
        signal.signal(signal.SIGALRM, self.handle_timeout)
        signal.alarm(self.seconds)

    def __exit__(self, type, value, traceback):
        signal.alarm(0)

    # def __init__(self, seconds=1, error_message='Timeout'):
    #     """
    #     :param seconds:
    #     :type seconds: :py:class:int
    #     :param error_message:
    #     :type error_message: :py:class:str
    #     """
    #     self.seconds = seconds
    #     self.error_message = error_message
    #
    # def handle_timeout(self):
    #     raise TimeoutError(self.error_message)
    #
    # def __enter__(self):
    #     self.timer = threading.Timer(self.seconds, self.handle_timeout)
    #     self.timer.start()
    #
    # def __exit__(self, type, value, traceback):
    #     self.timer.cancel()


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


def save_to_pickle(data, filename: str) ->None:
    """
    To save the data into pickle.

    :param data: data to be saved
    :param filename: the file to save the data
    :return: None.
    """
    with open(filename, 'wb') as outfile:
        pickle.dump(data, outfile)


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


def open_pickle(filename: str):
    """
    To load data via pickle.

    :param filename: the file to be loaded.
    :return: the decoded data from the file.
    """
    with open(filename, 'rb') as infile:
        data = pickle.load(infile)
    return data



