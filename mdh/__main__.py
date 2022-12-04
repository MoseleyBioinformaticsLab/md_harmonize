#!/usr/bin/python3

import docopt
from . import cli
from . import __version__


def main():

    cli.cli(args)


if __name__ == "__main__":

    args = docopt.docopt(cli.__doc__, version=__version__)
    import cProfile
    profiler = cProfile.Profile()
    profiler.enable()
    main()
    profiler.disable()
    profiler.print_stats()
