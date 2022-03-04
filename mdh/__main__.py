#!/usr/bin/python3

import docopt
from . import cli
from . import __version__

if __name__ == "__main__":

    args = docopt.docopt(cli.__doc__, version=__version__)
    import cProfile
    profiler = cProfile.Profile()
    profiler.enable()
    cli.cli(args)
    profiler.disable()
    profiler.print_stats()
