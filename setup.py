#!/usr/bin/python3

import os
import sys
import re
from setuptools import setup, find_packages, Extension
import Cython.Build
import numpy


def find_version():
    with open('mdh/__init__.py', 'r') as fd:
        version = re.search(r'^__version__\s*=\s*[\'"]([^\'"]*)[\'"]',
                            fd.read(), re.MULTILINE).group(1)
        if not version:
            raise RuntimeError('Cannot find version information')
        return version

def readme():
    with open('README.rst') as readme_file:
        return readme_file.read()


REQUIRES = [
        "docopt",
        "jsonpickle",
        "numpy",
        "cython",
        "epam.indigo",
        "ctfile",
        "timeout-decorator",
        "cython"
]


EXTENSIONS = [
    Extension("mdh.BASS",sources=["mdh/BASS.pyx"], extra_compile_args=['-O3'], include_dirs=[numpy.get_include()])
]

setup(
        name='MDH',
        version=find_version(),
        packages=find_packages(),
        author='Huan Jin, Hunter N.B. Moseley',
        author_email='jinhuan0905@gmail.com, hunter.moseley@gmail.com',
        description='Methods for harmonizing public metabolic databases',
        keywords='metabolite, metabolic reaction',
        license='Modified Clear BSD License',
        url='',
        ext_modules=EXTENSIONS,
        cmdclass={'build_ext': Cython.Build.build_ext},
        install_requires=REQUIRES,
        long_description=readme(),
        platforms='any',
        classifiers=[
            'Development Status :: 4 - Beta',
            'Environment :: Console',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: BSD License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            ],
)
