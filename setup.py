# -*- coding: utf-8 -*-
#
# Copyright (c) 2021 Ravindra Shinde (TREX CoE)
#
# This file is part of TREX (http://trex-coe.eu) and is distributed under
# the terms of the BSD 3-Clause License.

"""qc2champ: Convert a quantum chmical output file to a CHAMP input file.

qc2champ is a Python library that connverts outputs from various
quantum chemical calculations to CHAMP input files. It relies on cclib
for parsers.
"""

import setuptools


# Chosen from http://www.python.org/pypi?:action=list_classifiers
classifiers = """Development Status :: 5 - Production/Stable
Environment :: Console
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: BSD License
Natural Language :: English
Operating System :: OS Independent
Programming Language :: Python
Topic :: Scientific/Engineering :: Chemistry
Topic :: Scientific/Engineering :: Physics
Topic :: Software Development :: Libraries :: Python Modules"""


def setup_qc2champ():

    doclines = __doc__.split("\n")

    setuptools.setup(
        name="qc2champ",
        version="1.0.0",
        url="http://qc2champ.github.io/",
        author="Ravindra Shinde",
        author_email="neelravi@gmail.com",
        maintainer="Ravindra Shinde",
        maintainer_email="neelravi@gmail.com",
        license="BSD 3-Clause License",
        description=doclines[0],
        long_description="\n".join(doclines[2:]),
        classifiers=classifiers.split("\n"),
        platforms=["Any."],
        packages=setuptools.find_packages(exclude=['*test*']),
        entry_points={
            'console_scripts': [
                'qcget=qc2champ.scripts.qcget:qcget',
                'qcwrite=qc2champ.scripts.qcwrite:main',
            ]
        },
        install_requires=[
            "packaging>=19.0",
            "numpy",
            "periodictable",
            "scipy>=1.2.0",
            "cclib>=1.7.0",
        ],

    )


if __name__ == '__main__':
    setup_qc2champ()
