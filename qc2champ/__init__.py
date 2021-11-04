# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, Ravindra Shinde
#
# This file is part of TREX (http://trex-coe.eu) and is distributed under
# the terms of the BSD 3-Clause License.

"""
A converter for parsing and interpreting results from computational chemistry packages,
and generating input files for CHAMP code.

The qc2champ relies on cclib library to parse and extract data from the output files
of computational chemistry packages. New write functions are added to get
the input files for CHAMP code.

Currently supported programs whose output files are parsed by cclib are:
    ADF, Firefly, GAMESS(US), GAMESS-UK, Gaussian,
    Jaguar, Molpro, MOPAC, NWChem, ORCA, Psi, Q-Chem, MOLCAS

"""

__version__ = "1.0.0"

from qc2champ import io
import cclib

# The test module can be imported if it was installed with cclib.
try:
    from qc2champ import test
except ImportError:
    pass

# The objects below constitute our public API. These names will not change
# over time. Names in the sub-modules will typically also be backwards
# compatible, but may sometimes change when code is moved around.
qcopen = cclib.io.ccopen
qcopen = cclib.io.ccopen
qcwrite = cclib.io.ccwrite
