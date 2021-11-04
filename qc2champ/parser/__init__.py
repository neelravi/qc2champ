# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, Ravindra Shinde
#
"""Contains parsers for all supported programs"""


# These import statements are added for the convenience of users...
# Rather than having to type:
#         from cclib.parser.gaussianparser import Gaussian
# they can use:
#         from cclib.parser import Gaussian

from cclib.parser.adfparser import ADF
from cclib.parser.daltonparser import DALTON
from cclib.parser.fchkparser import FChk
from cclib.parser.gamessparser import GAMESS
from cclib.parser.gamessukparser import GAMESSUK
from cclib.parser.gaussianparser import Gaussian
from cclib.parser.jaguarparser import Jaguar
# from cclib.parser.molcasparser import Molcas
from qc2champ.parser.molcas84parser import Molcas84
from cclib.parser.molproparser import Molpro
from cclib.parser.mopacparser import MOPAC
from cclib.parser.nwchemparser import NWChem
from cclib.parser.orcaparser import ORCA
from cclib.parser.psi3parser import Psi3
from cclib.parser.psi4parser import Psi4
from cclib.parser.qchemparser import QChem
from cclib.parser.turbomoleparser import Turbomole
from qc2champ.parser import logfileparser
from cclib.parser.data import ccData

# This allows users to type:
#         from cclib.parser import ccopen
from cclib.io.ccio import ccopen
