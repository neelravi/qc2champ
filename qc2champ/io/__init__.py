# -*- coding: utf-8 -*-
#
"""Contains all writers for standard chemical representations."""


from cclib.io.xyzreader import XYZ as XYZReader
from cclib.io.xyzwriter import XYZ as XYZWriter

from cclib.io.ccio import ccopen
from cclib.io.ccio import ccread
from cclib.io.ccio import ccwrite
from qc2champ.io.qcio import write_champ_sym
from qc2champ.io.qcio import write_champ_eigenvalues
from qc2champ.io.qcio import write_champ_geometry
from qc2champ.io.qcio import write_champ_v2_lcao
from qc2champ.io.qcio import write_champ_v2_det

from qc2champ.io.trexio import write_trexio
from cclib.io.ccio import URL_PATTERN
