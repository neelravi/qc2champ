# -*- coding: utf-8 -*-
#
# Copyright (c) 2021, Ravindra Shinde
#
# This file is part of TREX (http://trex-coe.eu) and is distributed under
# the terms of the BSD 3-Clause License.
"""Tools for identifying, reading and writing files"""

import sys
import numpy as np
from cclib.parser.utils import PeriodicTable
from cclib.parser.utils import convertor



try:
    import trexio
except:
    print("Error: The TREXIO Python library is not installed")
    sys.exit(1)


def write_trexio(ccobj, outputdest=None, motype=None):
    """Writes the parsed information from the quantum
    chemistry calculation to trexio hdf5 file.

    Inputs:
        ccobj - Either a job (from ccopen) or a data (from job.parse()) object
        outputdest - A filename or file object for writing. Example, "rhf.geo", "cn3.geo"

    Returns:
        None as a function value
    """

    # If the output filename is mentioned, then write to that file
    # This will write in the new format that CHAMP recognizes.

    ## trexio part
    if outputdest is not None:
        if isinstance(outputdest, str):
            trexio_file = trexio.File(outputdest + ".hdf5", mode='w', back_end=trexio.TREXIO_HDF5)

            ## Metadata -------------------------
            trexio.write_metadata_description(trexio_file, "Trexio file written using qc2trexio package authored by Ravindra Shinde")

            if 'package' in ccobj.metadata.keys():
                trexio.write_metadata_code_num(trexio_file, 1)
                trexio.write_metadata_code(trexio_file, [ccobj.metadata['package']])

            if 'author' in ccobj.metadata.keys():
                trexio.write_metadata_author_num(trexio_file, 1)
                trexio.write_metadata_author(trexio_file, [ccobj.metadata['author']])

            ## Electron group -------------------
            # ## Revisit this group. Don't make it package specific
            # trexio.write_electron_up_num(trexio_file, int(ccobj.wfn_info['nelec_closed_shell']))
            # trexio.write_electron_dn_num(trexio_file, int(ccobj.wfn_info['nelec_closed_shell']))


            ## Nucleus group ---------------------
            #  Number of nuclei
            if ccobj.natom is not None:
                trexio.write_nucleus_num(trexio_file, ccobj.natom)

            #  nucleus coordinates
            if ccobj.atomcoords is not None:
                trexio.write_nucleus_coord(trexio_file, convertor(ccobj.atomcoords.flatten(), "Angstrom", "bohr"))

            #  nucleus labels
            if ccobj.atomnos is not None:
                trexio.write_nucleus_charge(trexio_file, ccobj.atomnos)
                pt = PeriodicTable()
                element_list = [pt.element[Z] for Z in ccobj.atomnos]
                trexio.write_nucleus_label(trexio_file, element_list)

            ## Molecular Group -------------------
            trexio.write_mo_type(trexio_file, motype)  # parse this from the file.

            # number of AO, number of MO
            if ccobj.mocoeffs[0] is not None:
                trexio.write_ao_num(trexio_file, ccobj.mocoeffs[0].shape[0])
                trexio.write_mo_num(trexio_file, ccobj.nmo)

            # MO symmetry irreps
            if ccobj.mosyms[0] is not None:
                trexio.write_mo_symmetry(trexio_file, ccobj.mosyms[0])

            # MO eigenvalues
            # if ccobj.moenergies[0] is not None:
            #     trexio.write_mo_energy(trexio_file, ccobj.moenergies[0])

            # MO coefficients
            if ccobj.mocoeffs[0] is not None:
                trexio.write_mo_coefficient(trexio_file, ccobj.mocoeffs[0])

        else:
            raise ValueError
    # If outputdest is None, return a string representation of the output.
    else:
        return None

