# -*- coding: utf-8 -*-
#
"""Tools for writing files into champ format."""


def write_champ_v2_sym(ccobj, outputdest=None):
    """Writes the parsed geometry, symmetry, determinants, MO coefficients data from the quantum
    chemistry calculation to old format of champ .sym, .geom, .det, and .lcao file.

    Inputs:
        ccobj - Either a job (from ccopen) or a data (from job.parse()) object
        outputdest - A filename or file object for writing. Example, "rhf.sym", "cn3.sym"

    Returns:
        None as a function value
    """
    # If the output filename is mentioned, then write to that file
    # This will write in the old format that CHAMP recognizes.


    if outputdest is not None:
        if isinstance(outputdest, str):
            ## Write down a symmetry file in old champ format
            with open(outputdest + ".sym", 'w') as file:

                values, counts = np.unique(ccobj.mosyms, return_counts=True)
                # point group symmetry independent line printed below
                file.write("sym_labels " + str(len(counts)) + " " + str(len(ccobj.mosyms[0]))+"\n")

                irrep_string = ""
                irrep_correspondence = {}
                for i, val in enumerate(values):
                    irrep_correspondence[val] = i+1
                    irrep_string += " " + str(i+1) + " " + str(val)

                if all(irreps in ccobj.mosyms[0] for irreps in values):
                    file.write(f"{irrep_string} \n")   # This defines the rule

                    for item in ccobj.mosyms[0]:
                        for key, val in irrep_correspondence.items():
                            if item == key:
                                file.write(str(val)+" ")
                    file.write("\n")
                file.write("end\n")
            file.close()

        elif isinstance(outputdest, io.IOBase):
            outputdest.write(output)
        else:
            raise ValueError
    # If outputdest is None, return a string representation of the output.
    else:
        return None


def write_champ_v2_geometry(ccobj, outputdest=None):
    """Writes the parsed geometry data from the quantum
    chemistry calculation to old format of champ .geo  file.

    Inputs:
        ccobj - Either a job (from ccopen) or a data (from job.parse()) object
        outputdest - A filename or file object for writing. Example, "rhf.geo", "cn3.geo"

    Returns:
        None as a function value
    """


    # If the output filename is mentioned, then write to that file
    # This will write in the old format that CHAMP recognizes.


    if outputdest is not None:
        if isinstance(outputdest, str):
            ## Write down a symmetry file in old champ format
            with open(outputdest + ".xyz", 'w') as file:

                # header line printed below
                file.write(str(ccobj.natom) + "\n" )
                file.write("# Coordinates are in Bohr units \n")

                coords = [[ccobj.atomcoords[0][i][j] for j in range(3)] for i in range(len(ccobj.atomnos))]
                coords = np.array(coords)/0.5291772109 #angstrom_to_bohr conversion

               for element in range(len(ccobj.atomnos)):
                   file.write("{: 0.6f} {: 0.6f} {: 0.6f} {} \n".format(coords[element][0], coords[element][1], coords[element][2], indices[element]+1))



            file.close()

        elif isinstance(outputdest, io.IOBase):
            outputdest.write(output)
        else:
            raise ValueError
    # If outputdest is None, return a string representation of the output.
    else:
        return None


def write_champ_v2_lcao(ccobj, outputdest=None):
    """Writes the parsed geometry data from the quantum
    chemistry calculation to old format of champ .lcao  file.

    Inputs:
        ccobj - Either a job (from ccopen) or a data (from job.parse()) object
        outputdest - A filename or file object for writing. Example, "rhf.lcao", "cn3.lcao"

    Returns:
        None as a function value
    """


    # If the output filename is mentioned, then write to that file
    # This will write in the old format that CHAMP recognizes.


    if outputdest is not None:
        if isinstance(outputdest, str):
            ## Write down a symmetry file in old champ format
            with open(outputdest + ".lcao", 'w') as file:

                # header line printed below
                file.write("# Comments about the system being studied \n")
                file.write("lcao " + str(len(ccobj.mocoeffs[0][0])) + " " + str(len(ccobj.mocoeffs[0][0])) + "\n" )
                np.savetxt(file, ccobj.mocoeffs[0][:])
                file.write("end\n")
            file.close()

        elif isinstance(outputdest, io.IOBase):
            outputdest.write(output)
        else:
            raise ValueError
    # If outputdest is None, return a string representation of the output.
    else:
        return None


def write_champ_v2_det(ccobj, outputdest=None):
    """Writes the parsed determinants data from the quantum
    chemistry calculation to old format of champ .det file.

    Inputs:
        ccobj - Either a job (from ccopen) or a data (from job.parse()) object
        outputdest - A filename or file object for writing. Example, "rhf.det", "cn3.det"

    Returns:
        None as a function value
    """

    # If the output filename is mentioned, then write to that file
    # This will write in the old format that CHAMP recognizes.


    if outputdest is not None:
        if isinstance(outputdest, str):
            ## Write down a symmetry file in old champ format
            with open(outputdest + ".det", 'w') as file:

                if ccobj.scftype in ["RHF", "UHF", "ROHF"]:

                    file.write(f"&electrons  nelec {ccobj.number_alpha_valence+ccobj.number_beta_valence} nup {ccobj.number_alpha_valence} \n" )
                    file.write(f"\n" )
                    file.write(f"determinants {1} {1} \n")
                    file.write(f"      {1:.6f} \n")

                    alpha_occupation = np.arange(ccobj.number_alpha_valence) + 1
                    beta_occupation  = np.arange(ccobj.number_alpha_valence) + 1
                    np.savetxt(file, np.row_stack((alpha_occupation, beta_occupation)), fmt='  %i', delimiter='  ', newline='')

                    file.write("\n")
                    file.write("end\n")
                    file.close()

                elif ccobj.scftype in ["MCSCF"]:
                    raise Warning("being implemented")

        elif isinstance(outputdest, io.IOBase):
            outputdest.write(output)
        else:
            raise ValueError
    # If outputdest is None, return a string representation of the output.
    else:
        return None



