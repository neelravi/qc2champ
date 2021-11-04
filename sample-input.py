import cclib
import qc2champ
import numpy as np
import sys
import os



for filename in [ "sample_acrolein_cas.out"]:
    basename = (filename.split(".", 1)[0])



    data = cclib.io.qcread(filename)

    print(dir(data))
    print(f"There are {data.natom} atoms and {data.nbasis} number of cartesian basis")
    # print(f"There are core electrons {data.coreelectrons} ")

    # print(f"There are {data.nelectrons} ")

    # print(f"There are {data.number_alpha} ")
    # print(f"There are {data.number_beta} ")

    # print(f"There are {data.number_alpha_valence} ")
    # print(f"There are {data.number_beta_valence} ")

    # print(f"There are {data.scftype} ")



#    for element in range(len(data.atomnos)):
#        print(data.atomnos[element], data.atomcoords[0][element])


    # print(f"The MO coefficients {data.mocoeffs[0][0]}")

    # print(f"The MO energies {data.moenergies}")

    # print(f"The MO symmetries {data.nmo}")
#    print(f"The MO symmetries {data.mosyms}")

    # values, counts = numpy.unique(data.mosyms, return_counts=True)

    # print (len(counts))

    # print(f"unique irreducible representations {dict(zip(values,counts))}")

    # print(f"The basis set {data.gbasis[0]}")

    # print(f"The SCF energies {data.scfenergies}")

    # print(f"The metadata {data.metadata}")


#    qc2champ.io.qcwrite(data, outputtype="xyz", outputdest="output.xyz")



    qc2champ.io.write_champ_old_sym(data, outputdest="CN3_" + basename)
    qc2champ.io.write_champ_old_geo(data, outputdest="CN3_" + basename)
    qc2champ.io.write_champ_old_lcao(data, outputdest="CN3_" + basename)
    qc2champ.io.write_champ_old_det(data, outputdest="CN3_" + basename)

