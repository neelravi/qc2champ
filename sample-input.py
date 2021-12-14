import cclib
import qc2champ
import numpy as np
import sys
import os

# TREXIO interface
try:
    import trexio
except:
    print("Error: The TREXIO Python library is not installed")
    sys.exit(1)



for filename in [sys.argv[1]]:
    basename = (filename.split(".", 1)[0])

    data = cclib.io.ccread(filename)

    print(dir(data))
    print(f"There are {data.natom} atoms and {data.nbasis} number of cartesian basis")



    print(f"The SCF energies {data.scfenergies}")
    print(f"The metadata {data.metadata}")
    print(f"The wfn information {data.wfn_info}")
    # print(f"The CI information {data.ci}")

    # print(f"The molecular coefficients {data.mocoeffs}")
    # print(f"The molecular symmetry mosyms  {data.mosyms}")
    # print(f"The molecular symmetry symm info {data.symm_info}")
    # print(f"The molecular coefficients {np.shape(data.mocoeffs)}")

    # qc2champ.io.ccwrite(data, outputtype="xyz", outputdest="output.xyz")
    # qc2champ.io.write_champ_v2_sym(data, outputdest="MOLCAS_" + basename)
    # qc2champ.io.write_champ_v2_geometry(data, outputdest="MOLCAS_" + basename)

    #if trexio in sys.modules:
        ## TREXIO All in One
#    qc2champ.io.write_trexio(data, outputdest="MOLCAS_" + basename)

    # qc2champ.io.write_champ_v2_lcao(data, outputdest="CN3_" + basename)
    # qc2champ.io.write_champ_v2_det(data, outputdest="MOLCAS_" + basename)

