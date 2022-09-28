import pytest
import cclib
import qc2champ
import numpy as np
import sys
import os



def test_cclib():
	assert data.natom == 2
	assert data.nbasis == 4
	np.testing.assert_almost_equal(data.scfenergies[0], -30.552003252900377)


def test_qc2champ():
	assert data.ci['Number of CSFs'] == '2'
	assert data.ci['Number of determinants'] == '2'
	assert data.ci['Number of root(s) required'] == '1'
	assert data.wfn_info['nelec_active_shell'] == '2'
	assert data.metadata['package_version'] == '8.4'


for filename in ["sample_h2_rasscf.out"]:
	basename = (filename.split(".", 1)[0])
	data = cclib.io.ccread(filename)
	test_cclib()
	test_qc2champ()
	# qc2champ.io.ccwrite(data, outputtype="xyz", outputdest="output.xyz")
	# qc2champ.io.write_champ_v2_geometry(data, outputdest="CN3_" + basename)
	# qc2champ.io.write_champ_v2_sym(data, outputdest="CN3_" + basename)
	# qc2champ.io.write_trexio(data, outputdest="CN3_" + basename)
	# qc2champ.io.write_champ_v2_det(data, outputdest="CN3_" + basename)