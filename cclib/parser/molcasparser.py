# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Molcas output files"""

import re
import string
import math
import numpy
from numpy.core.defchararray import center
from collections import Counter, OrderedDict

from cclib.parser import logfileparser
from cclib.parser import utils
from cclib.parser.utils import uniquelist
import sys


class Molcas(logfileparser.Logfile):
    """A Molcas log file."""

    def __init__(self, *args, **kwargs):

        # Call the __init__ method of the superclass
        super(Molcas, self).__init__(logname="Molcas", *args, **kwargs)

    def __str__(self):
        """Return a string repeesentation of the object."""
        return f"Molcas log file {self.filename}"

    def __repr__(self):
        """Return a representation of the object."""
        return f'Molcas("{self.filename}")'

    def normalisesym(self, label):
        """Normalise the symmetries used by Molcas.

        The labels are standardized except for the first character being lowercase.
        """
        return label[0].upper() + label[1:]

    def after_parsing(self):
        for element, ncore in self.core_array:
            self._assign_coreelectrons_to_element(element, ncore)

        if "package_version" in self.metadata:
            # Use the short version as the legacy version.
            self.metadata["legacy_package_version"] = self.metadata["package_version"]
            # If there is both a tag and the full hash, place the tag
            # first. Both are chosen to be local, since there isn't a
            # distinction between development and release builds in their
            # version cycle.
            if "tag" in self.metadata and "revision" in self.metadata:
                self.metadata["package_version"] = f"{self.metadata['package_version']}+{self.metadata['tag']}.{self.metadata['revision']}"
            elif "tag" in self.metadata:
                self.metadata["package_version"] = f"{self.metadata['package_version']}+{self.metadata['tag']}"
            elif "revision" in self.metadata:
                self.metadata["package_version"] = f"{self.metadata['package_version']}+{self.metadata['revision']}"

    def before_parsing(self):
        # Compile the regex for extracting the element symbol from the
        # atom label in the "Molecular structure info" block.
        self.re_atomelement = re.compile(r'([a-zA-Z]+)\d?')

        # Compile the dashes-and-or-spaces-only regex.
        self.re_dashes_and_spaces = re.compile(r'^[\s-]+$')

        # Molcas can do multiple calculations in one job, and each one
        # starts from the gateway module. Onle parse the first.
        # TODO: It would be best to parse each calculation as a separate
        # ccData object and return an iterator - something for 2.x
        self.gateway_module_count = 0

    def extract(self, inputfile, line):
        """Extract information from the file object inputfile."""

        if "Start Module: gateway" in line:
            self.gateway_module_count += 1

        if self.gateway_module_count > 1:
            return

        self.metadata["author"] = "Ravindra Shinde"
        # Extract the version number and optionally the Git tag and hash.

        if "Version" in line:
            match = re.search(r"\s{2,}Version:?\s(\d*\.\d*)", line)
            if match:
                self.metadata["package_version"] = match.groups()[0]
        if "version" in line:
            match = re.search(r"\s{2,}version:?\s(\d*\.\d*)", line)
            if match:
                self.metadata["package_version"] = match.groups()[0]
        if "tag" in line:
            self.metadata["tag"] = line.split()[-1]
        if "build" in line:
            match = re.search(r"\*\s*build\s(\S*)\s*\*", line)
            if match:
                self.metadata["revision"] = match.groups()[0]
        if "is licensed to" in line:
            self.metadata["owner"] = " ".join(line.split()[-2:])

        ## This section is present when executing &GATEWAY.
        # ++    Molecular structure info:
        #       -------------------------

        #                     ************************************************
        #                     **** Cartesian Coordinates / Bohr, Angstrom ****
        #                     ************************************************

        #      Center  Label                x              y              z                     x              y              z
        #         1      C1               0.526628      -2.582937       0.000000              0.278679      -1.366832       0.000000
        #         2      C2               2.500165      -0.834760       0.000000              1.323030      -0.441736       0.000000
        if line[25:63] == 'Cartesian Coordinates / Bohr, Angstrom':
            if not hasattr(self, 'atomnos'):
                self.atomnos = []

            self.skip_lines(inputfile, ['stars', 'blank', 'header'])

            line = next(inputfile)

            atomelements = []
            atomcoords = []

            while line.strip() not in ('', '--'):
                sline = line.split()
                atomelement = sline[1].rstrip(string.digits).title()
                atomelements.append(atomelement)
                atomcoords.append(list(map(float, sline[5:])))
                line = next(inputfile)

            self.append_attribute('atomcoords', atomcoords)

            if self.atomnos == []:
                self.atomnos = [self.table.number[ae.title()] for ae in atomelements]

            if not hasattr(self, 'natom'):
                self.set_attribute('natom', len(self.atomnos))

        # ## This section is present when executing &SCF.
        # # ++    Orbital specifications:
        # #       -----------------------

        # #       Symmetry species                           1   2   3   4
        # #                                                 a1  b1  a2  b2
        # #       Frozen orbitals                            0   0   0   0
        # #       Inactive orbitals                          6   5   0   0
        # #       Active orbitals                            0   0   4   6
        # #       Secondary orbitals                        45  39  10  12
        # #       Deleted orbitals                           0   0   0   0
        # #       Number of basis functions                 51  44  14  18
        # # --

        # if not hasattr(self, 'symm_info'):
        #     self.symm_info = {}

        # if line[:29] == '++    Orbital specifications:':

        #     self.skip_lines(inputfile, ['dashes', 'blank'])
        #     line = next(inputfile)

        #     symmetry_count = 1
        #     while not line.startswith('--'):
        #         if line.strip().startswith('Symmetry species'):
        #             symmetry_count = int(line.split()[-1])
        #             self.symm_info["symmetry_count"] = int(symmetry_count)
        #             line = next(inputfile)
        #             irreps = []
        #             for i in range(symmetry_count):
        #                 irreps.append(line.split()[i])
        #             self.symm_info["irreps"] = irreps

        #         if line.strip().startswith('Frozen orbitals'):
        #             frozen_orbitals = []
        #             for i in range(2,symmetry_count+2):
        #                 frozen_orbitals.append(line.split()[i])
        #             self.symm_info["frozen_orbitals"] = list(map(int, frozen_orbitals))

        #         if line.strip().startswith('Inactive orbitals'):
        #             inactive_orbitals = []
        #             for i in range(2,symmetry_count+2):
        #                 inactive_orbitals.append(line.split()[i])
        #             self.symm_info["inactive_orbitals"] = list(map(int, inactive_orbitals))

        #         if line.strip().startswith('Active orbitals'):
        #             active_orbitals = []
        #             for i in range(2,symmetry_count+2):
        #                 active_orbitals.append(line.split()[i])
        #             self.symm_info["active_orbitals"] = list(map(int, active_orbitals))

        #         if line.strip().startswith('Total number of orbitals'):
        #             nmos = line.split()[-symmetry_count:]
        #             self.set_attribute('nmo', sum(map(int, nmos)))
        #             orbitals_per_irrep = []
        #             for i in range(4,symmetry_count+4):
        #                 orbitals_per_irrep.append(line.split()[i])
        #             self.symm_info["orbitals_per_irrep"] = list(map(int, orbitals_per_irrep))

        #         if line.strip().startswith('Number of basis functions'):
        #             nbasis = line.split()[-symmetry_count:]
        #             self.set_attribute('nbasis', sum(map(int, nbasis)))
        #             basis_per_irrep = []
        #             for i in range(4,symmetry_count+4):
        #                 basis_per_irrep.append(line.split()[i])
        #             self.symm_info["basis_per_irrep"] = list(map(int, basis_per_irrep))
        #         line = next(inputfile)
        #     self.set_attribute('symm_info', self.symm_info)

        # print ("after symm info block ", self.symm_info)

        if line.strip().startswith(('Molecular charge', 'Total molecular charge')):
            self.set_attribute('charge', int(float(line.split()[-1])))

        #  ++    Molecular charges:
        #  ------------------

        #  Mulliken charges per centre and basis function type
        #  ---------------------------------------------------

        #         C1
        #  1s     2.0005
        #  2s     2.0207
        #  2px    0.0253
        #  2pz    0.1147
        #  2py    1.8198
        #  *s    -0.0215
        #  *px    0.0005
        #  *pz    0.0023
        #  *py    0.0368
        #  *d2+   0.0002
        #  *d1+   0.0000
        #  *d0    0.0000
        #  *d1-   0.0000
        #  *d2-   0.0000
        #  *f3+   0.0000
        #  *f2+   0.0001
        #  *f1+   0.0000
        #  *f0    0.0001
        #  *f1-   0.0001
        #  *f2-   0.0000
        #  *f3-   0.0003
        #  *g4+   0.0000
        #  *g3+   0.0000
        #  *g2+   0.0000
        #  *g1+   0.0000
        #  *g0    0.0000
        #  *g1-   0.0000
        #  *g2-   0.0000
        #  *g3-   0.0000
        #  *g4-   0.0000
        #  Total  6.0000

        #  N-E    0.0000

        #  Total electronic charge=    6.000000

        #  Total            charge=    0.000000
        #--
        if line[:24] == '++    Molecular charges:':

            atomcharges = []

            while line[6:29] != 'Total electronic charge':
                line = next(inputfile)
                if line[6:9] == 'N-E':
                    atomcharges.extend(map(float, line.split()[1:]))

            # Molcas only performs Mulliken population analysis.
            self.set_attribute('atomcharges', {'mulliken': atomcharges})

            # Ensure the charge printed here is identical to the
            # charge printed before entering the SCF.
            self.skip_line(inputfile, 'blank')
            line = next(inputfile)
            assert line[6:30] == 'Total            charge='
            if hasattr(self, 'charge'):
                assert int(float(line.split()[2])) == self.charge

        # This section is present when executing &SCF
        # This section parses the total SCF Energy.
        # *****************************************************************************************************************************
        # *                                                                                                                           *
        # *                                             SCF/KS-DFT Program, Final results                                             *
        # *                                                                                                                           *
        # *                                                                                                                           *
        # *                                                                                                                           *
        # *                                                       Final Results                                                       *
        # *                                                                                                                           *
        # *****************************************************************************************************************************

        # ::    Total SCF energy                                -37.6045426484
        if line[:22] == '::    Total SCF energy' or line[:25] == '::    Total KS-DFT energy':
            if not hasattr(self, 'scfenergies'):
                self.scfenergies = []
            scfenergy = float(line.split()[-1])
            self.scfenergies.append(utils.convertor(scfenergy, 'hartree', 'eV'))

        ## Parsing the scftargets in this section
        #  ++    Optimization specifications:
        #  ----------------------------

        #  SCF Algorithm: Conventional
        #  Minimized density differences are used
        #  Number of density matrices in core                9
        #  Maximum number of NDDO SCF iterations           400
        #  Maximum number of HF  SCF iterations            400
        #  Threshold for SCF energy change            0.10E-08
        #  Threshold for density matrix               0.10E-03
        #  Threshold for Fock matrix                  0.15E-03
        #  Threshold for linear dependence            0.10E-08
        #  Threshold at which DIIS is turned on       0.15E+00
        #  Threshold at which QNR/C2DIIS is turned on 0.75E-01
        #  Threshold for Norm(delta) (QNR/C2DIIS)     0.20E-04
        if line[:34] == '++    Optimization specifications:':
            self.skip_lines(inputfile, ['d', 'b'])
            line = next(inputfile)
            if line.strip().startswith('SCF'):
                scftargets = []
                self.skip_lines(inputfile,
                                ['Minimized', 'Number', 'Maximum', 'Maximum'])
                lines = [next(inputfile) for i in range(7)]
                targets = [
                    'Threshold for SCF energy change',
                    'Threshold for density matrix',
                    'Threshold for Fock matrix',
                    'Threshold for Norm(delta)',
                ]
                for y in targets:
                    scftargets.extend([float(x.split()[-1]) for x in lines if y in x])

                self.append_attribute('scftargets', scftargets)

        #  ++ Convergence information
        #                                     SCF        iterations: Energy and convergence statistics
        #
        #  Iter     Tot. SCF       One-electron     Two-electron   Energy   Max Dij or  Max Fij    DNorm      TNorm     AccCon    Time
        #             Energy          Energy          Energy       Change   Delta Norm                                          in Sec.
        #     1    -36.83817703    -50.43096166     13.59278464  0.00E+00   0.16E+00*  0.27E+01*   0.30E+01   0.33E+02   NoneDa    0.
        #     2    -36.03405202    -45.74525152      9.71119950  0.80E+00*  0.14E+00*  0.93E-02*   0.26E+01   0.43E+01   Damp      0.
        #     3    -37.08936118    -48.41536598     11.32600480 -0.11E+01*  0.12E+00*  0.91E-01*   0.97E+00   0.16E+01   Damp      0.
        #     4    -37.31610460    -50.54103969     13.22493509 -0.23E+00*  0.11E+00*  0.96E-01*   0.72E+00   0.27E+01   Damp      0.
        #     5    -37.33596239    -49.47021484     12.13425245 -0.20E-01*  0.59E-01*  0.59E-01*   0.37E+00   0.16E+01   Damp      0.
        # ...
        #           Convergence after 26 Macro Iterations
        # --
        if line[46:91] == 'iterations: Energy and convergence statistics':

            self.skip_line(inputfile, 'blank')

            while line.split() != ['Energy', 'Energy', 'Energy', 'Change', 'Delta', 'Norm', 'in', 'Sec.']:
                line = next(inputfile)

            iteration_regex = (r"^([0-9]+)"                                  # Iter
                               r"( [ \-0-9]*\.[0-9]{6,9})"                   # Tot. SCF Energy
                               r"( [ \-0-9]*\.[0-9]{6,9})"                   # One-electron Energy
                               r"( [ \-0-9]*\.[0-9]{6,9})"                   # Two-electron Energy
                               r"( [ \-0-9]*\.[0-9]{2}E[\-\+][0-9]{2}\*?)"   # Energy Change
                               r"( [ \-0-9]*\.[0-9]{2}E[\-\+][0-9]{2}\*?)"   # Max Dij or Delta Norm
                               r"( [ \-0-9]*\.[0-9]{2}E[\-\+][0-9]{2}\*?)"   # Max Fij
                               r"( [ \-0-9]*\.[0-9]{2}E[\-\+][0-9]{2}\*?)"   # DNorm
                               r"( [ \-0-9]*\.[0-9]{2}E[\-\+][0-9]{2}\*?)"   # TNorm
                               r"( [ A-Za-z0-9]*)"                           # AccCon
                               r"( [ \.0-9]*)$")                             # Time in Sec.

            scfvalues = []
            line = next(inputfile)
            while not line.strip().startswith("Convergence"):

                match = re.match(iteration_regex, line.strip())
                if match:
                    groups = match.groups()
                    cols = [g.strip() for g in match.groups()]
                    cols = [c.replace('*', '') for c in cols]

                    energy = float(cols[4])
                    density = float(cols[5])
                    fock = float(cols[6])
                    dnorm = float(cols[7])
                    scfvalues.append([energy, density, fock, dnorm])

                if line.strip() == "--":
                    self.logger.warning('File terminated before end of last SCF!')
                    break

                line = next(inputfile)

            self.append_attribute('scfvalues', scfvalues)

        #  Harmonic frequencies in cm-1
        #
        #  IR Intensities in km/mol
        #
        #                         1         2         3         4         5         6
        #
        #      Frequency:       i60.14    i57.39    128.18    210.06    298.24    309.65
        #
        #      Intensity:    3.177E-03 2.129E-06 4.767E-01 2.056E-01 6.983E-07 1.753E-07
        #      Red. mass:      2.42030   2.34024   2.68044   3.66414   2.61721   3.34904
        #
        #      C1         x   -0.00000   0.00000   0.00000  -0.05921   0.00000  -0.06807
        #      C1         y    0.00001  -0.00001  -0.00001   0.00889   0.00001  -0.02479
        #      C1         z   -0.03190   0.04096  -0.03872   0.00001  -0.12398  -0.00002
        #      C2         x   -0.00000   0.00001   0.00000  -0.06504   0.00000  -0.03487
        #      C2         y    0.00000  -0.00000  -0.00000   0.01045   0.00001  -0.05659
        #      C2         z   -0.03703  -0.03449  -0.07269   0.00000  -0.07416  -0.00001
        #      C3         x   -0.00000   0.00001   0.00000  -0.06409  -0.00001   0.05110
        #      C3         y   -0.00000   0.00001   0.00000   0.00152   0.00000  -0.03263
        #      C3         z   -0.03808  -0.08037  -0.07267  -0.00001   0.07305   0.00000
        # ...
        #      H20        y    0.00245  -0.00394   0.03215   0.03444  -0.10424  -0.10517
        #      H20        z    0.00002  -0.00001   0.00000  -0.00000  -0.00000   0.00000
        #
        #
        #
        # ++ Thermochemistry
        if line[1:29] == 'Harmonic frequencies in cm-1':

            self.skip_line(inputfile, 'blank')
            line = next(inputfile)

            while 'Thermochemistry' not in line:

                if 'Frequency:' in line:
                    if not hasattr(self, 'vibfreqs'):
                        self.vibfreqs = []
                    vibfreqs = [float(i.replace('i', '-')) for i in line.split()[1:]]
                    self.vibfreqs.extend(vibfreqs)

                if 'Intensity:' in line:
                    if not hasattr(self, 'vibirs'):
                        self.vibirs = []
                    vibirs = map(float, line.split()[1:])
                    self.vibirs.extend(vibirs)

                if 'Red.' in line:
                    if not hasattr(self, 'vibrmasses'):
                        self.vibrmasses = []
                    vibrmasses = map(float, line.split()[2:])
                    self.vibrmasses.extend(vibrmasses)

                    self.skip_line(inputfile, 'blank')
                    line = next(inputfile)
                    if not hasattr(self, 'vibdisps'):
                        self.vibdisps = []
                    disps = []
                    for n in range(3*self.natom):
                        numbers = [float(s) for s in line[17:].split()]
                        # The atomindex should start at 0 instead of 1.
                        atomindex = int(re.search(r'\d+$', line.split()[0]).group()) - 1
                        numbermodes = len(numbers)
                        if len(disps) == 0:
                            # Appends empty array of the following
                            # dimensions (numbermodes, natom, 0) to disps.
                            for mode in range(numbermodes):
                                disps.append([[] for x in range(0, self.natom)])
                        for mode in range(numbermodes):
                            disps[mode][atomindex].append(numbers[mode])
                        line = next(inputfile)
                    self.vibdisps.extend(disps)

                line = next(inputfile)

        ## Parsing thermochemistry attributes here
        #  ++ Thermochemistry
        #
        #   *********************
        #   *                   *
        #   *  THERMOCHEMISTRY  *
        #   *                   *
        #   *********************
        #
        #   Mass-centered Coordinates (Angstrom):
        #   ***********************************************************
        # ...
        #   *****************************************************
        #   Temperature =     0.00 Kelvin, Pressure =   1.00 atm
        #   -----------------------------------------------------
        #   Molecular Partition Function and Molar Entropy:
        #                          q/V (M**-3)    S(kcal/mol*K)
        #   Electronic            0.100000D+01        0.000
        #   Translational         0.100000D+01        0.000
        #   Rotational            0.100000D+01        2.981
        #   Vibrational           0.100000D+01        0.000
        #   TOTAL                 0.100000D+01        2.981
        #
        #   Thermal contributions to INTERNAL ENERGY:
        #   Electronic           0.000 kcal/mol      0.000000 au.
        #   Translational        0.000 kcal/mol      0.000000 au.
        #   Rotational           0.000 kcal/mol      0.000000 au.
        #   Vibrational        111.885 kcal/mol      0.178300 au.
        #   TOTAL              111.885 kcal/mol      0.178300 au.
        #
        #   Thermal contributions to
        #   ENTHALPY           111.885 kcal/mol      0.178300 au.
        #   GIBBS FREE ENERGY  111.885 kcal/mol      0.178300 au.
        #
        #   Sum of energy and thermal contributions
        #   INTERNAL ENERGY                       -382.121931 au.
        #   ENTHALPY                              -382.121931 au.
        #   GIBBS FREE ENERGY                     -382.121931 au.
        #   -----------------------------------------------------
        # ...
        #   ENTHALPY                              -382.102619 au.
        #   GIBBS FREE ENERGY                     -382.179819 au.
        #   -----------------------------------------------------
        #  --
        #
        #  ++    Isotopic shifts:
        if line[4:19] == 'THERMOCHEMISTRY':

            while "ZPVE" not in line:
                line = next(inputfile)
            self.set_attribute("zpve", float(line.split()[3]))

            temperature_values = []
            pressure_values = []
            entropy_values = []
            internal_energy_values = []
            enthalpy_values = []
            free_energy_values = []

            while 'Isotopic' not in line:

                if line[1:12] == 'Temperature':
                    temperature_values.append(float(line.split()[2]))
                    pressure_values.append(float(line.split()[6]))

                if line[1:48] == 'Molecular Partition Function and Molar Entropy:':
                    while 'TOTAL' not in line:
                        line = next(inputfile)
                    # Molcas reports entropy values in kcal/mol*K but actually writes them in cal/mol*K
                    entropy_values.append(utils.convertor(float(line.split()[2]), 'kcal/mol', 'hartree') / 1000)

                if line[1:40] == 'Sum of energy and thermal contributions':
                    internal_energy_values.append(float(next(inputfile).split()[2]))
                    enthalpy_values.append(float(next(inputfile).split()[1]))
                    free_energy_values.append(float(next(inputfile).split()[3]))

                line = next(inputfile)
            # When calculations for more than one temperature value are
            # performed, the values corresponding to room temperature (298.15 K)
            # are returned and if no calculations are performed for 298.15 K, then
            # the values corresponding last temperature value are returned.
            index = -1
            if 298.15 in temperature_values:
                index = temperature_values.index(298.15)

            self.set_attribute('temperature', temperature_values[index])
            if len(temperature_values) > 1:
                self.logger.warning('More than 1 values of temperature found')

            self.set_attribute('pressure', pressure_values[index])
            if len(pressure_values) > 1:
                self.logger.warning('More than 1 values of pressure found')

            self.set_attribute('entropy', entropy_values[index])
            if len(entropy_values) > 1:
                self.logger.warning('More than 1 values of entropy found')

            self.set_attribute('enthalpy', enthalpy_values[index])
            if len(enthalpy_values) > 1:
                self.logger.warning('More than 1 values of enthalpy found')

            self.set_attribute('freeenergy', free_energy_values[index])
            if len(free_energy_values) > 1:
                self.logger.warning('More than 1 values of freeenergy found')

        ## Parsing Geometrical Optimization attributes in this section.
        #  ++       Slapaf input parameters:
        #  ------------------------
        #
        # Max iterations:                            2000
        # Convergence test a la Schlegel.
        # Convergence criterion on gradient/para.<=: 0.3E-03
        # Convergence criterion on step/parameter<=: 0.3E-03
        # Convergence criterion on energy change <=: 0.0E+00
        # Max change of an internal coordinate:     0.30E+00
        # ...
        # ...
        #  **********************************************************************************************************************
        #  *                                    Energy Statistics for Geometry Optimization                                     *
        #  **********************************************************************************************************************
        #                          Energy     Grad      Grad              Step                 Estimated   Geom       Hessian
        #  Iter      Energy       Change     Norm      Max    Element    Max     Element     Final Energy Update Update   Index
        #    1   -382.30023222  0.00000000 0.107221  0.039531 nrc047   0.085726  nrc047     -382.30533799 RS-RFO  None      0
        #    2   -382.30702964 -0.00679742 0.043573  0.014908 nrc001   0.068195  nrc001     -382.30871333 RS-RFO  BFGS      0
        #    3   -382.30805348 -0.00102384 0.014883  0.005458 nrc010  -0.020973  nrc001     -382.30822089 RS-RFO  BFGS      0
        # ...
        # ...
        #   18   -382.30823419 -0.00000136 0.001032  0.000100 nrc053   0.012319  nrc053     -382.30823452 RS-RFO  BFGS      0
        #   19   -382.30823198  0.00000221 0.001051 -0.000092 nrc054   0.066565  nrc053     -382.30823822 RS-RFO  BFGS      0
        #   20   -382.30820252  0.00002946 0.001132 -0.000167 nrc021  -0.064003  nrc053     -382.30823244 RS-RFO  BFGS      0
        #
        #         +----------------------------------+----------------------------------+
        #         +    Cartesian Displacements       +    Gradient in internals         +
        #         +  Value      Threshold Converged? +  Value      Threshold Converged? +
        #   +-----+----------------------------------+----------------------------------+
        #   + RMS + 5.7330E-02  1.2000E-03     No    + 1.6508E-04  3.0000E-04     Yes   +
        #   +-----+----------------------------------+----------------------------------+
        #   + Max + 1.2039E-01  1.8000E-03     No    + 1.6711E-04  4.5000E-04     Yes   +
        #   +-----+----------------------------------+----------------------------------+
        if 'Convergence criterion on energy change' in line:
            self.energy_threshold = float(line.split()[6])
            # If energy change threshold equals zero,
            # then energy change is not a criteria for convergence.
            if self.energy_threshold == 0:
                self.energy_threshold = numpy.inf

        if 'Energy Statistics for Geometry Optimization' in line:
            if not hasattr(self, 'geovalues'):
                self.geovalues = []

            self.skip_lines(inputfile, ['stars', 'header'])
            line = next(inputfile)
            assert 'Iter      Energy       Change     Norm' in line
            # A variable keeping track of ongoing iteration.
            iter_number = len(self.geovalues) + 1
            # Iterate till blank line.
            while line.split() != []:
                for i in range(iter_number):
                    line = next(inputfile)
                self.geovalues.append([float(line.split()[2])])
                line = next(inputfile)
            # Along with energy change, RMS and Max values of change in
            # Cartesian Diaplacement and Gradients are used as optimization
            # criteria.
            self.skip_lines(inputfile, ['border', 'header', 'header', 'border'])
            line = next(inputfile)
            assert '+ RMS +' in line
            line_rms = line.split()
            line = next(inputfile)
            line_max = next(inputfile).split()
            if not hasattr(self, 'geotargets'):
                # The attribute geotargets is an array consisting of the following
                # values: [Energy threshold, Max Gradient threshold, RMS Gradient threshold, \
                #          Max Displacements threshold, RMS Displacements threshold].
                max_gradient_threshold = float(line_max[8])
                rms_gradient_threshold = float(line_rms[8])
                max_displacement_threshold = float(line_max[4])
                rms_displacement_threshold = float(line_rms[4])
                self.geotargets = [self.energy_threshold, max_gradient_threshold, rms_gradient_threshold, max_displacement_threshold, rms_displacement_threshold]

            max_gradient_change = float(line_max[7])
            rms_gradient_change = float(line_rms[7])
            max_displacement_change = float(line_max[3])
            rms_displacement_change = float(line_rms[3])
            self.geovalues[iter_number - 1].extend([max_gradient_change, rms_gradient_change, max_displacement_change, rms_displacement_change])

        #   *********************************************************
        #   * Nuclear coordinates for the next iteration / Angstrom *
        #   *********************************************************
        #    ATOM              X               Y               Z
        #    C1               0.235560       -1.415847        0.012012
        #    C2               1.313797       -0.488199        0.015149
        #    C3               1.087050        0.895510        0.014200
        # ...
        # ...
        #    H19             -0.021327       -4.934915       -0.029355
        #    H20             -1.432030       -3.721047       -0.039835
        #
        #  --
        if 'Nuclear coordinates for the next iteration / Angstrom' in line:
            self.skip_lines(inputfile, ['s', 'header'])
            line = next(inputfile)

            atomcoords = []
            while line.split() != []:
                atomcoords.append([float(c) for c in line.split()[1:]])
                line = next(inputfile)

            if len(atomcoords) == self.natom:
                self.atomcoords.append(atomcoords)
            else:
                self.logger.warning(
                        f"Parsed coordinates not consistent with previous, skipping. This could be due to symmetry being turned on during the job. Length was {len(self.atomcoords[-1])}, now found {len(atomcoords)}. New coordinates: {str(atomcoords)}")

        #  **********************************************************************************************************************
        #  *                                    Energy Statistics for Geometry Optimization                                     *
        #  **********************************************************************************************************************
        #                         Energy     Grad      Grad              Step                 Estimated   Geom       Hessian
        #  Iter      Energy       Change     Norm      Max    Element    Max     Element     Final Energy Update Update   Index
        #    1   -382.30023222  0.00000000 0.107221  0.039531 nrc047   0.085726  nrc047     -382.30533799 RS-RFO  None      0
        # ...
        # ...
        #   23   -382.30823115 -0.00000089 0.001030  0.000088 nrc053   0.000955  nrc053     -382.30823118 RS-RFO  BFGS      0
        #
        #         +----------------------------------+----------------------------------+
        #         +    Cartesian Displacements       +    Gradient in internals         +
        #         +  Value      Threshold Converged? +  Value      Threshold Converged? +
        #   +-----+----------------------------------+----------------------------------+
        #   + RMS + 7.2395E-04  1.2000E-03     Yes   + 2.7516E-04  3.0000E-04     Yes   +
        #   +-----+----------------------------------+----------------------------------+
        #   + Max + 1.6918E-03  1.8000E-03     Yes   + 8.7768E-05  4.5000E-04     Yes   +
        #   +-----+----------------------------------+----------------------------------+
        #
        #   Geometry is converged in  23 iterations to a Minimum Structure
        if 'Geometry is converged' in line:
            if not hasattr(self, 'optdone'):
                self.optdone = []
            self.optdone.append(len(self.atomcoords))

        #   *********************************************************
        #   * Nuclear coordinates of the final structure / Angstrom *
        #   *********************************************************
        #    ATOM              X               Y               Z
        #    C1               0.235547       -1.415838        0.012193
        #    C2               1.313784       -0.488201        0.015297
        #    C3               1.087036        0.895508        0.014333
        # ...
        # ...
        #    H19             -0.021315       -4.934913       -0.029666
        #    H20             -1.431994       -3.721026       -0.041078
        if 'Nuclear coordinates of the final structure / Angstrom' in line:
            self.skip_lines(inputfile, ['s', 'header'])
            line = next(inputfile)

            atomcoords = []

            while line.split() != []:
                atomcoords.append([float(c) for c in line.split()[1:]])
                line = next(inputfile)

            if len(atomcoords) == self.natom:
                self.atomcoords.append(atomcoords)
            else:
                self.logger.error(
                        f'Number of atoms ({len(atomcoords)}) in parsed atom coordinates is smaller than previously ({int(self.natom)}), possibly due to symmetry. Ignoring these coordinates.')

        ## Parsing Molecular Gradients attributes in this section.
        # ()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
        #
        #                                               &ALASKA
        #
        #                                    only a single process is used
        #                        available to each process: 2.0 GB of memory, 1 thread
        # ()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
        # ...
        # ...
        #  **************************************************
        #  *                                                *
        #  *              Molecular gradients               *
        #  *                                                *
        #  **************************************************
        #
        #   Irreducible representation: a
        #  ---------------------------------------------------------
        #                      X             Y             Z
        #  ---------------------------------------------------------
        #   C1               -0.00009983   -0.00003043    0.00001004
        #   ...
        #   H20              -0.00027629    0.00010546    0.00003317
        #  ---------------------------------------------------
        # WARNING: "Molecular gradients, after ESPF" is found for ESPF QM/MM calculations
        if "Molecular gradients " in line:

            if not hasattr(self, "grads"):
                self.grads = []

            self.skip_lines(inputfile, ['stars', 'stars', 'blank', 'header',
                                        'dashes', 'header', 'dashes'])

            grads = []
            line = next(inputfile)
            while len(line.split()) == 4:
                tmpgrads = list(map(float, line.split()[1:]))
                grads.append(tmpgrads)
                line = next(inputfile)

            self.append_attribute('grads', grads)

        # This code here works, but QM/MM gradients are printed after QM ones.
        # Maybe another attribute is needed to store them to have both.
        if "Molecular gradients, after ESPF" in line:

            self.skip_lines(inputfile, ['stars', 'stars', 'blank', 'header',
                                        'dashes', 'header', 'dashes'])

            grads = []
            line = next(inputfile)
            while len(line.split()) == 4:
                tmpgrads = list(map(float, line.split()[1:]))
                grads.append(tmpgrads)
                line = next(inputfile)

            self.grads[-1] = grads

        ## This section is present when executing &SCF.
        # ++    Orbital specifications:
        #       -----------------------

        #       Symmetry species                           1   2   3   4
        #                                                 a1  b1  a2  b2
        #       Frozen orbitals                            0   0   0   0
        #       Inactive orbitals                          6   5   0   0
        #       Active orbitals                            0   0   4   6
        #       Secondary orbitals                        45  39  10  12
        #       Deleted orbitals                           0   0   0   0
        #       Number of basis functions                 51  44  14  18
        # --

        if not hasattr(self, 'symm_info'):
            self.symm_info = {}

        if line[:29] == '++    Orbital specifications:':

            self.skip_lines(inputfile, ['dashes', 'blank'])
            line = next(inputfile)

            symmetry_count = 1
            while not line.startswith('--'):
                if line.strip().startswith('Symmetry species'):
                    symmetry_count = int(line.split()[-1])
                    self.symm_info["symmetry_count"] = int(symmetry_count)
                    line = next(inputfile)
                    irreps = []
                    for i in range(symmetry_count):
                        irreps.append(line.split()[i])
                    self.symm_info["irreps"] = irreps

                if line.strip().startswith('Frozen orbitals'):
                    frozen_orbitals = []
                    for i in range(2,symmetry_count+2):
                        frozen_orbitals.append(line.split()[i])
                    self.symm_info["frozen_orbitals"] = list(map(int, frozen_orbitals))

                if line.strip().startswith('Inactive orbitals'):
                    inactive_orbitals = []
                    for i in range(2,symmetry_count+2):
                        inactive_orbitals.append(line.split()[i])
                    self.symm_info["inactive_orbitals"] = list(map(int, inactive_orbitals))

                if line.strip().startswith('Active orbitals'):
                    active_orbitals = []
                    for i in range(2,symmetry_count+2):
                        active_orbitals.append(line.split()[i])
                    self.symm_info["active_orbitals"] = list(map(int, active_orbitals))

                if line.strip().startswith('Total number of orbitals'):
                    nmos = line.split()[-symmetry_count:]
                    self.set_attribute('nmo', sum(map(int, nmos)))
                    orbitals_per_irrep = []
                    for i in range(4,symmetry_count+4):
                        orbitals_per_irrep.append(line.split()[i])
                    self.symm_info["orbitals_per_irrep"] = list(map(int, orbitals_per_irrep))

                if line.strip().startswith('Number of basis functions'):
                    nbasis = line.split()[-symmetry_count:]
                    self.set_attribute('nbasis', sum(map(int, nbasis)))
                    basis_per_irrep = []
                    for i in range(4,symmetry_count+4):
                        basis_per_irrep.append(line.split()[i])
                    self.symm_info["basis_per_irrep"] = list(map(int, basis_per_irrep))
                line = next(inputfile)
            self.set_attribute('symm_info', self.symm_info)

##      Some more information about the symmetry and AOs:
        if 'Nuclear Potential Energy' in line:

            self.skip_lines(inputfile, ['d', 'b'])
            line = next(inputfile)
            line = next(inputfile)

            if line.strip().startswith('Symmetry species'):
                tokens = line.split()
                symmetry_species = tokens[2:]
                symmetry_count = len(symmetry_species)
                line = next(inputfile)

            if line.strip().startswith('Basis functions'):
                basis_per_irrep = []
                for i in range(2,symmetry_count+2):
                    basis_per_irrep.append(line.split()[i])
                basis_per_irrep = list(map(int, basis_per_irrep))
                line = next(inputfile)

        line_counter = 0; species_counter = 0
        if '++    SO/AO info:' in line:

            self.skip_lines(inputfile, ['d', 'b','s'])
            line = next(inputfile)

            num_irrep = 8
            symm_adapted_basis = {}
            list_irrep = []
            basis_function_irrep = []
            self.mosyms = []

            self.skip_lines(inputfile, ['d', 'b', 's'])
            line = next(inputfile)
            for irrep in range(num_irrep):
                tokens = line.split()
                while tokens and tokens[0] != '--':
                    tokens = line.split()

                    if line.strip().startswith('Irreducible representation'):
                        species_counter += 1
                        tokens = line.split()
                        list_irrep.append(tokens[3])
                        line = next(inputfile)

                    if line.strip().startswith('Basis function(s) of irrep'):
                        tokens = line.split()
                        basis_function_irrep.append(tokens[4:])
                        line = next(inputfile)

                    if line.strip().startswith('Basis Label'):
                        line = next(inputfile)

                    tokens = line.split()
                    if len(tokens) != 0 and tokens[0] != '--':
                        self.mosyms.append(list_irrep[species_counter-1])
                        line_counter += 1
                        symm_adapted_basis.setdefault("iter", []).append(tokens[0])
                        symm_adapted_basis.setdefault("atom", []).append(tokens[1])
                        symm_adapted_basis.setdefault("type", []).append(tokens[2])
                        symm_adapted_basis.setdefault("center", []).append([tokens[i] for i in range(3, len(tokens), 2)])
                        symm_adapted_basis.setdefault("phase", []).append([tokens[i] for i in range(4, len(tokens), 2)])
                        # symm_adapted_basis.append("{atom}_{orbital}".format(atom=tokens[1], orbital=tokens[2]))
                    line = next(inputfile)
            self.set_attribute('mosyms', self.mosyms)
            print ("mo symms ", self.mosyms)


            line=next(inputfile)
            print ("basis_function_irrep ", basis_function_irrep)
            print ("symmetry adapted basis data ", symm_adapted_basis)         # uncomment later
            print ("equivalent centers")
            # Finding the equivalent centers needed for mo coefficient reconstruction
            equiv_center = []; sublist = []
            for l in symm_adapted_basis["center"]:
                # if atom is not equivalent (only one pair of center-phase)
                if len(l) == 1:
                    [equiv_center.append(l[0]) if l[0] not in equiv_center else None]
                else:
                    # if atom is equivalent (multiple pairs of center-phase)
                    [sublist.append(l) if l not in sublist else None]
            equiv_center.append(sublist)
            print ("equivalent list ", equiv_center)

            # Finding and counting types of basis functions
            types_of_basis_functions = Counter([item[1:2] for item in symm_adapted_basis["type"]])
            print ("types of basis functions", types_of_basis_functions)

            # testorb=[]; m = 0
            # orbfile = open('cn5-c2v-BFD-Da-cas610_sym2.ScfOrb','r')
            # line2 = orbfile.readline()
            # new_coeffs=[]; m = 0
            # while ( (len(line2.split())<3) or (line2.split()[1]!='ORBITAL') or (int(line2.split()[2])!=1+1) ):
            #     line2=orbfile.readline();  m=m+1
            # line2=orbfile.readline()
            # print ("line2 ", line2)
            num_irrep =  len(list_irrep)
            if (num_irrep == 1):
                self.mosyms = ["a"]
                self.set_attribute('mosyms', self.mosyms)
            irreps = OrderedDict(Counter(self.mosyms).items())
            # print (irreps, len(irreps))
            for _, bas in irreps.items():
                print (bas)
            ## AO information ends here

        ###
        #        All orbitals with orbital energies smaller than  E(LUMO)+0.5 are printed
        #
        #  ++    Molecular orbitals:
        #        -------------------
        #
        #        Title: RKS-DFT orbitals
        #
        #        Molecular orbitals for symmetry species 1: a
        #
        #            Orbital        1         2         3         4         5         6         7         8         9        10
        #            Energy      -10.0179  -10.0179  -10.0075  -10.0075  -10.0066  -10.0066  -10.0056  -10.0055   -9.9919   -9.9919
        #            Occ. No.      2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    2.0000    2.0000
        #
        #          1 C1    1s     -0.6990    0.6989    0.0342    0.0346    0.0264   -0.0145   -0.0124   -0.0275   -0.0004   -0.0004
        #          2 C1    2s     -0.0319    0.0317   -0.0034   -0.0033   -0.0078    0.0034    0.0041    0.0073   -0.0002   -0.0002
        # ...
        # ...
        #         58 H18   1s      0.2678
        #         59 H19   1s     -0.2473
        #         60 H20   1s      0.1835
        #  --
        if '++    Molecular orbitals:' in line:

            self.skip_lines(inputfile, ['d', 'b'])
            line = next(inputfile)

            # We don't currently support parsing natural orbitals or active space orbitals.
            # if 'Natural orbitals' not in line and "Pseudonatural" not in line and 'Quasi-canonical orbitals' not in line:
            if 'Natural orbitals' not in line and "Pseudonatural" not in line:
                self.skip_line(inputfile, 'b')

                homos = 0
                num_irrep =  self.symm_info["symmetry_count"]
                basis_per_irrep =  self.symm_info["basis_per_irrep"]
                mocoeffs_per_irrep = [[] for i in range(num_irrep)]
                moenergies_per_irrep = [[] for i in range(num_irrep)]
                # aonames_per_irrep = [[] for i in range(num_irrep)]
                orbital_index_per_irrep = [[] for i in range(num_irrep)]
                mocoeffs = [[]]

                if (num_irrep == 1):
                    self.mosyms = [item for item in ["a"] for i in range(basis_per_irrep[0])]
                    self.set_attribute('mosyms', self.mosyms)
                irreps = OrderedDict(Counter(self.mosyms).items())


                line = next(inputfile)
                tokens = line.split()
                irrep = 0
                while not line.strip() == "--":
                    line = next(inputfile)
                    tokens = line.split()
                    if line.strip().startswith('Molecular orbitals for symmetry species'):
                        irrep += 1
                        line = next(inputfile)

                    if line.strip().startswith('Orbital'):
                        orbital_index_per_irrep[irrep].extend(line.split()[1:])
                        line = next(inputfile)

                    if line.strip().startswith('Energy'):
                        energies = [utils.convertor(float(x), 'hartree', 'eV') for x in line.split()[1:]]
                        moenergies_per_irrep[irrep].extend(energies)
                        line = next(inputfile)

                    if 'Occ. No.' in line:
                        for i in line.split()[2:]:
                            if float(i) != 0:
                                homos += 1
                        self.skip_line(inputfile, 'b')
                        line = next(inputfile)

                    tokens = line.split()
                    # aonames_per_irrep[irrep].append("{atom}_{orbital}".format(atom=tokens[1], orbital=tokens[2]))

                    info = tokens[3:]
                    j = 0
                    mocoeffs_per_irrep[irrep].append([ float(x) for x in tokens[3:]])
                    print ( "some float ", [float(x) for x in tokens[3:]])
                    # self.set_attribute('aonames', aonames)
                orbitals_per_irrep = [len(i) for i in orbital_index_per_irrep]
                print("orbital index per irrep   ", orbital_index_per_irrep)
                print("orbitals per irrep   ", orbitals_per_irrep)



        #         ## Get the intermediate list of indices needed for reshaping the coefficients array
        #         intermediate = []
        #         for ind, i in enumerate(orbital_index_per_irrep):
        #             print (ind, i[9::10], int(i[-1])%10 )
        #             intermediate.append( [ j for j in i[9::10] ])
        #             intermediate[ind].append(i[-1])
        #             intermediate[ind] = uniquelist(intermediate[ind])
        #             intermediate[ind] = list(map(int, intermediate[ind]))

        #         # print ("intermediate ", intermediate)
        #         ## Convert the parsed mocoeffs data into an ordered arrays.

        #         ### Step 1: Remove empty lists from the list of lists
        #         for irrep in range(num_irrep):
        #             mocoeffs_per_irrep[irrep] = [ele for ele in mocoeffs_per_irrep[irrep] if ele != []]


        #         ### Step 2: Get the basis set information
        #         irreps = OrderedDict(Counter(self.mosyms).items())
        #         basis_per_irrep = []
        #         for _, bas in irreps.items():
        #             basis_per_irrep.append(bas)

        #         # print ("basis_per_irrep ", basis_per_irrep)

        #         # Number of splitted blocks of 10 orbitals per irrep
        #         blocks = []
        #         for blk in intermediate:
        #             blocks.append(len(blk))

        #         # print ("blocks ", blocks)

        #         # print ("mocoeffs_per_irrep ", mocoeffs_per_irrep)


        #         npmocoeff = numpy.zeros([num_irrep, max(basis_per_irrep), max(orbitals_per_irrep)], dtype=float)

        #         numpy.set_printoptions(threshold=sys.maxsize)

        #         if num_irrep == 1:
        #             for irrep in range(num_irrep):                          # Run loop over all irreps
        #                 for ind in range(blocks[irrep]):                    # Run loop over all blocks of 10 orbitals
        #                     blockoften =  intermediate[0][ind] - 10*ind     # Get the block number
        #                     for i in range(basis_per_irrep[irrep]):         # Run loop over all basis functions
        #                         for j in range(blockoften):                 # Run loop over all orbitals in block
        #                             npmocoeff[irrep, i, 10*ind + j] = mocoeffs_per_irrep[0][basis_per_irrep[irrep]*ind+i][j]
        #         else:
        #             for irrep in range(num_irrep):                          # Run loop over all irreps
        #                 for ind in range(blocks[irrep]):                    # Run loop over all blocks of 10 orbitals
        #                     blockoften =  intermediate[0][ind] - 10*ind     # Get the block number
        #                     for i in range(basis_per_irrep[irrep]):         # Run loop over all basis functions
        #                         for j in range(blockoften):                 # Run loop over all orbitals in block
        #                             npmocoeff[irrep, i, 10*ind + j] = mocoeffs_per_irrep[0][basis_per_irrep[irrep]*ind+i][j]

        #         if num_irrep == 1:
        #             self.append_attribute('mocoeffs', (numpy.transpose(npmocoeff[0])).tolist())
                # print ("self mocoeffs 3rd orbs ", self.mocoeffs[0][2])

                # i = 0; ik=0
                # for irrep in range(num_irrep):
                #     for iB in range(bas(irrep)):
                #         i += 1
                #         for j in range(nbas):
                #             if label[j] == symbol[i] + cNumber[ik]:
                #                 for k in range(8):
                #                     if coord[k,i] == Cent3[j]:
                #                         for ii in range(nbas):
                #                             if weight == 0:
                #                                 V[j,ii] = V[j,ii]+(C2[i,ii]*r_Norm[j])*phase(k,i)/Cent2(i)
                #                             else:
                #                                 V[j,ii] = V[j,ii]+(C2[i,ii]*r_Norm[j])*Phase(k,i)/numpy.sqrt(Cent2[i])




                # for ind, i in enumerate(orbitals_per_irrep):
                    # print ("mocoeffs per irrep 3rd ", i, ind, [mocoeffs_per_irrep[ind][j::i] for j in range(i)])

                        # if len(moenergies_per_irrep[irrep]) != self.symm_info["orbitals_per_irrep"][irrep]:
                        #     moenergies_per_irrep[irrep].extend([numpy.nan for x in range(self.symm_info["orbitals_per_irrep"][irrep] - len(moenergies_per_irrep[irrep]))])

                        # self.append_attribute('moenergies', moenergies_per_irrep[irrep])

                        # if not hasattr(self, 'homos'):
                        #     self.homos = []
                        # self.homos.extend([homos-1])

                        # while len(mocoeffs_per_irrep[irrep]) < self.symm_info["orbitals_per_irrep"][irrep]:
                        #     nan_array = [numpy.nan for i in range(self.symm_info["basis_per_irrep"][irrep])]
                        #     mocoeffs_per_irrep[irrep].append(nan_array)

                        # self.append_attribute('mocoeffs', mocoeffs)
                # print ("aonames ", aonames_per_irrep)         # uncomment later


        # This part is to read the symmetry species and basis functions per symmetry species. It is required by the
        # following block of AO/SO.



        ## Parsing MP energy from the &MBPT2 module.
        #  Conventional algorithm used...
        #
        #         SCF energy                           =      -74.9644564043 a.u.
        #         Second-order correlation energy      =       -0.0364237923 a.u.
        #
        #         Total energy                         =      -75.0008801966 a.u.
        #         Reference weight ( Cref**2 )         =        0.98652
        #
        #  ::    Total MBPT2 energy                              -75.0008801966
        #
        #
        #         Zeroth-order energy (E0)             =      -36.8202538520 a.u.
        #
        #         Shanks-type energy S1(E)             =      -75.0009150108 a.u.
        if 'Total MBPT2 energy' in line:
            mpenergies = []
            mpenergies.append(utils.convertor(utils.float(line.split()[4]), 'hartree', 'eV'))
            if not hasattr(self, 'mpenergies'):
                self.mpenergies = []
            self.mpenergies.append(mpenergies)

        # Parsing data ccenergies from &CCSDT module.
        #  --- Start Module: ccsdt at Thu Jul 26 14:03:23 2018 ---
        #
        #  ()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()()
        #
        #                                                 &CCSDT
        # ...
        # ...
        #          14          -75.01515915      -0.05070274      -0.00000029
        #          15          -75.01515929      -0.05070289      -0.00000014
        #          16          -75.01515936      -0.05070296      -0.00000007
        #       Convergence after                    17  Iterations
        #
        #
        #      Total energy (diff) :     -75.01515936      -0.00000007
        #      Correlation energy  :        -0.0507029554992
        if 'Start Module: ccsdt' in line:
            self.skip_lines(inputfile, ['b', '()', 'b'])
            line = next(inputfile)
            if '&CCSDT' in line:
                while not line.strip().startswith('Total energy (diff)'):
                    line = next(inputfile)

                ccenergies = utils.convertor(utils.float(line.split()[4]), 'hartree', 'eV')
                if not hasattr(self, 'ccenergies'):
                    self.ccenergies= []
                self.ccenergies.append(ccenergies)

        #  ++    Primitive basis info:
        #        ---------------------
        #
        #
        #                      *****************************************************
        #                      ******** Primitive Basis Functions (Valence) ********
        #                      *****************************************************
        #
        #
        #   Basis set:C.AUG-CC-PVQZ.........
        #
        #                    Type
        #                     s
        #             No.      Exponent    Contraction Coefficients
        #             1  0.339800000D+05   0.000091  -0.000019   0.000000   0.000000   0.000000   0.000000
        #             2  0.508900000D+04   0.000704  -0.000151   0.000000   0.000000   0.000000   0.000000
        # ...
        # ...
        #             29  0.424000000D+00   0.000000   1.000000
        #
        #   Number of primitives                                   93
        #   Number of basis functions                              80
        #
        #  --
        if line.startswith('++    Primitive basis info:'):
            self.skip_lines(inputfile, ['d', 'b', 'b', 's', 'header', 's', 'b'])
            line = next(inputfile)
            gbasis_array = []
            while '--' not in line and '****' not in line:
                if 'Basis set:' in line:
                    basis_element_patterns = re.findall(r'Basis set:([A-Za-z]{1,2})\.', line)
                    assert len(basis_element_patterns) == 1
                    basis_element = basis_element_patterns[0].title()
                    gbasis_array.append((basis_element, []))

                if 'Type' in line:
                    line = next(inputfile)
                    shell_type = line.split()[0].upper()

                    self.skip_line(inputfile, 'headers')
                    line = next(inputfile)

                    exponents = []
                    coefficients = []
                    func_array = []
                    while line.split():
                        exponents.append(utils.float(line.split()[1]))
                        coefficients.append([utils.float(i) for i in line.split()[2:]])
                        line = next(inputfile)

                    for i in range(len(coefficients[0])):
                        func_tuple = (shell_type, [])
                        for iexp, exp in enumerate(exponents):
                            coeff = coefficients[iexp][i]
                            if coeff != 0:
                                func_tuple[1].append((exp, coeff))
                        gbasis_array[-1][1].append(func_tuple)

                line = next(inputfile)

            atomsymbols = [self.table.element[atomno] for atomno in self.atomnos]
            self.gbasis = [[] for i in range(self.natom)]
            for element, gbasis in gbasis_array:
                mask = [element == possible_element for possible_element in atomsymbols]
                indices = [i for (i, x) in enumerate(mask) if x]
                for index in indices:
                    self.gbasis[index] = gbasis

        #  ++    Basis set information:
        #        ----------------------
        # ...
        #        Basis set label: MO.ECP.HAY-WADT.5S6P4D.3S3P2D.14E-LANL2DZ.....
        #
        #        Electronic valence basis set:
        #        ------------------
        #        Associated Effective Charge  14.000000 au
        #        Associated Actual Charge     42.000000 au
        #        Nuclear Model: Point charge
        # ...
        #
        #        Effective Core Potential specification:
        #        =======================================
        #
        #         Label   Cartesian Coordinates / Bohr
        #
        #   MO                 0.0006141610       -0.0006141610        0.0979067106
        #  --
        if '++    Basis set information:' in line:
            self.core_array = []
            basis_element = None
            ncore = 0

            while line[:2] != '--':
                if 'Basis set label' in line:
                    try:
                        basis_element = line.split()[3].split('.')[0]
                        basis_element = basis_element[0] + basis_element[1:].lower()
                    except:
                        self.logger.warning('Basis set label is missing!')
                        basis_element = ''
                if 'valence basis set:' in line.lower():
                    self.skip_line(inputfile, 'd')
                    line = next(inputfile)
                    if 'Associated Effective Charge' in line:
                        effective_charge = float(line.split()[3])
                        actual_charge = float(next(inputfile).split()[3])
                        element = self.table.element[int(actual_charge)]
                        ncore = int(actual_charge - effective_charge)
                        if basis_element:
                            assert basis_element == element
                        else:
                            basis_element = element

                if basis_element and ncore:
                    self.core_array.append((basis_element, ncore))
                    basis_element = ''
                    ncore = 0

                line = next(inputfile)

        #                               &RASSCF
        # ++    Wave function specifications:
        #       -----------------------------
        #
        #       Number of closed shell electrons          20
        #       Number of electrons in active shells       6
        #       Max number of holes in RAS1 space          0
        #       Max nr of electrons in RAS3 space          0
        #       Number of inactive orbitals               10
        #       Number of active orbitals                  5
        #       Number of secondary orbitals             320
        #       Spin quantum number                      0.0
        #       State symmetry                             1
        # --

        if '++    Wave function specifications:' in line:
            self.wfn_info = {}

            while line[:2] != '--':
                if 'Number of closed shell electrons' in line:
                    try:
                        self.wfn_info["nelec_closed_shell"] = line.split()[5]
                    except:
                        self.logger.warning('Number of closed shell electrons label is missing!')
                        self.wfn_info["nelec_closed_shell"] = ""
                if 'Number of electrons in active shells' in line:
                    try:
                        self.wfn_info["nelec_active_shell"] = line.split()[6]
                    except:
                        self.logger.warning('Number of active shell electrons label is missing!')
                        self.wfn_info["nelec_active_shell"] = ""
                if 'Number of inactive orbitals' in line:
                    try:
                        self.wfn_info["norb_inactive"] = line.split()[4]
                    except:
                        self.logger.warning('Number of inactive orbitals label is missing!')
                        self.wfn_info["norb_inactive"] = ""
                if 'Number of active orbitals' in line:
                    try:
                        self.wfn_info["norb_active"] = line.split()[4]
                    except:
                        self.logger.warning('Number of active orbitals label is missing!')
                        self.wfn_info["norb_active"] = ""
                if 'Number of secondary orbitals' in line:
                    try:
                        self.wfn_info["norb_secondary"] = line.split()[4]
                    except:
                        self.logger.warning('Number of secondary orbitals label is missing!')
                        self.wfn_info["norb_secondary"] = ""
                if 'State symmetry' in line:
                    try:
                        self.wfn_info["state_symmetry"] = line.split()[2]
                    except:
                        self.logger.warning('State symmetry label is missing!')
                        self.wfn_info["state_symmetry"] = ""

                line = next(inputfile)
            self.set_attribute('wfn_info', self.wfn_info)


        # ++    CI expansion specifications:
        #       ----------------------------
        #
        #       Number of CSFs                            28
        #       Number of determinants                    31
        #       Number of root(s) required                 3
        #       CI roots used                              1     2     3
        #       weights                                0.333 0.333 0.333
        #       highest root included in the CI            3
        #       Root passed to geometry opt.               3
        # --

        if '++    CI expansion specifications:' in line:
            # Read the first instance of the CI expansion specifications line only
            if not hasattr(self, 'ci'):
                self.ci = {}
            while line[:2] != '--':
                if 'Number of CSFs' in line:
                    try:
                        self.ci["Number of CSFs"] = line.split()[3]
                    except:
                        self.logger.warning('Number of CSFs label is missing!')
                        self.ci["Number of CSFs"] = ""
                if 'Number of determinants' in line:
                    try:
                        self.ci["Number of determinants"] = line.split()[3]
                    except:
                        self.logger.warning('Number of determinants label is missing!')
                        self.ci["Number of determinants"] = ""
                if 'Number of root(s) required' in line:
                    try:
                        self.ci["Number of root(s) required"] = line.split()[4]
                    except:
                        self.logger.warning('Number of root(s) required label is missing!')
                        self.ci["Number of root(s) required"] = ""
                if 'Root chosen for geometry opt' in line:
                    try:
                        self.ci["Root chosen for geometry opt"] = line.split()[5]
                    except:
                        self.logger.warning('Root chosen for geometry opt label is missing!')
                        self.ci["Root chosen for geometry opt"] = ""
                if 'CI root' in line:
                    try:
                        self.ci["CI root"] = line.split()[3:]
                    except:
                        self.logger.warning('CI roots used label is missing!')
                        self.ci["CI root"] = []
                if 'weights' in line:
                    try:
                        self.ci["weights"] = line.split()[1:]
                    except:
                        self.logger.warning('weights label is missing!')
                        self.ci["weights"] = []

                line = next(inputfile)


        #       ************************************************************************************************************************
        #                                                       Wave function printout:
        #                        occupation of active orbitals, and spin coupling of open shells (u,d: Spin up or down)
        #       ************************************************************************************************************************
        #
        #       Note: transformation to natural orbitals
        #       has been made, which may change the order of the CSFs.
        #
        #       printout of CI-coefficients larger than  0.00 for root  1
        #       energy=    -188.396249
        #       conf/sym  11111     Coeff  Weight
        #              1  22200  -0.84226 0.70941
        #
        #    + sqrt(     1/1     )  |22200|
        #
        #              2  22ud0   0.41173 0.16952
        #
        #    + sqrt(     1/2     )  |22ab0|
        #    - sqrt(     1/2     )  |22ba0|
        #
        #              3  22020   0.07400 0.00548

        if 'Wave function printout:' in line:
            from fractions import Fraction
            number_of_roots = int(self.ci["Number of root(s) required"])
            number_of_determinants = int(self.ci["Number of determinants"])
            number_of_csfs = int(self.ci["Number of CSFs"])

            ci_energy = numpy.ndarray(shape=(number_of_roots), dtype=float)
            csf_coeff = numpy.ndarray(shape=(number_of_roots, number_of_csfs), dtype=float)
            dets_per_csf = numpy.zeros(shape=(number_of_roots, number_of_csfs), dtype=int)
            csf_occupations = numpy.ndarray(shape=(number_of_roots, number_of_csfs), dtype=object)
            ci_occupations = [numpy.array([],dtype=object) for i in range(number_of_roots)]
            ci_coeff       = [numpy.array([],dtype=object) for i in range(number_of_roots)]


            while line[6:54] != 'Natural orbitals and occupation numbers for root':

                if 'printout of CI-coefficients larger than' in line:
                    try:
                        root_number = int(line.split()[8])
                        print ("root number ", root_number )
                    except:
                        self.logger.warning('CI coefficient root label is missing!')
                        root_number = 0

                if 'energy=' in line:
                    print (line)
                    try:
                        ci_energy[root_number-1] = line.split()[1]
                    except:
                        self.logger.warning('ci energy label is missing!')
                        ci_energy[root_number-1] = 0.0

                if 'conf/sym' in line:
                    line = next(inputfile)

                    icsf_counter = 0; csfmap_counter = 0; det_per_csf_counter = 0
                    while icsf_counter < number_of_csfs:

                        if "sqrt" in line:
                            dets_per_csf[root_number-1,icsf_counter-1] += 1
                            ci_occupations[root_number-1] = numpy.append( ci_occupations[root_number-1], line.split()[4])
                            # Read the phase factor + or -
                            if line.split()[0] == "-":
                                coeff = -numpy.sqrt(float(Fraction(line.split()[2])))
                            else:
                                coeff = numpy.sqrt(float(Fraction(line.split()[2])))

                            ci_coeff[root_number-1]       = numpy.append( ci_coeff[root_number-1], coeff)
                            csfmap_counter += 1
                            line = next(inputfile)
                        elif not line.strip():
                            line = next(inputfile)
                        else:
                            icsf_counter += 1
                            icsf = int(line.split()[0])
                            csf_occupations[root_number-1,icsf-1]  = line.split()[1]
                            csf_coeff[root_number-1,icsf-1] = line.split()[2]
                            line = next(inputfile)
                line = next(inputfile)
            self.ci["CSF_Mappings"] = csfmap_counter
            self.ci["CI_Energy"] = ci_energy

            # Replace the occupation strings with champ formatted numbers
            # print ("ci occupations original ", ci_occupations[0].shape)
            # ci_occupations = numpy.vectorize(utils.molcas_occup_strings_to_numbers)(ci_occupations)
            #
            self.ci["Dets_Per_CSF"] = dets_per_csf
            # self.ci["CI_Occupations"] = ci_occupations
            self.ci["CI_Coefficients"] = ci_coeff
            self.ci["CSF_Occupations"] = csf_occupations
            self.ci["CSF_Coefficients"] = csf_coeff
            self.set_attribute('ci', self.ci)


