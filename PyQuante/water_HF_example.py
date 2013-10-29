#!/usr/bin/python


# Mark J. Williamson Sept. 2012
#
# An example HF calculation on Water using PyQuante, to extract
# the resulting MO coefficients from the converged wavefunction.
# This is to aid debugging of an issue with Jon and NWChem


# Useful Refs:
#  General MO theory
#    http://www.chm.bris.ac.uk/pt/harvey/elstruct/hf_method.html
#  Water MOs
#    http://www.lsbu.ac.uk/water/h2oorb.html
#  Expt. for Water
#    http://dx.doi.org/10.1016/j.chemphys.2007.09.030



from PyQuante.Ints import getbasis,getints
from PyQuante.Molecule import Molecule
from PyQuante.hartree_fock import *
from PyQuante.Basis.sto3g import basis_data
#from PyQuante.Basis.p631ss import basis_data
from PyQuante.Constants import hartree2ev

#import logging
#LOG = "log"
#logging.basicConfig(filename=LOG, level=logging.DEBUG)

# Printing formatting for matrix output
import numpy as np
np.set_printoptions(threshold='nan')
np.set_printoptions(linewidth=100)
np.set_printoptions(precision=3)








h2o = Molecule('H2O',
                   [('O',  ( 0.00000000,     0.000000,   0.119748)),
                    ('H',  ( 0.00000000,     0.761561,  -0.478993)),
                    ('H',  ( 0.00000000,    -0.761561,  -0.478993))],
                   units='Angstrom')


atoms = h2o
bfs = getbasis(atoms,basis_data)
S,h,Ints = getints(bfs,atoms)
en,orbe,orbs = rhf(atoms,integrals=(S,h,Ints))



print "Orbital Energy:  (eV)"
print "1a_1, 2a_1, 1b_2, 3a_1, 1b_1, || 4a_1, 2b_2, 3b_2"
print  orbe / hartree2ev
print ""

print ""
print "HF Energy  ",en / hartree2ev, " (eV)"

print ""
print "Overlap matrix: S"
print "(to what degree each AO overlaps with another AO)"
print "Overlap matrix: S"
print np.array(S)

print ""
print "Eigenvector matrix: "
print "Col = AO coefficient"
print "Row = MO level"
print np.array(orbs)


print ""
print "Density Matrix" 
# How many closed MOs and open MOs
nclosed,nopen =  atoms.get_closedopen()
nocc = nclosed
# Form a density matrix C*Ct given eigenvectors C,nstart,nstop]
D = mkdens(orbs,0,nocc)
print D


print ""
print "Number of electrons"
print h2o.get_nel()
# Calculate HOMO and LUMO MO index (offset for array indexing starting from zero)
# h2o.get_nel is the number of electrons in our molecule
homo_mo_index =  (h2o.get_nel()/2)-1
lumo_mo_index =  (h2o.get_nel()/2)+1-1

homo_mo_index = 4
lumo_mo_index = 5


print "HOMO (4a_1) energy: ", orbe[homo_mo_index] / hartree2ev, "eV"
print "LUMO (1b_1) energy: ", orbe[lumo_mo_index] / hartree2ev, "eV"
print "HOMO - LUMO energy diff: ", (orbe[homo_mo_index] - orbe[lumo_mo_index]) / hartree2ev, "eV"





