#!/usr/bin/python


# An example MINDO3 calculation on Water using PyQuante

# Links
# http://en.wikipedia.org/wiki/MINDO



from PyQuante import SCF
from PyQuante.Molecule import Molecule
from PyQuante.Basis.sto3g import basis_data


# Used to prettify matrix printing
import numpy as np
np.set_printoptions(threshold='nan')
np.set_printoptions(linewidth=100)
np.set_printoptions(precision=3)


# Build our molecule
h2o = Molecule('H2O',
                   [(8,  ( 0.00000000,     0.00000000,     0.00000000)),
                    (1,  ( 0.00000000,     1.43042809,    -1.10715266)),
                    (1,  ( 0.00000000,    -1.43042809,    -1.10715266))],
                   units='Angstrom')

# Carry out the calculation
h2o_mindo3 = SCF(h2o,method="MINDO3")
h2o_mindo3.iterate()


# Display the results
print "MINDO3 Energy: =", 
print ""
print h2o_mindo3.energy
print ""
print "MINDO3 Orbital energies:", 
print ""
print h2o_mindo3.orbe
print ""
print "MINDO3 Orbitals:", 
print ""
print "Col = AO coefficient"
print "Row = MO level"
print ""
print h2o_mindo3.orbs
print h2o_mindo3.nel


# Calculate HOMO and LUMO MO index (offset for array indexing starting from zero)
# h2o_mindo3.nel is the number of electrons in our molecule
homo_mo_index =  (h2o_mindo3.nel/2)-1
lumo_mo_index =  (h2o_mindo3.nel/2)+1-1

print "HOMO energy: ", h2o_mindo3.orbe[homo_mo_index]
print "LUMO energy: ", h2o_mindo3.orbe[lumo_mo_index]
print "HOMO - LUMO energy diff: ", h2o_mindo3.orbe[homo_mo_index] - h2o_mindo3.orbe[lumo_mo_index]






