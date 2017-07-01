#!/usr/bin/env python

from pyscf import gto
from pyscf import dft


h2o = gto.M(atom='O 0.00000000,  0.000000,  0.119748; H 0.00000000,  0.761561, -0.478993; H 0.00000000, -0.761561, -0.478993', basis='cc-pvtz')

h2o.build()

#h2o.verbose = 5

method = dft.RKS(h2o)

method.grids.prune = dft.gen_grid.treutler_prune
method.grids.atom_grid = {"H": (50, 194), "O": (50, 194),}

method.xc = 'b3lyp'
method.scf()


#print h2o.atom_coords()
#print h2o.nao_nr()
#print method.make_rdm1()
