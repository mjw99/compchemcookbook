#!/usr/bin/env python

from pyscf import gto
from pyscf import dft
from pyscf.geomopt import berny_solver

from pyscf.tools import cubegen


h2o = gto.M(atom='	O 0.00000000,  0.000000,  0.119748;\
			H 0.00000000,  0.761561, -0.478993;\
			H 0.00000000, -0.761561, -0.478993',\
			basis='sto-3g')
			#basis='6-31G*')
			#basis='ccpvdz')
			#basis='aug-cc-pvdz')


h2o.build()

h2o.verbose = 5

method = dft.RKS(h2o)
method.grids.prune = dft.gen_grid.treutler_prune
method.grids.atom_grid = {"H": (50, 194), "O": (50, 194),}
method.xc = 'b3lypg'
method.scf()

mol1 = berny_solver.optimize(method, verbose=5)
#mol1 = berny_solver.optimize(method)

print(mol1)



