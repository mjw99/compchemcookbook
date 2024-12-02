#!/usr/bin/env python

from pyscf import gto
from pyscf import dft
from pyscf.geomopt.geometric_solver import optimize

h2o = gto.M(atom='O 0.00000000  0.000000  0.119748;\
			H 0.00000000  0.761561 -0.478993;\
			H 0.00000000 -0.761561 -0.478993',\
			#basis='sto-3g')
			basis='6-31G*')
			#basis='ccpvdz')
			#basis='aug-cc-pvdz')


h2o.build()

h2o.verbose = 3

method = dft.RKS(h2o)
method.grids.prune = dft.gen_grid.treutler_prune
#method.grids.atom_grid = {"H": (50, 194), "O": (50, 194),}
#method.grids.level = 9
method.xc = 'b3lypg'
method.scf()

# Create callback to run at every geometry opt. iteration
# https://github.com/pyscf/pyscf/issues/798
energies = []
def cb(envs):
  mf = envs["g_scanner"].base
  energies.append(mf.e_tot)

mol1 = optimize(method, callback=cb, convergence_set="GAU_TIGHT")

print(energies[-1])

mol1.tofile("h2o.xyz", format="xyz")



