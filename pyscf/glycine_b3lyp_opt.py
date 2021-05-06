#!/usr/bin/env python

from pyscf import gto
from pyscf import lib
from pyscf import dft

from pyscf.geomopt import berny_solver



mol = gto.Mole()

mol.atom = open('glycine.xyz').read()

mol.basis = {"H": '6-31g', "O": '6-31g', "N": '6-31g', "C": '6-31g'}
#mol.basis = 'cc-pvtz'
mol.build()

#h2o.verbose = 5

method = dft.RKS(mol)
method.grids.prune = dft.gen_grid.treutler_prune
method.grids.atom_grid = {"H": (50, 194), "O": (50, 194),}
method.xc = 'b3lyp'


mol1 = berny_solver.optimize(method, verbose=5)
method.scf()
method.analyze()
