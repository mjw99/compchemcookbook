#!/usr/bin/env python

from pyscf import gto, scf, lo, tools

mol = gto.M(atom='O 0.00000000,  0.000000,  0.119748; H 0.00000000,  0.761561, -0.478993; H 0.00000000, -0.761561, -0.478993', basis='cc-pvtz')

mol.verbose = 5


#mf = scf.RHF(mol)
#mf.kernel()


mf = scf.RHF(mol).run()
orb = lo.Boys(mol).kernel(mf.mo_coeff[:,:180])
tools.molden.from_mo(mol, 'h2o.molden', orb)

