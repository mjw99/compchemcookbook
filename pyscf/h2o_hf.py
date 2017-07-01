#!/usr/bin/env python

from pyscf import gto
from pyscf import scf


h2o = gto.M(atom='O 0.00000000,  0.000000,  0.119748; H 0.00000000,  0.761561, -0.478993; H 0.00000000, -0.761561, -0.478993', basis='ccpvtz')

h2o.build()

#h2o.verbose = 5

mf = scf.RHF(h2o)
mf.kernel()


