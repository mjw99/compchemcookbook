#!/usr/bin/env python

from pyscf import gto, scf, cc
from pyscf.tools import cubegen


h2o = gto.M(atom='	O 0.00000000,  0.000000,  0.119748;\
			H 0.00000000,  0.761561, -0.478993;\
			H 0.00000000, -0.761561, -0.478993',\
			basis='ccpvdz')

h2o.build()

#h2o.verbose = 5

method = scf.RHF(h2o).run()
cc = cc.CCSD(method).run()
et = cc.ccsd_t()
print('CCSD(T) correlation energy', cc.e_corr + et)

dm = cc.make_rdm1()

#cubegen.density(h2o, 'water_den.cube', dm, nx=grid_points, ny=grid_points, nz=grid_points)
#cubegen.density(h2o, 'water_den.cube', dm, step=0.088, pad=2.0)
#cubegen.density(h2o, 'water_den.cube', dm, pad=4.0)
#cubegen.density(h2o, 'water_den.cube', dm, step=0.088, pad=2.0)

# Coarse (3 points/Bohr)
#cubegen.density(h2o, 'water_den.cube', dm, gridspacing=0.1763)
#cubegen.mep(h2o, 'water_pot.cube', dm, gridspacing=0.1763)

# Medium (6 points/Bohr)
#cubegen.density(h2o, 'water_den.cube', dm, gridspacing=0.088)
#cubegen.mep(h2o, 'water_pot.cube', dm, gridspacing=0.088)

#cubegen.isomep(h2o, 'water_den.cube', dm, electronic_iso=0.001, gridspacing=0.088)
#cubegen.isomep(h2o, 'water_den.cube', dm, electronic_iso=0.002, gridspacing=0.088)
#cubegen.isomep(h2o, 'water_den_cc.cube', dm, electronic_iso=0.002, gridspacing=0.088)
#cubegen.isomep(h2o, 'water_den.cube', dm, electronic_iso=0.002, gridspacing=0.0441)

# Fine (12 points/Bohr)
#cubegen.density(h2o, 'water_den.cube', dm, pad=2.0, gridspacing=0.0441)
#cubegen.mep(h2o, 'water_pot.cube', dm, gridspacing=0.0441)



#print h2o.atom_coords()
#print h2o.nao_nr()
#print method.make_rdm1()
