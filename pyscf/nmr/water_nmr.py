import numpy
from pyscf import gto
from pyscf import dft
from pyscf.prop import nmr

mol = gto.Mole()
mol.verbose = 5

mol.atom = '''
     O      0.   0.       0.
     H      0.  -0.757    0.587
     H      0.   0.757    0.587'''
mol.basis = 'ccpvdz'
mol.build()

def finger(mat):
    w = numpy.cos(numpy.arange(mat.size))
    return numpy.dot(w, mat.ravel())


mf = dft.RKS(mol)
mf.conv_tol_grad = 1e-6
mf.grids.prune = False
mf.xc = 'b3lypg'
mf.scf()

m = nmr.RKS(mf)
msc = m.kernel()

#print(finger(msc))

print(msc)

