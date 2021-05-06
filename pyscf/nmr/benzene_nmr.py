import numpy
from pyscf import gto
from pyscf import dft
from pyscf.prop import nmr
from numpy import linalg as LA

mol = gto.Mole()
#mol.verbose = 5

mol.atom = '''
    C                   -1.03071033     0.94315337     0.00000000
    C                   -1.33671671    -0.41962694     0.00000000
    C                   -0.30953506    -1.36596640     0.00000000
    C                    1.02358258    -0.94962752     0.00000000
    C                    1.32958443     0.41308298     0.00000000
    C                    0.30244216     1.35945863     0.00000000
    H                    1.82296797    -1.68619635     0.00000000
    H                    2.36713226     0.73723910     0.00000000
    H                   -2.37429908    -0.74372191     0.00000000
    H                   -0.54764449    -2.42657221     0.00000000
    H                   -1.83020449     1.67964543     0.00000000
    H                    0.54060403     2.42005736     0.00000000'''

mol.basis = 'ccpvdz'
mol.build()

def finger(mat):
    w = numpy.cos(numpy.arange(mat.size))
    return numpy.dot(w, mat.ravel())

mf = dft.RKS(mol)
mf.conv_tol_grad = 1e-6
mf.grids.prune = False
mf.xc = 'b3lypg'
#mf.xc = 'b3lyp'
#mf.xc = 'm062x'
mf.scf()

m = nmr.RKS(mf)
msc = m.kernel()


#TMS shielding values of 189.78 (C) and 32.18 (H).
for i in range(11):
  print(numpy.sum(LA.eigvalsh(msc[i]))/3.0)
