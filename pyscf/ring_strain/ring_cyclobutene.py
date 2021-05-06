#!/usr/bin/env python

# Attempt at reproducing the computation experiment in
# http://dx.doi.org/10.1021/ma0501287

from pyscf import gto
from pyscf import dft

from pyscf.geomopt import berny_solver

#basis='sto-3g'
#basis='ccpvdz'
#basis='6-31G**'
basis='6-31G*'
#basis='aug-cc-pvdz'

ethene = gto.Mole()
ethene.atom = '''
   C        0.00000000    -0.66547425     0.00000000
   C        0.00000000     0.66547425     0.00000000
   H        0.92358611    -1.23957529     0.00000000
   H       -0.92358611    -1.23957529     0.00000000
   H        0.92358611     1.23957529     0.00000000
   H       -0.92358611     1.23957529     0.00000000'''
ethene.basis = basis
ethene.build()

cyclobutene = gto.Mole()
cyclobutene.atom = '''
   C       -0.00170400    -0.78623897     0.70103078
   C        0.00170400     0.78623897     0.70103078
   C       -0.00171504     0.67017663    -0.81309212
   H        0.89554248     1.24325386     1.14335578
   H       -0.88299298     1.25157302     1.15270599
   C        0.00171504    -0.67017663    -0.81309212
   H       -0.89554248    -1.24325386     1.14335578
   H        0.88299298    -1.25157302     1.15270599
   H       -0.00353388     1.42022143    -1.60005736
   H        0.00353388    -1.42022143    -1.60005736'''
cyclobutene.basis = basis
cyclobutene.build()

one_five_hexadiene = gto.Mole()
one_five_hexadiene.atom = '''
   C       -0.09665729     1.47211137    -0.35142135
   C        0.90360762     2.20197994     0.14523780
   C        1.33977048    -1.54596319    -0.13124177
   C        0.04924829    -1.81698446     0.07479350
   H        1.68506669    -0.57527410    -0.47353830
   H        1.17646948     2.15525023     1.19766675
   C       -0.98398503     0.54834348     0.43918355
   H       -0.33178570     1.55262660    -1.41521903
   H        1.48662596     2.87487005    -0.47820737
   H        2.10298422    -2.29948153     0.04762345
   C       -1.12607765    -0.88282516    -0.13584173
   H       -1.99347521     0.98742529     0.46236607
   H       -0.63920925     0.49660963     1.47898925
   H       -1.35126151    -0.80829493    -1.21092133
   H       -2.01467215    -1.34047238     0.31530207
   H       -0.21275391    -2.81723932     0.42274517'''
one_five_hexadiene.basis = basis
one_five_hexadiene.build()






def opt_molecule(mol):

   method = dft.RKS(mol)
   method.grids.prune = dft.gen_grid.treutler_prune
   method.grids.atom_grid = {"H": (50, 194), "C": (50, 194),}
   method.xc = 'b3lyp'
   return method.scf()



a = opt_molecule(one_five_hexadiene)
b = opt_molecule(cyclobutene)
c = opt_molecule(ethene)

#print "b+c: " + str((b+c))
#print "a: " + str(a)
print(( (b + c) - a )*627.509)
