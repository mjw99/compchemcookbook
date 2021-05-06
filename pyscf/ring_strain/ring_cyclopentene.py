#!/usr/bin/env python

# Attempt at reproducing the computation experiment in
# http://dx.doi.org/10.1021/ma0501287

from pyscf import gto
from pyscf import dft

from pyscf.geomopt import berny_solver

basis='6-31G*'
#basis='STO-3G'

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

cyclopentene = gto.Mole()
cyclopentene.atom = '''
   C        0.27164387    -1.25067230    -0.05566139
   C        1.21303136    -0.03417642     0.15773275
   C        0.34701699     1.22103960    -0.13581856
   H        2.11332754    -0.08180598    -0.46278396
   H        1.54161395    -0.01026984     1.20285149
   C       -1.10729297    -0.63571896     0.03513673
   H        0.43816314    -2.04248466     0.68598713
   H        0.42173950    -1.71786017    -1.04129813
   C       -1.06666032     0.69831867    -0.00874188
   H        0.52375014     1.61392189    -1.14904079
   H        0.56235035     2.04734185     0.55381760
   H       -2.01453535    -1.23213185     0.08698041
   H       -1.93564448     1.35114061     0.00406174'''
cyclopentene.basis = basis
cyclopentene.build()

one_six_heptadiene = gto.Mole()
one_six_heptadiene.atom = '''
   C        0.65948863    -3.49253887    -0.13792007
   C       -0.12027685    -2.53665231     0.36849820
   C       -0.53939132    -1.28340735    -0.35044157
   C       -0.06097215     0.00009111     0.35436324
   H       -1.63800864    -1.25327892    -0.42261030
   H       -0.15934964    -1.30356420    -1.38052231
   C       -0.53141783     1.28346130    -0.35598341
   H        1.03468561    -0.00322085     0.41427222
   H       -0.42571460     0.00342475     1.39176403
   C       -0.10459023     2.53716899     0.35748025
   H       -1.63019680     1.25978449    -0.42808341
   H       -0.15130761     1.29675273    -1.38615578
   C        0.68097911     3.48604374    -0.15313747
   H       -0.47193243     2.64854036     1.37972924
   H        1.07356327     3.41868959    -1.16613960
   H        0.96049699     4.36805661     0.41707300
   H        0.93356567    -4.37370437     0.43619221
   H        1.05248972    -3.43202141    -1.15117035
   H       -0.48830186    -2.64123827     1.39123350'''
one_six_heptadiene.basis = basis
one_six_heptadiene.build()



def energy_molecule(mol):
   method = dft.RKS(mol)
   method.grids.prune = dft.gen_grid.treutler_prune
   method.grids.atom_grid = {"H": (50, 194), "C": (50, 194),}
   method.xc = 'b3lyp'
   return method.scf()


a = energy_molecule(one_six_heptadiene)
b = energy_molecule(cyclopentene)
c = energy_molecule(ethene)

print ( ( b + c ) - a )*627.509 


#print ( a  )*627.509 
#print (  ( b + c )  )*627.509 

