# A script to create N-acetyl-alanine-N′-methylamide

source leaprc.ff99SB
set default PBradii mbondi2
ACEALANME = {sequence ACE ALA NME }
savepdb ACEALANME snap.pdb
saveamberparm ACEALANME prmtop inpcrd
quit
