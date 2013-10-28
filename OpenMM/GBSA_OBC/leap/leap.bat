# A script to create N-acetyl-alanine-N′-methylamide (i.e. blocked alanine dipeptide)

source leaprc.ff99SB
set default PBradii mbondi2
ACEALANME = {sequence ACE ALA NME }
savepdb ACEALANME snap.pdb
saveamberparm ACEALANME prmtop inpcrd
quit
