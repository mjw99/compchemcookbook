# Load GAFF since some of Tom
source leaprc.gaff

# Load Tom's heme the related parameters.
HEM = loadMol2 HEM.mol2
loadAmberParams CPDI.frcmod

# Needed because mol2 does not contain head/tail info
set CYP head CYP.1.N
set CYP tail CYP.1.C
desc CYP

# Load the preparsed PDB
CYP3A4 = loadpdb heme.pdb

# Save a pdb for convience of debugging
savepdb CYP3A4 snap.pdb

# Generate the AMBER parameter and coordinate files
saveamberparm CYP3A4 prmtop inpcrd

quit

