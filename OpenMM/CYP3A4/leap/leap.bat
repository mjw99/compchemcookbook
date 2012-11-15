# We will use FF99SB for the protein
source leaprc.ff99SB

# Load GAFF since some of Tom
source leaprc.gaff

# Load Tom's heme the related parameters.
HEM = loadMol2 HEM.mol2
CYP = loadMol2 CYP.mol2
loadAmberParams CPDI.frcmod

# Needed because mol2 does not contain head/tail info
set CYP head CYP.1.N
set CYP tail CYP.1.C
desc CYP

# Load the preparsed PDB
CYP3A4 = loadpdb 1TQN_modded.pdb

# Solvate Oct
solvateBox CYP3A4 TIP3PBOX 9.0

# Add counter ions
addions CYP3A4 Cl- 0

# Internal sanity checks
#check CYP3A4

# Save a pdb for convience of debugging
savepdb CYP3A4 snap.pdb

# Generate the AMBER parameter and coordinate files
saveamberparm CYP3A4 prmtop inpcrd

quit

