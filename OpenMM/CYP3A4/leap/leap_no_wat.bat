# We will use FF99SB for the protein
source leaprc.ff99SB

# Load GAFF since some of Tom
source leaprc.gaff

# Load Tom's heme the related parameters.
loadOff ./tom_params/HEM.off
loadOff ./tom_params/CYP.off
loadAmberParams ./tom_params/ic6.fh.frcmod

# Generate a mapping between the PBD's heme atoms
# and the ones in Tom's parameters

# cat heme.pdb | awk '{print " { "$3 " TODO}"}'

addPdbAtomMap {
 { FE FE1}
 { CHA C13}
 { CHB C20}
 { CHC C26}
 { CHD C6}
 { NA N18}
 { C1A C15}
 { C2A C16}
 { C3A C19}
 { C4A C17}
 { CMA C52}
 { CAA C56}
 { CBA C62}
 { CGA C69}
 { O1A O70}
 { O2A O72}
 { NB N24}
 { C1B C21}
 { C2B C22}
 { C3B C29}
 { C4B C25}
 { CMB C40}
 { CAB C30}
 { CBB C32}
 { NC N1}
 { C1C C2}
 { C2C C4}
 { C3C C5}
 { C4C C3}
 { CMC C44}
 { CAC C35}
 { CBC C37}
 { ND N9}
 { C1D C8}
 { C2D C12}
 { C3D C11}
 { C4D C10}
 { CMD C48}
 { CAD C59}
 { CBD C65}
 { CGD C68}
 { O1D O71}
 { O2D O73}
}


# Load the preparsed PDB
CYP3A4 = loadpdb 1TQN_modded.pdb


# Internal sanity checks
check CYP3A4

# Save a pdb for convience of debugging
savepdb CYP3A4 snap_no_wat.pdb

# Generate the AMBER parameter and coordinate files
saveamberparm CYP3A4 prmtop_no_wat inpcrd_no_wat

quit

