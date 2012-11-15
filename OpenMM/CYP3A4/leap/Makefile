include ../config.h

all: prmtop


1TQN_with_missing_residues.pdb: 1TQN.pdb missing_4_residues.pdb
	# Add in the 4 missing residues
	# grep -e " N " -e CA -e "O " ../mmmd/1/aligned.pdb
	#
	# Begin
	#ATOM   2035  NE2 GLN A 279     -13.310 -26.229  11.067  1.00 78.32           N
	grep "ATOM   2035" -B 10000 1TQN.pdb > start.pdb
	#
	# End
	#ATOM   2066  N   LYS A 288      -5.848 -16.607   8.893  1.00 78.27           N
	grep "ATOM   2066" -A 10000 1TQN.pdb > end.pdb

	# Splice in missing 4 residues (plus two buffer residues either end):
	cat missing_4_residues.pdb >> start.pdb
	# Add on end
	cat end.pdb >> start.pdb
	mv start.pdb 1TQN_with_missing_residues.pdb
	rm end.pdb

prmtop:	leap.bat 1TQN_with_missing_residues.pdb  ./pqr/propka_2_with_4_missing_residues_added/1TQN_with_missing_residues.pqr
	rm -f leap.log ; \

	# Change the proximal cysteine to type CYP (deprotonated)
	cat 1TQN_with_missing_residues.pdb | sed s/"CYS A 442"/"CYP A 442"/g > 1TQN_mod2.pdb

	# Mutate the distal water into the ferryl oxygen.
 	# This is a workaround to fix a tleap bug:
	# http://bugzilla.ambermd.org/show_bug.cgi?id=151
	# O   HOH A 509 --> O2D HEM A 508
	cat 1TQN_mod2.pdb | sed s/"O   HOH A 509"/"O1  HEM A 508"/ > 1TQN_mod3.pdb

	# Assign correct protonation states
	# These are determined by ./pqr/propka_2_with_4_missing_residues_added/1TQN_with_missing_residues.pqr

	#HIP 30
	#HIP 54
	#HIE 65
	#HID 267
	#HID 287
	#HIE 324
	#HID 402

	cat 1TQN_mod3.pdb | \
	sed s/"HIS A  30"/"HIP A  30"/g |\
	sed s/"HIS A  54"/"HIP A  54"/g |\
	sed s/"HIS A  65"/"HIE A  65"/g |\
	sed s/"HIS A 267"/"HID A 267"/g |\
	sed s/"HIS A 287"/"HID A 287"/g |\
	sed s/"HIS A 324"/"HIE A 324"/g |\
	sed s/"HIS A 402"/"HID A 402"/g \
	> 1TQN_mod4.pdb

	# Update the HELIX records as well since gleap will fail
	# Note, this is not needed, since we are back to using tleap
	# and has been left in since it is technically the right thing to 
	# do, but tleap does not care about the secondary structure
	# sections in a PDB.
	cat 1TQN_mod4.pdb | \
	sed s/"HIS A  324"/"HIE A  324"/g \
	> 1TQN_modded.pdb

	$(AMBERHOME)/bin/tleap -f leap.bat ; \
	# Use specific incarnation of MAKEPDB since we are not working 
	# with a restart here
	$(AMBERHOME)/bin/ambpdb -p prmtop < inpcrd > snap.pdb
clean:
	rm -f inpcrd prmtop snap.pdb *.log 1TQN_mod* 1TQN_with_missing_residues.pdb
