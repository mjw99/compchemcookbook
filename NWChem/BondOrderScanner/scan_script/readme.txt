# The following NWChem script makes use of python to enumerate
# all C-H bonds in a molecule and then carry out a relaxed
# scan along each one of them. This specific example uses
# methane.
#
# The python script will:
# a) Enumerate all C-H bonds
# b) Optimise the structure to a minimum
# c) For each, peform a relaxed scan compressing and stretching the
#    bond in question by a defined amount away from the bond length's
#    optimised values (currently +/- 0.2A)
# d) The energy profile of each scanned bond will be written to a file
#    "bond_scan_i_j.dat", where i and j are the indices of the C-H atoms.
# e) These can be post-processed using the FindDe.py script to fit a 
     Morse potential and determin the dissociation energy D_e


Preamble
*) Ensure you have a python enabled build of NWchem; see googledoc on how to do this.
	source setup_nwchem_environment.sh

*) SCipy installed
	sudo apt-get install python-scipy


1) Run the methane scan script
	nwchem methane_scan.nw &> methane_scan.nwlog &

2) Post process the scanned energies. This will fit a Morse potential
   to the scanned energies and work out the dissociation energy D_e
	ls bond_scan_1_* | xargs -I {} ./wrapperToFindDe.py {}

3) An animation of the scan can be generated with the following command:
	jmol -s generate_movie.jmol
	   
  
