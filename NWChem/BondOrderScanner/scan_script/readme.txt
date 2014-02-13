# The following NWChem script makes use of python to enumerate
# all C-H bonds in a molecule and then carry out a relaxed
# scan along each one of them.
#
# For each C-H bond, a file "bond_scan_1_2.dat" for the index
# of each scanned bond will be generated. 
#


# input_formatter.py is in nwchem.py which is in methane_scan.nw


*) Ensure you have a python enabled build of NWchem; see googledoc on how to do this.
	source setup_nwchem_environment.sh

*) SCipy installed
	sudo apt-get install python-scipy


1) run the methane scan script
	nwchem methane_scan.nw &> methane_scan.log &

2) Post process the scanned energies. This will fit a Morse potential
   to the scanned energies and work out the dissociation energy D_e

	./wrapperToFindDe.py bond_scan_1_2.dat
   
  
