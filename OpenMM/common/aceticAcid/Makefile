include config.h

mdout: prmtop
	$(AMBERHOME)/bin/sander -O


prmtop: acetic_acid.mol2
	rm -f leap.log ; \
	$(AMBERHOME)/bin/tleap -f leap.bat 

acetic_acid.mol2: ./raw/acetic_acid.mol2
	$(AMBERHOME)/bin/antechamber -i ./raw/acetic_acid.mol2 -fi mol2 -c bcc -o acetic_acid.mol2 -fo mol2 -rn ACE





clean:
	rm -f ANTECHAMBER* sqm* acetic_acid.mol2 ATOMTYPE.INF
	rm -f prmtop inpcrd leap.log
	rm -f mdinfo mdout restrt snap.pdb
