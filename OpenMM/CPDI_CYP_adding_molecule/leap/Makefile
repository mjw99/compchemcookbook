AMBERHOME=/home/mw529/code/AMBER/amber11

all: ketoconazole.mol2

ketoconazole.mol2: 
	$(AMBERHOME)/bin/antechamber -i ketoconazole.pdb -fi pdb -c bcc -o ketoconazole.mol2 -fo mol2 -rn KTN

clean:
	rm -f inpcrd prmtop ANTECHAMBER_* sqm.in sqm.out ATOMTYPE.INF snap.pdb octanol.mol2
