Example AMBER system with bespoke residues, HEM and CYP, 
that can be built from a PDB using Modeller.

HEM is the heme moiety in the compound I state (CPDI) and 
CYP is a modified cysteine residue that is covalently 
bonded to the iron on the HEM.

These parameters are from http://dx.doi.org/10.1002/jcc.21922
and have been converted using JAMBER:
 https://bitbucket.org/mjw99/jamber

CPDI_CYP.xml contains the parameters for HEM and CYP.

CPDI_CYP_hydrogens.xml contains information about missing hydrogen 
atoms on both these residues and is used by modeller to add missing
hydrogens to these residues when loading a PDB. This needs to be
added to simtk/openmm/app/data/hydrogens.xml

CPDI_CYP_residues.xml contains information about the CPDI and
CYP residues. This needs to be added to 
simtk/openmm/app/data/residues.xml

A system object is created by Modeller. This is used in 
conjuction with an integrator to create a simulation 
object.

An example drug, ketoconazole is added to the system using modeller's add 
method. Parameters for ketoconazole have been generated using antechamber
and converted using JAMBER.
