Testing to see if the CPDI parameters by Cheatham et al. can
be successfully converted to OpenMM XML (using JAMBER) and then 
applied to a PDB using OpenMM's modeller tool.

Note, <ExternalBond from="27"/> in CPDI_CYP_ff.xml has been removed
for this example



1) Benchmark; single point energy using SANDER
==============================================

   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
      1       3.2344E+02     5.3098E+01     5.0919E+02     FE         28

 BOND    =      141.5828  ANGLE   =       88.8884  DIHED      =       65.4066
 VDWAALS =       -8.8151  EEL     =       12.3517  HBOND      =        0.0000
 1-4 VDW =        5.1883  1-4 EEL =       18.8373  RESTRAINT  =        0.0000



2) Single point energy using OpenMM via its prmtop+inpcrd readers
=================================================================
(PotentialEnergyDecomp.py)

Total potential energy is 323.44091565 kcal/mol

<class 'simtk.openmm.openmm.HarmonicBondForce'>
Found 81 HarmonicBondForce terms
<class 'simtk.openmm.openmm.HarmonicAngleForce'>
Found 148 HarmonicAngleForce terms
<class 'simtk.openmm.openmm.PeriodicTorsionForce'>
Found 236 PeriodicTorsionForce terms
<class 'simtk.openmm.openmm.NonbondedForce'>
Found 74 NonbondedForce terms )
Found 423 NonbondedForce exception terms (i.e. 1-4)
<class 'simtk.openmm.openmm.CMMotionRemover'>

Harmonic bond:                  141.582807548 kcal/mol
Harmonic angle:                 88.888382108 kcal/mol
Proper and improper torsions:   65.4065605046 kcal/mol
eE, VDW, 14EE and 14VDW:        27.5631654896 kcal/mol


=> Hence, this is all good. OpenMM's energies match SANDER
for the heme residue that has been converted using JAMBER,
to generate CPDI_CYP.xml .


3) Single point energy using OpenMM via modeller, reading
the output pdb (with hydrogens) from AMBER, but not adding
any hydrogens.

ModellerPotentialEnergyDecomp.py
================================

Total potential energy is 323.447925604 kcal/mol

<class 'simtk.openmm.openmm.HarmonicBondForce'>
Found 81 HarmonicBondForce terms
<class 'simtk.openmm.openmm.HarmonicAngleForce'>
Found 146 HarmonicAngleForce terms
<class 'simtk.openmm.openmm.PeriodicTorsionForce'>
Found 236 PeriodicTorsionForce terms
<class 'simtk.openmm.openmm.NonbondedForce'>
Found 74 NonbondedForce terms )
Found 423 NonbondedForce exception terms (i.e. 1-4) 
<class 'simtk.openmm.openmm.CMMotionRemover'>

Harmonic bond:                  141.583758575 kcal/mol
Harmonic angle:                 88.8821495736 kcal/mol
Proper and improper torsions:   65.4099409076 kcal/mol
EE, VDW, 14EE and 14VDW:        27.5720765474 kcal/mol


This is close enough for the moment; differences in the 3rd decimal place
will be due to rounding errors in parameters in the XML and will be 
addressed soon.

==> Hence modeller can parse the information in CPDI_CYP_ff.xml.

However, one quirk has appeared; if add addHydrogens() is used on the
complete ./leap_example/out.pdb , differences are seen:

<class 'simtk.openmm.openmm.HarmonicBondForce'>
Found 81 HarmonicBondForce terms
<class 'simtk.openmm.openmm.HarmonicAngleForce'>
Found 146 HarmonicAngleForce terms
<class 'simtk.openmm.openmm.PeriodicTorsionForce'>
Found 236 PeriodicTorsionForce terms
<class 'simtk.openmm.openmm.NonbondedForce'>
Found 74 NonbondedForce terms )
Found 423 NonbondedForce exception terms (i.e. 1-4) 
<class 'simtk.openmm.openmm.CMMotionRemover'>

Harmonic bond:                  141.357906796 kcal/mol
Harmonic angle:                 88.5768810531 kcal/mol
Proper and improper torsions:   40.7086009112 kcal/mol
EE, VDW, 14EE and 14VDW:        28.203629494 kcal/mol

There is around 20 kcal/mol of energy missing from the proper and improper torsions.

THIS IS FINE; addHydrogens rebuilds where the hydrogens are, hence the topology is the
same, but the coordinates of the hydrogens are slightly different.

==> Hence, modeller.addHydrogens is able to process the information in
CPDI_CYP_hydrogens.xml correctly for the HEM residue.


Using the ./1TQN_HEM.pdb structure (i.e. without hydrogens and adding them)
the same is seen:


4)  Single point energy using OpenMM via modeller, reading
the /1TQN_HEM.pdb , and reconstructing i.e. add hydrogens
(ModellerPotentialEnergyDecomp.py)

Total potential energy is 298.907093835 kcal/mol

<class 'simtk.openmm.openmm.HarmonicBondForce'>
Found 81 HarmonicBondForce terms
<class 'simtk.openmm.openmm.HarmonicAngleForce'>
Found 146 HarmonicAngleForce terms
<class 'simtk.openmm.openmm.PeriodicTorsionForce'>
Found 236 PeriodicTorsionForce terms
<class 'simtk.openmm.openmm.NonbondedForce'>
Found 74 NonbondedForce terms )
Found 423 NonbondedForce exception terms (i.e. 1-4) 
<class 'simtk.openmm.openmm.CMMotionRemover'>

Harmonic bond:                  141.359370846 kcal/mol
Harmonic angle:                 88.6099751581 kcal/mol
Proper and improper torsions:   40.7164024948 kcal/mol
EE, VDW, 14EE and 14VDW:        28.2213453359 kcal/mol


Again, this is fine because the coordinates of the hydrogens are rebuilt, but 
the same topology remains.
