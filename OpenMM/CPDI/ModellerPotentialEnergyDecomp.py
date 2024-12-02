# Obtaining per term energy contributions in OpenMM is hard.
# This takes an original system and then copies out each class
# of force to their own system.

from openmm.app import *
from openmm import *
from simtk.unit import *
from sys import stdout
import time
from copy import deepcopy 
import decimal 

platform = openmm.Platform.getPlatformByName("Reference")


# OpenMM 5.1 only
Topology.loadBondDefinitions('CPDI_CYP_residues.xml')
Modeller.loadHydrogenDefinitions('CPDI_CYP_hydrogens.xml')


# PDB
# Note, leap miscreates this file. The line:
# ATOM     28 FE   HEM     1     -15.846 -23.032 -11.293  1.00  0.00           FE
# should be
# ATOM     28 FE   HEM     1     -15.846 -23.032 -11.293  1.00  0.00          FE
pdb = PDBFile('./leap_example/out.pdb')
forceField = ForceField('CPDI_CYP_ff.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forceField)

# Dump modeller structure
PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(),open('modeller.pdb', 'w'))


system = forceField.createSystem(modeller.topology, nonbondedMethod=NoCutoff)




# Remember, this is being run NVE
integrator = VerletIntegrator(1*femtoseconds)

simulation = Simulation(pdb.topology, system, integrator, platform)

print("Platform: ", (simulation.context.getPlatform().getName()))

simulation.context.setPositions(pdb.positions)

# Entire system
system = simulation.context.getSystem()
state = simulation.context.getState( getEnergy=True)
print("Total potential energy is ",  str(state.getPotentialEnergy().in_units_of(kilocalorie/mole)))
print("")



# Create a map between index and name
indexToAtomNameDict = {}
atoms =  pdb.topology.atoms()
print("Creating map")
for atom in atoms:
  indexToAtomNameDict[atom.index] = atom.name
  print( atom.index, atom.name)



# Create new state with only Harmonic Bonds in
HarmonicBondSystem = System()
for j in range (system.getNumParticles()):
  HarmonicBondSystem.addParticle(system.getParticleMass(j) )
HarmonicBondIntegrator = VerletIntegrator(1*femtoseconds)


# Creat state with only Harmonic Angles in
HarmonicAngleSystem = System()
for j in range (system.getNumParticles()):
  HarmonicAngleSystem.addParticle(system.getParticleMass(j) )
HarmonicAngleIntegrator = VerletIntegrator(1*femtoseconds)


# Create state with only Torsion terms in
PeriodicTorsionSystem = System()
for j in range (system.getNumParticles()):
  PeriodicTorsionSystem.addParticle(system.getParticleMass(j) )
PeriodicTorsionIntegrator = VerletIntegrator(1*femtoseconds)


# Create state with only Nonbonded terms in
NonbondedSystem = System()
for j in range (system.getNumParticles()):
  NonbondedSystem.addParticle(system.getParticleMass(j) )
NonbondedIntegrator = VerletIntegrator(1*femtoseconds)



# Harmonic angle debug
#force = system.getForce(1)
#for i in range(force.getNumAngles()):
#   a = force.getAngleParameters(i)[0]
#   b = force.getAngleParameters(i)[1]
#   c = force.getAngleParameters(i)[2]
#   print indexToAtomNameDict[a], indexToAtomNameDict[b], indexToAtomNameDict[c] + " " + str(force.getAngleParameters(i))



#force = system.getForce(2)
#for i in range(force.getNumTorsions()):
#   a = force.getTorsionParameters(i)[0]
#   b = force.getTorsionParameters(i)[1]
#   c = force.getTorsionParameters(i)[2]
#   d = force.getTorsionParameters(i)[3]
#   print indexToAtomNameDict[a], indexToAtomNameDict[b], indexToAtomNameDict[c], indexToAtomNameDict[d] + " " + str(force.getTorsionParameters(i))

#Non-bonded
#force = system.getForce(3)
#for i in range(force.getNumParticles()):
#   print indexToAtomNameDict[i] + " " + str(force.getParticleParameters(i))



# Decompose and copy out respective force terms
for i in range(system.getNumForces()):

   force = system.getForce(i)
   print(type(force))


   if isinstance( force, openmm.HarmonicBondForce ):
     print( "Found ", str(force.getNumBonds()), " HarmonicBondForce terms")

     # Deep copy these forces into out new state
     copyOfForce = deepcopy(force)
     # Add this to our HarmonicBondSystem
     HarmonicBondSystem.addForce(copyOfForce)

   if isinstance( force, openmm.HarmonicAngleForce ):
     print("Found ", str(force.getNumAngles()), " HarmonicAngleForce terms")

     copyOfForce = deepcopy(force)
     HarmonicAngleSystem.addForce(copyOfForce)

   if isinstance( force, openmm.PeriodicTorsionForce ):
     print("Found ", str(force.getNumTorsions()), " PeriodicTorsionForce terms")

     copyOfForce = deepcopy(force)
     PeriodicTorsionSystem.addForce(copyOfForce)

   if isinstance( force, openmm.NonbondedForce ):
     print("Found ", str(force.getNumParticles()), " NonbondedForce terms )")
     print("Found ", str(force.getNumExceptions()), " NonbondedForce exception terms (i.e. 1-4) ")

    
     copyOfForce = deepcopy(force)
     NonbondedSystem.addForce(copyOfForce)





## Now evaluate each system

print ("")
# Harmonic Bond
HarmonicBondContext = Context(HarmonicBondSystem, HarmonicBondIntegrator, platform)
HarmonicBondContext.setPositions(pdb.positions)
HarmonicBondState = HarmonicBondContext.getState( getEnergy=True)
print ("Harmonic bond: " + "\t\t\t" + str( HarmonicBondState.getPotentialEnergy().in_units_of(kilocalorie/mole)))

# Harmonic Angle
HarmonicAngleContext = Context(HarmonicAngleSystem, HarmonicAngleIntegrator, platform)
HarmonicAngleContext.setPositions(pdb.positions)
HarmonicAngleState = HarmonicAngleContext.getState( getEnergy=True)
print ("Harmonic angle: " + "\t\t" + str(HarmonicAngleState.getPotentialEnergy().in_units_of(kilocalorie/mole)))

# Torsions
PeriodicTorsionContext = Context(PeriodicTorsionSystem, PeriodicTorsionIntegrator, platform)
PeriodicTorsionContext.setPositions(pdb.positions)
PeriodicTorsionState = PeriodicTorsionContext.getState( getEnergy=True)
print ("Proper and improper torsions: " + "\t" + str(PeriodicTorsionState.getPotentialEnergy().in_units_of(kilocalorie/mole)))

# NB
NonbondedContext = Context(NonbondedSystem, NonbondedIntegrator, platform)
NonbondedContext.setPositions(pdb.positions)
NonbondedState = NonbondedContext.getState( getEnergy=True)

print ("EE, VDW, 14EE and 14VDW: " + "\t" + str(NonbondedState.getPotentialEnergy().in_units_of(kilocalorie/mole)))



# Output forces
print("Forces")
forces = [None]
forces = simulation.context.getState(getForces=True).getForces()


for force in forces:
     #print("{0}".format( force.in_units_of(kilocalorie/mole/angstrom) ))
     print(" {0:+1.16E} {1:+1.16E} {2:+1.16E}".format(  force.value_in_unit(kilocalorie/mole/angstrom)[0] , force.value_in_unit(kilocalorie/mole/angstrom)[1], force.value_in_unit(kilocalorie/mole/angstrom)[2]) )


