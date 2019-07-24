# Obtaining per term energy contributions in OpenMM is hard.
# This takes an original system and then copies out each class
# of force to their own system.

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time
from copy import deepcopy 

#platform = openmm.Platform_getPlatformByName("OpenCL")
#platform = openmm.Platform_getPlatformByName("Cuda")
platform = openmm.Platform_getPlatformByName("Reference")



# Run on multiple cards
# 0  Tesla M2090
# 1  Tesla C2075
# 2  Tesla C2075
#platformProperties = {"OpenCLDeviceIndex":"0"}
#platformProperties = {"OpenCLDeviceIndex":"1"}
#platformProperties = {"OpenCLDeviceIndex":"0,1,2"}
#platformProperties = {"OpenCLDeviceIndex":"1,2"}
platformProperties = {}


prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd')

#system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)
#system = prmtop.createSystem(nonbondedMethod=NoCutoff, implicitSolvent=OBC2)
#system = prmtop.createSystem(nonbondedMethod=NoCutoff, removeCMMotion=False)
system = prmtop.createSystem(nonbondedMethod=NoCutoff)

# Remember, this is being run NVE
integrator = VerletIntegrator(1*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)
print("Platform: %s" % (simulation.context.getPlatform().getName()))
print("Number of atoms %i" , len(inpcrd.positions))

simulation.context.setPositions(inpcrd.positions)

# Entire system
system = simulation.context.getSystem()
state = simulation.context.getState( getEnergy=True)
print("Total potential energy is %s", str(state.getPotentialEnergy().in_units_of(kilocalorie/mole)))
print("")



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






# Decompose and copy out respective force terms
for i in range(system.getNumForces()):

   force = system.getForce(i)
   print(type(force))


   if isinstance( force, openmm.HarmonicBondForce ):
     print("Found %i", str(force.getNumBonds()), " HarmonicBondForce terms")

     # Deep copy these forces into out new state
     copyOfForce = deepcopy(force)
     # Add this to our HarmonicBondSystem
     HarmonicBondSystem.addForce(copyOfForce)

   if isinstance( force, openmm.HarmonicAngleForce ):
     print("Found %i", str(force.getNumAngles()),  " HarmonicAngleForce terms")

     copyOfForce = deepcopy(force)
     HarmonicAngleSystem.addForce(copyOfForce)

   if isinstance( force, openmm.PeriodicTorsionForce ):
     print("Found %i", str(force.getNumTorsions()), " PeriodicTorsionForce terms")

     copyOfForce = deepcopy(force)
     PeriodicTorsionSystem.addForce(copyOfForce)

   if isinstance( force, openmm.NonbondedForce ):
     print("Found %i" , str(force.getNumParticles()) , " NonbondedForce terms ")
     print("Found %i" , str(force.getNumExceptions()) , " NonbondedForce exception terms (i.e. 1-4) ")

    
     copyOfForce = deepcopy(force)
     NonbondedSystem.addForce(copyOfForce)





## Now evaluate each system

print ("")
# Harmonic Bond
HarmonicBondContext = Context(HarmonicBondSystem, HarmonicBondIntegrator, platform)
HarmonicBondContext.setPositions(inpcrd.positions)
HarmonicBondState = HarmonicBondContext.getState( getEnergy=True)
print("Harmonic bond: " + "\t\t\t" + str( HarmonicBondState.getPotentialEnergy().in_units_of(kilocalorie/mole)))

# Harmonic Angle
HarmonicAngleContext = Context(HarmonicAngleSystem, HarmonicAngleIntegrator, platform)
HarmonicAngleContext.setPositions(inpcrd.positions)
HarmonicAngleState = HarmonicAngleContext.getState( getEnergy=True)
print("Harmonic angle: " + "\t\t" + str(HarmonicAngleState.getPotentialEnergy().in_units_of(kilocalorie/mole)))

# Torsions
PeriodicTorsionContext = Context(PeriodicTorsionSystem, PeriodicTorsionIntegrator, platform)
PeriodicTorsionContext.setPositions(inpcrd.positions)
PeriodicTorsionState = PeriodicTorsionContext.getState( getEnergy=True)
print("Proper and improper torsions: " + "\t" + str(PeriodicTorsionState.getPotentialEnergy().in_units_of(kilocalorie/mole)))

# NB
NonbondedContext = Context(NonbondedSystem, NonbondedIntegrator, platform)
NonbondedContext.setPositions(inpcrd.positions)
NonbondedState = NonbondedContext.getState( getEnergy=True)

print("EE, VDW, 14EE and 14VDW: " + "\t" + str(NonbondedState.getPotentialEnergy().in_units_of(kilocalorie/mole)))



# Output forces
print("")
forces = [None]
forces = simulation.context.getState(getForces=True).getForces()


for force in forces:
     print(force)
