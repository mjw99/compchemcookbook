# DHRF Benchmark using OpenMM

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
platformProperties = {"OpenCLDeviceIndex":"0"}
#platformProperties = {"OpenCLDeviceIndex":"1"}
#platformProperties = {"OpenCLDeviceIndex":"0,1,2"}
#platformProperties = {"OpenCLDeviceIndex":"1,2"}



prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd')

#system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)
#system = prmtop.createSystem(nonbondedMethod=NoCutoff, implicitSolvent=OBC2)
system = prmtop.createSystem(nonbondedMethod=NoCutoff)

# Remember, this is being run NVE
integrator = VerletIntegrator(1*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)
print "Platform: %s" % (simulation.context.getPlatform().getName())
print "Number of atoms %i"      % len(inpcrd.positions)

simulation.context.setPositions(inpcrd.positions)

#simulation.reporters.append(PDBReporter('output.pdb', 1))
#simulation.reporters.append(StateDataReporter(stdout, 1, step=True,  totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True))


system = simulation.context.getSystem()
state = simulation.context.getState( getEnergy=True)

print "Total potential energy is " +  str(state.getPotentialEnergy().in_units_of(kilocalorie/mole))


for i in range(system.getNumForces()):

   force = system.getForce(i)
   print type(force)
   if isinstance( force, openmm.HarmonicBondForce ):
     print "Found " + str(force.getNumBonds()) + " HarmonicBondForce terms"

     # New state with only Harmonic Bonds in
     HarmonicBondSystem = System()
     # Deep copy these forces into out new state
     copyOfForce = deepcopy(force)
     # Add this to our HarmonicBondSystem
     HarmonicBondSystem.addForce(copyOfForce)

     # Add the particles?
     for j in range (system.getNumParticles()):
       HarmonicBondSystem.addParticle(system.getParticleMass(j) )

     HarmonicBondIntegrator = VerletIntegrator(1*femtoseconds)
     HarmonicBondContext = Context(HarmonicBondSystem, HarmonicBondIntegrator, platform)
     HarmonicBondContext.setPositions(inpcrd.positions)

     HarmonicBondState = HarmonicBondContext.getState( getEnergy=True)
    
     print (HarmonicBondState.getPotentialEnergy().in_units_of(kilocalorie/mole))
     print str(5.6366273691689504) + " (JAMBER)"

   if isinstance( force, openmm.HarmonicAngleForce ):
     print "Found " + str(force.getNumAngles()) + " HarmonicAngleForce terms"

     # New state with only Harmonic Bonds in
     HarmonicAngleSystem = System()
     # Deep copy these forces into out new state
     copyOfForce = deepcopy(force)
     # Add this to our HarmonicBondSystem
     HarmonicAngleSystem.addForce(copyOfForce)

     # Add the particles?
     for j in range (system.getNumParticles()):
       HarmonicAngleSystem.addParticle(system.getParticleMass(j) )

     HarmonicAngleIntegrator = VerletIntegrator(1*femtoseconds)
     HarmonicAngleContext = Context(HarmonicAngleSystem, HarmonicAngleIntegrator, platform)
     HarmonicAngleContext.setPositions(inpcrd.positions)

     HarmonicAngleState = HarmonicAngleContext.getState( getEnergy=True)

     print (HarmonicAngleState.getPotentialEnergy().in_units_of(kilocalorie/mole))
     print str(8.2046173826501736) + " (JAMBER)"

   if isinstance( force, openmm.PeriodicTorsionForce ):
     print "Found " + str(force.getNumTorsions()) + " PeriodicTorsionForce terms"

     # New state with only Harmonic Bonds in
     PeriodicTorsionSystem = System()
     # Deep copy these forces into out new state
     copyOfForce = deepcopy(force)
     # Add this to our HarmonicBondSystem
     PeriodicTorsionSystem.addForce(copyOfForce)

     # Add the particles?
     for j in range (system.getNumParticles()):
       PeriodicTorsionSystem.addParticle(system.getParticleMass(j) )

     PeriodicTorsionIntegrator = VerletIntegrator(1*femtoseconds)
     PeriodicTorsionContext = Context(PeriodicTorsionSystem, PeriodicTorsionIntegrator, platform)
     PeriodicTorsionContext.setPositions(inpcrd.positions)

     PeriodicTorsionState = PeriodicTorsionContext.getState( getEnergy=True)

     print (PeriodicTorsionState.getPotentialEnergy().in_units_of(kilocalorie/mole))
     print str(11.365686464834990) + " (JAMBER)"

   if isinstance( force, openmm.NonbondedForce ):
     print "Found " + str(force.getNumParticles()) + " NonbondedForce terms"
     print "Found " + str(force.getNumExceptions()) + " NonbondedForce exception terms"


     # New state with only Harmonic Bonds in
     NonbondedSystem = System()
     # Deep copy these forces into out new state
     copyOfForce = deepcopy(force)
     # Add this to our HarmonicBondSystem
     NonbondedSystem.addForce(copyOfForce)

     # Add the particles?
     for j in range (system.getNumParticles()):
       NonbondedSystem.addParticle(system.getParticleMass(j) )

     NonbondedIntegrator = VerletIntegrator(1*femtoseconds)
     NonbondedContext = Context(NonbondedSystem, NonbondedIntegrator, platform)
     NonbondedContext.setPositions(inpcrd.positions)

     NonbondedState = NonbondedContext.getState( getEnergy=True)



     print "Normal pairs"
     for k in range (force.getNumParticles()):
       print force.getParticleParameters(k)


     print "Exception pairs"
     for k in range (force.getNumExceptions()):
        print force.getExceptionParameters(k)

     print (NonbondedState.getPotentialEnergy().in_units_of(kilocalorie/mole))
     
     JAMBER_Result = -79.994516268914055 + -1.0455001424914672 + 4.4614567778860472 + 48.744345933282325
     print str( JAMBER_Result  ) + " (JAMBER)"










#simulation.step(1)


# https://bitbucket.org/mjw99/jamber/src/c3b915fd3a56/src/test/resources/name/mjw/jamber/IO/AMBER/blockedAlanineDipeptide/mdinfo?at=master
# -72.165632 kcal/mol

# Refs
# Python API docs
# https://simtk.org/api_docs/openmm/api4_1/python/

