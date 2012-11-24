# DHRF Benchmark using OpenMM

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time

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
print "Speed relative to reference is : " + str(platform.getSpeed())



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

print "Potential energy is " +  str(state.getPotentialEnergy().in_units_of(kilocalorie/mole))

for i in range(system.getNumForces()):
   force = system.getForce(i)
   if isinstance( force, openmm.HarmonicBondForce ):
     print ""
     print "Found " + str(force.getNumBonds()) + " HarmonicBondForce terms"
     for j in range(force.getNumBonds()):
       print force.getBondParameters(j)

   if isinstance( force, openmm.HarmonicAngleForce ):
     print ""
     print "Found " + str(force.getNumAngles()) + " HarmonicAngleForce terms"
     for j in range(force.getNumAngles()):
       print force.getAngleParameters(j)


   if isinstance( force, openmm.PeriodicTorsionForce ):
     print ""
     print "Found " + str(force.getNumTorsions()) + " PeriodicTorsionForce terms"
     for j in range(force.getNumTorsions()):
       print force.getTorsionParameters(j)


   if isinstance( force, openmm.NonbondedForce ):
     print ""
     print "Found " + str(force.getNumParticles()) + " NonbondedForce terms"
     for j in range(force.getNumParticles()):
       print force.getParticleParameters(j)

   if isinstance( force, openmm.GBSAOBCForce ):
     print force 



#simulation.step(1)


# https://bitbucket.org/mjw99/jamber/src/c3b915fd3a56/src/test/resources/name/mjw/jamber/IO/AMBER/blockedAlanineDipeptide/mdinfo?at=master
# -72.165632 kcal/mol

# Refs
# Python API docs
# https://simtk.org/api_docs/openmm/api4_1/python/

