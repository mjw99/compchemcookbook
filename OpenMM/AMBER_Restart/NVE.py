# DHRF Benchmark using OpenMM

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time
from AmberRestrtFile import AmberRestrtFile

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
inpcrd = AmberInpcrdFile('inpcrd', loadVelocities=True)

#system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)
#system = prmtop.createSystem(nonbondedMethod=NoCutoff, implicitSolvent=OBC2)
system = prmtop.createSystem(nonbondedMethod=NoCutoff)
# Remember, this is being run NVE
integrator = VerletIntegrator(1*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)
print "Platform: %s" % (simulation.context.getPlatform().getName())

print "Number of atoms %i"      % len(inpcrd.positions)
simulation.context.setPositions(inpcrd.positions)
simulation.context.setVelocities(inpcrd.velocities)

simulation.reporters.append(PDBReporter('output.pdb', 1))
simulation.reporters.append(StateDataReporter(stdout, 1, step=True,  totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True))


#simulation.step(1)

# Write the AMBER restart
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
# TODO box vectors
# TODO time? Needs to be parsed out to picoseconds
newInpcrd = AmberRestrtFile('restrt', positions, velocities)



# Refs
# Python API docs
# https://simtk.org/api_docs/openmm/api4_1/python/

