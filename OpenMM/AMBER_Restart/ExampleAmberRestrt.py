# An example use of saving a state with MJW's
# AmberRestrtFile object

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time
from amberrestrtfile import AmberRestrtFile


## Select which platform to use:
#platform = openmm.Platform_getPlatformByName("OpenCL")
#platform = openmm.Platform_getPlatformByName("Cuda")
platform = openmm.Platform_getPlatformByName("Reference")




## Select which card(s) to run on if using OpenCL
# Run on multiple cards
# 0  Tesla M2090
# 1  Tesla C2075
# 2  Tesla C2075

platformProperties = {"OpenCLDeviceIndex":"0"}
#platformProperties = {"OpenCLDeviceIndex":"1"}
#platformProperties = {"OpenCLDeviceIndex":"0,1,2"}
#platformProperties = {"OpenCLDeviceIndex":"1,2"}



# Load the AMBER topology and coordinate files
# Also, read the velocities from the inpcrd file as well
prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd', loadVelocities=True)


# Typical PBC system with PME
#system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)

# Explicit solvent model
#system = prmtop.createSystem(nonbondedMethod=NoCutoff, implicitSolvent=OBC2)

# No PBC, infinite cut off
system = prmtop.createSystem(nonbondedMethod=NoCutoff)

# Remember, this is being run NVE
integrator = VerletIntegrator(1*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)

print("Platform: ", (simulation.context.getPlatform().getName()))

print("Number of atoms ",   len(inpcrd.positions))

# Set the positions and velocities from the AMBER files
simulation.context.setPositions(inpcrd.positions)
simulation.context.setVelocities(inpcrd.velocities)


# We will skip this since, I want to check that I get the same
# result as I read in
#simulation.step(1)

# Write the AMBER restart
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
boxVectors = simulation.context.getState().getPeriodicBoxVectors()
time = simulation.context.getState().getTime()

restrt = AmberRestrtFile('restrt', positions, velocities, boxVectors, time)



# Refs
# Python API docs
# https://simtk.org/api_docs/openmm/api4_1/python/

