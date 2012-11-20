# DHRF Benchmark using OpenMM

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time


# Run on multiple cards
# 0  Tesla M2090
# 1  Tesla C2075
# 2  Tesla C2075
platformProperties = {"OpenCLDeviceIndex":"0"}
#platformProperties = {"OpenCLDeviceIndex":"1,2"}
#platformProperties = {"OpenCLDeviceIndex":"0,1,2"}
platform = openmm.Platform_getPlatformByName("OpenCL")
print "Speed relative to reference is : " + str(platform.getSpeed())




prmtop = AmberPrmtopFile('prmtop')
#inpcrd = AmberInpcrdFile('inpcrd',  loadVelocities=True, loadBoxVectors=True)
inpcrd = AmberInpcrdFile('inpcrd',  loadBoxVectors=True)

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=8*angstrom, constraints=HBonds )

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties )
print "Platform: %s" % (simulation.context.getPlatform().getName())

print "Number of atoms %i"      % len(inpcrd.positions)
#print inpcrd.positions
#print "Number of velocities %i" % len(inpcrd.velocities)


simulation.context.setPositions(inpcrd.positions)
#simulation.context.setVelocities(inpcrd.velocities)

LocalEnergyMinimizer.minimize(simulation.context)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('heating.pdb', 1000))

start_time = time.time()
simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# 100 seconds to run 0.02 ns 
# Hence it will take 1/0.02  * 100s to run one ns.
totalDynmaicRunTimeInSeconds = time.time() - start_time

timeNeedToRunOneNsinSeconds = (1/0.02) * totalDynmaicRunTimeInSeconds

NsPerDay = 86400 / timeNeedToRunOneNsinSeconds


print str(totalDynmaicRunTimeInSeconds) + " seconds"
print str(timeNeedToRunOneNsinSeconds) + " is the time needed to run 1 ns"
print str(NsPerDay)  + " ns/day"

# Refs
# Python API docs
# https://simtk.org/api_docs/openmm/api4_1/python/main.html



