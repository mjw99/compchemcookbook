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
#platformProperties = {"OpenCLDeviceIndex":"0"}
#platformProperties = {"OpenCLDeviceIndex":"1"}
platformProperties = {"OpenCLDeviceIndex":"0,1,2"}
#platformProperties = {"OpenCLDeviceIndex":"1,2"}
platform = openmm.Platform_getPlatformByName("OpenCL")
print "Speed relative to reference is : " + str(platform.getSpeed())





prmtop = AmberPrmtopFile('prmtop7')
inpcrd = AmberInpcrdFile('inpcrd.equil.openmm',  loadVelocities=True, loadBoxVectors=True)

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)
integrator = VerletIntegrator(2*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)
print "Platform: %s" % (simulation.context.getPlatform().getName())

print "Number of atoms %i"      % len(inpcrd.positions)
print "Number of velocities %i" % len(inpcrd.velocities)


simulation.context.setPositions(inpcrd.positions)
simulation.context.setVelocities(inpcrd.velocities)

simulation.reporters.append(PDBReporter('output.pdb', 1000))
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

start_time = time.time()
simulation.step(10000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# 100 seconds to run 0.02 ns 
# Hence it will take 1/0.02  * 100s to run one nS.
totalDynmaicRunTimeInSeconds = time.time() - start_time

timeNeedToRunOneNsinSeconds = (1/0.02) * totalDynmaicRunTimeInSeconds

NsPerDay = 86400 / timeNeedToRunOneNsinSeconds


print str(totalDynmaicRunTimeInSeconds) + " seconds"
print str(timeNeedToRunOneNsinSeconds) + " is the time needed to run 1 ns"
print str(NsPerDay)  + " nS/day"

# Refs
# Python API docs
# https://simtk.org/api_docs/openmm/api4_1/python/main.html



# http://wiki.simtk.org/openmm/BenchmarkOpenMMDHRF
#  1xC2070 = 30.9 ns/day
#
# http://ambermd.org/gpus/benchmarks.htm
#  1xM2090 = 43.74 ns/day
#  




# Results
# AMBER11, M2090,  ecc on 			38.12 ns/day  (45.28s run time)
# AMBER11, M2090,  ecc off 			42.48 ns/day  (40.67s run time)


# OpenMM 4.1.1/OpenCL, M2090,  ecc on 		17.38 ns/day  (99.41 run time)
# OpenMM 4.1.1/OpenCL, M2090,  ecc off 		18.40 ns/day  (93.91 run time)

# OpenMM 4.1.1/OpenCL, C2075,  ecc on 		12.31 ns/day  (140.32 run time)
# OpenMM 4.1.1/OpenCL, 2xC2075,  ecc on 	19.24 ns/day  (89.81 run time)
# OpenMM 4.1.1/OpenCL, 1XM2090+2xC2075, ecc on	25.55 ns/day  (67.60 run time)

# OpenMM 4.1.1/CUDA, M2090,  ecc off		11.77 ns/day (146.69 run time)

