# DHRF Benchmark using OpenMM

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time


# Reference
#platform=openmm.Platform_getPlatform(0)

# CUDA
#platform=openmm.Platform_getPlatform(1)

# OpenCL
platform=openmm.Platform_getPlatform(2)




prmtop = AmberPrmtopFile('prmtop')
#inpcrd = AmberInpcrdFile('inpcrd',  loadVelocities=True, loadBoxVectors=True)
inpcrd = AmberInpcrdFile('inpcrd',  loadBoxVectors=True)

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=8*angstrom, constraints=HBonds )

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform)
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
# Hence it will take 1/0.02  * 100s to run one nS.
totalDynmaicRunTimeInSeconds = time.time() - start_time

timeNeedToRunOneNsinSeconds = (1/0.02) * totalDynmaicRunTimeInSeconds

NsPerDay = 86400 / timeNeedToRunOneNsinSeconds


print str(totalDynmaicRunTimeInSeconds) + " seconds"
print str(timeNeedToRunOneNsinSeconds) + " is the time needed to run 1 ns"
print str(NsPerDay)  + " ns/day"

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
# AMBER11 ecc on 			38.12 ns/day  (45.28s run time)
# AMBER11 ecc off 			42.48 ns/day  (40.67s run time)

# OpenMM 4.1.1/OpenCL ecc on 		17.38 ns/day  (99.41)
# OpenMM 4.1.1/OpenCL ecc off 		18.40 ns/day  (93.91 run time)

# OpenMM 4.1.1/CUDA ecc off		11.77 ns/day (146.69 run time)

