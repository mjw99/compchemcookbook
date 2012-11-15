#
# Typical AMBER MD protocol 

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




###################
# 1) Minimisation #
###################
prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd',  loadBoxVectors=True)

print "Number of atoms %i"      % len(inpcrd.positions)

integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=8*angstrom )

simulation = Simulation(prmtop.topology, system, integrator, platform)

print "Platform: %s" % (simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
#simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True)


simulation.reporters.append(PDBReporter('min.pdb', 1))
# Minimise system with L-BFGS algorithm
LocalEnergyMinimizer.minimize(simulation.context,)

# Save context?
#for x in state.getPostions().value_in_unit(angstrom):
#  print x 

##############
# 2) Heating #
##############
#prmtop = AmberPrmtopFile('prmtop')
#inpcrd = AmberInpcrdFile('inpcrd.minized',  loadVelocities=True, loadBoxVectors=True)

#print "Number of atoms %i"      % len(inpcrd.positions)
#print "Number of velocities %i" % len(inpcrd.velocities)

#integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
#simulation = Simulation(prmtop.topology, system, integrator, platform)

#simulation.context.setPositions(inpcrd.positions)
#simulation.context.setVelocities(inpcrd.velocities)

#print "Platform: %s" % (simulation.context.getPlatform().getName())

#simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
#simulation.reporters.append(PDBReporter('heating.pdb', 1000))

#simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns


#Save context?


##########################
# NTP Density correction #
##########################

#prmtop = AmberPrmtopFile('prmtop')
#inpcrd = AmberInpcrdFile('inpcrd.minized',  loadVelocities=True, loadBoxVectors=True)

#print "Number of atoms %i"      % len(inpcrd.positions)
#print "Number of velocities %i" % len(inpcrd.velocities)

#integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
# Barostat


















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

