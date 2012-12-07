# DHRF Benchmark using OpenMM

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time
from ambernetcdfreporter import AmberNetCDFReporter

platform = openmm.Platform_getPlatformByName("OpenCL")
#platform = openmm.Platform_getPlatformByName("Cuda")
#platform = openmm.Platform_getPlatformByName("Reference")



# Run on multiple cards
# 0  Tesla M2090
# 1  Tesla C2075
# 2  Tesla C2075
platformProperties = {"OpenCLDeviceIndex":"0"}
#platformProperties = {"OpenCLDeviceIndex":"1"}
#platformProperties = {"OpenCLDeviceIndex":"0,1,2"}
#platformProperties = {"OpenCLDeviceIndex":"1,2"}
print "Speed relative to reference is : " + str(platform.getSpeed())





prmtop = AmberPrmtopFile('prmtop7')
inpcrd = AmberInpcrdFile('inpcrd.equil.openmm',  loadVelocities=True, loadBoxVectors=True)

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.8*nanometer, constraints=HBonds)
# Remember, this is being run NVE
integrator = VerletIntegrator(2*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)
print "Platform: %s" % (simulation.context.getPlatform().getName())

print "Number of atoms %i"      % len(inpcrd.positions)
print "Number of velocities %i" % len(inpcrd.velocities)

simulation.context.setPositions(inpcrd.positions)
simulation.context.setVelocities(inpcrd.velocities)

simulation.reporters.append(AmberNetCDFReporter('output.nc', 1))
simulation.reporters.append(StateDataReporter(stdout, 1, step=True,  totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True))

start_time = time.time()
simulation.step(3) # i.e. 20,000 fs == 20 ps == 0.02 ns
