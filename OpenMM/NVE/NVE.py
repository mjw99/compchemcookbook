# ALAALA NVE example

from openmm.app import *
from openmm import *
from simtk.unit import *
from sys import stdout
import time

#platform = openmm.Platform.getPlatformByName("Reference")
platform = openmm.Platform.getPlatformByName("CUDA")
#platform = openmm.Platform.getPlatformByName("OpenCL")

# OpenCL precision
#platformProperties = {"OpenCLPrecision":"mixed"}
# CUDA precision
platformProperties = {"CudaPrecision":"mixed"}


prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd')

system = prmtop.createSystem(nonbondedMethod=NoCutoff, removeCMMotion=False)

#integrator = VerletIntegrator(1*femtoseconds)
integrator = VerletIntegrator(0.01*femtoseconds)
simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)

# Set the positions
simulation.context.setPositions(inpcrd.positions)

simulation.reporters.append(PDBReporter('output.pdb', 1))
#simulation.reporters.append(StateDataReporter(stdout, 1, step=True,  totalEnergy=True, kineticEnergy=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(StateDataReporter(stdout, 1,  totalEnergy=True))

simulation.step(10000)


