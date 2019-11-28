# Example force printing of an AMBER sytem using OpenMM
# Mark J. Williamson

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time

platform = openmm.Platform_getPlatformByName("Reference")




prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd')

system = prmtop.createSystem(nonbondedMethod=NoCutoff)
integrator = VerletIntegrator(1*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(inpcrd.positions)



print("Platform: %s" % (simulation.context.getPlatform().getName()))
print("Number of atoms %i"      % len(inpcrd.positions))




system = simulation.context.getSystem()
state = simulation.context.getState( getEnergy=True)
forces = [None]
forces = simulation.context.getState(getForces=True).getForces()


print("Potential energy is " +  str(state.getPotentialEnergy().in_units_of(kilocalorie/mole)))

for force in forces:
     print(force.in_units_of(kilocalorie/mole/angstrom))


