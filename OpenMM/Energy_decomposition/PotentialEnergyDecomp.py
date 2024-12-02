# Obtain per term energy contributions in OpenMM.

from openmm.app import *
from openmm import *
from simtk.unit import *
from sys import stdout

platform = openmm.Platform.getPlatformByName("Reference")

prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd')

system = prmtop.createSystem(nonbondedMethod=NoCutoff)
integrator = VerletIntegrator(1*femtoseconds)
simulation = Simulation(prmtop.topology, system, integrator, platform)
simulation.context.setPositions(inpcrd.positions)

print("Platform: ", (simulation.context.getPlatform().getName()))
print("Number of atoms ", len(inpcrd.positions))

print("")
# Entire system
state = simulation.context.getState( getEnergy=True)
print("Total potential energy is ", str(state.getPotentialEnergy().in_units_of(kilocalorie/mole)))
print("")

for i in range(system.getNumForces()):
    system.getForce(i).setForceGroup(i)

for i in range(system.getNumForces()):
    print(type(system.getForce(i)))
    state = simulation.context.getState(getEnergy=True, groups=1<<i)
    print(state.getPotentialEnergy().in_units_of(kilocalorie/mole))


print("")
# Harmonic Bond
print( str(5.6366273691689504), " kcal/mol (JAMBER)")
# Harmonic Angle
print( str(8.2046173826501736), " kcal/mol (JAMBER)")
# Torsions
print( str(11.365686464834990), " kcal/mol (JAMBER)")
# NB
JAMBER_Result = -79.994516268914055 + -1.0455001424914672 + 4.4614567778860472 + 48.744345933282325
print(str( JAMBER_Result  ) , " kcal/mol (JAMBER)")
