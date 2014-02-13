# Example force printing for a dihedral term
# Mark J. Williamson

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time

platform = openmm.Platform_getPlatformByName("Reference")



system = System()
topology = Topology()
newChain = topology.addChain()
newResidue = topology.addResidue("UNK", newChain)

system.addParticle(12.0)
topology.addAtom("Carbon", element.carbon, newResidue)
system.addParticle(12.0)
topology.addAtom("Carbon", element.carbon, newResidue)
system.addParticle(12.0)
topology.addAtom("Carbon", element.carbon, newResidue)
system.addParticle(12.0)
topology.addAtom("Carbon", element.carbon, newResidue)


dihedral = PeriodicTorsionForce()

#addTorsion(PeriodicTorsionForce self, int particle1, int particle2, int particle3, int particle4, int periodicity, double phase, double k) 
dihedral.addTorsion(0, 1, 2, 3, 4, 2, 10*kilocalorie/mole)
system.addForce(dihedral)


integrator = VerletIntegrator(1*femtoseconds)

simulation = Simulation(topology, system, integrator, platform)
#simulation.context.setPositions( [Vec3(0,0,0)*angstrom, Vec3(1,1,0)*angstrom, Vec3(2,1,0)*angstrom, Vec3(3,2,0)*angstrom]  )
simulation.context.setPositions( [Vec3(30.958572, 73.709137, 38.031029)*angstrom, Vec3(31.924915, 72.577698, 37.985279)*angstrom, Vec3(31.643818, 71.563255, 36.873047)*angstrom, Vec3(31.253489, 70.131218, 37.371773)*angstrom]  )




print "Platform: %s" % (simulation.context.getPlatform().getName())



system = simulation.context.getSystem()
state = simulation.context.getState( getEnergy=True)
forces = [None]
forces = simulation.context.getState(getForces=True).getForces()


print "Potential energy is " +  str(state.getPotentialEnergy().in_units_of(kilocalorie/mole))

for i in range(system.getNumForces()):
    system.getForce(i).setForceGroup(i)

for i in range(system.getNumForces()):
    print type(system.getForce(i))
    state = simulation.context.getState(getEnergy=True, groups=1<<i)
    print state.getPotentialEnergy().in_units_of(kilocalorie/mole)

for force in forces:
     print force.in_units_of(kilocalorie/mole/angstrom)


