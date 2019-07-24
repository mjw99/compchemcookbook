# Example VDW calculation

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from array import *


# =============================================================================
# Specify simulation parameters
# =============================================================================

nparticles = 2 # number of particles

mass = 39.9 * amu # mass
sigma = 3.4 * angstrom # Lennard-Jones sigma
epsilon = 0.238 * kilocalories_per_mole # Lennard-Jones well-depth
timestep = 1 * femtosecond # integrator timestep

cutoff = 999 * angstrom
print("sigma = ", sigma)
print("cutoff = ", cutoff)

# =============================================================================
# Build system
# =============================================================================

# Create argon system where first particle is alchemically modified by lambda_value.
lambda_value = 1.0

system = System()

topology = Topology()
newChain = topology.addChain()



nonbondedForce = NonbondedForce()
nonbondedForce.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
nonbondedForce.setCutoffDistance(cutoff)

# Assign force types
for particle_index in range(nparticles):
  system.addParticle(mass)
  newResidue = topology.addResidue("UNK", newChain)

  if (particle_index == 0 ):
     # Add alchemically-modified particle.
     topology.addAtom("MODD", element.argon, newResidue)
     nonbondedForce.addParticle(0, sigma, epsilon)
     
  else:
     # Add normal particle
     topology.addAtom("Argo", element.argon, newResidue)
     nonbondedForce.addParticle(0, sigma, epsilon)

system.addForce(nonbondedForce)

print("System.getNumParticles() is ", system.getNumParticles())
print("system.getNumForces() is ", system.getNumForces())



# Create Integrator
integrator = VerletIntegrator(1*femtosecond)
simulation = Simulation(topology, system, integrator)

# Create initial positions.
positions = [  Vec3(0, 0, 0), Vec3(0, 0, 0.6) ]
simulation.context.setPositions(positions)

print("simulation.system.getNumForces() is ", simulation.system.getNumForces())

print(simulation.context.getState(getEnergy=True).getPotentialEnergy())



