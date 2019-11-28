# Example illustrating use of CustomNonbondedForce to create VDW term

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

#print sigma
#print epsilon
#print sigma.in_units_of(nanometer)
#print epsilon.in_units_of(kilojoule/mole)

timestep = 1 * femtosecond # integrator timestep

cutoff = 999 * angstrom
#print "sigma = %s" % sigma
#print "cutoff = %s" % cutoff

# =============================================================================
# Build system
# =============================================================================


system = System()

topology = Topology()
newChain = topology.addChain()


###################################################################################

# 12-6 LJ
# Note, the Lorentz-Bertelot rules are being invoked here....
customNonbondedForceVDW = CustomNonbondedForce("4 * epsilon * (  (sigma/r)^12 - (sigma/r)^6 ) ;"
"sigma=0.5*(sigma1 + sigma2);"
"epsilon=sqrt(epsilon1*epsilon2)");


customNonbondedForceVDW.addPerParticleParameter("sigma")
customNonbondedForceVDW.addPerParticleParameter("epsilon")

customNonbondedForceVDW.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
customNonbondedForceVDW.setCutoffDistance(cutoff)





# Assign force types
for particle_index in range(nparticles):
  system.addParticle(mass)
  newResidue = topology.addResidue("UNK", newChain)

  if (particle_index == 0 ):
     # Add alchemically-modified particle.
     topology.addAtom("MODD", element.argon, newResidue)
     customNonbondedForceVDW.addParticle( [3.4*angstrom, 0.238*kilocalories_per_mole ]  )
     
  else:
     # Add normal particle
     topology.addAtom("Argo", element.argon, newResidue)
     customNonbondedForceVDW.addParticle( [3.4*angstrom, 0.238*kilocalories_per_mole ]  )

system.addForce(customNonbondedForceVDW)

print("System.getNumParticles() is ", system.getNumParticles())
print("system.getNumForces() is " , system.getNumForces())



# Create Integrator
integrator = VerletIntegrator(1*femtosecond)
simulation = Simulation(topology, system, integrator)

# Create initial positions.
positions = [  Vec3(0, 0, 0), Vec3(0, 0, 0.6) ]
simulation.context.setPositions(positions)

print("simulation.system.getNumForces() is ", simulation.system.getNumForces())

for i in range(system.getNumForces()):
    system.getForce(i).setForceGroup(i)

for i in range(system.getNumForces()):
    print(type(system.getForce(i)))
    state = simulation.context.getState(getEnergy=True, groups=1<<i)
    print(state.getPotentialEnergy())


print(simulation.context.getState(getEnergy=True).getPotentialEnergy())




