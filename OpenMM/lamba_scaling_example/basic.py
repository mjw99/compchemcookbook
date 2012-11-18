# Example illustrating use of CustomNonbondedForce to lambda map in an Atom
# Based upon argon-chemical-potential.py

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import numpy
from array import *


# =============================================================================
# Specify simulation parameters
# =============================================================================

nparticles = 216 # number of particles

mass = 39.9 * amu # mass
sigma = 3.4 * angstrom # Lennard-Jones sigma
epsilon = 0.238 * kilocalories_per_mole # Lennard-Jones well-depth
charge = 0.0 * elementary_charge # argon model has no charge

nequil_steps = 5000 # number of dynamics steps for equilibration
nprod_steps = 50000 # number of dynamics steps per production iteration

reduced_density = 0.85 # reduced density rho*sigma^3
temperature = 300 * kelvin # temperature
pressure = 1 * atmosphere # pressure
collision_rate = 5.0 / picosecond # collision rate for Langevin thermostat
timestep = 1 * femtosecond # integrator timestep


# =============================================================================
# Compute box size.
# =============================================================================

volume = nparticles*(sigma**3)/reduced_density
box_edge = volume**(1.0/3.0)
cutoff = min(box_edge*0.49, 2.5*sigma) # Compute cutoff
print "sigma = %s" % sigma
print "box_edge = %s" % box_edge
print "cutoff = %s" % cutoff

# =============================================================================
# Build system
# =============================================================================

# Create argon system where first particle is alchemically modified by lambda_value.
# It is currently zero
lambda_value = 0.0


system = System()
system.setDefaultPeriodicBoxVectors(Vec3(box_edge, 0, 0), Vec3(0, box_edge, 0), Vec3(0, 0, box_edge))

topology = Topology()
newChain = topology.addChain()


# Normal nonbonded
#nonBondedForce = NonbondedForce()
#nonBondedForce.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
#nonBondedForce.setCutoffDistance(cutoff)


# Modified VdW
customNonBondedForce = CustomNonbondedForce("4*epsilon*l12*( 1/( (alphaLJ*(1-l12) + (r/sigma)^6)^2) - 1/( alphaLJ*(1-l12) + (r/sigma)^6) ) ;sigma=0.5*(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2); alphaLJ=0.5; l12=1-(1-lambda)*step(useLambda1+useLambda2-0.5)");

customNonBondedForce.addPerParticleParameter("sigma")
customNonBondedForce.addPerParticleParameter("epsilon")

customNonBondedForce.addGlobalParameter("lambda", 1.0)
customNonBondedForce.addPerParticleParameter("useLambda")

customNonBondedForce.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
customNonBondedForce.setCutoffDistance(cutoff)


# Assign force types
for particle_index in range(nparticles):
  system.addParticle(mass)
  newResidue = topology.addResidue("UNK", newChain)

  if (particle_index == 0 ):
     # Add alchemically-modified particle.
     topology.addAtom("MODD", "Ar", newResidue)
     # Is there a better way to do this?
     customNonBondedForce.addParticle( array('d',[3.4*NmPerAngstrom, 0.238*KJPerKcal, 1 ])  )
     
  else:
     # Add normal particle.
     topology.addAtom("Argo", "ar", newResidue)
     # Is there a better way to do this?
     customNonBondedForce.addParticle( array('d',[3.4*NmPerAngstrom, 0.238*KJPerKcal, 0 ])  )

system.addForce(customNonBondedForce)

print "System.getNumParticles() is %i"  % system.getNumParticles()
print "system.getNumForces() is %i " % system.getNumForces()


# Create random initial positions.
import numpy.random
positions = Quantity(numpy.random.uniform(high=box_edge/angstroms, size=[nparticles,3]), angstrom)


# Create Integrator
integrator = LangevinIntegrator(temperature, collision_rate, timestep)
simulation = Simulation(topology, system, integrator)

# Initiate from Random set of positions
simulation.context.setPositions(positions)

# Minimize energy from coordinates.
simulation.minimizeEnergy()

print "simulation.system.getNumForces() is %i " % simulation.system.getNumForces()

# Equilibrate.
simulation.step(nequil_steps)

# Run dynamics.
simulation.reporters.append( PDBReporter('output.pdb', 10000) )
simulation.reporters.append( StateDataReporter(stdout, 10000, step=True, potentialEnergy=True, temperature=True) )



# =============================================================================
# Lambda map in the Argon atom
# =============================================================================

for i in range(10):
  print i
  simulation.context.setParameter("lambda",  i*0.1  )
  print "Current lambda value is " + str(simulation.context.getParameter("lambda"))
  simulation.step(10000)


