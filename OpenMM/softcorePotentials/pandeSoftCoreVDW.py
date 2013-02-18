# Example illustrating use of CustomNonbondedForce to lambda map out an Atom
# Based upon argon-chemical-potential.py

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
print "sigma = %s" % sigma
print "cutoff = %s" % cutoff

# =============================================================================
# Build system
# =============================================================================

# Create argon system where first particle is alchemically modified by lambda_value.
lambda_value = 1.0

system = System()

topology = Topology()
newChain = topology.addChain()


###################################################################################

# Softcore VdW; see page E in http://dx.doi.org/10.1021/ct300857j
# which is actually the VDW from http://dx.doi.org/10.1063/1.1877132 Equ. 4
pandeSoftCoreVDW = CustomNonbondedForce("4*epsilon*l12*( 1/( (alphaLJ*(1-l12) + (r/sigma)^6)^2) - 1/( alphaLJ*(1-l12) + (r/sigma)^6) ) ;"
"sigma=0.5*(sigma1*sigma2);"
"epsilon=sqrt(epsilon1*epsilon2);"
"alphaLJ=0.5;"
"l12=1-(1-lambda)*step(useLambda1+useLambda2-0.5)");

# Note, the Lorentz-Bertelot rules are being invoked here....

pandeSoftCoreVDW.addPerParticleParameter("sigma")
pandeSoftCoreVDW.addPerParticleParameter("epsilon")
pandeSoftCoreVDW.addPerParticleParameter("useLambda")

# "useLamba" is a per particle parameter that is a function of the global parameter, lambda
# This enables a subset of particles to be affected by the lambda value.

pandeSoftCoreVDW.addGlobalParameter("lambda", 1.0)

pandeSoftCoreVDW.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
pandeSoftCoreVDW.setCutoffDistance(cutoff)



# Assign force types
for particle_index in range(nparticles):
  system.addParticle(mass)
  newResidue = topology.addResidue("UNK", newChain)

  if (particle_index == 0 ):
     # Add alchemically-modified particle.
     topology.addAtom("MODD", element.argon, newResidue)
     pandeSoftCoreVDW.addParticle( [3.4*angstrom, 0.238*kilocalories, 1 ]  )
     
  else:
     # Add normal particle
     topology.addAtom("Argo", element.argon, newResidue)
     pandeSoftCoreVDW.addParticle( [3.4*angstrom, 0.238*kilocalories, 0 ]  )

system.addForce(pandeSoftCoreVDW)

print "System.getNumParticles() is %i"  % system.getNumParticles()
print "system.getNumForces() is %i " % system.getNumForces()



# Create Integrator
integrator = VerletIntegrator(1*femtosecond)
simulation = Simulation(topology, system, integrator)

# Create initial positions.
positions = [  Vec3(0, 0, 0), Vec3(0, 0, 0.1) ]
simulation.context.setPositions(positions)

print "simulation.system.getNumForces() is %i " % simulation.system.getNumForces()

# Run dynamics.
simulation.reporters.append( PDBReporter('output.pdb', 100) )
simulation.reporters.append( StateDataReporter(stdout, 100, step=True, potentialEnergy=True) )


# =============================================================================
# Lambda map in the Argon atom
# =============================================================================

# This is NVE; kinetic energy has come from the slight perturbation of the interatomic distance
for i in range(10):
  #1, 0.9, 0.8
  simulation.context.setParameter("lambda",  1 - (i*0.1)  )
  #simulation.context.setParameter("lambda",  1 )
  print "Current lambda value is " + str(simulation.context.getParameter("lambda"))
  simulation.step(10000)


