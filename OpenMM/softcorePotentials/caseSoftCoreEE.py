# Example illustrating use of CustomNonbondedForce to lambda map out an Atom

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from array import *



platform = openmm.Platform_getPlatformByName("Reference")


# =============================================================================
# Specify simulation parameters
# =============================================================================

nparticles = 2 # number of particles
mass = 1 * amu # mass
timestep = 1 * femtosecond # integrator timestep
charge = 1 * elementary_charge


cutoff = 999 * angstrom
print "cutoff = %s" % cutoff

# =============================================================================
# Build system
# =============================================================================

# Create argon system where first particle is alchemically modified by lambda_value.
lambda_value = 1.0

system = System()

topology = Topology()
newChain = topology.addChain()


####################################################################################

# caseSoftCoreEE; see Eq. 5 in http://dx.doi.org/10.1002/jcc.21909
caseSoftCoreEE = CustomNonbondedForce("(1-l12) * ((q1 * q2) / ( 4 * PI * EPSILON0  * ((l12 * BETAC) + r^M)^(1/M) ))  ;"
"PI=3.141592653589;"
"EPSILON0=5.72765E-4;"
"M=6;"
"BETAC=2.5;"
"l12=1-(1-lambda)*step(useLambda1+useLambda2-0.5)");

#For e0, see EPSILON0 in ./SimTKUtilities/SimTKOpenMMRealType.h

# For values of M and B, see page 3260 of 10.1002/jcc.21909


caseSoftCoreEE.addPerParticleParameter("q")
caseSoftCoreEE.addPerParticleParameter("useLambda")

caseSoftCoreEE.addGlobalParameter("lambda", 1.0)

caseSoftCoreEE.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
caseSoftCoreEE.setCutoffDistance(cutoff)


######################################################################################

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
     caseSoftCoreEE.addParticle([charge, 1])
     #nonbondedForce.addParticle(charge, 0 , 0)
     
  else:
     # Add normal particle
     topology.addAtom("Argo", element.argon, newResidue)
     caseSoftCoreEE.addParticle([charge, 0])
     #nonbondedForce.addParticle(charge, 0 , 0)

system.addForce(caseSoftCoreEE)
#system.addForce(nonbondedForce)

print "system.getNumParticles() is %i"  % system.getNumParticles()
print "system.getNumForces() is %i " % system.getNumForces()



# Create Integrator
integrator = VerletIntegrator(1*femtosecond)
simulation = Simulation(topology, system, integrator)

# Create initial positions; 1 A apart
positions = [  Vec3(0, 0, 0), Vec3(0, 0, 0.1) ]
simulation.context.setPositions(positions)

print "simulation.system.getNumForces() is %i " % simulation.system.getNumForces()



# Should be 1387.28173828 kJ/mol for:
#
# nonbondedForce = NonbondedForce()
# nonbondedForce.setNonbondedMethod(NonbondedForce.CutoffNonPeriodic)
# nonbondedForce.setCutoffDistance(999 * angstrom)

simulation.context.setParameter("lambda",  1 )
print "Current lambda value is " + str(simulation.context.getParameter("lambda"))
print simulation.context.getState(getEnergy=True).getPotentialEnergy()



