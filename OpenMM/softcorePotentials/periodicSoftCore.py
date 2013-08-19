# WIP Example use of OpenMM 5.1 to carry out a TI (like) calculation
# with bespoke softcore potentials.
# 
# This example uses water and copies a system over to a new system,
# editing the nonbonded forces during the process.
#
# It is testing that:
#   1) The NonbondedForce of the original system can be disected and 
#      replaced with a new NonbondedForce in the new system.
#
#   2) The NonbondedForce of the original system can be disected and 
#      replaced a CustomNonbondedForce that has the identical form
#      of NonbondedForce.
#
# It is currently demonstrating that the correct denstiy of water
# cannot be achieved without PME or RF.
#
# With the CustomNonbondedForce.CutoffPeriodic a decrease in density from 
# ~0.6 to 0.15 g/mL is seen in the "Density correction under NPT" phase.
# With NonbondedForce.CutoffPeriodic, the density converges to ~0.99 g/mL in
# the "Density correction under NPT" phase. Special care must be taken with
# the setting of rfDielectric (and enabling) when copying the system across.
# 
#
# This difference is because NonbondedForce.CutoffPeriodic has a hidden RF term, 
# whereas CustomNonbondedForce.CutoffPeriodic does not.
#
# With NonbondedForce.PME, density converges to 0.99 g/mL.
#
# CustomNonbondedForce.PME does not exist :(
#
# Careful; two things are happening here:
#       i) VDW has long range dispersion correction added i.e. 8*PI*N^2*V....etc
#		--> setUseDispersionCorrection(bool useCorrection)
#          This is defaults to true for a system from modeller
#
#      ii) Coulomb has a Reaction Field (RF) applied; everything beyond cutoff is 
#          treated with a uniform field.
#		--> setReactionFieldDielectric(double dielectric)
#          This defaults to 78.3 for a system from modeller
#
# Mark J. Williamson

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from sys import stdout
from copy import deepcopy


#useCustomNonbondedForce = 1
useCustomNonbondedForce = 0


platformProperties = {}

# Set up platform
#platform = openmm.Platform_getPlatformByName("Reference")

platform = openmm.Platform_getPlatformByName("CUDA")
#platform = openmm.Platform_getPlatformByName("OpenCL")
platformProperties = {}
# CUDA precision
platformProperties['CudaPrecision'] = 'mixed'
# CUDA parallel
platformProperties['CudaDeviceIndex'] = '0'





print "Building system"
forceField = ForceField('tip3p.xml');
watPdb = PDBFile('WAT.pdb')
modeller = Modeller(watPdb.topology, watPdb.positions)


modeller.addSolvent(forceField, model='tip3p', padding=8*angstrom)
# Dump modeller structure
PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(), open('modeller.pdb', 'w'))




# Note


if useCustomNonbondedForce:

  # Now, define a customNonbondedForce
  # Recall, it asssumes r = Units of nm
  combinedVDWEE = CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6)+138.935456*q/r;"
  "q=q1*q2;" 
  "sigma=0.5*(sig1+sig2);" 
  "eps=sqrt(eps1*eps2)");
  combinedVDWEE.addPerParticleParameter("eps")
  combinedVDWEE.addPerParticleParameter("q")
  combinedVDWEE.addPerParticleParameter("sig")
  #combinedVDWEE.setNonbondedMethod(CustomNonbondedForce.CutoffNonPeriodic)
  combinedVDWEE.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
  #combinedVDWEE.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
  combinedVDWEE.setCutoffDistance(8*angstrom)
  # 5.2 ish only, but only VDW
  #combinedVDWEE.setUseLongRangeCorrection(True)

else:
  nonbondedForce = NonbondedForce();
  #nonbondedForce.setNonbondedMethod(NonbondedForce.NoCutoff)
  nonbondedForce.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
  #nonbondedForce.setNonbondedMethod(NonbondedForce.PME)
  nonbondedForce.setCutoffDistance(8*angstrom)

  # These are important for nonbondedForce, since they are not copied across or set
  #nonbondedForce.setUseDispersionCorrection(True)
  nonbondedForce.setUseDispersionCorrection(False)
  #nonbondedForce.setReactionFieldDielectric(78.3)
  nonbondedForce.setReactionFieldDielectric(0)




# Create the original system
#system = forceField.createSystem(modeller.topology, nonbondedMethod=NoCutoff, removeCMMotion=False)
system = forceField.createSystem(modeller.topology, nonbondedMethod=CutoffPeriodic, nonbondedCutoff=8*angstrom, removeCMMotion=False)
#system = forceField.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=8*angstrom, removeCMMotion=False)


#Energy
integrator = VerletIntegrator(1*femtosecond)

simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(modeller.positions)

state = simulation.context.getState( getEnergy=True)
print "Original system Total potential energy is " +  str(state.getPotentialEnergy().in_units_of(kilocalorie/mole))

for i in range(system.getNumForces()):
    system.getForce(i).setForceGroup(i)

for i in range(system.getNumForces()):
    print type(system.getForce(i))
    state = simulation.context.getState(getEnergy=True, groups=1<<i)
    print state.getPotentialEnergy().in_units_of(kilocalorie/mole)




# Serialise it
f = open('original_system.xml','w')
f.write(XmlSerializer.serializeSystem(system))
f.close()


# Create new System which will have NonbondedForce terms replaced with CustomNonbondedForce
newSystem = System()

# Copy across particles
for j in range (system.getNumParticles()):
  newSystem.addParticle( system.getParticleMass(j) )

# Copy across constraints
for j in range (system.getNumConstraints()):
  particle1, particle2, distance = system.getConstraintParameters(j)
  newSystem.addConstraint(particle1, particle2, distance)

# Copy across Box vectors
a, b, c = system.getDefaultPeriodicBoxVectors()
newSystem.setDefaultPeriodicBoxVectors(a,b,c)


# Decompose and copy across respective force terms
for i in range(system.getNumForces()):

   force = system.getForce(i)

   # We want to substitute NonbondedForce for caseSoftCoreVDW & caseSoftCoreEE
   if isinstance( force, openmm.NonbondedForce ):

     for j in range(force.getNumParticles()):
       charge, sigma, epsilon =  force.getParticleParameters(j)
       if useCustomNonbondedForce:
         combinedVDWEE.addParticle( [epsilon,charge,sigma ]  )
       else:
         nonbondedForce.addParticle( charge,sigma,epsilon  )


     for k in range(force.getNumExceptions()):
       particle1, particle2, chargeProd, sigma, epsilon  = force.getExceptionParameters(k)

       # These are probably bonded; hence exclude
       if ( (chargeProd == 0.0*elementary_charge * elementary_charge) & (epsilon == 0.0*kilojoule/mole) ):
         if useCustomNonbondedForce:
           combinedVDWEE.addExclusion(particle1, particle2)
         else:
           nonbondedForce.addException(particle1, particle2, chargeProd, sigma, epsilon)

     if useCustomNonbondedForce:
       newSystem.addForce(combinedVDWEE)
     else:
       newSystem.addForce(nonbondedForce)



#Serialise it
f = open('new_system.xml','w')
f.write(XmlSerializer.serializeSystem(newSystem))
f.close()



######################
# 1) Minimisation    #
######################
integrator = VerletIntegrator(1*femtosecond)
simulation = Simulation(modeller.topology, newSystem, integrator, platform, platformProperties)
simulation.context.setPositions(modeller.positions)

state = simulation.context.getState( getEnergy=True)
print "Copied system Total potential energy is " +  str(state.getPotentialEnergy().in_units_of(kilocalorie/mole))
#print "Copied system r is " +  str(state.getPotentialEnergy());


print "Minimising system"
simulation.minimizeEnergy(maxIterations=1000)

# Saving minimised positions
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))



################################
# 2) Thermalisation under NVT  #
################################
print "Heating system under NVT"
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtoseconds)

simulation = Simulation(modeller.topology, newSystem, integrator, platform, platformProperties)
simulation.context.setPositions(positions)


simulation.reporters.append(StateDataReporter("heating.txt", 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('heating.pdb', 1000))

simulation.step(35000) 

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()


#clear reporters
simulation.reporters = []


####################################
# 4) Density correction under NPT  #
####################################
print "Density correction under NPT"

newSystem.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtoseconds)

simulation = Simulation(modeller.topology, newSystem, integrator, platform, platformProperties)

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)

simulation.reporters.append(StateDataReporter("density.txt", 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('density.pdb', 1000))

simulation.step(35000) 

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()

#clear reporters
simulation.reporters = []






