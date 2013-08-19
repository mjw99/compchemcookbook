# Example to demonstrate a difference in behaviour between 
# NonbondedForce.CutoffPeriodic and CustomNonbondedForce.CutoffPeriodic

# Using NoCutoff
# NonbondedForce Total potential energy is
# 18.4478996817 kcal/mol
# CustomNonbondedForce Total potential energy is
# 18.4478996817 kcal/mol
#
# Using CutoffPeriodic
# NonbondedForce Total potential energy is
# 5274.61198913 kcal/mol
# CustomNonbondedForce Total potential energy is
# 5277.20604162 kcal/mol
#
# Where has the 2.59405249 kcal/mol gone?


from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *





def CustomNonbondedForceCutoffPeriodicTest(useNoCutOff):
  # Standard system
  systemNonbondedForce = System()

  systemNonbondedForce.addParticle( 35.4532 * dalton )
  systemNonbondedForce.addParticle( 35.4532 * dalton )

  systemNonbondedForce.setDefaultPeriodicBoxVectors( Vec3(2*nanometer,0,0), Vec3(0,2*nanometer,0), Vec3(0,0,2*nanometer) )


  nonbondedForce = NonbondedForce();
  if (useNoCutOff):
    nonbondedForce.setNonbondedMethod(NonbondedForce.NoCutoff)
  else:
    nonbondedForce.setNonbondedMethod(NonbondedForce.CutoffPeriodic)
  nonbondedForce.setCutoffDistance(8*angstrom)

  nonbondedForce.setUseDispersionCorrection(False)
  nonbondedForce.setReactionFieldDielectric(0)

  #nonbondedForce.addParticle( charge,sigma,epsilon  )
  nonbondedForce.addParticle(-1.0 *elementary_charge, 0.440103966761*nanometer, 0.4184*kilojoule_per_mole)
  nonbondedForce.addParticle(-1.0 *elementary_charge, 0.440103966761*nanometer, 0.4184*kilojoule_per_mole)

  systemNonbondedForce.addForce(nonbondedForce)

  integrator = VerletIntegrator(1*femtosecond)

  contextNonbondedForce = Context(systemNonbondedForce, integrator)

  initPos = [ ]
  initPos.append( (1,0,0)*angstrom )
  initPos.append( (19,0,0)*angstrom )

  contextNonbondedForce.setPositions(initPos)

  stateNonbondedForce = contextNonbondedForce.getState( getEnergy=True)
  print "NonbondedForce Total potential energy is " 
  print  str(stateNonbondedForce.getPotentialEnergy().in_units_of(kilocalorie/mole))










  # Custom system
  systemCustomNonbondedForce = System()

  systemCustomNonbondedForce.addParticle( 35.4532 * dalton )
  systemCustomNonbondedForce.addParticle( 35.4532 * dalton )

  systemCustomNonbondedForce.setDefaultPeriodicBoxVectors( Vec3(2*nanometer,0,0), Vec3(0,2*nanometer,0), Vec3(0,0,2*nanometer) )



  combinedVDWEE = CustomNonbondedForce("4*eps*((sigma/r)^12-(sigma/r)^6)+138.935456*q/r;"
"q=q1*q2;"
"sigma=0.5*(sig1+sig2);"
"eps=sqrt(eps1*eps2)");
  #combinedVDWEE = CustomNonbondedForce("r")
  combinedVDWEE.addPerParticleParameter("eps")
  combinedVDWEE.addPerParticleParameter("q")
  combinedVDWEE.addPerParticleParameter("sig")
  if (useNoCutOff):
    combinedVDWEE.setNonbondedMethod(CustomNonbondedForce.NoCutoff)
  else:
    combinedVDWEE.setNonbondedMethod(CustomNonbondedForce.CutoffPeriodic)
  combinedVDWEE.setCutoffDistance(8*angstrom)


  #combinedVDWEE.addParticle( [epsilon,charge,sigma ]  )
  #  <Atom type="1954" charge="-1.0" sigma="0.440103966761" epsilon="0.4184"/>
  combinedVDWEE.addParticle( [0.4184*kilojoule_per_mole, -1.0 *elementary_charge , 0.440103966761*nanometer ]  )
  combinedVDWEE.addParticle( [0.4184*kilojoule_per_mole, -1.0 *elementary_charge , 0.440103966761*nanometer ]  )


  systemCustomNonbondedForce.addForce(combinedVDWEE)

  integrator = VerletIntegrator(1*femtosecond)



  contextCustomNonbondedForce = Context(systemCustomNonbondedForce, integrator)

  initPos = [ ]
  initPos.append( (1,0,0)*angstrom )
  initPos.append( (19,0,0)*angstrom )

  contextCustomNonbondedForce.setPositions(initPos)

  stateCustomNonbondedForce = contextCustomNonbondedForce.getState( getEnergy=True)
  print "CustomNonbondedForce Total potential energy is " 
  print str(stateCustomNonbondedForce.getPotentialEnergy().in_units_of(kilocalorie/mole))






## SETUP
platformProperties = {}
# Set up platform
platform = openmm.Platform_getPlatformByName("CUDA")

platformProperties = {}
# CUDA precision
platformProperties['CudaPrecision'] = 'mixed'
# CUDA parallel
platformProperties['CudaDeviceIndex'] = '0'





print "1) Using NoCutoff"
CustomNonbondedForceCutoffPeriodicTest(1)

print ""
print "2) Using CutoffPeriodic"
CustomNonbondedForceCutoffPeriodicTest(0)


