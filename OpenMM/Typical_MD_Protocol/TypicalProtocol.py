# Typical AMBER MD protocol 

import simtk.openmm.app as app
import simtk.openmm as mm
import simtk.unit as unit

from sys import stdout


#platform = mm.Platform_getPlatformByName("Reference")
#platform = mm.Platform_getPlatformByName("OpenCL")
platform = mm.Platform_getPlatformByName("CUDA")


platformProperties = {}
# OpenCL precision
#platformProperties['OpenCLPrecision'] = 'mixed'

# CUDA precision
platformProperties['CudaPrecision'] = 'mixed'


# Run on multiple cards
# 0  Tesla M2090
# 1  Tesla C2075
# 2  Tesla C2075

#OpenCL parallel
#platformProperties['OpenCLDeviceIndex'] = '0,1,2'
#platformProperties['OpenCLDeviceIndex'] = '1'
#platformProperties['OpenCLDeviceIndex'] = '1,2'

# CUDA parallel
#platformProperties['CudaDeviceIndex'] = '1,2'
platformProperties['CudaDeviceIndex'] = '1'





#####################
#####################
## Leap like stages##
#####################
#####################

#######################################################
# 1) Build system from PDB and assign AMBER FF99SB FF #
#######################################################
print "Building system"

# wget "http://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId=1UBQ" -O 1UBQ.pdb
#ff6700fb140ab9289134ba6555f87d0e  1UBQ.pdb

pdb = app.PDBFile('1UBQ.pdb')
forceField = app.ForceField('amber99sb.xml', 'tip3p.xml')


# Add missing hydrogens
modeller = app.Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forceField)

# Add TIP3P solvent
modeller.addSolvent(forceField, model='tip3p', padding=10*unit.angstrom)



print "Number of atoms %i"      % len(modeller.positions)

######################
# 2) Minimisation    #
######################
print "Minimising system"

system = forceField.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=8*unit.angstrom)

integrator = mm.VerletIntegrator(1*unit.femtosecond)

simulation = app.Simulation(modeller.topology, system, integrator, platform, platformProperties)

print "Platform: %s" % (simulation.context.getPlatform().getName())
#print platform.getPropertyValue(simulation.context, "CUDADeviceIndex")

simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=1000)

# Saving minimised positions
positions = simulation.context.getState(getPositions=True).getPositions()
app.PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))

#######################
#######################
## PMEMD like stages ##
#######################
#######################

dt = 2*unit.femtoseconds
friction = 1*(1/unit.picosecond)
constraints = app.HBonds
hydrogenMass = None

################################
# 3) Thermalisation under NVT  #
################################
print "Heating system under NVT"
integrator = mm.LangevinIntegrator(300*unit.kelvin, friction, dt)

# Note, new system, with SHAKE
system = forceField.createSystem(modeller.topology, nonbondedMethod=app.PME, nonbondedCutoff=8*unit.angstrom, constraints=constraints, hydrogenMass=hydrogenMass)

# Set the COM Removal to something sensible
for i in range(system.getNumForces()):
   if (type(system.getForce(i)) == mm.CMMotionRemover):
      system.getForce(i).setFrequency(1000)


simulation = app.Simulation(modeller.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)

simulation.reporters.append(app.StateDataReporter("heating.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=35000))
simulation.reporters.append(app.PDBReporter('heating.pdb', 1000))

simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()


#clear reporters
simulation.reporters = []


####################################
# 4) Density correction under NPT  #
####################################
print "Density correction under NPT"

system.addForce(mm.MonteCarloBarostat(1*unit.bar, 300*unit.kelvin))
integrator = mm.LangevinIntegrator(300*unit.kelvin, friction, dt)

simulation = app.Simulation(modeller.topology, system, integrator, platform, platformProperties)

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)

simulation.reporters.append(app.StateDataReporter("density.csv", 1000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=35000))
simulation.reporters.append(app.PDBReporter('density.pdb', 1000))

simulation.step(35000)

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()

#clear reporters
simulation.reporters = []


####################################
# 5) Production under NPT          #
####################################
print "Production under NPT"

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)


# Report every 0.1 ns / 100 ps
simulation.reporters.append(app.StateDataReporter("production.csv", 50000, time=True, potentialEnergy=True, temperature=True, density=True, remainingTime=True, speed=True, totalSteps=5000000))
simulation.reporters.append(app.PDBReporter('production.pdb', 50000))

# 10 ns
simulation.step(5000000) 

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()









