from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time
from AmberRestrtFile import AmberRestrtFile


# Run on multiple cards
# 0  Tesla M2090
# 1  Tesla C2075
# 2  Tesla C2075
#platformProperties = {"OpenCLDeviceIndex":"0"}
#platformProperties = {"OpenCLDeviceIndex":"1,2"}
platformProperties = {"OpenCLDeviceIndex":"0,1,2"}

platform = openmm.Platform_getPlatformByName("OpenCL")
print "Speed relative to reference is : " + str(platform.getSpeed())





# Prepared via LEAP; I'm hoping to replace this with modeller,
# but I need to parse the Heme frcmods to XML first.
prmtop = AmberPrmtopFile('prmtop')
inpcrd = AmberInpcrdFile('inpcrd',  loadBoxVectors=True)

system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=8*angstrom, constraints=HBonds )
integrator = VerletIntegrator(1*femtosecond)

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties )
print "Platform: %s" % (simulation.context.getPlatform().getName())
print "Number of atoms %i"      % len(inpcrd.positions)



######################
# 2) Minimisation    #
######################
print "Minimising system"
integrator = VerletIntegrator(1*femtosecond)

simulation.context.setPositions(inpcrd.positions)
LocalEnergyMinimizer.minimize(simulation.context)

# Saving minimised positions
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))

# Write the AMBER restart
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
boxVectors = simulation.context.getState().getPeriodicBoxVectors()
time = simulation.context.getState().getTime()

restrt = AmberRestrtFile('minimsation.restrt', positions, velocities, boxVectors, time)


#######################
#######################
## PMEMD like stages ##
#######################
#######################

################################
# 3) Thermalisation under NVT  #
################################

# Use all GPUs
platformProperties = {"OpenCLDeviceIndex":"0,1,2"}

print "Heating system under NVT"
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)

# Note, new system
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=8*angstrom, constraints=HBonds )

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('heating.pdb', 1000))

simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()

# Write the AMBER restart
boxVectors = simulation.context.getState().getPeriodicBoxVectors()
time = simulation.context.getState().getTime()

restrt = AmberRestrtFile('NVT.restrt', positions, velocities, boxVectors, time)



#clear reporters
simulation.reporters = []


####################################
# 4) Density correction under NPT  #
####################################
print "Density correction under NPT"

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)

simulation = Simulation(prmtop.topology, system, integrator, platform, platformProperties)

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('density.pdb', 1000))

simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()

# Write the AMBER restart
boxVectors = simulation.context.getState().getPeriodicBoxVectors()
time = simulation.context.getState().getTime()

restrt = AmberRestrtFile('NPT.restrt', positions, velocities, boxVectors, time)


#clear reporters
simulation.reporters = []


####################################
# 5) Production under NPT          #
####################################
print "Production under NPT"

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)


# Report every 0.1 ns / 100 ps
simulation.reporters.append(StateDataReporter(stdout, 50000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('production.pdb', 50000))

# 10 ns
simulation.step(5000000) 

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()

# Write the AMBER restart
boxVectors = simulation.context.getState().getPeriodicBoxVectors()
time = simulation.context.getState().getTime()

restrt = AmberRestrtFile('production.restrt', positions, velocities, boxVectors, time)






