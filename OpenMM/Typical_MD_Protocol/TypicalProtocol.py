# Typical AMBER MD protocol 

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from sys import stdout


# Reference
#platform=openmm.Platform_getPlatform(0)
# CUDA
#platform=openmm.Platform_getPlatform(1)
# OpenCL
platform=openmm.Platform_getPlatform(2)




#####################
#####################
## Leap like stages##
#####################
#####################

#######################################################
# 1) Build system from PDB and assign AMBER FF99SB FF #
#######################################################
print "Building system"

pdb = PDBFile('1UBQ.pdb')
forceField = ForceField('amber99sb.xml', 'tip3p.xml')


# Add missing hydrogens
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forceField)

# Add TIP3P solvent
modeller.addSolvent(forceField, model='tip3p', padding=10*angstrom)



print "Number of atoms %i"      % len(modeller.positions)

######################
# 2) Minimisation    #
######################
print "Minimising system"

system = forceField.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=8*angstrom)

integrator = VerletIntegrator(1*femtosecond)
simulation = Simulation(modeller.topology, system, integrator, platform)

print "Platform: %s" % (simulation.context.getPlatform().getName())

simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=1000)

# Saving minimised positions
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minisation.pdb', 'w'))

#######################
#######################
## PMEMD like stages ##
#######################
#######################

################################
# 3) Thermalisation under NVT  #
################################
print "Heating system under NVT"
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)

# Note, new system
system = forceField.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=8*angstrom, constraints=HBonds )

simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(positions)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('heating.pdb', 1000))

simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()



####################################
# 4) Density correction under NPT  #
####################################
print "Density correction under NPT"

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)

simulation = Simulation(modeller.topology, system, integrator, platform)

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('density.pdb', 1000))

simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()


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









