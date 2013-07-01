# Typical AMBER GB MD protocol with N-acetyl-alanine-N'-methylamide
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from sys import stdout

platform = openmm.Platform_getPlatformByName("CUDA")

platformProperties = {}

# CUDA precision
platformProperties['CudaPrecision'] = 'mixed'
platformProperties['CudaDeviceIndex'] = '1'


#######################################################
# 1) Build system from PDB and assign AMBER FF99SB FF #
#######################################################
print "Building system"

pdb = PDBFile('./leap/snap.pdb')

forceField = ForceField('amber99sb.xml', 'amber99_obc.xml')
modeller = Modeller(pdb.topology, pdb.positions)

# Dump modeller structure
PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(),open('modeller.pdb', 'w'))



######################
# 2) Minimisation    #
######################
print "Minimising system"

system = forceField.createSystem(modeller.topology, nonbondedCutoff=20*angstrom, nonbondedMethod=CutoffNonPeriodic)

f = open('system.xml','w')
f.write(XmlSerializer.serializeSystem(system))
f.close()

integrator = VerletIntegrator(1*femtosecond)

simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)

simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=1000)

# Saving minimised positions
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))

#######################
#######################
## PMEMD like stages ##
#######################
#######################

#####################
# 3) Thermalisation #
#####################
print "Heating system"
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtoseconds)

# Note, new system
system = forceField.createSystem(modeller.topology, nonbondedCutoff=20*angstrom, nonbondedMethod=CutoffNonPeriodic)

# The reaction field method is turned on for NonbondedForce objects even when 
# GBSAOBCForce in present in a system. This is wrong and causes problems
#  https://github.com/SimTk/openmm/issues/26 
#  https://simtk.org/tracker/index.php?func=detail&aid=1881&group_id=161&atid=435
# This turns the reaction field off.
for i in range(system.getNumForces()):
   print type(system.getForce(i))
   if (type(system.getForce(i)) == openmm.NonbondedForce):
      system.getForce(i).setReactionFieldDielectric(1)


simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, time=True, potentialEnergy=True, temperature=True))
simulation.reporters.append(PDBReporter('heating.pdb', 1000))
simulation.step(35000)


# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()

#clear reporters
simulation.reporters = []

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()

#clear reporters
simulation.reporters = []


##################
# 4) Production  #
##################
print "Production"

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)

# Report every 0.1 ps
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, time=True, potentialEnergy=True, temperature=True) )
simulation.reporters.append(PDBReporter('production.pdb', 1000))

# 1 ns
simulation.step(1000000) 
