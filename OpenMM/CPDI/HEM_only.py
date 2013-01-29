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
#platform=openmm.Platform_getPlatform(2)


# Run on multiple cards
# 0  Tesla M2090
# 1  Tesla C2075
# 2  Tesla C2075
#platformProperties = {"OpenCLDeviceIndex":"0,1,2"}
#platformProperties = {"OpenCLDeviceIndex":"1,2"}
platform = openmm.Platform_getPlatformByName("OpenCL")
print "Speed relative to reference is : " + str(platform.getSpeed())




#####################
#####################
## Leap like stages##
#####################
#####################

#######################################################
# 1) Build system from PDB and assign AMBER FF99SB FF #
#######################################################
print "Building system"

pdb = PDBFile('1TQN_HEM.pdb')

forceField = ForceField('CPDI_CYP.xml')

# for templateSignatures in forceField._templateSignatures:
#  print templateSignatures

modeller = Modeller(pdb.topology, pdb.positions)
# Add missing hydrogens
modeller.addHydrogens(forceField)


# Write state after hydrogens added
PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(),open('postHydrogenAddition.pdb', 'w'))


print "Number of atoms %i"      % len(modeller.positions)

######################
# 2) Minimisation    #
######################
print "Minimising system"

system = forceField.createSystem(modeller.topology, nonbondedMethod=NoCutoff)

integrator = VerletIntegrator(1*femtosecond)

#simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)
simulation = Simulation(modeller.topology, system, integrator, platform)

print "Platform: %s" % (simulation.context.getPlatform().getName())
#print platform.getPropertyValue(simulation.context, "OpenCLDeviceIndex")

simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=1000)

# Saving minimised positions
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))


################################
# 3) Thermalisation under NVT  #
################################
print "Heating system under NVT"
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 1*femtoseconds)

# Note, new system, with SHAKE
system = forceField.createSystem(modeller.topology, nonbondedMethod=NoCutoff  )

# Set the COM Removal to something sensible
#for i in range(system.getNumForces()):
#   if (type(system.getForce(i)) == openmm.CMMotionRemover):
#      system.getForce(i).setFrequency(1000)


simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(positions)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('heating.pdb', 1000))

simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()


#clear reporters
simulation.reporters = []

