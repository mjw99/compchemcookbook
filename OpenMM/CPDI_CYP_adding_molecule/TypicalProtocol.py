# Typical AMBER MD protocol 

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from sys import stdout


#platform = openmm.Platform_getPlatformByName("Reference")
platform = openmm.Platform_getPlatformByName("CUDA")
#platform = openmm.Platform_getPlatformByName("OpenCL")


# OpenCL precision
#platformProperties = {"OpenCLPrecision":"mixed"}
# CUDA precision
platformProperties = {"CudaPrecision":"mixed"}


# Run on multiple cards
# 0  Tesla M2090
# 1  Tesla C2075
# 2  Tesla C2075

#OpenCL parallel
#platformProperties = {"OpenCLDeviceIndex":"0,1,2"}
#platformProperties = {"OpenCLDeviceIndex":"1,2"}

# CUDA parallel
#platformProperties = {"CudaDeviceIndex":"0,1,2"}
platformProperties = {"CudaDeviceIndex":"1"}



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

pdb = PDBFile('1TQN_modded.pdb')


# Be very careful with order here, if amber99sb is
# first, CYP will be matched against CYX in amber99sb.xml
# and not CYP in CPDI_CYP.xml
forceField = ForceField('CPDI_CYP.xml', 'amber99sb.xml', 'tip3p.xml', 'ketoconazole_ff.xml')


modeller = Modeller(pdb.topology, pdb.positions)


# Add missing hydrogens
modeller.addHydrogens(forceField)

# Add TIP3P solvent
modeller.addSolvent(forceField, model='tip3p', padding=14*angstrom)


# Add ketoconazole
ketoconazole = PDBFile('./leap/ketoconazole.pdb')
modeller.add(ketoconazole.topology, ketoconazole.positions)


# Dump modeller structure
PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(), open('modeller.pdb', 'w'))


######################
# 2) Minimisation    #
######################
print "Minimising system"

system = forceField.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=8*angstrom)
integrator = VerletIntegrator(1*femtosecond)

simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)

simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=1000)

# Saving minimised positions
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))



################################
# 3) Thermalisation under NVT  #
################################
print "Heating system under NVT"
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)

# Note, new system, with SHAKE
system = forceField.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=8*angstrom, constraints=HBonds )

# Set the COM Removal to something sensible
for i in range(system.getNumForces()):
   if (type(system.getForce(i)) == openmm.CMMotionRemover):
      system.getForce(i).setFrequency(1000)


simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('heating.pdb', 1000))

simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()


#clear reporters
simulation.reporters = []


