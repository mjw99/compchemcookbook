# Typical AMBER MD protocol with CYP3A4

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from sys import stdout


#platform = openmm.Platform_getPlatformByName("Reference")
platform = openmm.Platform_getPlatformByName("CUDA")
#platform = openmm.Platform_getPlatformByName("OpenCL")


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
#platformProperties['OpenCLDeviceIndex'] = '0'

# CUDA parallel
platformProperties['CudaDeviceIndex'] = '1,2'
#platformProperties['CudaDeviceIndex'] = '0'




#####################
#####################
## Leap like stages##
#####################
#####################

#######################################################
# 1) Build system from PDB and assign AMBER FF99SB FF #
#######################################################
print "Building system"


# OpenMM 5.1 only
Topology.loadBondDefinitions('CPDI_CYP_residues.xml')
Modeller.loadHydrogenDefinitions('CPDI_CYP_hydrogens.xml')

pdb = PDBFile('1TQN_modded.pdb')

# if amber99sb.xml is first, CYP will be matched against CYX in amber99sb.xml
# and not CYP in CPDI_CYP_ff.xml
forceField = ForceField('CPDI_CYP_ff.xml','amber99sb.xml','tip3p.xml')

#
#for templateSignatures in forceField._templateSignatures:
#  print templateSignatures

modeller = Modeller(pdb.topology, pdb.positions)

# Add missing hydrogens
modeller.addHydrogens(forceField)

# Add TIP3P solvent

modeller.addSolvent(forceField, model='tip3p', padding=8*angstrom)

# Dump modeller structure
PDBFile.writeFile(modeller.getTopology(), modeller.getPositions(),open('modeller.pdb', 'w'))



print "Number of atoms %i"      % len(modeller.positions)

######################
# 2) Minimisation    #
######################
print "Minimising system"

system = forceField.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=8*angstrom)

integrator = VerletIntegrator(1*femtosecond)

simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)

print "Platform: %s" % (simulation.context.getPlatform().getName())

simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=1000)

# Saving minimised positions
positions = simulation.context.getState(getPositions=True, enforcePeriodicBox=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('minimisation.pdb', 'w'))

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


# Note, new system, with SHAKE
system = forceField.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=8*angstrom, constraints=HBonds, rigidWater=True )


# Set the COM Removal to something sensible
for i in range(system.getNumForces()):
   if (type(system.getForce(i)) == openmm.CMMotionRemover):
      system.getForce(i).setFrequency(1000)

simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)
simulation.context.setPositions(positions)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
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

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)

simulation = Simulation(modeller.topology, system, integrator, platform, platformProperties)

simulation.context.setPositions(positions)
simulation.context.setVelocities(velocities)

simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('density.pdb', 1000))

simulation.step(35000) # i.e. 20,000 fs == 20 ps == 0.02 ns

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
simulation.reporters.append(StateDataReporter(stdout, 50000, step=True, potentialEnergy=True, temperature=True, density=True))
simulation.reporters.append(PDBReporter('production.pdb', 50000))

# 10 ns
simulation.step(5000000) 

# Save the positions and velocities
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()









