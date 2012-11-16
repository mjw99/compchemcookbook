# Typical AMBER MD protocol 

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


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

####################################################
# Build system from PDB and assign AMBER FF99SB FF #
####################################################

pdb = PDBFile('1UBQ.pdb')
forceField = ForceField('amber99sb.xml', 'tip3p.xml')


# Add missing hydrogens
modeller = Modeller(pdb.topology, pdb.positions)
modeller.addHydrogens(forceField)

# Add TIP3P solvent
modeller.addSolvent(forceField, model='tip3p', padding=10*angstrom)



print "Number of atoms %i"      % len(modeller.positions)

###################
# Minimisation    #
###################

system = forceField.createSystem(modeller.topology, nonbondedMethod=PME)

integrator = VerletIntegrator(1*femtosecond)
simulation = Simulation(modeller.topology, system, integrator, platform)
simulation.context.setPositions(modeller.positions)
simulation.minimizeEnergy(maxIterations=100)

# Saving minimised positions
positions = simulation.context.getState(getPositions=True).getPositions()
PDBFile.writeFile(simulation.topology, positions, open('output.pdb', 'w'))

#######################
#######################
## PMEMD like stages ##
#######################
#######################

####################
# Thermalisation   #
####################
