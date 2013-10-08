# Obtaining per term energy contributions in OpenMM is hard.
# This takes an original system and then copies out each class
# of force to their own system.

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import time
from copy import deepcopy 

platform = openmm.Platform_getPlatformByName("Reference")

platformProperties = {"OpenCLDeviceIndex":"0"}

print "Building system"
forceField = ForceField('amber99sb.xml');
alaala = PDBFile('alaala.pdb')
modeller = Modeller(alaala.topology, alaala.positions)

system = forceField.createSystem(modeller.topology, nonbondedMethod=NoCutoff)
# Serialise it
f = open('system.xml','w')
f.write(XmlSerializer.serializeSystem(system))
f.close()
del system


system = forceField.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HBonds)
# Serialise it
f = open('system_HBonds.xml','w')
f.write(XmlSerializer.serializeSystem(system))
f.close()
del system

system = forceField.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=HAngles)

# Serialise it
f = open('system_HAngles.xml','w')
f.write(XmlSerializer.serializeSystem(system))
f.close()
del system


