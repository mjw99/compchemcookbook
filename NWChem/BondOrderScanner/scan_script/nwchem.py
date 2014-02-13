
import math
from nwgeom import *

# Variables
atomZValueI = 6
atomZValueJ = 1
cutOffDistance = 1.2


indexToZ = {}
indexToCoordinates = {}
bondsToConsider = []


# How many atoms are there?
numberOfAtoms = len(rtdb_get("geometry:geometry:charges"))

# Array of Z values per atom
ZArray = rtdb_get("geometry:geometry:charges")


rawGeometryInAu = rtdb_get("geometry:geometry:coords")
# Get scaling factor
factor = rtdb_get('geometry:geometry:angstrom_to_au')
# Carry out conversion
rawGeometryInAng = map(lambda x: x/factor, rawGeometryInAu)


# Note this runs 0 --> numberOfAtoms-1
for i in range(numberOfAtoms):
  # Add Atom to index
  indexToZ[i+1] = int(ZArray[i])

  # Add Cartesians to index
  # This pretty mess is due to the rawGeometryInAng starting from 0
  x = (3*(i)) + 0
  y = (3*(i)) + 1
  z = (3*(i)) + 2


  indexToCoordinates[i+1] = rawGeometryInAng[x], rawGeometryInAng[y], rawGeometryInAng[z]


# Show atom index, Zvalue and coordinates
#print "Index, Zvalue, Coordinates"
#for key in indexToZ.keys():
#  print key, indexToZ[key], indexToCoordinates[key]



# Step through index
for keyI in indexToZ.keys():
  # Select current Atom C
  if indexToZ[keyI] == atomZValueI:
    # Step through index
    for keyJ in indexToZ.keys():
      # Select atom H
      if indexToZ[keyJ] == atomZValueJ:
        # Calculate distance from C_i to H_j
        x = float(indexToCoordinates[keyI][0]) - float(indexToCoordinates[keyJ][0])
        y = float(indexToCoordinates[keyI][1]) - float(indexToCoordinates[keyJ][1])
        z = float(indexToCoordinates[keyI][2]) - float(indexToCoordinates[keyJ][2])

        dist =  math.sqrt(x*x + y*y + z*z)

        # Typical C-H bond length
        if dist < cutOffDistance :
          #print keyI, indexToZ[keyI], indexToCoordinates[keyI], keyJ, indexToZ[keyJ], indexToCoordinates[keyJ], dist
          bondsToConsider.append( (keyI,keyJ,dist) )

print "Bonds to consider for this molecule"
for item in bondsToConsider:
  print item


# Save Optimised geom for the moment
optimised_list = geom_get_coords('geometry')

for currentBond in bondsToConsider:
  # Restore optimised geom
  geom_set_coords('geometry',optimised_list)

  # Set our bond parameters
  atomI  = currentBond[0]
  atomJ  = currentBond[1]
  eq = currentBond[2]
  print "Atom i is " + str(atomI) + " atom j is " + str(atomJ) + " eq is " + str(eq)



  # Scan

  geomWithOutBondInformation = '''
    geometry adjust
      zcoord
        bond  %i %i    %s cccc constant
      end
    end
  '''

  # This is truly dirty; I will go to Hell for this...
  geomWithBondInformation =  geomWithOutBondInformation % (atomI,atomJ,"%f")


  # How far aside of the equilibrium bond length should we go?
  min = eq - 0.2
  max = eq + 0.2
  results=scan_input(geomWithBondInformation,[min],[max],25,'dft',task_optimize)

  # Temp arrays for r and E
  x = []
  y = []

  for i in range(0, len(results) ):
     x.append(results[i][0][0])
     y.append(results[i][1][0])

  #Write data
  fileName = "bond_scan_" + str(atomI) + "_" + str(atomJ) + ".dat"
  f = open(fileName,'w')

  for i in range(0, len(x) ):
    f.write("%s, %s \n" % (float(x[i]),float(y[i])) )

  f.close()




