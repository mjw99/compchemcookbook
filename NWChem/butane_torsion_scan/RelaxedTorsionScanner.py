from nwgeom import *



def scanTorsion(atomI, atomJ, atomK, atomL):
  """
  Carries out a relaxed scan about a torsion.

  Refs:
      http://www.nwchem-sw.org/index.php/Release65:Python
      http://verahill.blogspot.com.au/2013/08/503-relaxed-pes-scanning-in-nwchem.html
  
  """

  geomWithOutTorsionInformation = '''
    geometry noprint adjust
      zcoord
        torsion %i %i %i %i %s cccc constant
      end
    end
  '''

  geomWithTorsionInformation =  geomWithOutTorsionInformation % (atomI, atomJ, atomK, atomL, "%f")

  results = scan_input(geomWithTorsionInformation, [180.0], [540.0], 72, 'dft', task_optimize)


  # Temp arrays for angle and energy
  angle = []
  energy = []

  for i in range(0, len(results) ):
     # Angle in degrees
     angle.append(results[i][0][0])

     # Energy in Hartrees
     energy.append(results[i][1][0])


  # Form filename
  fileName = "dihedral_scan_" + str(atomI) + "_" + str(atomJ) + "_" + str(atomK) + "_" + str(atomL) + ".dat"

  # Write data
  f = open(fileName,'w')

  for i in range(0, len(angle) ):
    f.write("%s, %s \n" % (float(angle[i]), float(energy[i])) )

  f.close()

