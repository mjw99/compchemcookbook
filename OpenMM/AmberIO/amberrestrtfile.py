"""
amberrestrtfile.py: Used for writing AMBER restrt files.
Please see http://ambermd.org/formats.html#restart

"""
__author__ = "Mark J. Williamson"
__version__ = "0.1"

from simtk.unit import Quantity, angstroms, picoseconds
import datetime
try:
    import numpy
except:
    pass


class AmberRestrtFile():
    """AmberRestrtFile saves an AMBER restrt file."""
    
    def __init__(self, filename, coordinates, velocities, boxVectors, time):
        """Save a restrt file.
        
        Parameters:
         - file (string) the name of the file to save
        """
        self.writeAmberCoordinates(filename, coordinates, velocities, boxVectors, time)


    def writeAmberCoordinates(self, filename, coordinates, velocities, boxVectors, time):
       # Open coordinate file for writing.
       outfile = open(filename, 'w')

       # Write title
       now = datetime.datetime.now()
       outfile.write('{0}{1:%Y-%m-%d %H:%M:%S}\n'.format(
                   "AMBER restart file written by AmberRestrtFile.py at ",now))

       # Coordinates is populated with Vec3's hence its length will be the same
       # as the number of atoms in the system
       natoms = len(coordinates)

       # Write number of atoms + time
       #FORMAT(I5,5E15.7) NATOM,TIME
       time = time.value_in_unit(picoseconds)
       outfile.write( '{0:5}{1:15.7f}\n'.format(natoms, time ) )


       # coordinates and velocities are of type simtk.unit.quantity.Quantity

       ## Coordinates first

       # Translate units for coordinates
       coordinates = coordinates.value_in_unit(angstroms)

       # Local counter to check if one should be adding a newline
       val = 0
       for item in coordinates:
          if ( val%2 == 0 ):
            outfile.write( '{0:12.7f}{1:12.7f}{2:12.7f}'.format(item[0], item[1], item[2]) )
            val += 1
          else:
            outfile.write( '{0:12.7f}{1:12.7f}{2:12.7f}\n'.format(item[0], item[1], item[2]) )
            val += 1

       ## Velocities next

       # Translate units for velocities
       velocities = velocities.value_in_unit(angstroms/picoseconds)

       # See here for 20.455 scaling
       # http://ambermd.org/Questions/units.html
       def divideByFactor(x): return x/20.455
       velocities = map(divideByFactor, velocities)


       # Velocities next
       val = 0
       for item in velocities:
          if ( val%2 == 0 ):
            outfile.write( '{0:12.7f}{1:12.7f}{2:12.7f}'.format(item[0], item[1], item[2]) )
            val += 1
          else:
            outfile.write( '{0:12.7f}{1:12.7f}{2:12.7f}\n'.format(item[0], item[1], item[2]) )
            val += 1

       # Box vectors
       boxVectors = boxVectors.value_in_unit(angstroms)

       outfile.write( '{0:12.7f}{1:12.7f}{2:12.7f}{3:12.7f}{3:12.7f}{3:12.7f}\n'.format(
                           boxVectors[(0)][0], 
                           boxVectors[(1)][1], 
                           boxVectors[(2)][2], 
                           90.0) )


       outfile.close()

