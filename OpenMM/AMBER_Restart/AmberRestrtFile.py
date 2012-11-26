"""
amberrestrtfile.py: Used for writing AMBER restrt files.

"""
__author__ = "Mark J. Williamson"
__version__ = "1.0"

from simtk.unit import Quantity, angstroms, picoseconds
try:
    import numpy
except:
    pass


class AmberRestrtFile():
    """AmberRestrtFile saves an AMBER restrt file."""
    
    def __init__(self, filename, coordinates, velocities):
        """Save a restrt file.
        
        Parameters:
         - file (string) the name of the file to save
        """
        self.writeAmberCoordinates(filename, coordinates, velocities)


    def writeAmberCoordinates(self, filename, coordinates, velocities):
       # Open coordinate file for writing.
       outfile = open(filename, 'w')
       # Write title
       outfile.write("foo" + "\n")

       # Coordinates is populated with Vec3's hence its length will be the same
       # as the number of atoms in the system
       natoms = len(coordinates)

       # Write number of atoms + time
       #FORMAT(I5,5E15.7) NATOM,TIME
       outfile.write( '{0:5}{1:15.7f}\n'.format(natoms, 0.0 ) )


       # coordinates and velocities are of type simtk.unit.quantity.Quantity

       # Translate units for coordinates
       coordinates = coordinates.value_in_unit(angstroms)
       # Translate units for velocities
       velocities = velocities.value_in_unit(angstroms/picoseconds) 

       def divideByFactor(x): return x/20.455
       velocities = map(divideByFactor, velocities)


       # Coordinates first#

       # Local counter to check if one should be adding a newline
       val = 0
       for item in coordinates:
          if ( val%2 == 0 ):
            outfile.write( '{0:12.7f}{1:12.7f}{2:12.7f}'.format(item[0], item[1], item[2]) )
            val += 1
          else:
            outfile.write( '{0:12.7f}{1:12.7f}{2:12.7f}\n'.format(item[0], item[1], item[2]) )
            val += 1

       # Velocities next
       val = 0
       for item in velocities:
          if ( val%2 == 0 ):
            outfile.write( '{0:12.7f}{1:12.7f}{2:12.7f}'.format(item[0], item[1], item[2]) )
            val += 1
          else:
            outfile.write( '{0:12.7f}{1:12.7f}{2:12.7f}\n'.format(item[0], item[1], item[2]) )
            val += 1

       # Box vectors?
       outfile.close()

