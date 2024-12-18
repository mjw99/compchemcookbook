"""
ambernetcdffile.py: Used for writing Amber NetCDF files.
Please see http://ambermd.org/netcdf/nctraj.html

Based upon dcdfile.py by Peter Eastman
"""
__author__ = "Mark J. Williamson"
__version__ = "0.1"

import numpy as np
import math
from Scientific.IO.NetCDF import NetCDFFile
from simtk.unit import angstroms, is_quantity, picoseconds, norm



class AmberNetCDFFile(object):
  """AmberNetCDFFile provides a method for wrting AMBER NetCDF files.

  To use this class, create a AmberNetCDFFile object, then call writeModel() once for each model in the file."""

  def __init__(self, filename, topology):
    """Create an AmberNetCDFFile, and write the initial header data.

    Parameters:
      - file (file) A file to write to
      - topology (Topology) The Topology defining the molecular system being written
    """

    self._filename = filename
    self._topology = topology
    self._modelCount = 0

    self.nc =  NetCDFFile(self._filename, 'w')
    
    # Global attributes
    setattr(self.nc, 'title', 'default_name')
    setattr(self.nc, 'application', 'AMBER')
    setattr(self.nc, 'program', 'OpenMM')
    # TODO, make this automatic
    setattr(self.nc, 'programVersion', '4.1.1')
    setattr(self.nc, 'Conventions', 'AMBER')
    setattr(self.nc, 'ConventionVersion', '1.0')

    # dimensions
    self.nc.createDimension('frame', None)
    self.nc.createDimension('spatial', 3)
    self.nc.createDimension('atom', len(list(self._topology.atoms())) )   
    self.nc.createDimension('label', 5)
    self.nc.createDimension('cell_spatial', 3)
    self.nc.createDimension('cell_angular', 3)

    # variable
    self.spatial = self.nc.createVariable('spatial', 'c', ('spatial',))

    self.time = self.nc.createVariable('time', 'd', ('frame',))
    self.time.units = 'picosecond'
    self.time.long_name = 'Time point of the simulations'

    self.coordinates = self.nc.createVariable('coordinates', 'f', ('frame','atom','spatial'))
    self.coordinates.units = 'angstrom'
    self.coordinates.long_name = 'Coordinates of all the atoms in the simulation'

    self.cell_spatial = self.nc.createVariable('cell_spatial', 'c', ('cell_spatial',))
    self.cell_angular = self.nc.createVariable('cell_angular', 'c', ('cell_angular','label',))

    self.cell_lengths = self.nc.createVariable('cell_lengths', 'd', ('frame', 'cell_spatial'))
    self.cell_lengths.units = 'angstrom'

    self.cell_angles = self.nc.createVariable('cell_angles', 'd', ('frame', 'cell_angular'))
    self.cell_angles.units = 'degrees'


    # Variables
    self.spatial[:] = 'xyz'
    self.cell_spatial[:] = 'abc'
    self.cell_angular[:] = 'alpha', 'beta ', 'gamma'

  def writeModel(self, positions, time):
    """Write out a model to the AMBER NetCDF file.

    Parameters:
      - positions (list) The list of atomic positions to write.
      - time The current time of the frame being written

    """


    # Basic sanity checking
    if len(list(self._topology.atoms())) != len(positions):
            raise ValueError('The number of positions must match the number of atoms')
    if is_quantity(positions):
          positions = positions.value_in_unit(angstroms)
    if is_quantity(time):
          time = time.value_in_unit(picoseconds)
    if any(math.isnan(norm(pos)) for pos in positions):
           raise ValueError('Particle position is NaN')
    if any(math.isinf(norm(pos)) for pos in positions):
            raise ValueError('Particle position is infinite')


    self.time[self._modelCount] = time
    
    self.coordinates[self._modelCount] = positions

    boxSize = self._topology.getUnitCellDimensions()
    if boxSize is not None:
      size = boxSize.value_in_unit(angstroms)
      self.cell_lengths[self._modelCount] =  size[0], size[1], size[2]
      self.cell_angles[self._modelCount] = 90.0, 90.0, 90.0


    self._modelCount += 1
    self.nc.sync()

  def close(self):
    self.nc.close()
