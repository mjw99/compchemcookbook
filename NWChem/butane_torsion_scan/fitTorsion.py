#!/usr/bin/python

from pylab import *
from scipy import *
from scipy.optimize import curve_fit

import numpy as np
import sys


def plot(x,y,yfit):
  plt.figure()
  plt.title('Torsion')

  plt.xlabel(r' Torsion angle (degrees)')
  plt.ylabel(r'Energy (kcal/mol)')

  plt.plot(x, y, marker='x')
  plt.plot(x, yfit, marker='x')

  plt.show()




def fitTorsion(fileName):
  """Fits a torsion energy profile to a cosine series of terms

  """

  HartreeInKcalmol = 627.509

  rawDataX = []
  rawDataY = []

  # Read
  f = open(fileName,'r')

  for line in iter(f):
    wordList = line.rstrip('\n').split(',')
    rawDataX.append(float(wordList[0]))
    rawDataY.append(float(wordList[1]))

  f.close()

  # scipy.optimize needs the arrays in this form
  dataX = np.array(rawDataX)
  dataY = np.array(rawDataY)


  # Offset the data for fitting
  zeroEnergyOffset = dataY.min()
  offsetDataY = dataY - zeroEnergyOffset

  # X is fine...
  offsetDataX = dataX


  # Define Torsion
  def Torsion(x, A, n, d):
     return A*(1.0 + cos( n * np.radians(x) - np.radians(d) ))

  # Define Torsion
  def Torsion3(x,
A, An, Ad,
B, Bn, Bd,
C, Cn, Cd,
):
     return Torsion(x, A, An, Ad) + Torsion(x, B, Bn, Bd) + Torsion(x, C, Cn, Cd)


  # Actually do the fit
  popt, pcov = curve_fit(Torsion3, offsetDataX, offsetDataY, maxfev=100000)

  print popt

  plot( offsetDataX, offsetDataY,Torsion3(offsetDataX,
popt[0], popt[1], popt[2], 
popt[3], popt[4], popt[5], 
popt[6], popt[7], popt[8]) )







fitTorsion('dihedral_scan_1_2_6_9.dat')
