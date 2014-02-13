#!/usr/bin/python

from pylab import *
from scipy import *
from scipy.optimize import curve_fit

import numpy as np
import sys



def findDe(fileName):
  """Calculates a bond's dissociation energy.

  Give a list of bond lengths (in Angstroms) and corresponding energies
  (in Hartrees), this will fit a Morse function to the data and calculate
  the bond's dissociation energy, De.

  """

  HartreeInKcalmol = 627.509

  rawDataX = []
  rawDataY = []

  # Read
  #fileName = "bond_scan_" + str(atomI) + "_" + str(atomJ) + ".dat"
  #fileName = "bond_scan_1_12.dat"
  #fileName = sys.argv[1]
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
  zeroBondOffset =  dataX[dataY.argmin()]

  offsetDataX = dataX - zeroBondOffset
  offsetDataY = dataY - zeroEnergyOffset


  # Define Morse potential
  def Morse(x, De, alpha):
     return De*(1-exp( -alpha*( x ) ))**2

  # Actually do the fit
  popt, pcov = curve_fit(Morse, offsetDataX, offsetDataY)


  # Tidy up the naming
  words = fileName.split('_');


  De = popt[0]
  print "De is " + str(De * HartreeInKcalmol) + " kcal/mole of bond " + str(words[2]) + " " + str(words[3]) 


