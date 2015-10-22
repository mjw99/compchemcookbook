import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def plotEnergyFromTorsionDataFile(atomI, atomJ, atomK, atomL):
  HartreeToKJMol = 2625.5 

  torsionStringName = str(atomI) + "_" + str(atomJ) + "_" + str(atomK) + "_" + str(atomL)

  # Temp arrays for angle and E
  angle = []
  energy = []

  fileName = "dihedral_scan_" + torsionStringName + ".dat"
  f = open(fileName,'r')

  for columns in ( raw.strip().split(",") for raw in f ):  
    angle.append(float(columns[0]))
    energy.append(float(columns[1]))

  f.close()

  energyInKJMol = [x*HartreeToKJMol for x in energy]
  minValue = min(energyInKJMol)
  energyInKJMol = [x-minValue for x in energyInKJMol]


  plt.figure()
  plt.title('Butane Total Energy B3LYP/6-31G*')

  plt.xlabel(r'' + torsionStringName + ' Torsion angle (degrees)')
  plt.ylabel(r'Energy (kJ/mol)')

  plt.plot(angle, energyInKJMol, marker='x')

  fileName = "dihedral_scan_" +  torsionStringName + ".png"
  plt.savefig(fileName)





plotEnergyFromTorsionDataFile(1,2,6,9)

