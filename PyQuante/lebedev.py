import numpy as np
from PyQuante import Lebedev
from PyQuante.Atom import *
from PyQuante.AtomicGrid import *

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


at1 = Atom(6,0,0,0)

print at1
print at1.atuple()

ag = AtomicGrid(at1, nrad=2)
print ag


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


for point in ag.points:
    print point.xyzw()
    ax.scatter(point._x, point._y, point._z)


plt.show()
