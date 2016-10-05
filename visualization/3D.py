# НЕ ПОЛЬЗОВАТЬСЯ! НА СТАДИИ РАЗРАБОТКИ!

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pylab
import matplotlib.pyplot as plt
import numpy

Z = numpy.loadtxt('OUTPUT.txt', unpack=True)
Y = numpy.arange (-1, 1, 0.004)
X = numpy.arange (0, 2, 1)
X, Y = numpy.meshgrid(X, Y)

fig = pylab.figure()
axes = Axes3D(fig)
axes.plot_surface(X, Y, Z)

plt.show()