from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pylab
import matplotlib.pyplot as plt
import numpy as np
import os
import re

path = os.path.abspath(os.path.join(os.path.dirname(__file__),
									"..", "initial", "INPUT.txt"))

with open(path, 'r') as file:
    pattern = re.compile('[A-Za-z]+=-?\d+')
    setting = { line.split('=')[0] : float(line.split('=')[1])
                for line in file if pattern.match(line) }

xStart = setting['XSTART']
xFinish = setting['XEND']
NX = setting['NX']
tStart = setting['TSTART']
tEnd = setting['TFINISH']
dt = setting['dt']

step = abs(xFinish-xStart)/NX

Y = np.arange(xStart, xFinish, step)
X = np.arange(0, 11, 1)

X, Y = np.meshgrid(X, Y)


Zim = np.loadtxt('./../result/complex/Im.txt', unpack=True)
Zre = np.loadtxt('./../result/complex/Re.txt', unpack=True)


fig = pylab.figure()
ax1 = fig.add_subplot(122, projection = '3d')
ax1.plot_surface(X, Y, Zim,
                 color = 'y',
                 rstride=100,
                 cstride=100,
                 cmap='coolwarm',
                 linewidth=0.4)


ax1.set_xlabel('TIME')
ax1.set_ylabel('X')
ax1.set_zlabel('U(X)')

ax2 = fig.add_subplot(121, projection = '3d')
ax2.plot_surface(X, Y, Zre,
                 color = 'y',
                 rstride=100,
                 cstride=100,
                 cmap='coolwarm',
                 linewidth=0.4)

ax2.set_xlabel('TIME')
ax2.set_ylabel('X')
ax2.set_zlabel('U(X)')

plt.show()