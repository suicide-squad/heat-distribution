import numpy
import matplotlib.pyplot as plt

import os
import re

path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "initial", "INPUT.txt"))
file = open(path, 'r')

pattern = re.compile('[A-Za-z]+=-?\d+')
setting = { line.split('=')[0] : float(line.split('=')[1]) 
			for line in file if pattern.match(line) }

file.close()

xStart = setting['XSTART']
xFinish = setting['XEND']
NX = setting['NX']

step = abs(xFinish-xStart)/NX

x = numpy.arange (xStart, xFinish, step)

file = open(path, 'r')
pattern = re.compile('\d+\.?\d*e?[+-]?\d*')
yStart = numpy.array([float(line) for line in file if pattern.match(line)])
file.close()

path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "result", "PetrovResult.txt"))
yFinish = numpy.loadtxt(path)

# Рисование графиков
plt.plot(x, yStart, label ='start time')
plt.plot(x, yFinish, label ='end time')
plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
plt.xlabel('x', fontsize=14)
plt.ylabel('U(x)', fontsize=14)
plt.grid(True) 
plt.show() 
