import numpy as np
from math import cos, pi

with open('INPUT.txt') as file:
	setting = {line.split('=')[0] :float(line.split('=')[1]) for line in file}
	file.close()

xStart = setting['XSTART']
xEnd = setting['XEND']
NX = setting['NX']

X = np.linspace(xStart, xEnd, NX, dtype = float)

U = [cos(x*pi) if -0.5 < x < 0.5 else 0 for x in X]

U = list(map(str,U))
print (len(U))

with open('INPUT.txt', 'a') as file:
	file.writelines('\n'.join(U))
	file.close()
	