from numpy import arange
from math import cos, pi

with open('INPUT.txt') as file:
	setting = {line.split('=')[0] :float(line.split('=')[1]) for line in file}
	file.close()

xStart = setting['XSTART']
xFinish = setting['XEND']
NX = setting['NX']

step = (abs(xStart)+abs(xFinish))/NX

massX = arange(xStart, xFinish, step, dtype = float)

massY = [cos(x*pi) if -0.5 < x < 0.5 else 0 for x in massX]

massY = list(map(str,massY))
print (len(massY))

with open('INPUT.txt', 'a') as file:
	file.writelines('\n'.join(massY))
	file.close()
	