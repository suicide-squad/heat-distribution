from numpy import arange
from math import cos, pi

with open('INPUT_COMPLEX.txt') as file:
	setting = {line.split('=')[0] : line.split('=')[1] for line in file}
	file.close()

for key, value in setting.items():
	if key == 'XSTART' or key == 'XEND':
		setting[key] = complex(value.replace('i','j'))
	else:
		setting[key] = float(value)


xStart = setting['XSTART']
xFinish = setting['XEND']
NX = setting['NX']

step = abs(xStart.real - xFinish.real)/NX

massX = arange(xStart.real, xFinish.real, step, dtype = float)

massY = [cos(x*pi) if -0.5 < x < 0.5 else 0 for x in massX]


massY = [str(i)+'+0.0i' for i in massY]
print (len(massY))

with open('INPUT_COMPLEX.txt', 'a') as file:
	file.writelines('\n'.join(massY))
	file.close()
	