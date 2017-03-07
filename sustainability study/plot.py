import os
import re
import sys
import numpy
import matplotlib.pyplot as plt

def fault(file1, file2):
    array1 = numpy.loadtxt(file1)
    array2 = numpy.loadtxt(file2)

    yAbs = [abs(xi-xj) for xi, xj in zip(array1, array2)]

    absoluteFault = max(yAbs)
    return absoluteFault

directory = os.path.abspath(os.path.join(os.path.dirname(__file__),
									"KEKUREKUS", "Runge–Kutt"))

fileDefault = os.path.join(directory, "OUTPUT_Runge.txt")

files = [os.path.join(directory, file) for file in os.listdir(directory)
         if file[-4:] == '.txt' and file != "OUTPUT.txt" and file != "OUTPUT_Runge.txt"]

print(files)
number = re.compile(r'\d+')
numbers = [int(number.findall(file)[0]) for file in files]

absoluteFaults = [fault(file, fileDefault) for file in files]

x = zip(numbers, absoluteFaults)
print(list(x))
#####################################################################
#                    Рисование графика                              #
#####################################################################

plt.figure(num = 'FAULT', facecolor = (1, 1, .54))

plt.plot(numbers, absoluteFaults, 'ro',  label ='absolut', color = 'green')
plt.yscale('symlog')

# plt.legend(loc = 2)
# plt.xlabel('x', fontsize = 14)
# plt.grid(True)
# ax = plt.gca()
# f = lambda x,y: '{:.15e}'.format(x)
# ax.yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{:.0e}'.format(x)))

plt.show()

# print(dic)