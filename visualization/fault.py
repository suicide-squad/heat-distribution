import os
import re
import sys
import numpy
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QFileDialog, QApplication

from matplotlib.ticker import FuncFormatter

path = os.path.abspath(os.path.join(os.path.dirname(__file__),
									"..", "initial", "INPUT.txt"))
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

app = QApplication(sys.argv)
dialog = QFileDialog()
file, ok = dialog.getOpenFileName(None, 'Первый файл',
                                  '/home/kirill/heat-distribution/result',
                                  'Text files (*.txt)')


array1 = numpy.loadtxt(file)

file, ok = dialog.getOpenFileName(None, 'Второй файл',
                                  '/home/kirill/heat-distribution/result',
                                  'Text files (*.txt)')

array2 = numpy.loadtxt(file)

yAbs = [xi-xj for xi, xj in zip(array1, array2)]
yRelat = [(xi-xj)/max(xi, xj) for xi, xj in zip(array1, array2)]
yResult1 = array1[:]
yResult2 = array2[:]

absoluteFault = max(map(abs,yAbs))
relativeFault = max(map(abs,yRelat))
print("абсолютная:\t%.15f" % absoluteFault)
print("относительная:\t%.15f" % relativeFault)

#####################################################################
#                    Рисование графиков                             #
#####################################################################

try:
    assert len(x) == len(yAbs)

    plt.figure(num = 'FAULT', facecolor = (1, 1, .54))

    plt.subplot(221)
    plt.plot(x, yAbs, label ='absolut', color = 'green')

    plt.legend(loc = 2)
    plt.xlabel('x', fontsize = 14)
    plt.grid(True)
    ax = plt.gca()
    f = lambda x,y: '{:.15e}'.format(x)
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{:.0e}'.format(x)))
    plt.yscale('log')

    plt.subplot(222)
    plt.plot(x, yRelat, label = 'relative', color = 'red')
    plt.legend(loc = 2)
    plt.xlabel('x', fontsize=14)
    plt.grid(True)

    plt.subplot(223)
    plt.plot(x, yResult1, label = 'result 1', color = 'blue')
    plt.legend(loc = 2)
    plt.xlabel('x', fontsize = 14)
    plt.grid(True)

    plt.subplot(224)
    plt.plot(x, yResult2, label = 'result 2', color = 'brown')
    plt.legend(loc = 2)
    plt.xlabel('x', fontsize = 14)
    plt.grid(True)

    plt.show()
except AssertionError:
    print ('ERROR! Не совпадают размерности!')

sys.exit(app.exec_())

