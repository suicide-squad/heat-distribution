import os
import re
import sys
import numpy
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QFileDialog, QApplication

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

absoluteFault = max(yAbs)
relativeFault = max(abs(xi-xj)/max(xi, xj) for xi, xj in zip(array1, array2))
print ("абсолютная:	%.15f" % absoluteFault)
print ("относительная:	%.15f" % relativeFault)

try:
    assert len(x) == len(yAbs)
    # Рисование графиков
    plt.plot(x, yAbs, label ='absolut')
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
    plt.xlabel('x', fontsize=14)
    plt.ylabel('abs', fontsize=14)
    plt.grid(True)
    # plt.yscale('log')
    plt.show()
except AssertionError:
	print ('ERROR! Не совпадают размерности!')

sys.exit(app.exec_())
