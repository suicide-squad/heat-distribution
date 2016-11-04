import sys
import numpy
from PyQt5.QtWidgets import QFileDialog, QApplication

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

absoluteFault = max(abs(xi-xj) for xi, xj in zip(array1, array2))
relativeFault = max(abs(xi-xj)/max(xi, xj) for xi, xj in zip(array1, array2))
print ("абсолютная:	%.15f" % absoluteFault)
print ("относительная:	%.15f" % relativeFault)
sys.exit(app.exec_())
