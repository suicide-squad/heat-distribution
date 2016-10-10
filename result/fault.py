import sys
import numpy
from PyQt5.QtWidgets import QFileDialog, QApplication

if __name__ == '__main__':
    app = QApplication(sys.argv)
    dialog = QFileDialog()
    print("Выберете первый файл для сравнения")
    file, ok = dialog.getOpenFileName(None, 'Выбрать файл',
                                      '/home/kirill/heat-distribution/result',
                                      'Text files (*.txt)')

    array1 = numpy.loadtxt(file)

    print ("Выберете второй файл для сравнения")
    file, ok = dialog.getOpenFileName (None, 'Выбрать файл', '/home/kirill/heat-distribution/result',
                                       'Text files (*.txt)')

    array2 = numpy.loadtxt(file)

    m = max(abs(xi-xj) for xi,xj in zip(array1,array2))
    print (m)
    sys.exit(app.exec_())