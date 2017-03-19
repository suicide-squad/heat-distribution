import os
import re
import sys
import numpy
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

def fault(file1, file2):
    array1 = numpy.loadtxt(file1)
    array2 = numpy.loadtxt(file2)

    yAbs = [abs(xi-xj) for xi, xj in zip(array1, array2)]

    absoluteFault = max(yAbs)
    return absoluteFault

dirMethods = os.path.join(os.path.dirname(__file__), "KEKUREKUS")
methods = os.listdir(dirMethods)

print(methods)
# Для каждой папки строяться по файлам графики!
for method in methods:
    directory = os.path.join(dirMethods, method)

    files = [os.path.join(directory, file) for file in os.listdir(directory)
             if file[-4:] == '.txt']

    number = re.compile(r'\d+\.?\d*')
    numbers = [float(number.findall(file)[0]) for file in files]
    # numbers = numbers[:-8]

    dic = dict(zip(numbers, files))

    fileDefault = dic.get(max(numbers))

    numbers.sort()
    files = [dic[i] for i in numbers]

    # absoluteFaults = [abs(xi-xj)/max(xi, xj) for xi, xj in zip(array1, array2)]

    absoluteFaults = [fault(file, fileDefault) for file in files]

    #####################################################################
    #                    Рисование графика                              #
    #####################################################################

    plt.figure(num = method, facecolor = (1, 1, .54))

    plt.plot(numbers, absoluteFaults,'bo',  label = method, color = 'black',)
    plt.plot(numbers, absoluteFaults,':k', color = 'green')
    # plt.yscale('log')

    plt.legend(loc = 1)
    plt.xlabel('timesize', fontsize = 12)
    plt.grid(True)
    ax = plt.gca()
    ax.yaxis.set_major_formatter(FuncFormatter(lambda x,y: '{:.0e}'.format(x)))
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x,y: '2^{:.0f}'.format(x)))

    plt.savefig(os.path.join('plots', method + '.png'))
    # plt.show()
