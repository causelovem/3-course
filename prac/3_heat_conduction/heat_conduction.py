import sys
import matplotlib.pyplot as plt
import numpy as np
#import scipy as sc


#def slau(matrixA, matrixB, res):
#    res = sc.solve(matrixA. matrixB)
#
#    return


def reverse_metod(data, num, dist, a, step, fin):
    return


def evident_metod(data, num, dist, a, step, fin):
    for i in range(1, data.shape[0]):
        data[i][0] = data[i - 1][0]
        data[i][num - 1] = data[i - 1][num - 1]
        for j in range(1, data.shape[1] - 1):
            point = (data[i - 1][j + 1] + data[i - 1][j - 1])
            point = point - 2 * data[i - 1][j]
            data[i][j] = step * a * point / (dist ** 2) + data[i - 1][j]
    return


def plot(data, step):
    plt.figure()
    plt.ion()
    for i in range(data.shape[0]):
        plt.plot(data[i])
        plt.title('Temperature_changes')
        plt.grid(True)
        plt.pause(0.1)
#        plt.clf()
# Uncomment code string above to draw different plots

    plt.show(block=True)
    return


if len(sys.argv) != 3:
    print ">Unexpected quantity of arguments, check your comand string.\n"
    exit(-1)

in_file = open(sys.argv[1])

file = in_file.readlines()

tmp = file[0].split()

num = int(tmp[0])
dist = float(tmp[1])
a = float(tmp[2])
step = float(tmp[3])
fin = int(tmp[4])

tmp = file[1].split()
data = np.zeros((fin, num))

for i in range(data.shape[1]):
    data[0][i] = float(tmp[i])

evident_metod(data, num, dist, a, step, fin)
plot(data, step)

in_file.close()
