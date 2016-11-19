# !!!FILE FORMAT!!!
# <number of points> <distance between points> <a ^ 2> <tau> <final time>
# <first point temperature> <...> <last point temperature>

import sys
import math
import matplotlib.pyplot as plt
import numpy as np


def reverse_metod(data, num, dist, a, t):
    tmp = dist / (t * a)
    time = 0
    for i in range(1, data.shape[0]):
        time += t
        al = np.zeros(num)
        be = np.zeros(num)
        al[0] = 0
        be[0] = data[i - 1][0]
        # data[i][0] = data[i - 1][0]
        # data[i][-1] = data[i - 1][num - 1]
        data[i][0] = 0
        data[i][-1] = math.sin(time * 180 / math.pi)

        for j in range(1, num):
            al[j] = 1 / (2 + tmp - al[j - 1])
            be[j] = (be[j - 1] + tmp * data[i - 1][j - 1]) / (2 + tmp - al[j - 1])

        # data[i][-1] = (be[-1] - math.sqrt(dist) * math.cos(time * 180 / math.pi)) / (1 - al[-1])
        for j in range(data.shape[1] - 2, 0, -1):
            data[i][j] = al[j + 1] * data[i][j + 1] + be[j + 1]

        # data[i][0] = data[i][1] - math.sqrt(dist) * 0

    return


def evident_metod(data, num, dist, a, t):
    tmp = (t * a) / dist
    time = 0
    for i in range(1, data.shape[0]):
        time += t
        # data[i][0] = data[i - 1][0]
        # data[i][num - 1] = data[i - 1][num - 1]
        data[i][0] = 0
        data[i][-1] = math.sin(time * 180 / math.pi)
        for j in range(1, data.shape[1] - 1):
            point = (data[i - 1][j + 1] + data[i - 1][j - 1])
            point = point - 2 * data[i - 1][j]
            data[i][j] = tmp * point + data[i - 1][j]

        # data[i][0] = data[i][1] - math.sqrt(dist) * 0
        # data[i][-1] = data[i][num - 2] + math.sqrt(dist) * math.cos(time * 180 / math.pi)

    return


def plot(data, xcord, ycord):
    plt.figure()
    plt.ion()
    for i in range(data.shape[0]):
        plt.plot(data[i])
        plt.xlim(0, xcord)
        plt.ylim(-ycord, ycord)
        plt.title('Temperature_changes')
        plt.grid(True)
        plt.pause(0.01)
        plt.clf()
    plt.show(block=True)
    return


if len(sys.argv) != 2:
    print ">Unexpected quantity of arguments, check your comand string.\n"
    exit(-1)

num = float(100)
# dist = float(0.1)
dist = float((50.0 - 0.0) / num)
print "DIST"
print dist
a = float(1000)
t = float(0.1)
fin = 10
dist = dist ** 2

if (int(sys.argv[1]) == 1) & (t >= dist / (2 * a)):
    print ">Your TAU >= (h ** 2) / (2 * a); It changed to (h ** 2) / (4 * a)"
    t = dist / (4 * a)
print "TAU"
print t
# Used to check convergence of evident metod

k = int(fin / t + 1)
data = np.zeros((k, int(num)))
data[0][0] = 0
data[0][-1] = 0
print "STEPS"
print data.shape[0]

if int(sys.argv[1]) == 1:
    evident_metod(data, int(num), dist, a, t)
else:
    reverse_metod(data, int(num), dist, a, t)

plot(data, int(num) + 10, 1)
