import math
import matplotlib.pyplot as plt
import numpy as np


def func(x):
    # return math.exp(x)
    # return math.sin(x)
    # return math.cos(x)
    # return math.sin(x ** 2 + 3)
    # return math.sin(3 * x) - 5 * math.sin(2 * x / 3)
    # return x
    # return math.sqrt(1 - x ** 2)
    # return (x ** 3 + 4 * (x ** 2) + 7)
    # return math.exp(- 1 / (x ** 2))
    # return math.sin(1 / x)
    # return math.sin(math.log(math.fabs(x)))
    # return math.exp(x) * (1 / (x - 5))

    if (x < 0):
        return -1
    else:
        return 1


def plot(data, eps):
    plt.figure()
    plt.ion()
    cords = np.arange(-eps, eps, 2 * eps / data.shape[1])
    for i in range(data.shape[0]):
        plt.title('Aproximation Fouire')
        plt.xlim(-eps, eps)
        plt.ylim(data.min(), data.max())
        # Asixs
        plt.plot([-eps, eps], [0, 0], 'k-', lw=1)
        plt.plot([0, 0], [data.min(), data.max()], 'k-', lw=1)
        # Func
        plt.plot(cords, data[0])
        # Aprox
        plt.plot(cords, data[i])
        plt.grid(True)
        plt.pause(0.5)
        plt.clf()
        # Uncomment/comment code-string above to draw different plots
    plt.show(block=True)
    return


def integralAn(eps, n):
    sum = 0.0
    step = (2 * eps) / point_num
    for x in range(0, xcord.shape[0]):
        sum += func(xcord[x]) * math.cos(math.pi * n * xcord[x] / eps) * step

    return sum / eps


def integralBn(eps, n):
    sum = 0.0
    step = (2 * eps) / point_num
    for x in range(0, xcord.shape[0]):
        sum += func(xcord[x]) * math.sin(math.pi * n * xcord[x] / eps) * step

    return sum / eps


def fouire(derivative, eps):
    for i in range(2, derivative.shape[0]):
        for x in range(0, derivative.shape[1]):
            an = integralAn(eps, (i - 1)) * math.cos(math.pi * (i - 1) * xcord[x] / eps)
            bn = integralBn(eps, (i - 1)) * math.sin(math.pi * (i - 1) * xcord[x] / eps)
            derivative[i][x] = derivative[i - 1][x] + an + bn


eps = 5.0
point_num = 100
der_num = 51
derivative = np.zeros((der_num, point_num))
xcord = np.zeros(point_num)

step = -eps

for x in range(0, point_num):
    xcord[x] = step
    step += (2 * eps) / point_num

for x in range(0, point_num):
    derivative[0][x] = func(xcord[x])
    derivative[1][x] = integralAn(eps, 0) / 2

fouire(derivative, eps)

plot(derivative, eps)
