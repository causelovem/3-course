import math
import matplotlib.pyplot as plt
import numpy as np


def func(x):
    return math.exp(x)
    # return 2 ** x
    # return math.sin(x)


def plot(data, a, eps):
    plt.figure()
    plt.ion()
    cords = np.arange(a - eps / 2, a + eps / 2, eps / data.shape[1])
    for i in range(data.shape[0]):
        plt.title('Aproximation Taylor')
        plt.xlim(a - eps / 2, a + eps / 2)
        plt.ylim(data.min(), data.max())
        # plt.ylim(-2, 2)
        # Asixs
        plt.plot([a - eps / 2, a + eps / 2], [0, 0], 'k-', lw=1)
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


def taylor(derivative, a):
    for i in range(2, derivative.shape[0]):
        factor = math.factorial(i - 1)
        deriv = math.exp(a) / factor
        # deriv = (2 ** a) * (math.log(2) ** (i - 1)) / factor
        # deriv = math.sin(a + (math.pi * (i - 1)) / 2) / factor
        for x in range(0, derivative.shape[1]):
            derivative[i][x] = derivative[i - 1][x] + deriv * ((xcord[x] - a) ** (i - 1))


a = 1.0
eps = 5.0
point_num = 100
der_num = 51
derivative = np.zeros((der_num, point_num))
xcord = np.zeros(point_num)

step = a - eps / 2

for x in range(0, point_num):
    xcord[x] = step
    step += eps / point_num

for x in range(0, point_num):
    derivative[0][x] = func(xcord[x])
    derivative[1][x] = func(a)

taylor(derivative, a)

plot(derivative, a, eps)
