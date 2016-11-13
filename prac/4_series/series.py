import sys
import math
import matplotlib.pyplot as plt
import numpy as np


def plot(data, xcord, ycord):
    plt.figure()
    plt.ion()
    for i in range(data.shape[0]):
        plt.plot(data[i])
        plt.xlim(-xcord, xcord)
        plt.ylim(-ycord, ycord)
        plt.title('Temperature_changes')
        plt.grid(True)
        plt.pause(0.5)
        plt.clf()
# Uncomment code-string above to draw different plots
    plt.show(block=True)
    return


def taylor(data):
    # function 2 ^ x

    a = 1.0

    for x in range(0, 100):
        derivative[0][x] = math.log(a)
        # derivative[0][x] = 2 ** a

    for i in range(1, 10):
        factor = math.factorial(i)
        deriv = ((1 / ((-a) ** i)) * (-math.factorial(i - 1))) / factor
        for x in range(0, 100):
            derivative[i][x] = derivative[i - 1][x] + deriv * ((x - a) ** i)
            # derivative[i][x] = derivative[i - 1][x] + (math.log(2) ** i) * ((x - a) ** i) / factor


derivative = np.zeros((10, 100))

taylor(derivative)

plot(derivative, 50, 50)
