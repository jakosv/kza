from os import path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

import sys
sys.path.append(path.dirname(path.dirname(__file__)))

import libKZ_py

def kz2d(x, m, k = 3):
    xp = x
    dim = 1
    window = m
    res = np.array(libKZ_py.kz2d(xp, window, k))
    return res

def plot_surface(x, y, z):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    X, Y = np.meshgrid(x, y)
    ax.plot_surface(X, Y, z, cmap=cm.coolwarm)
    ax.set_xlabel('X Label')
    ax.set_ylabel('Y Label')
    ax.set_zlabel('Z Label')

    plt.show()


# 2 dimensions
np.random.seed(1)
a = np.zeros((100, 100))
a[35:71:, 35:71:] = 1
noise = np.random.normal(loc=0, scale=1.0, size=(100, 100))
z = a + noise
x = np.arange(0, 100, 1);
y = np.arange(0, 100, 1);

plot_surface(x, y, a)

plot_surface(x, y, z)

filtered = kz2d(z, m=[20, 5])

plot_surface(x, y, filtered)
