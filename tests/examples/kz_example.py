from os import path
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append(path.dirname(path.dirname(__file__)))

import libKZ_py

def kz1d(x, m, k = 3):
    xp = x
    dim = 1
    size = len(x)
    window = m
    res = libKZ_py.kz(xp, dim, size, window, k)
    ans = res[:len(x)]
    return ans

def kza1d(x, m, y = None, k = 3, min_size = None, tol = 1.0e-5,
        impute_tails = False):
    xp = x
    if y != None:
        yp = y
    else:
        yp = kz1d(x,m)
    if min_size == None:
        min_size = round(0.05*m)
    dim = 1
    size = len(x)
    window = m
    res = libKZ_py.kza(xp, dim, size, yp, window, k,
                     min_size, tol)
    ans = res[:len(x)]
    return ans

# Example
np.random.seed(1)
t = np.linspace(0, 20, 20*365)
e = np.random.normal(loc=0, scale=2.0, size=len(t))
y = np.sin(3*np.pi*t) + e
z = kz1d(y, 30)

plt.subplot(211)
plt.title("y = sin(3*pi*t) + noise")
plt.plot(t, y)
plt.ylim([-5, 5])

plt.subplot(212)
plt.title("KZ filter")
plt.plot(t, z)
plt.plot(t, np.sin(3*np.pi*t))
plt.ylim([-5, 5])
plt.show()
