from os import path
import numpy as np
import matplotlib.pyplot as plt

import sys
sys.path.append(path.dirname(path.dirname(__file__))+"/build")

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
#np.random.seed(1)
yrs = 20
m = 365
t = np.linspace(0, yrs, yrs*m)
e = np.random.normal(loc=0, scale=1.0, size=len(t))
trend = np.linspace(0, -1, len(t))

bkpt = 3452
brk = np.empty(len(t))
brk[:bkpt] = 0
brk[bkpt:] = 0.5

signal = trend + brk

y = np.sin(2*np.pi*t) + signal + e

kz = kz1d(y, m)

z = kza1d(y, m, y=kz, min_size=10)

plt.subplot(211)
plt.title("y = sin(2*pi*t) + trend + break point + noise")
xs = np.arange(0, len(t), 1)
plt.plot(xs, y)
plt.ylim([-3, 3])

plt.subplot(212)
plt.title("Signal and KZA Reconstruction")
plt.plot(xs, signal)
plt.plot(xs, z)
#plt.plot(t, np.sin(2*np.pi*t))
plt.ylim([-3, 3])
plt.show()
