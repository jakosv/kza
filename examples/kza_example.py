import os
import ctypes as ct
import numpy as np
import matplotlib.pyplot as plt

script_path = os.path.realpath(__file__)
libkza_path = os.path.dirname(script_path) + "/../src/libkza.so"

libkza = ct.CDLL(libkza_path)

libkza.kz.argtypes = [
    ct.POINTER(ct.c_double), 
    ct.c_int, 
    ct.POINTER(ct.c_int), 
    ct.POINTER(ct.c_int), 
    ct.c_int
]
libkza.kz.restype = ct.POINTER(ct.c_double)

libkza.kza.argtypes = [
    ct.POINTER(ct.c_double), 
    ct.c_int, 
    ct.POINTER(ct.c_int), 
    ct.POINTER(ct.c_double), 
    ct.POINTER(ct.c_int), 
    ct.c_int,
    ct.c_double,
    ct.c_double
]
libkza.kza.restype = ct.POINTER(ct.c_double)

libkza.kza_free.argtypes = [ct.POINTER(ct.c_double)]
libkza.kza_free.restype = None


def kz1d(x, m, k = 3):
    xp = (ct.c_double * len(x))(*x)
    dim = 1
    size = (ct.c_int)(len(x))
    window = (ct.c_int)(m)
    res = libkza.kz(xp, dim, ct.byref(size), ct.byref(window), k)
    ans = res[:len(x)]
    libkza.kza_free(res)
    return ans

def kza1d(x, m, y = None, k = 3, min_size = None, tol = 1.0e-5,
        impute_tails = False):
    xp = (ct.c_double * len(x))(*x)
    if y != None:
        yp = (ct.c_double * len(y))(*y)
    if min_size == None:
        min_size = round(0.05*m)
    dim = 1
    size = (ct.c_int)(len(x))
    window = (ct.c_int)(m)
    res = libkza.kza(xp, dim, ct.byref(size), yp, ct.byref(window), k,
                     min_size, tol)
    ans = res[:len(x)]
    libkza.kza_free(res)
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
