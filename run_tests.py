from os import listdir, getcwd
from os.path import isfile, join
import ctypes as ct
import numpy as np

def get_dynamic_lib(path):
    libkza_path = path + "/libkza.so"
    libkza = ct.CDLL(libkza_path)
    return libkza

def get_kz_func(libkza):
    libkza.kz.argtypes = [
        ct.POINTER(ct.c_double), 
        ct.c_int, 
        ct.POINTER(ct.c_int), 
        ct.POINTER(ct.c_int), 
        ct.c_int
    ]
    libkza.kz.restype = ct.POINTER(ct.c_double)
    return libkza.kz

def get_kza_func(libkza):
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
    return libkza.kza

def get_kza_free_func(libkza):
    libkza.kza_free.argtypes = [ct.POINTER(ct.c_double)]
    libkza.kza_free.restype = None
    return libkza.kza_free

def get_input_files(path):
    input_files = []
    for f in listdir(path):
        if isfile(join(path, f)):
            file_suffixes = f.split('.')
            if len(file_suffixes) > 0:
                data_type = file_suffixes[1]
                if data_type == "in":
                    input_files.append(f)
    return input_files

def kz1d(x, m, k, kz, kza_free):
    xp = (ct.c_double * len(x))(*x)
    dim = 1
    size = (ct.c_int)(len(x))
    window = (ct.c_int)(m)
    res = kz(xp, dim, ct.byref(size), ct.byref(window), k)
    ans = res[:len(x)]
    kza_free(res)
    return ans

def kza1d(x, m, k, kza, kz, kza_free):
    min_size = 10
    tol = 1.0e-5
    y = kz1d(x, m, k, kz, kza_free) 
    xp = (ct.c_double * len(x))(*x)
    yp = (ct.c_double * len(y))(*y)
    dim = 1
    size = (ct.c_int)(len(x))
    window = (ct.c_int)(m)
    res = kza(xp, dim, ct.byref(size), yp, ct.byref(window), k, min_size, tol)
    ans = res[:len(x)]
    kza_free(res)
    return ans

def run_func_tests(func_name, path, funcs):
    tests_path = join(path, func_name + '/')
    input_files = get_input_files(tests_path)

    kz = funcs["kz"]
    kza = funcs["kza"]
    kza_free = funcs["kza_free"]

    for f in input_files:
        file_name = f.split('.')[0]

        x = np.loadtxt(join(tests_path, f))
        ans = np.loadtxt(join(tests_path, file_name + '.out'))
        res = None 

        if func_name == "kz1d":
            win_size = 30
            iterations = 3
            res = kz1d(x, win_size, iterations, kz, kza_free)
        elif func_name == "kza1d":
            win_size = 365
            iterations = 3
            res = kza1d(x, win_size, iterations, kza, kz, kza_free)

        if not np.allclose(res, ans):
            diff = res - ans
            print(diff)
            print(f"test [{func_name}/{file_name}] failed")
            

def main():
    script_path = getcwd()

    funcs = {
        "kz": None,
        "kza": None,
        "kza_free": None
    }

    libkza = get_dynamic_lib(script_path)
    funcs["kz"] = get_kz_func(libkza)
    funcs["kza"] = get_kza_func(libkza)
    funcs["kza_free"] = get_kza_free_func(libkza)

    tests_path = join(script_path, 'tests/')
    run_func_tests("kz1d", tests_path, funcs)
    run_func_tests("kza1d", tests_path, funcs)

if __name__ == "__main__":
    main()
