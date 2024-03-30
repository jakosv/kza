from os import listdir, path 
from os.path import isfile, join
import numpy as np

import sys
sys.path.append(path.dirname(path.dirname(__file__))+"/build")

import KZ_Lib

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

def kz1d(x, m, k = 3):
    xp = x
    dim = 1
    size = len(x)
    window = m
    res = KZ_Lib.kz(xp, dim, size, window, k)
    ans = res[:len(x)]
    #KZ_Lib.kza_free(res) res is copy
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
    res = KZ_Lib.kza(xp, dim, size, yp, window, k,
                     min_size, tol)
    ans = res[:len(x)]
    #KZ_Lib.kza_free(res) res is copy
    return ans

def run_func_tests(func_name, path):
    tests_path = join(path, func_name + '/')
    input_files = get_input_files(tests_path)

    for f in input_files:
        file_name = f.split('.')[0]

        x = np.loadtxt(join(tests_path, f))
        ans = np.loadtxt(join(tests_path, file_name + '.out'))
        res = None 

        if func_name == "kz1d":
            win_size = 30
            iterations = 3
            res = kz1d(x, win_size, k=iterations)
        elif func_name == "kza1d":
            win_size = 365
            iterations = 3
            res = kza1d(x, win_size, k=iterations, min_size=10)

        if not np.allclose(res, ans):
            diff = res - ans
            print(diff)
            print("test "+str(func_name)+"/"+str(file_name)+" failed")
            

def main():
    script_path = path.dirname(path.realpath(__file__))

    tests_path = script_path
    run_func_tests("kz1d", tests_path)
    run_func_tests("kza1d", tests_path)

if __name__ == "__main__":
    main()