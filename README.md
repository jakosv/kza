# Kolmogorov-Zurbenko Adaptive Filter

## Pybind11 installation
To build this project pybind11 is needed

If you dont want to provide pybind11 path in CMakeLists.txt, you can install it globaly with:
```console
$ pip install "pybind11[global]"
```

## Build 
CMake has multiple options:
- prefix_sum - prefix sum optimization;
- kza_timer - time measurement;

From _kza/src_ directory:

```console
$ cmake . -Dkza_timer=true -Dprefix_sum=true
$ make
```

If you get the following error: "Imported target "pybind11::module" includes non-existent path", run:
```console
$ sudo apt install python[version]-dev 
```
minimum version 3.6 ( recommended 3.11 ) 

Current target directory for the shared library is set to _kza/tests_

## Usage 
Makefile builds shared library. For example, we can use it in python.

```console
$ python kza_example.py
$ python kz_example.py 
```

## Testing 

```console
$ python run_tests.py
```
