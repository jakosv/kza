# Kolmogorov-Zurbenko Adaptive Filter

## Build 
CMake has multiple options:
- prefix_sum - prefix sum optimization;
- kza_timer - time measurement;
- python_library - build pybind11 module

Build shared library:
```console
$ cmake src
$ make
```
The target libKZ.so file is located in the build directory

To build python module pybind11 is needed
If you dont want to provide pybind11 path in CMakeLists.txt, you can 
install it globaly with:
```console
$ pip install "pybind11[global]"
```
If you get the following error: "Imported target "pybind11::module" 
includes non-existent path", run:
```console
$ sudo apt install python[version]-dev 
```
minimum version 3.6 ( recommended 3.11 ) 

Build bybind11 module:
```console
$ cmake src -Dpython_library=true
$ make
```


## Usage examples
Makefile builds shared library. For example, we can use it in python.

```console
$ cmake src -Dpython_library=true -Dprefix_sum=true
$ make
$ python examples/kz_example.py 
$ python examples/kza_example.py
```

## Testing 

```console
$ cmake src -Dpython_library=true -Dkza_timer=true
$ make
$ python tests/run_tests.py
```
