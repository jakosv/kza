# Kolmogorov-Zurbenko Adaptive Filter

## Build 
CMake has multiple options:
- `-DPREFIX_SUM` - prefix sum optimization;
- `-DTIMER` - time measurement;



## Usage examples
Makefile builds shared library. For example, we can use it in python.
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
$ cd tests 
$ cmake .
$ make
```
Run examples:
```console
$ python examples/kz1d_example.py 
$ python examples/kz2d_example.py 
$ python examples/kza1d_example.py
$ python examples/kza2d_example.py
```

## Testing 

```console
$ cd tests 
$ cmake . -Dkza_timer=true
$ make
$ python run_tests.py
```
