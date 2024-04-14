# Kolmogorov-Zurbenko Adaptive Filter
Header-only library with implementation of
Adaptive Kolmogorov-Zurbenko filter for different dimensions.

[Documentation](
https://github.com/jakosv/kza/blob/main/doc/kza_project_documentation.pdf)

See also R [kza package](
https://cran.r-project.org/web/packages/kza/index.html)


## Build options
CMake has multiple options:
- `-DPREFIX_SUM` - prefix sum optimization;
- `-DTIMER` - time measurement;



## Usage
To use the library, you need to include the header file kza.hpp in 
your project. After that, you can use one of the algorithms, for 
example one-dimensional Kolmogorov-Zurbenko Adaptive filter:
```
KZ1D<double, size_t, unsigned> kz1d(window, iterations);
std::vector<double> kz_ans = kz1d(x);
KZA1D<double, size_t, unsigned> kz1d(window, min_win, tolerance, iterations);
std::vector<double> kza_ans = kza1d(x, kz_ans);
```

You can also run examples of using the library in Python.
To build a Python module, you need pybind11.
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

## Authors:
[Vadim Marchenko](https://github.com/jakosv) - worked on kza.hpp

[Polina Zolotareva](https://github.com/polin-drom) - 
worked on doc/kza_project_documentation.pdf

[Kirill Avdeichuk](https://github.com/DotaSlaer) - worked on tests

Scientific supervisor: Denis Vasilyevich Parfenov (promasterden@yandex.ru)

