# Kolmogorov-Zurbenko Adaptive Filter
Header-only library with implementation of
Adaptive Kolmogorov-Zurbenko filter for different dimensions.

## KZ and KZA filters
The Kolmogorov-Zurbenko (KZ) filter is a low-pass filter designed to remove noise from time series and two- or three-dimensional data. It is based on a simple moving average filter and relies on two parameters: the smoothing window size and the number of iterations of the MA filter.

The Kolmogorov-Zurbenko adaptive (KZA) filter is an extension of the KZ filter. It was designed to detect breaks in both one-dimensional and high-dimensional data by adjusting the size of the KZ smoothing window when a break occurs.

[Documentation](https://jakosv.github.io/kza/html/md_README.html)

Theory: docs/kza_project.pdf

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
[Vadim Marchenko](https://github.com/jakosv) 
<jakosvadim at gmail dot com> - worked on kza.hpp

[Polina Zolotareva](https://github.com/polin-drom) / 
[Nataliia Shalimova](https://github.com/LostOwlNata) - 
worked on docs/kza_project.pdf

[Kirill Avdeichuk](https://github.com/DotaSlaer) - worked on tests

Scientific supervisor: Denis Vasilyevich Parfenov (promasterden@yandex.ru)

