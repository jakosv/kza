#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "kza.h"
#include <vector>

namespace py = pybind11;

std::vector<double> kz_wrap(std::vector<double> x, int dim, const int *size, const int *window, int iterations) {
    auto res = kz(&x[0], dim, size, window, iterations);
    return std::vector<double>(res,res + *size);
}

std::vector<double> kza_wrap(std::vector<double> x, int dim, const int *size, std::vector<double> y, 
         const int *window, int iterations, double min_window, double tol) {
    auto res = kza(&x[0], dim, size, &y[0], window, iterations, min_window, tol);
    return std::vector<double>(res,res + *size);
}

PYBIND11_MODULE(libKZ_py, m) {
    m.doc() = "optional module docstring";

    m.def("kz", &kz_wrap, "A kz wrapping function");
    m.def("kza", &kza_wrap, "A kza wrapping function");
}