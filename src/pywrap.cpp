#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "kza.h"
#include <vector>

namespace py = pybind11;
constexpr auto byref = py::return_value_policy::reference_internal;

std::vector<double> kz_wrap(std::vector<double> x, int dim, const int *size, const int *window, int iterations) {
    auto res = kz(&x[0], dim, size, window, iterations);
    return std::vector<double>(res,res + *size);
}

std::vector<double> kza_wrap(std::vector<double> x, int dim, const int *size, std::vector<double> y, 
         const int *window, int iterations, double min_window, double tol) {
    auto res = kza(&x[0], dim, size, &y[0], window, iterations, min_window, tol);
    return std::vector<double>(res,res + *size);
}

PYBIND11_MODULE(KZ_Lib, m) {
    m.doc() = "optional module docstring";

    /*py::class_<MyClass>(m, "MyClass")
    .def(py::init<double, double, int>())  
    .def("run", &MyClass::run, py::call_guard<py::gil_scoped_release>())
    .def_readonly("v_data", &MyClass::v_data, byref)
    .def_readonly("v_gamma", &MyClass::v_gamma, byref)
    ;*/
    m.def("kz", &kz_wrap, "A function ...");
    m.def("kza", &kza_wrap, "A function ...");
    //m.def("kza_free", &kza_free, "A function ...");
}