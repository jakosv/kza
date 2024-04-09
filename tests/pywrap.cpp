#include "../kza.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

std::vector<double> kz_wrap(std::vector<double> x, size_t dim,
                            const size_t *size, const unsigned *window, 
                            size_t iterations)
{
    auto res = kz<double, size_t, unsigned>(&x[0], dim, size, window, iterations);
    return std::vector<double>(res, res + size[0]);
}

std::vector<double> kza_wrap(std::vector<double> x, unsigned dim,
                             const size_t *size, std::vector<double> y, 
                             const unsigned *window, unsigned iterations,
                             unsigned min_window, double tol)
{
    auto res = kza<double, size_t, unsigned>(&x[0], dim, size, &y[0],
                                               window, min_window,
                                               iterations, tol);
    return std::vector<double>(res, res + size[0]);
}

PYBIND11_MODULE(libKZ_py, m) {
    m.doc() = "Kolmogorov-Zurbenko Adaptive Filter";

    m.def("kz", &kz_wrap, "A kz wrapping function");
    m.def("kza", &kza_wrap, "A kza wrapping function");
}
