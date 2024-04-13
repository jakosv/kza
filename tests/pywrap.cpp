#include "../kza.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>

namespace py = pybind11;

std::vector<double> kz1d_wrap(const std::vector<double> &x,
                              unsigned window, size_t iterations)
{
    KZ1D<double, size_t, unsigned> kz1d_algo(window, iterations);

    return kz1d_algo(x);
}

/*
std::vector<std::vector<double>> kz2d_wrap(
                            const std::vector<std::vector<double>> &x,
                            const std::vector<unsigned> &window, 
                            size_t iterations)
{
    KZ2D<double, size_t, unsigned> kz2d_algo(window, iterations);

    return kz2d_algo(x);
}
*/

std::vector<double> kza_wrap(const std::vector<double> &x,
                             const std::vector<double> &y, 
                             unsigned window, unsigned iterations,
                             unsigned min_window, double tol)
{
    KZA1D<double, size_t, unsigned> kza1d_algo(window, min_window,
                                              tol, iterations);

    return kza1d_algo(x, y);
}

PYBIND11_MODULE(libKZ_py, m) {
    m.doc() = "Kolmogorov-Zurbenko Adaptive Filter";

    m.def("kz1d", &kz1d_wrap, "A kz1d wrapping function");
    //m.def("kz2d", &kz2d_wrap, "A kz2d wrapping function");
    m.def("kza", &kza_wrap, "A kza wrapping function");
}
