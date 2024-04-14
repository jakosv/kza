#include "../kza.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/numpy.h>
#include <vector>

namespace py = pybind11;

std::vector<double> kz1d_wrap(const std::vector<double> &x,
                              unsigned window, size_t iterations)
{
    KZ1D<double, size_t, unsigned> kz1d_algo(window, iterations);

    return kz1d_algo(x);
}

py::array_t<double> to_matrix(std::vector<std::vector<double>> &vec2d)
{
    size_t N = vec2d.size();
    size_t M = vec2d[0].size();

    py::array_t<float, py::array::c_style> arr({N, M});

    auto ra = arr.mutable_unchecked();

    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            ra(i, j) = vec2d[i][j];
        }
    }

    return arr;
}

std::vector<std::vector<double>> kz2d_wrap(
                            const std::vector<std::vector<double>> &x,
                            const std::vector<unsigned> &window, 
                            size_t iterations)
{
    KZ2D<double, size_t, unsigned> kz2d_algo(window, iterations);

    return kz2d_algo(x);
}

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
    m.def("kz2d", &kz2d_wrap, "A kz2d wrapping function");
    m.def("kza", &kza_wrap, "A kza wrapping function");
    m.def("to_matrix", &to_matrix);
}
