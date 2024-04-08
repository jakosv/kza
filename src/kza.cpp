#include "kza.h"
#include "KZGeneric.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <functional>
#include <cstring>
#include <concepts>

#ifdef TIMER
#include "timer.h"
#endif

#include <thread>
#include <barrier>

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZA: public KZGeneric<ValueT, SizeT, WinSizeT> {
public:

    KZA(WinSizeT window_size, WinSizeT min_window_size,
        ValueT tolerance, const ValueT *data, SizeT data_size, 
        const ValueT *kz_result) : 
            KZGeneric<ValueT, SizeT, WinSizeT>(window_size, data, data_size),
            min_window_size(min_window_size), tolerance(tolerance),
            kz_diff(data_size), kz_derivative(data_size)
    {
        this->differenced(kz_result, data_size, window_size);
    }

private:
    WinSizeT min_window_size;
    ValueT tolerance;
    std::vector<ValueT> kz_diff;
    ValueT max_diff;
    std::vector<ValueT> kz_derivative;

    using KZGeneric<ValueT, SizeT, WinSizeT>::window_size;
    using KZGeneric<ValueT, SizeT, WinSizeT>::data;
    using KZGeneric<ValueT, SizeT, WinSizeT>::ans;

    inline WinSizeT normalize_left_win(SizeT left_win, SizeT window_center)
    {
        SizeT new_size = (left_win < min_window_size)
                            ? min_window_size
                            : left_win;

        new_size = (new_size > window_center)
                    ? window_center
                    : new_size;

        return new_size;
    }

    inline WinSizeT normalize_right_win(SizeT right_win,
                                        SizeT window_center)
    {
        SizeT new_size = (right_win < min_window_size)
                            ? min_window_size
                            : right_win;

        SizeT max_right_win_size = data.size() - window_center - 1;
        new_size = (new_size > max_right_win_size)
                    ? max_right_win_size
                    : new_size;

        return new_size;
    }

    inline void calc_adaptive_window_sizes(WinSizeT &left_win_size, 
                                           WinSizeT &right_win_size, 
                                           SizeT t)  
    {
        SizeT adaptive_win_size = 
            std::floor(window_size * (1 - kz_diff[t]/max_diff));

        // derivative[t] = 0
        if (std::fabs(kz_derivative[t]) < tolerance) {
            left_win_size = 
                this->normalize_left_win(adaptive_win_size, t);
            right_win_size =
                this->normalize_right_win(adaptive_win_size, t);
        } else if (kz_derivative[t] < 0) {
            left_win_size = 
                this->normalize_left_win(adaptive_win_size, t);
            right_win_size = this->normalize_right_win(window_size, t);
        } else {
            right_win_size =
                this->normalize_right_win(adaptive_win_size, t);
            left_win_size = this->normalize_left_win(window_size, t);
        }
    }

    void calc_kz_diff(const ValueT *y, SizeT n, WinSizeT q)
    {
        // calculate d = |y(i+q) - y(i-q)|
        for (SizeT i = 0; i < q; ++i) {
            kz_diff[i] = std::fabs(y[i+q] - y[0]);
            max_diff = std::max(max_diff, kz_diff[i]);
        }
        for (SizeT i = q; i < n-q; ++i) {
            kz_diff[i] = std::fabs(y[i+q] - y[i-q]);
            max_diff = std::max(max_diff, kz_diff[i]);
        }
        for (SizeT i = n-q; i < n; ++i) {
            kz_diff[i] = std::fabs(y[n-1] - y[i-q]);
            max_diff = std::max(max_diff, kz_diff[i]);
        }
    }

    void differenced(const ValueT *kz, SizeT n, WinSizeT q)
    {
        this->calc_kz_diff(kz, n, q);

        /* d'(t) = d(i+1)-d(i) */
        for (SizeT i = 0; i < n-1; i++)
            kz_derivative[i] = kz_diff[i+1] - kz_diff[i];
        kz_derivative[n-1] = kz_derivative[n-2];
    }

    void perform_single_iteration(SizeT start_idx, SizeT end_idx)
    {
        WinSizeT left_win_size, right_win_size;

        for (SizeT time = start_idx; time <= end_idx; ++time) {
            this->calc_adaptive_window_sizes(left_win_size,
                                             right_win_size,
                                             time);

            ans[time] = this->average(time - left_win_size,
                                      time + right_win_size);
        }
    }

};

static double *kza1d(const double *data, int data_size, const double *kz_res, 
                     int win_size, int iterations, int min_win_size, 
                     double tolerance)
{
    KZA<double, size_t, unsigned short> kza_data(win_size, min_win_size,
                                                 tolerance,
                                                 data, data_size,
                                                 kz_res);

    kza_data.perform_iterations(iterations);

    return kza_data.get_ans();
}

double *kza(const double *x, int dim, const int *size, const double *y, 
            const int *window, int iterations, double min_window, 
            double tolerance)
{
    double *kz_ans = NULL;
    double *kza_ans = NULL;

    if (!x || !size || !window) {
        std::cerr << "kza: Incorrect params\n";
        return NULL;
    }

    if (dim > 1) {
        std::cerr << "kza: Not yet implemented\n";
        return NULL;
    }

    int data_size = size[0];

    if (!y) {
        kz_ans = kz(x, dim, size, window, iterations);
    } else {
        try {
            kz_ans = new double[data_size];
        } 
        catch (const std::bad_alloc&) {
            return nullptr; 
        }
        std::copy(y, y + data_size, kz_ans);
    }

#ifdef TIMER
    Timer timer;
    timer.start();
#endif
    kza_ans = kza1d(x, data_size, kz_ans, window[0], iterations, min_window,
                    tolerance);
#ifdef TIMER
    timer.stop();
    timer.print_elapsed("kza(): ");
#endif

    kz_free(kz_ans);

    return kza_ans;
}

void kza_free(double *data)
{
    kz_free(data);
}
