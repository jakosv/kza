#include "kza.h"
#include "KZGeneric.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <functional>
#include <cstring>

#ifdef TIMER
#include "timer.h"
#endif

#include <thread>
#include <barrier>

template<class ValueType, class SizeType>
class KZA: public KZGeneric<ValueType, SizeType> {
public:

    KZA(SizeType window_size, SizeType min_window_size,
        ValueType tolerance, const ValueType *data, SizeType data_size, 
        const ValueType *kz_result) : 
            KZGeneric<ValueType, SizeType>(window_size, data, data_size),
            min_window_size(min_window_size), tolerance(tolerance),
            kz_diff(data_size), kz_derivative(data_size)
    {
        this->differenced(kz_result, data_size, window_size);
    }

private:
    SizeType min_window_size;
    ValueType tolerance;
    std::vector<ValueType> kz_diff;
    ValueType max_diff;
    std::vector<ValueType> kz_derivative;

    using KZGeneric<ValueType, SizeType>::window_size;
    using KZGeneric<ValueType, SizeType>::data;
    using KZGeneric<ValueType, SizeType>::ans;


    inline void normalize_left_window_size(SizeType &left_win, 
                                           SizeType window_center)
    {
        left_win = (left_win < min_window_size) ? min_window_size : left_win;
        left_win = (left_win > window_center) ? window_center : left_win;
    }

    inline void normalize_right_window_size(SizeType &right_win,
                                            SizeType window_center)
    {
        right_win = (right_win < min_window_size) ? 
                                    min_window_size : 
                                    right_win;

        /* check bounds */
        SizeType max_right_bound = data.size() - window_center - 1;
        right_win = 
            (right_win > max_right_bound) ? max_right_bound : right_win;

    }

    inline void calc_adaptive_window_sizes(SizeType &left_win_size, 
                                           SizeType &right_win_size, 
                                           SizeType t)  
    {
        SizeType adaptive_win_size = std::floor(window_size * 
                                     (1 - kz_diff[t]/max_diff));
        // derivative[t] = 0
        if (fabs(kz_derivative[t]) < tolerance) {
            left_win_size = adaptive_win_size;
            right_win_size = adaptive_win_size;
        } else if (kz_derivative[t] < 0) {
            right_win_size = window_size;
            left_win_size = adaptive_win_size;
        } else {
            right_win_size = adaptive_win_size;
            left_win_size = window_size;
        }
    }

    void differenced(const ValueType *y, SizeType n, SizeType q)
    {
        /* calculate d = |Z(i+q) - Z(i-q)| */
        for (SizeType i = 0; i < q; ++i) {
            kz_diff[i] = std::fabs(y[i+q] - y[0]);
            max_diff = std::max(max_diff, kz_diff[i]);
        }
        for (SizeType i = q; i < n-q; ++i) {
            kz_diff[i] = std::fabs(y[i+q] - y[i-q]);
            max_diff = std::max(max_diff, kz_diff[i]);
        }
        for (SizeType i = n-q; i < n; ++i) {
            kz_diff[i] = std::fabs(y[n-1] - y[i-q]);
            max_diff = std::max(max_diff, kz_diff[i]);
        }

        /* d'(t) = d(i+1)-d(i) */
        for (SizeType i = 0; i < n-1; i++)
            kz_derivative[i] = kz_diff[i+1] - kz_diff[i];
        kz_derivative[n-1] = kz_derivative[n-2];
    }

    void perform_single_iteration(SizeType start_idx, SizeType end_idx)
    {
        SizeType left_win_size, right_win_size;

        for (SizeType t = start_idx; t <= end_idx; ++t) {
            this->calc_adaptive_window_sizes(left_win_size, right_win_size, t);
            this->normalize_left_window_size(left_win_size, t);
            this->normalize_right_window_size(right_win_size, t);

            ans[t] = this->average(t - left_win_size, t + right_win_size);
        }
    }

};

static double *kza1d(const double *data, int data_size, const double *kz_res, 
                     int win_size, int iterations, int min_win_size, 
                     double tolerance)
{
    KZA<double, int> kza_data(win_size, min_win_size, tolerance, 
                              data, data_size, kz_res);

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
