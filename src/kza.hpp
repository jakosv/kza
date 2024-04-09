#pragma once

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


typedef std::barrier<std::function<void()>> barrier_t;

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZGeneric {
public:

#ifdef PREFIX_SUM
    KZGeneric(WinSizeT window_size, const ValueT *data, SizeT data_size)
                : window_size(window_size),
                  ans(data_size), data(data, data + data_size),
                  pref_sum(data_size + 1), pref_finite_cnt(data_size + 1)
    {
        threads_cnt = std::thread::hardware_concurrency();
  #ifdef DEBUG
        printf("Number of cores: %d\n", threads_cnt);
  #endif

        task_size = (data_size + threads_cnt-1) / threads_cnt;

        update_prefix_sum(this->data); 
    }

#else

    KZGeneric(WinSizeT window_size, const ValueT *data, SizeT data_size)
                    : window_size(window_size),
                      ans(data_size), data(data, data + data_size)
    {
        threads_cnt = std::thread::hardware_concurrency();
#ifdef DEBUG
        printf("Number of cores: %d\n", threads_cnt);
#endif

        task_size = (data_size + threads_cnt-1) / threads_cnt;
    }

#endif

    virtual ~KZGeneric() = default;

    void perform_iterations(SizeT iterations)
    {

        if (iterations <= 0)
            return;

        auto on_iteration_complete = [this]() {
#ifdef PREFIX_SUM
            this->update_prefix_sum(ans);
#else
            std::swap(data, ans);
#endif
        };

        barrier_t sync_iteration(threads_cnt, on_iteration_complete);
        std::vector<ThreadData> th(threads_cnt);

        this->start_threads(th, iterations, sync_iteration);

        for (auto &thread_data : th)
            thread_data.thread.join();
    }

    ValueT *get_ans()
    {
        ValueT *answer;
        try {
            answer = new ValueT[data.size()];
        } 
        catch (const std::bad_alloc&) {
            return nullptr; 
        }

        std::copy(ans.begin(), ans.end(), answer);

        return answer;
    }

protected:
    struct ThreadData {
        std::thread thread;
        SizeT start_idx, end_idx;
        SizeT iterations;
    };

    WinSizeT window_size;
    SizeT task_size;
    SizeT threads_cnt;
    std::vector<ValueT> ans;
    std::vector<ValueT> data;

#ifdef PREFIX_SUM
    std::vector<ValueT> pref_sum;
    std::vector<SizeT> pref_finite_cnt;

    void update_prefix_sum(const std::vector<ValueT> &data)
    {
        pref_sum[0] = 0;
        pref_finite_cnt[0] = 0;

        for (SizeT i = 1; size_t(i) <= data.size(); ++i) {
            bool is_finite_flag = std::isfinite(data[i-1]);
            pref_sum[i] = pref_sum[i-1] + is_finite_flag * data[i-1];
            pref_finite_cnt[i] = pref_finite_cnt[i-1] + is_finite_flag;
        }
    }
#endif

    virtual void perform_single_iteration(SizeT start_idx, SizeT end_idx) {}

    void worker(ThreadData &t_data,
                SizeT iterations,
                barrier_t &sync_iteration)
    {
        for (SizeT i = 0; i < iterations - 1; ++i) {
            this->perform_single_iteration(t_data.start_idx, t_data.end_idx);

            // waiting for other threads to complete iteration
            sync_iteration.arrive_and_wait();
        }
        // perform last iteration and finish
        this->perform_single_iteration(t_data.start_idx, t_data.end_idx);
    }

    void start_threads(std::vector<ThreadData> &th,
                       SizeT iterations,
                       barrier_t &sync_iteration)
    {
        SizeT start_idx = 0;

        for (auto &t_data : th) {
            t_data.start_idx = start_idx;
            t_data.end_idx = (size_t(start_idx + task_size) >= data.size())
                                ? data.size() - 1
                                : start_idx + task_size;

            start_idx = t_data.end_idx;

            t_data.thread = std::thread(&KZGeneric::worker, this, 
                                        std::ref(t_data), iterations,
                                        std::ref(sync_iteration));
        }
    }

#ifdef PREFIX_SUM

    ValueT average(SizeT start_idx, SizeT end_idx)
    {
        /* (window sum) = (sum containig window) - (sum before window) */
        SizeT z =
            (pref_finite_cnt[end_idx + 1] - pref_finite_cnt[start_idx]);

        if (z == 0)
            return std::numeric_limits<ValueT>::quiet_NaN();

        return (pref_sum[end_idx+1] - pref_sum[start_idx]) / z;
    }

#else

    ValueT average(SizeT start_idx, SizeT end_idx)
    {
        SizeT z = 0;
        ValueT s = 0;

        for (SizeT i = start_idx; i <= end_idx; ++i) {
            bool is_finite_flag = std::isfinite(data[i]);
            z += is_finite_flag;
            s += is_finite_flag * data[i];
        }

        if (z == 0)
            return std::numeric_limits<ValueT>::quiet_NaN();

        return s / z;
    }

#endif
};


template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZ: public KZGeneric<ValueT, SizeT, WinSizeT> {
public:
    KZ(WinSizeT window_size, const ValueT *data, SizeT data_size) : 
        KZGeneric<ValueT, SizeT, WinSizeT>(window_size, data, data_size)
    {}

private:
    using KZGeneric<ValueT, SizeT, WinSizeT>::window_size;
    using KZGeneric<ValueT, SizeT, WinSizeT>::data;
    using KZGeneric<ValueT, SizeT, WinSizeT>::ans;

    // window length is 2*window_size+1
    inline SizeT window_left_bound(SizeT win_center)
    {
        return (win_center >= SizeT(window_size))
                ? (win_center - window_size)
                : 0;
    }

    inline SizeT window_right_bound(SizeT win_center)
    {
        return (size_t(win_center + window_size) >= data.size())
                ? data.size() - 1
                : win_center + window_size;
    }

    inline void perform_single_iteration(SizeT start_idx,
                                         SizeT end_idx)
    {
        for (SizeT time = start_idx; time <= end_idx; ++time)
            ans[time] = this->average(this->window_left_bound(time),
                                      this->window_right_bound(time));
    }
};

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
ValueT *kz1d(const ValueT *data, SizeT data_size, SizeT win_size,
             SizeT iterations)
{
    KZ<ValueT, SizeT, WinSizeT> kz_data(win_size, data, data_size);

    kz_data.perform_iterations(iterations);

    return kz_data.get_ans();
}

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
ValueT *kz(const ValueT *data, SizeT dimention, const SizeT *data_sizes,
           const WinSizeT *window_sizes, SizeT iterations)
{
    ValueT *ans = NULL;

#ifdef TIMER
    Timer timer;
#endif

    switch (dimention) {
    case 1:
#ifdef TIMER
        timer.start();
#endif
        ans = kz1d<ValueT, SizeT, WinSizeT>(data, data_sizes[0],
                                            window_sizes[0]>>1,
                                            iterations); 
#ifdef TIMER
        timer.stop();
        timer.print_elapsed("kz(): ");
#endif
        break;
    default:
        fprintf(stderr, "kza: Too many dimensions\n");
        break;
    }

    return ans;
}


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


template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
static ValueT *kza1d(const ValueT *data, SizeT data_size, const ValueT *kz_res, 
                     SizeT win_size, SizeT min_win_size, SizeT iterations, 
                     ValueT tolerance)
{
    KZA<ValueT, SizeT, WinSizeT> kza_data(win_size, min_win_size, tolerance,
                                          data, data_size, kz_res);

    kza_data.perform_iterations(iterations);

    return kza_data.get_ans();
}

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
ValueT *kza(const ValueT *data, SizeT dimention, const SizeT *data_sizes,
            const ValueT *kz_filter_data,
            const WinSizeT *window_sizes, WinSizeT min_window, 
            SizeT iterations, ValueT tolerance)
{
    ValueT *kz_ans = NULL;
    ValueT *kza_ans = NULL;

    if (!data || !data_sizes || !window_sizes) {
        std::cerr << "kza: Incorrect params\n";
        return NULL;
    }

    if (dimention > 1) {
        std::cerr << "kza: Not yet implemented\n";
        return NULL;
    }

    SizeT data_size = data_sizes[0];

    if (!kz_filter_data) {
        kz_ans = kz<ValueT, SizeT, WinSizeT>(data, dimention, data_sizes, 
                                             window_sizes, iterations);
    } else {
        try {
            kz_ans = new ValueT[data_size];
        } 
        catch (const std::bad_alloc&) {
            return nullptr; 
        }
        std::copy(kz_filter_data, kz_filter_data + data_size, kz_ans);
    }

#ifdef TIMER
    Timer timer;
    timer.start();
#endif
    kza_ans = kza1d<ValueT, SizeT, WinSizeT>(data, data_size, kz_ans,
                                             window_sizes[0], min_window,
                                             iterations, tolerance);
#ifdef TIMER
    timer.stop();
    timer.print_elapsed("kza(): ");
#endif

    delete[] kz_ans;

    return kza_ans;
}
