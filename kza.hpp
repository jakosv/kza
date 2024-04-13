#ifndef KZA_HPP_SENTRY
#define KZA_HPP_SENTRY

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <functional>
#include <cstring>
#include <concepts>
#include <thread>
#include <barrier>


typedef std::barrier<std::function<void()>> barrier_t;

#ifdef TIMER
#include <chrono>

class Timer
{
public:
    inline void start()
    {
        begin = std::chrono::steady_clock::now();
    }

    inline void stop()
    {
        end = std::chrono::steady_clock::now();
    }

    inline long elapsed_time()
    {
        auto diff = end - begin;
        return std::chrono::duration_cast<std::chrono::milliseconds>(diff).count();
    }

    inline void print_elapsed(const char *msg)
    {
        long elapsed = this->elapsed_time();
        std::cout << msg << elapsed << "ms" << std::endl;
    }

private:
    std::chrono::steady_clock::time_point begin, end;
};
#endif


template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZ
{
public:
    KZ(SizeT iterations): iterations(iterations)
    {
        threads_cnt = std::thread::hardware_concurrency();
#ifdef DEBUG
        printf("Number of cores: %d\n", threads_cnt);
#endif
    }

    virtual ~KZ() = default;

    inline void set_iterations(SizeT iterations)
    {
        this->iterations = iterations;
    }

protected:
    struct ThreadData {
        std::thread thread;
        SizeT start_idx, end_idx;
    };

    SizeT iterations;
    SizeT task_size;
    SizeT threads_cnt;

    virtual void perform_single_iteration(SizeT start_idx, SizeT end_idx) = 0;

#ifdef PREFIX_SUM
    virtual void update_prefix_sum(const std::vector<ValueT> &data) = 0;
#endif

    void perform_iterations(std::function<void()> on_iteration_complete)
    {
        barrier_t sync_iteration(threads_cnt, on_iteration_complete);
        std::vector<ThreadData> th(threads_cnt);

#ifdef TIMER
    Timer timer;
    timer.start();
#endif
        this->start_threads(th, sync_iteration);

        for (auto &thread_data : th)
            thread_data.thread.join();
#ifdef TIMER
    timer.stop();
    timer.print_elapsed("time: ");
#endif
    }

    void worker(ThreadData &t_data, barrier_t &sync_iteration)
    {
        for (SizeT i = 0; i < iterations - 1; ++i) {
            this->perform_single_iteration(t_data.start_idx, t_data.end_idx);

            // waiting for other threads to complete iteration
            sync_iteration.arrive_and_wait();
        }
        // perform last iteration and finish
        this->perform_single_iteration(t_data.start_idx, t_data.end_idx);
    }

    void start_threads(std::vector<ThreadData> &th, barrier_t &sync_iteration)
    {
        SizeT thread_task_size = (task_size + threads_cnt-1) / threads_cnt;
        SizeT start_idx = 0;

        for (auto &t_data : th) {
            t_data.start_idx = start_idx;
            t_data.end_idx = (start_idx + thread_task_size >= task_size)
                                ? task_size - 1
                                : start_idx + thread_task_size;

            start_idx = t_data.end_idx;

            t_data.thread = std::thread(&KZ::worker, this, std::ref(t_data),
                                        std::ref(sync_iteration));
        }
    }
};


template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZ1D: public KZ<ValueT, SizeT, WinSizeT>
{
public:
    KZ1D(WinSizeT window_size, SizeT iterations) 
        : KZ<ValueT, SizeT, WinSizeT>(iterations) 
    {
        half_window = window_size>>1; // window half (window is 2*m+1)
    }

    std::vector<ValueT> operator()(const std::vector<ValueT> &vec)
    {
        task_size = vec.size();

        data.assign(vec.cbegin(), vec.cend());
        ans.resize(vec.size());

#ifdef PREFIX_SUM
        pref_sum.resize(data.size() + 1);
        pref_finite_cnt.resize(data.size() + 1);
        update_prefix_sum(data); 
#endif

        if (iterations <= 0)
            return std::move(data);

        this->perform_iterations([this]() {
#ifdef PREFIX_SUM
            this->update_prefix_sum(ans);
#else
            std::swap(data, ans);
#endif
        });

        return std::move(ans);
    }

    inline void set_window(WinSizeT window)
    {
        half_window = window >> 1;
    }

protected:

#ifdef PREFIX_SUM
    std::vector<ValueT> pref_sum;
    std::vector<SizeT> pref_finite_cnt;

    void update_prefix_sum(const std::vector<ValueT> &data) override
    {
        pref_sum[0] = 0;
        pref_finite_cnt[0] = 0;

        for (SizeT i = 1; size_t(i) <= data.size(); ++i) {
            bool is_finite_flag = std::isfinite(data[i-1]);
            pref_sum[i] = pref_sum[i-1] + is_finite_flag * data[i-1];
            pref_finite_cnt[i] = pref_finite_cnt[i-1] + is_finite_flag;
        }
    }

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

private:

    // window length is 2*half_window+1
    inline SizeT window_left_bound(SizeT win_center)
    {
        return (win_center >= SizeT(half_window))
                ? (win_center - half_window)
                : 0;
    }

    inline SizeT window_right_bound(SizeT win_center)
    {
        return (size_t(win_center + half_window) >= data.size())
                ? data.size() - 1
                : win_center + half_window;
    }

    inline void perform_single_iteration(SizeT start_idx,
                                         SizeT end_idx) override
    {
        for (SizeT time = start_idx; time <= end_idx; ++time)
            ans[time] = this->average(this->window_left_bound(time),
                                      this->window_right_bound(time));
    }

protected:
    WinSizeT half_window;
    std::vector<ValueT> data;
    std::vector<ValueT> ans;

    using KZ<ValueT, SizeT, WinSizeT>::task_size;
    using KZ<ValueT, SizeT, WinSizeT>::iterations;
    using KZ<ValueT, SizeT, WinSizeT>::threads_cnt;
};

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZ2D: public KZ<ValueT, SizeT, WinSizeT>
{
public:
    KZ2D(const std::vector<WinSizeT> window_sizes, SizeT iterations) 
        : KZ<ValueT, SizeT, WinSizeT>(iterations)
    {
        /*
        half_window.resize(window_sizes.size());
        for (int i = 0; i < 2; i++)
            half_window[i] >>= 1;
            */
        this->set_window(window_sizes);
    }

    std::vector<std::vector<ValueT>> operator()(
            const std::vector<std::vector<ValueT>> &vec)
    {
        task_size = vec.size();

        data = vec;
        rows = data.size();
        cols = data[0].size();
        ans.resize(rows, std::vector<ValueT>(cols));

#ifdef PREFIX_SUM
        pref_sum.resize(rows + 1, std::vector<ValueT>(cols + 1));
        pref_finite_cnt.resize(rows + 1, std::vector<ValueT>(cols + 1));
        update_prefix_sum(data); 
#endif

        if (iterations <= 0)
            return std::move(data);

        this->perform_iterations([this]() {
#ifdef PREFIX_SUM
            this->update_prefix_sum(ans);
#else
            std::swap(data, ans);
#endif
        });

        return std::move(ans);
    }

    inline void set_window(const std::vector<WinSizeT> &window)
    {
        half_window.resize(window.size());
        for (int i = 0; i < 2; i++)
            half_window[i] >>= 1;
    }

protected:

#ifdef PREFIX_SUM
    std::vector<std::vector<ValueT>> pref_sum;
    std::vector<std::vector<SizeT>> pref_finite_cnt;

    void update_prefix_sum(
            const std::vector<std::vector<ValueT>> &data) override
    {
        // calc 2D prefix sum implemetation
    }

    ValueT average(SizeT start_idx, SizeT end_idx)
    {
        // 2D prefix sum average implemetation
    }

#else

    ValueT average(SizeT start_row, SizeT end_row, 
                   SizeT start_col, SizeT end_col)
    {
        SizeT z = 0;
        ValueT s = 0;

        for (SizeT i = start_row; i <= end_row; ++i) {
            for (SizeT j = start_col; i <= end_col; ++j) {
                bool is_finite_flag = std::isfinite(data[i][j]);
                z += is_finite_flag;
                s += is_finite_flag * data[i][j];
            }
        }

        if (z == 0)
            return std::numeric_limits<ValueT>::quiet_NaN();

        return s / z;
    }

#endif

private:

    // window size is 2*half_window+1
    inline SizeT window_left_bound(SizeT win_center, WinSizeT half_win)
    {
        return (win_center >= SizeT(half_win))
                ? (win_center - half_win)
                : 0;
    }

    inline SizeT window_right_bound(SizeT win_center, WinSizeT half_win,
                                    SizeT data_size)
    {
        return (size_t(win_center + half_win) >= data_size)
                ? data_size - 1
                : win_center + half_win;
    }

    inline void perform_single_iteration(SizeT start_idx,
                                         SizeT end_idx) override
    {
        for (SizeT i = 0; i < rows; ++i)
            for (SizeT j = 0; j < cols; ++j)
                ans[i][j] = this->average(
                        this->window_left_bound(i, half_window[0]),
                        this->window_right_bound(i, half_window[0], rows),
                        this->window_left_bound(j, half_window[1]),
                        this->window_right_bound(j, half_window[1], cols));
    }

protected:
    std::vector<WinSizeT> half_window;
    std::vector<std::vector<ValueT>> data;
    std::vector<std::vector<ValueT>> ans;
    SizeT rows, cols;

    using KZ<ValueT, SizeT, WinSizeT>::task_size;
    using KZ<ValueT, SizeT, WinSizeT>::iterations;
    using KZ<ValueT, SizeT, WinSizeT>::threads_cnt;
};


template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZAExtension
{
public:
    inline void set_min_window(WinSizeT min_window)
    {
        min_window_size = min_window; 
    }

    inline void set_tolerance(ValueT tol)
    {
        tolerance = tol; 
    }

protected:
    KZAExtension(WinSizeT min_window_size, ValueT tolerance)
        : min_window_size(min_window_size), tolerance(tolerance) {}

    WinSizeT min_window_size;
    ValueT tolerance;
    std::vector<ValueT> kz_diff;
    ValueT max_diff;
    std::vector<ValueT> kz_derivative;
};


template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZA1D: public KZ1D<ValueT, SizeT, WinSizeT>,
             public KZAExtension<ValueT, SizeT, WinSizeT>
{
public:
    KZA1D(WinSizeT window_size, WinSizeT min_window_size,
          ValueT tolerance, SizeT iterations)
            : KZ1D<ValueT, SizeT, WinSizeT>(window_size, iterations),
              KZAExtension<ValueT, SizeT, WinSizeT>(min_window_size,
                                                    tolerance)
    {
        // unlike of KZ algorithm, KZA uses (2*window_size+1) window 
        half_window = window_size;           
    }

    std::vector<ValueT> operator()(const std::vector<ValueT> &data,
                                   const std::vector<ValueT> &kz_result)
    {
        kz_diff.resize(data.size());
        kz_derivative.resize(data.size());

        calc_kz_difference(kz_result);
        calc_kz_dirivative();

        // call KZ algorithm with overrided perform_single_iteration() method
        return KZ1D<ValueT, SizeT, WinSizeT>::operator()(data); 
    }

    inline void set_window(WinSizeT window)
    {
        half_window = window;           
    }

private:

    using KZ<ValueT, SizeT, WinSizeT>::iterations;

    using KZ1D<ValueT, SizeT, WinSizeT>::half_window;
    using KZ1D<ValueT, SizeT, WinSizeT>::data;
    using KZ1D<ValueT, SizeT, WinSizeT>::ans;

    using KZAExtension<ValueT, SizeT, WinSizeT>::min_window_size;
    using KZAExtension<ValueT, SizeT, WinSizeT>::tolerance;
    using KZAExtension<ValueT, SizeT, WinSizeT>::kz_diff;
    using KZAExtension<ValueT, SizeT, WinSizeT>::max_diff;
    using KZAExtension<ValueT, SizeT, WinSizeT>::kz_derivative;


    void perform_single_iteration(SizeT start_idx, SizeT end_idx) override
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
            std::floor(half_window * (1 - kz_diff[t]/max_diff));

        // derivative[t] = 0
        if (std::fabs(kz_derivative[t]) < tolerance) {
            left_win_size = 
                this->normalize_left_win(adaptive_win_size, t);
            right_win_size =
                this->normalize_right_win(adaptive_win_size, t);
        } else if (kz_derivative[t] < 0) {
            left_win_size = 
                this->normalize_left_win(adaptive_win_size, t);
            right_win_size = this->normalize_right_win(half_window, t);
        } else {
            right_win_size =
                this->normalize_right_win(adaptive_win_size, t);
            left_win_size = this->normalize_left_win(half_window, t);
        }
    }

    inline void calc_kz_difference(const std::vector<ValueT> &y)
    {
        // calculate d = |y(i+q) - y(i-q)|
        SizeT n = y.size();
        WinSizeT q = half_window;

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

    inline void calc_kz_dirivative()
    {
        /* d'(t) = d(i+1)-d(i) */
        SizeT n = kz_derivative.size();

        for (SizeT i = 0; i < n-1; i++)
            kz_derivative[i] = kz_diff[i+1] - kz_diff[i];
        kz_derivative[n-1] = kz_derivative[n-2];
    }
};


/*
template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
ValueT *kz(const ValueT *data, SizeT dim, const SizeT *data_sizes,
           const WinSizeT *window_sizes, SizeT iterations)
{
    ValueT *ans = NULL;

    if (dim == 1) {
        KZ1D<ValueT, SizeT, WinSizeT> kz1d_algo(win_size, iterations);

        std::vector<ValueT> vec(data, data + data_sizes[0]);
        vec = kz1d_algo(vec);
        ans = new ValueT[vec.size()];
        std::copy(vec.cbegin(), vec.cend(), ans);
    } else if (dim = 2) {
        KZ2D<ValueT, SizeT, WinSizeT> kz2d_algo(win_size, iterations);

        std::vector<std::vector<ValueT>> vec(
            std::vector<ValueT>(begin(data[0]), end(data[0])),
            std::vector<ValueT>(begin(data[1]), end(data[1])));
        vec = kz1d_algo(vec);
        ans = new ValueT[vec.size() * vec[0].size()];
        for (i = 0; i < vec.size(); ++i)
            std::copy(vec[i].cbegin(), vec[i].cend(), ans + i*vec[0].size());
    } else {
        std::cerr << "kz: Too many dimensions\n";
    }

    return ans;
}


template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
ValueT *kza(const ValueT *data, SizeT dimension, const SizeT *data_sizes,
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

    if (dimension > 1) {
        std::cerr << "kza: Not yet implemented\n";
        return NULL;
    }

    SizeT data_size = data_sizes[0];

    if (!kz_filter_data) {
        kz_ans = kz<ValueT, SizeT, WinSizeT>(data, dimension, data_sizes, 
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

    return kza_ans;
}
*/

#endif
