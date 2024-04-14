/*

    Kolmogorov-Zurbenko Adaptive Filter
    Copyright (C) 2024 Vadim Marchenko

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.


*/


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

    void update_prefix_sum(const std::vector<ValueT> &data)
    {
        pref_sum[0] = 0;
        pref_finite_cnt[0] = 0;
        bool is_finite_flag;

        for (SizeT i = 1; size_t(i) <= data.size(); ++i) {
            is_finite_flag = std::isfinite(data[i-1]);
            pref_sum[i] = pref_sum[i-1] + (is_finite_flag ? data[i-1] : 0);
            pref_finite_cnt[i] = pref_finite_cnt[i-1] + is_finite_flag;
        }
    }

    ValueT average(SizeT start_idx, SizeT end_idx)
    {
        /* (window sum) = (sum containig window) - (sum before window) */
        SizeT z =
            (pref_finite_cnt[end_idx + 1] - pref_finite_cnt[start_idx]);

        return (z != 0) ? (pref_sum[end_idx + 1] - pref_sum[start_idx]) / z
                        : std::numeric_limits<ValueT>::quiet_NaN();
    }

#else

    ValueT average(SizeT start_idx, SizeT end_idx)
    {
        SizeT z = 0;
        ValueT s = 0;
        bool is_finite_flag;

        for (SizeT i = start_idx; i <= end_idx; ++i) {
            is_finite_flag = std::isfinite(data[i]);
            z += is_finite_flag;
            s += (is_finite_flag ? data[i] : 0);
        }

        return (z != 0) ? s/z : std::numeric_limits<ValueT>::quiet_NaN();
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
        pref_sum.resize(rows + 1, std::vector<ValueT>(cols + 1, 0));
        pref_finite_cnt.resize(rows + 1, std::vector<SizeT>(cols + 1, 0));
        update_prefix_sum(data); 
#endif

        if (iterations <= 0)
            return std::move(data);

        this->perform_iterations([this]() {
#ifdef PREFIX_SUM
            this->update_prefix_sum(ans);
#else
            data.swap(ans);
#endif
        });

        return std::move(ans);
    }

    inline void set_window(const std::vector<WinSizeT> &window)
    {
        half_win_rows = window[0] >> 1;
        half_win_cols = window.size() > 1 ? window[1] >> 1 : half_win_rows;
    }

protected:

#ifdef PREFIX_SUM
    std::vector<std::vector<ValueT>> pref_sum;
    std::vector<std::vector<SizeT>> pref_finite_cnt;

    void update_prefix_sum(const std::vector<std::vector<ValueT>> &data)
    {
        bool is_finite_flag;

        for (SizeT i = 1; i <= rows; ++i) {
            for (SizeT j = 1; j <= cols; ++j) {
                is_finite_flag = std::isfinite(data[i-1][j-1]);
                pref_sum[i][j] = 
                    pref_sum[i-1][j] + pref_sum[i][j-1] - pref_sum[i-1][j-1];
                pref_sum[i][j] += (is_finite_flag ? data[i-1][j-1] : 0);
                pref_finite_cnt[i][j] = pref_finite_cnt[i-1][j] + 
                                        pref_finite_cnt[i][j-1] -
                                        pref_finite_cnt[i-1][j-1] +
                                        is_finite_flag;
            }
        }
    }

    ValueT average(SizeT start_row, SizeT end_row, 
                   SizeT start_col, SizeT end_col)
    {
        ValueT s = pref_sum[end_row+1][end_col+1] -
                  pref_sum[end_row+1][start_col] -
                  pref_sum[start_row][end_col+1] +
                  pref_sum[start_row][start_col];

        SizeT z = pref_finite_cnt[end_row+1][end_col+1] -
                  pref_finite_cnt[end_row+1][start_col] -
                  pref_finite_cnt[start_row][end_col+1] +
                  pref_finite_cnt[start_row][start_col];

        return (z != 0) ? s/z : std::numeric_limits<ValueT>::quiet_NaN();
    }

#else

    ValueT average(SizeT start_row, SizeT end_row, 
                   SizeT start_col, SizeT end_col)
    {
        SizeT z = 0;
        ValueT s = 0;

        for (SizeT i = start_row; i <= end_row; ++i) {
            for (SizeT j = start_col; j <= end_col; ++j) {
                bool is_finite_flag = std::isfinite(data[i][j]);
                z += is_finite_flag;
                s += (is_finite_flag ? data[i][j] : 0);
            }
        }

        return (z != 0) ? s/z : std::numeric_limits<ValueT>::quiet_NaN();
    }

#endif

    // window size is 2*half_window+1
    inline WinSizeT window_left_bound(SizeT win_center, WinSizeT half_win)
    {
        return (win_center >= SizeT(half_win))
                ? (win_center - half_win)
                : 0;
    }

    inline WinSizeT window_right_bound(SizeT win_center, WinSizeT half_win,
                                       SizeT data_size)
    {
        return (win_center + half_win >= data_size)
                ? data_size - 1
                : win_center + half_win;
    }


private:

    inline void perform_single_iteration(SizeT start_idx,
                                         SizeT end_idx) override
    {
        for (SizeT row = start_idx; row <= end_idx; ++row) {
            for (SizeT col = 0; col < cols; ++col) {
                ans[row][col] = this->average(
                    window_left_bound(row, half_win_rows),
                    window_right_bound(row, half_win_rows, rows),
                    window_left_bound(col, half_win_cols),
                    window_right_bound(col, half_win_cols, cols));
            }
        }
    }

protected:
    WinSizeT half_win_rows, half_win_cols;
    std::vector<std::vector<ValueT>> data, ans;
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

    inline SizeT adaptive_window(WinSizeT window, ValueT kz_d, ValueT max_d)
    {
        return (std::fabs(max_d) > 0) 
               ? std::floor(window * (1 - kz_d/max_d))
               : window;
    }

protected:
    KZAExtension(WinSizeT min_window_size, ValueT tolerance)
        : min_window_size(min_window_size), tolerance(tolerance) {}

    WinSizeT min_window_size;
    ValueT tolerance;
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
    std::vector<ValueT> kz_diff;
    ValueT max_diff;
    std::vector<ValueT> kz_derivative;

    using KZ<ValueT, SizeT, WinSizeT>::iterations;

    using KZ1D<ValueT, SizeT, WinSizeT>::half_window;
    using KZ1D<ValueT, SizeT, WinSizeT>::data;
    using KZ1D<ValueT, SizeT, WinSizeT>::ans;

    using KZAExtension<ValueT, SizeT, WinSizeT>::min_window_size;
    using KZAExtension<ValueT, SizeT, WinSizeT>::tolerance;


    void perform_single_iteration(SizeT start_idx, SizeT end_idx) override
    {
        WinSizeT left_win_size, right_win_size;

        for (SizeT time = start_idx; time <= end_idx; ++time) {
            this->calc_adaptive_half_windows(left_win_size, right_win_size,
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

    inline void calc_adaptive_half_windows(WinSizeT &left_win_size, 
                                           WinSizeT &right_win_size, 
                                           SizeT t)  
    {
        SizeT adaptive_size = this->adaptive_window(half_window, 
                                                    kz_diff[t], max_diff);

        // derivative[t] = 0
        if (std::fabs(kz_derivative[t]) < tolerance) {
            left_win_size = normalize_left_win(adaptive_size, t);
            right_win_size = normalize_right_win(adaptive_size, t);
        } else if (kz_derivative[t] < 0) {
            left_win_size = normalize_left_win(adaptive_size, t);
            right_win_size = normalize_right_win(half_window, t);
        } else {
            right_win_size = normalize_right_win(adaptive_size, t);
            left_win_size = normalize_left_win(half_window, t);
        }
    }

    inline void calc_kz_difference(const std::vector<ValueT> &y)
    {
        // calculate d = |y(i+q) - y(i-q)|
        SizeT n = y.size();
        WinSizeT q = half_window;
        max_diff = 0;

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

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZA2D: public KZ2D<ValueT, SizeT, WinSizeT>,
             public KZAExtension<ValueT, SizeT, WinSizeT>
{
public:
    KZA2D(const std::vector<WinSizeT> &window, WinSizeT min_window_size,
          ValueT tolerance, SizeT iterations)
            : KZ2D<ValueT, SizeT, WinSizeT>(window, iterations),
              KZAExtension<ValueT, SizeT, WinSizeT>(min_window_size,
                                                    tolerance)
    {
        // use KZA specific window size: (2*window+1) 
        this->set_window(window);
    }

    std::vector<std::vector<ValueT>> operator()(
                const std::vector<std::vector<ValueT>> &data,
                const std::vector<std::vector<ValueT>> &kz_result)
    {
        rows = data.size();
        cols = data[0].size();

        kz_dx.resize(rows, std::vector<ValueT> (cols));
        kz_dy.resize(rows, std::vector<ValueT> (cols));
        x_derivative.resize(rows, std::vector<ValueT> (cols));
        y_derivative.resize(rows, std::vector<ValueT> (cols));

        calc_kz_difference(kz_result);
        calc_kz_dirivative();

        // call KZ algorithm with overrided perform_single_iteration() method
        return KZ2D<ValueT, SizeT, WinSizeT>::operator()(data); 
    }

    inline void set_window(const std::vector<WinSizeT> &window)
    {
        // unlike of KZ algorithm, KZA uses (2*window+1) window 
        half_win_rows = window[0];           
        half_win_cols = window.size() > 1 ? window[1] : half_win_rows;           
    }

private:
    std::vector<std::vector<ValueT>> kz_dx, kz_dy, x_derivative, y_derivative;
    ValueT max_dx, max_dy;

    using KZ<ValueT, SizeT, WinSizeT>::iterations;

    using KZ2D<ValueT, SizeT, WinSizeT>::half_win_rows;
    using KZ2D<ValueT, SizeT, WinSizeT>::half_win_cols;
    using KZ2D<ValueT, SizeT, WinSizeT>::data;
    using KZ2D<ValueT, SizeT, WinSizeT>::ans;
    using KZ2D<ValueT, SizeT, WinSizeT>::cols;
    using KZ2D<ValueT, SizeT, WinSizeT>::rows;

    using KZAExtension<ValueT, SizeT, WinSizeT>::min_window_size;
    using KZAExtension<ValueT, SizeT, WinSizeT>::tolerance;


    void perform_single_iteration(SizeT start_idx, SizeT end_idx) override
    {
        std::pair<WinSizeT, WinSizeT> rows_win;
        std::pair<WinSizeT, WinSizeT> cols_win;

        for (SizeT row = start_idx; row <= end_idx; ++row) {
            for (SizeT col = 0; col < cols; ++col) {
                rows_win = adaptive_half_windows(row, half_win_rows, rows,
                                                 kz_dy[row][col], max_dy, 
                                                 y_derivative[row][col]);
                cols_win = adaptive_half_windows(col, half_win_cols, cols,
                                                 kz_dx[row][col], max_dx, 
                                                 x_derivative[row][col]);

                // check that area is not zero
                ans[row][col] = (rows_win.second + rows_win.first +
                                 cols_win.second + cols_win.first)

                                ? this->average(row - rows_win.first,
                                                row + rows_win.second,
                                                col - cols_win.first,
                                                col + cols_win.second)

                                : data[row][col];
            }
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
                                        SizeT window_center,
                                        SizeT data_size)
    {
        SizeT new_size = (right_win < min_window_size)
                            ? min_window_size
                            : right_win;

        SizeT max_right_win_size = data_size - window_center - 1;
        new_size = (new_size > max_right_win_size)
                    ? max_right_win_size
                    : new_size;

        return new_size;
    }

    std::pair<WinSizeT, WinSizeT> adaptive_half_windows(
                                    SizeT win_center, WinSizeT window,
                                    SizeT data_size,
                                    ValueT kz_d, ValueT max_d,
                                    ValueT derivative)
    {
        std::pair<WinSizeT, WinSizeT> adaptive_win;
        SizeT adaptive = this->adaptive_window(window, kz_d, max_d);

        // derivative = 0
        if (std::fabs(derivative) < tolerance) {
            adaptive_win.first = normalize_left_win(adaptive, win_center);
            adaptive_win.second = 
                normalize_right_win(adaptive, win_center, data_size);
        } else
        if (derivative < 0) {
            adaptive_win.first = normalize_left_win(adaptive, win_center);
            adaptive_win.second =
                normalize_right_win(window, win_center, data_size);
        } else {
            adaptive_win.second = 
                normalize_right_win(adaptive, win_center, data_size);
            adaptive_win.first = normalize_left_win(window, win_center);
        }

        return adaptive_win;
    }

    inline void calc_kz_difference(
            const std::vector<std::vector<ValueT>> &kz)
    {
        WinSizeT q1 = half_win_rows;
        WinSizeT q2 = half_win_cols;

        max_dx = 0;
        max_dy = 0;

        // calculate d1 = |Z(x+q2,y) - Z(x-q2,y)|
        // d2 = |Z(x,y+q1) - Z(x,y-q1)|
        for (SizeT i = 0; i < rows; ++i) { // y coord
            for (SizeT j = 0; j < cols; ++j) { // x coord
                SizeT row_tail = this->window_left_bound(i, q1);
                SizeT col_tail = this->window_left_bound(j, q2);
                SizeT row_head = this->window_right_bound(i, q1, rows);
                SizeT col_head = this->window_right_bound(j, q2, cols);

                kz_dx[i][j] = std::fabs(kz[i][col_head] - kz[i][col_tail]);
                kz_dy[i][j] = std::fabs(kz[row_head][j] - kz[row_tail][j]);

                max_dx = std::max(max_dx, kz_dx[i][j]);
                max_dy = std::max(max_dy, kz_dy[i][j]);
            }
        }
    }

    inline void calc_kz_dirivative()
    {
        /* d'(t) = d(i+1)-d(i) */

        for (SizeT i = 0; i < rows-1; ++i) {
            for (SizeT j = 0; j < cols-1; ++j) {
                x_derivative[i][j] = kz_dx[i][j+1] - kz_dx[i][j];
                y_derivative[i][j] = kz_dy[i+1][j] - kz_dy[i][j];
            }
        }
        for (SizeT i = 0; i < rows-1; ++i) {
            x_derivative[i][cols-1] = x_derivative[i][cols-2];
            y_derivative[i][cols-1] = kz_dy[i+1][cols-1] - kz_dy[i][cols-1];
        }
        for (SizeT j = 0; j < cols-1; ++j) {
            y_derivative[rows-1][j] = y_derivative[rows-2][j];
            x_derivative[rows-1][j] = kz_dx[rows-1][j+1] - kz_dx[rows-1][j];
        }
        /*
        for (SizeT i = 0; i < rows; ++i) {
            SizeT j;
            for (j = 0; j < cols-1; ++j) {
                x_derivative[i][j] = kz_dx[i][j+1] - kz_dx[i][j];
            }
            x_derivative[i][j] = x_derivative[i][j-1];
        }
        for (SizeT j = 0; j < cols; ++j) {
            for (SizeT i = 0; i < rows-1; ++i) {
                y_derivative[i][j] = kz_dy[i+1][j] - kz_dy[i][j];
            }
            y_derivative[rows-1][j] = y_derivative[rows-2][j];
        }
        */
    }
};

#endif
