/*

    Kolmogorov-Zurbenko Adaptive Filter
    Copyright (C) 2023, 2024 Vadim V. Marchenko <jakosvadim at gmail dot com>
    Copyright (C) 2005, 2015 Brian D. Close

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


/* -DKZA_TIMER compiler option for time measurement */
#ifdef KZA_TIMER
#include <chrono>
class Timer
{
public:
    Timer(const char *name = "time"): name(name) {}

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
        return std::chrono::duration_cast
                                <std::chrono::milliseconds>(diff).count();
    }

    inline void print_elapsed()
    {
        long elapsed = this->elapsed_time();
        std::cout << name << ": " << elapsed << "ms" << std::endl;
    }

    inline void set_name(const char *new_name)
    {
        name = new_name;
    }

private:
    std::chrono::steady_clock::time_point begin, end;
    const char *name;
};
#endif


/**
 * \brief Base template class of Kolmogorov-Zurbenko filter with general 
          parallel algorithm.
    * \param ValueT type of data values
    * \param SizeT type of data size
    * \param WinSizeT type of window size
*/
template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZGeneric
{
    typedef std::barrier<std::function<void()>> barrier_t;

public:
    /**
     * \brief Filter iterations setter. 
        \param iterations number of filter iterations
    */
    inline void set_iterations(SizeT iterations)
    {
        this->iterations = iterations;
    }

protected:
    /**
     * \param thread Thread identificator.
     * \param start_idx Starting data index of thread.
     * \param end_idx Ending data index of thread.
    */
    struct ThreadData {
        std::thread thread;
        SizeT start_idx, end_idx;
    };

    /** \brief Iterations number */
    SizeT iterations;
    /** \brief Data size to be processed */
    SizeT task_size;
    /** \brief Number of CPU cores available for multithreading */
    SizeT threads_cnt;
#ifdef KZA_TIMER
    Timer timer;
#endif


    /**
        \param iterations number of filter iterations
    */
    KZGeneric(SizeT iterations): iterations(iterations)
    {
        threads_cnt = std::thread::hardware_concurrency();
#ifdef KZA_DEBUG
        std::cout << "Number of cores: " << threads_cnt << std::endl;
#endif
    }

    virtual ~KZGeneric() = default;

    /**
     * \brief Single iteration step of the filter algorithm for each thred.
     * This method must be overridden by derived classes
         * \param start_idx Starting data index.
         * \param end_idx Ending (inclusive) data index.
    */
    virtual void perform_single_iteration(SizeT start_idx, SizeT end_idx) = 0;

    /**
     * \brief This method runs iterations of the filter algorithm. 
         At each iteration, parallel threads call 
         perform_single_iteration() on their piece of data. 
         The threads are synchronized using a barrier that calls 
         the on_iteration_complete() lambda function to synchronize the
         result of each thread.

         * \param on_iteration_complete Lambda function for threads results
                  synchronization.
    */
    void perform_iterations(std::function<void()> on_iteration_complete)
    {
        barrier_t sync_iteration(threads_cnt, on_iteration_complete);
        std::vector<ThreadData> th(threads_cnt);

#ifdef KZA_TIMER
    timer.start();
#endif
        this->start_threads(th, sync_iteration);

        for (auto &thread_data : th)
            thread_data.thread.join();

#ifdef KZA_TIMER
    timer.stop();
    timer.print_elapsed();
#endif
    }

    /**
     * \brief Calculates the left bound of a window with a given 
                center location. If the left bound of the window is 
                less than 0, then it is cut to 0. Whole window size 
                is 2*half_window+1.
         \param win_center Window center index.
         \param half_win Half-window size.
    */
    inline WinSizeT window_left_bound(SizeT win_center, WinSizeT half_win)
    {
        return (win_center >= half_win)
                ? (win_center - half_win)
                : 0;
    }

    /**
     * \brief Calculates the right bound of a window with a given 
                center location. If the right bound of the window 
                goes beyond the boundaries of the data array, then the 
                right border is trimmed to the maximum index.
                Whole window size is 2*half_window+1.

         \param win_center Window center index.
         \param half_win Half-window size.
         \param data_size Data size.
    */
    inline WinSizeT window_right_bound(SizeT win_center, WinSizeT half_win,
                                       SizeT data_size)
    {
        return (win_center + half_win >= data_size)
                ? data_size - 1
                : win_center + half_win;
    }

private:
    /**
     * \brief Thread main function.
         * \param t_data Thread data.
         * \param sync_iteration barrier to synchronize threads after each 
                iteration.
    */
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

    /**
     * \brief This function distributes the task between threads and 
            starts them.
         * \param th Threads array.
         * \param sync_iteration barrier to synchronize threads after each 
                iteration.
    */
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

            t_data.thread = std::thread(&KZGeneric::worker, this, 
                                        std::ref(t_data),
                                        std::ref(sync_iteration));
        }
    }
};


/**
 * \brief Kolmogorov-Zurbenko 1D filter.
*/
template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZ1D: public KZGeneric<ValueT, SizeT, WinSizeT>
{
public:
    /**
     * \param window_size Window size.
     * \param iterations Number of filter iterations. 
    */
    KZ1D(WinSizeT window_size, SizeT iterations) 
        : KZGeneric<ValueT, SizeT, WinSizeT>(iterations) 
    {
        half_window = window_size >> 1; // window half (window is 2*m+1)
#ifdef KZA_TIMER
        KZGeneric<ValueT, SizeT, WinSizeT>::timer.set_name("kz1d");
#endif
    }

    /**
     * \param vec Data to be processed by filter.
     * \return Result vector.
    */
    std::vector<ValueT> operator()(const std::vector<ValueT> &vec)
    {
        task_size = vec.size();

        data.assign(vec.cbegin(), vec.cend());
        ans.resize(vec.size());

#ifdef KZA_PREFIX_SUM
        pref_sum.resize(data.size() + 1);
        pref_finite_cnt.resize(data.size() + 1);
        update_prefix_sum(data); 
#endif

        if (iterations <= 0)
            return std::move(data);

        this->perform_iterations([this]() {
#ifdef KZA_PREFIX_SUM
            this->update_prefix_sum(ans);
#else
            std::swap(data, ans);
#endif
        });

        return std::move(ans);
    }

    /**
     * \brief Filter window size setter. Adaptive filter algorithm override
            this function since it uses different half-window sizes.
         \param window Window size.
    */
    inline virtual void set_window(WinSizeT window)
    {
        half_window = window >> 1;
    }

protected:

#ifdef KZA_PREFIX_SUM
    /** \brief -DKZA_PREFIX_SUM preprocessor option required. */
    std::vector<ValueT> pref_sum;
    std::vector<SizeT> pref_finite_cnt;

    /**
     * \brief -DKZA_PREFIX_SUM preprocessor option required.
            Calculate data prefix sums. The function also 
            calculates the number of finite elements on data array.
         \param data Initial data array.
    */
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

    /**
     * \brief -DKZA_PREFIX_SUM preprocessor option required.
            The method calculates the average value on the interval in
            data array usign prefix sum.
         \param start_idx Index of the first element.
         \param end_idx Index of the last element (inclusive).
    */
    ValueT average(SizeT start_idx, SizeT end_idx)
    {
        /* (window sum) = (sum containig window) - (sum before window) */
        SizeT z =
            (pref_finite_cnt[end_idx + 1] - pref_finite_cnt[start_idx]);

        return (z != 0) ? (pref_sum[end_idx + 1] - pref_sum[start_idx]) / z
                        : std::numeric_limits<ValueT>::quiet_NaN();
    }

#else

    /**
     * \brief The method calculates the average value on the interval in
            data array.
         \param start_idx Index of the first element.
         \param end_idx Index of the last element (inclusive).
    */
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

    inline void perform_single_iteration(SizeT start_idx, 
                                         SizeT end_idx) override
    {
        for (SizeT time = start_idx; time <= end_idx; ++time)
            ans[time] = this->average(
                this->window_left_bound(time, half_window),
                this->window_right_bound(time, half_window, data.size()));
    }

protected:
    /** \brief Half window size. Whole window size is 2*half_window+1. */
    WinSizeT half_window;
    /** \brief Data array. */
    std::vector<ValueT> data;
    /** \brief Algorithm result after iteration. */
    std::vector<ValueT> ans;

    using KZGeneric<ValueT, SizeT, WinSizeT>::task_size;
    using KZGeneric<ValueT, SizeT, WinSizeT>::iterations;
    using KZGeneric<ValueT, SizeT, WinSizeT>::threads_cnt;
};

/**
 * \brief Kolmogorov-Zurbenko 2D filter.
*/
template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZ2D: public KZGeneric<ValueT, SizeT, WinSizeT>
{
public:
    /**
     * \param window_sizes Window sizes in row and column spaces.
     * \param iterations Number of filter iterations. 
    */
    KZ2D(const std::vector<WinSizeT> window_sizes, SizeT iterations) 
        : KZGeneric<ValueT, SizeT, WinSizeT>(iterations)
    {
        this->set_window(window_sizes);
#ifdef KZA_TIMER
        KZGeneric<ValueT, SizeT, WinSizeT>::timer.set_name("kz2d");
#endif
    }

    /**
     * \param vec 2D vector of initial data to be processed by filter.
     * \return Result 2D vector.
    */
    std::vector<std::vector<ValueT>> operator()(
            const std::vector<std::vector<ValueT>> &vec)
    {
        task_size = vec.size();

        data = vec;
        rows = data.size();
        cols = data[0].size();
        ans.resize(rows, std::vector<ValueT>(cols));

#ifdef KZA_PREFIX_SUM
        pref_sum.resize(rows + 1, std::vector<ValueT>(cols + 1, 0));
        pref_finite_cnt.resize(rows + 1, std::vector<SizeT>(cols + 1, 0));
        update_prefix_sum(data); 
#endif

        if (iterations <= 0)
            return std::move(data);

        this->perform_iterations([this]() {
#ifdef KZA_PREFIX_SUM
            this->update_prefix_sum(ans);
#else
            data.swap(ans);
#endif
        });

        return std::move(ans);
    }

    inline virtual void set_window(const std::vector<WinSizeT> &window)
    {
        half_win_rows = window[0] >> 1;
        half_win_cols = window.size() > 1 ? window[1] >> 1 : half_win_rows;
    }

protected:

#ifdef KZA_PREFIX_SUM
    /** \brief -DKZA_PREFIX_SUM preprocessor option required. */
    std::vector<std::vector<ValueT>> pref_sum;
    std::vector<std::vector<SizeT>> pref_finite_cnt;

    /**
     * \brief -DKZA_PREFIX_SUM preprocessor option required.
            Calculate 2D data prefix sums. The function also 
            calculates the number of finite elements on data array.
         \param data 2D data array.
    */
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

    /**
     * \brief -DKZA_PREFIX_SUM preprocessor option required.
         The method calculates the average in a rectangle with given 
         corners coords using 2D prefix sum.
         \param start_row Upper corner row number.
         \param end_row Lower corner row number (inclusive).
         \param start_col Upper corner column number.
         \param end_col Lower corner column number (inclusive).
    */
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

    /**
     * \brief The method calculates the average in a rectangle with given 
             cornes coords.
         \param start_row Upper corner row number.
         \param end_row Lower corner row number (inclusive).
         \param start_col Upper corner column number.
         \param end_col Lower corner column number (inclusive).
    */
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

private:

    inline void perform_single_iteration(SizeT start_idx,
                                         SizeT end_idx) override
    {
        for (SizeT row = start_idx; row <= end_idx; ++row) {
            for (SizeT col = 0; col < cols; ++col) {
                ans[row][col] = this->average(
                    this->window_left_bound(row, half_win_rows),
                    this->window_right_bound(row, half_win_rows, rows),
                    this->window_left_bound(col, half_win_cols),
                    this->window_right_bound(col, half_win_cols, cols));
            }
        }
    }

protected:
    /** \brief Half window size in the row space. */
    WinSizeT half_win_rows;

    /** \brief Half window size in the column space. */
    WinSizeT half_win_cols;

    std::vector<std::vector<ValueT>> data, ans;
    SizeT rows, cols;

    using KZGeneric<ValueT, SizeT, WinSizeT>::task_size;
    using KZGeneric<ValueT, SizeT, WinSizeT>::iterations;
    using KZGeneric<ValueT, SizeT, WinSizeT>::threads_cnt;
};


/**
 * \brief Kolmogorov-Zurbenko filter adaptive extension.
*/
template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZAExtension
{
public:
    /**
     * \brief Adaptive filter minimal window setter. 
     * \param min_window The lower bound of the window size is 
        used during the calculation of the adaptive window size.
    */
    inline void set_min_window(WinSizeT min_window)
    {
        min_window_size = min_window; 
    }

    /**
     * \brief Adaptive filter tolerance setter. 
     * \param tol Accuracy when comparing the derivative to zero 
                        when calculating the adaptive window.
    */
    inline void set_tolerance(ValueT tol)
    {
        tolerance = tol; 
    }

    /**
     * \brief Calculates adaptive window size.
     * \param window Initial window size.
     * \param kz_d KZ filter result difference value.
     * \param max_d Maximum if KZ filter result difference values.
     * \returns Adaptive window size.
    */
    inline SizeT adaptive_window(WinSizeT window, ValueT kz_d, ValueT max_d)
    {
        return (std::fabs(max_d) > 0) 
               ? std::floor(window * (1 - kz_d/max_d))
               : window;
    }

protected:
    WinSizeT min_window_size;
    ValueT tolerance;

    /**
     * \param min_window_size The lower bound of the window size is 
        used during the calculation of the adaptive window size.
     * \param tolerance Accuracy when comparing the derivative to zero 
            when calculating the adaptive window.
    */
    KZAExtension(WinSizeT min_window_size, ValueT tolerance)
        : min_window_size(min_window_size), tolerance(tolerance) {}

    /**
     * \brief Normalize left window size. If the window is smaller than 
            the minimum window length, it extends it to this minimum. 
            If the window goes beyond the boundaries of the array, 
            it cuts it off.
     * \param left_win Initial left window size.
     * \param window_center Index of window center.
     * \return New left window size.
    */
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

    /**
     * \brief Normalize left window size. If the window is smaller than 
            the minimum window length, it extends it to this minimum. 
            If the window goes beyond the boundaries of the array, 
            it cuts it off.
     * \param right_win Initial right window size.
     * \param window_center Index of window center.
     * \param data_size Data array size.
     * \return New right window size.
    */
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
};


/**
 * \brief Kolmogorov-Zurbenko 1D adaptive filter.
*/
template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZA1D: public KZ1D<ValueT, SizeT, WinSizeT>,
             public KZAExtension<ValueT, SizeT, WinSizeT>
{
public:
    /**
     * \param window_size Window size.
     * \param min_window_size The lower bound of the window size is 
        used during the calculation of the adaptive window size.
     * \param tolerance Accuracy when comparing the derivative to zero 
            when calculating the adaptive window.
     * \param iterations Number of filter iterations. 
    */
    KZA1D(WinSizeT window_size, WinSizeT min_window_size,
          ValueT tolerance, SizeT iterations)
            : KZ1D<ValueT, SizeT, WinSizeT>(window_size, iterations),
              KZAExtension<ValueT, SizeT, WinSizeT>(min_window_size,
                                                    tolerance)
    {
        // unlike of KZ algorithm, KZA uses (2*window_size+1) window 
        half_window = window_size;           
#ifdef KZA_TIMER
        KZGeneric<ValueT, SizeT, WinSizeT>::timer.set_name("kza1d");
#endif
    }

    /**
     * \param data Initial data to be processed by filter.
     * \param kz_result The result of processing data with a 
                KZ1D filter.
     * \return Result vector.
    */
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

    inline void set_window(WinSizeT window) override
    {
        half_window = window;           
    }

private:
    /** \brief Array of KZ1D filter result differences |y(i+q) - y(i-q)| */
    std::vector<ValueT> kz_diff;
    ValueT max_diff;

    /** \brief Array of KZ1D filter result discrete derivatives */
    std::vector<ValueT> kz_derivative;

    using KZGeneric<ValueT, SizeT, WinSizeT>::iterations;

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

    /**
     * \brief Calculates the adaptive size of the window halves.
     * \param left_half_win A reference to the variable of the 
     *                      left half window.
     * \param right_half_win A reference to the variable of the 
     *                       right half window.
     * \param t Index of window center. Window size is 
                left_half_win + right_half_win + 1.
    */
    void calc_adaptive_half_windows(WinSizeT &left_half_win, 
                                    WinSizeT &right_half_win, 
                                    SizeT t)  
    {
        SizeT data_size = data.size();
        SizeT adaptive_size = this->adaptive_window(half_window, 
                                                    kz_diff[t], max_diff);

        // derivative[t] = 0
        if (std::fabs(kz_derivative[t]) < tolerance) {
            left_half_win = this->normalize_left_win(adaptive_size, t);
            right_half_win = 
                this->normalize_right_win(adaptive_size, t, data_size);
        } else if (kz_derivative[t] < 0) {
            left_half_win = this->normalize_left_win(adaptive_size, t);
            right_half_win = 
                this->normalize_right_win(half_window, t, data_size);
        } else {
            right_half_win = 
                this->normalize_right_win(adaptive_size, t, data_size);
            left_half_win = this->normalize_left_win(half_window, t);
        }
    }

    /** 
     * \brief Calculates differences of KZ1D filter result:
                 |y[i+q] - y[i-q]| when q is a half window size.
        \param y KZ1D filter result array.
    */
    void calc_kz_difference(const std::vector<ValueT> &y)
    {
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

    /** 
     * \brief Calculates discrete derivative of KZ1D filter result:
              d'[t] = d[t+1] - d[i].
    */
    inline void calc_kz_dirivative()
    {
        SizeT n = kz_diff.size();

        if (n < 2)
            return;

        for (SizeT i = 0; i < n-1; i++)
            kz_derivative[i] = kz_diff[i+1] - kz_diff[i];
        kz_derivative[n-1] = kz_derivative[n-2];
    }
};

/**
 * \brief Kolmogorov-Zurbenko 2D adaptive filter.
*/
template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZA2D: public KZ2D<ValueT, SizeT, WinSizeT>,
             public KZAExtension<ValueT, SizeT, WinSizeT>
{
public:
    /**
     * \param window Window size.
     * \param min_window_size The lower bound of the window size is 
        used during the calculation of the adaptive window size.
     * \param tolerance Accuracy when comparing the derivative to zero 
            when calculating the adaptive window.
     * \param iterations Number of filter iterations. 
    */
    KZA2D(const std::vector<WinSizeT> &window, WinSizeT min_window_size,
          ValueT tolerance, SizeT iterations)
            : KZ2D<ValueT, SizeT, WinSizeT>(window, iterations),
              KZAExtension<ValueT, SizeT, WinSizeT>(min_window_size,
                                                    tolerance)
    {
        // use KZA specific window size: (2*window+1) 
        this->set_window(window);
#ifdef KZA_TIMER
        KZGeneric<ValueT, SizeT, WinSizeT>::timer.set_name("kza2d");
#endif
    }

    /**
     * \param data Initial 2D data array to be processed by filter.
     * \param kz_result The result of processing data with a 
                        KZ2D filter.
     * \return Result 2D vector.
    */
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

    inline void set_window(const std::vector<WinSizeT> &window) override
    {
        // unlike of KZ algorithm, KZA uses (2*window+1) window 
        half_win_rows = window[0];           
        half_win_cols = window.size() > 1 ? window[1] : half_win_rows;           
    }

private:
    /** \brief Array of KZ2D filter result differences in column space */
    std::vector<std::vector<ValueT>> kz_dx;

    /** \brief Array of KZ2D filter result differences in row space */
    std::vector<std::vector<ValueT>> kz_dy;

    ValueT max_dx, max_dy;

    /** 
     \brief Array of KZ1D filter result discrete derivatives in column space 
    */
    std::vector<std::vector<ValueT>> x_derivative;

    /** 
     \brief Array of KZ1D filter result discrete derivatives in row space 
    */
    std::vector<std::vector<ValueT>> y_derivative;

    using KZGeneric<ValueT, SizeT, WinSizeT>::iterations;

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
                rows_win = this->adaptive_window_halves(
                    row, half_win_rows, rows,
                    kz_dy[row][col], max_dy, 
                    y_derivative[row][col]
                );
                cols_win = this->adaptive_window_halves(
                    col, half_win_cols, cols,
                    kz_dx[row][col], max_dx, 
                    x_derivative[row][col]
                );

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

    /**
     * \brief Calculates the adaptive size of the window halves in
     *          column or row space.
     * \param win_center Window center coord.
     * \param window Window size.
     * \param data_size Data size in row or column space.
     * \param kz_d KZ2D result difference in window_center coord.
     * \param max_d Maximum of KZ2D result differences.
     * \param derivative Discrete derivative if KZ2D result 
     *                  in window_center coord.
     * \returns Pair of adaptive window left and right halves.
    */
    std::pair<WinSizeT, WinSizeT> adaptive_window_halves(
                                    SizeT win_center, WinSizeT window,
                                    SizeT data_size,
                                    ValueT kz_d, ValueT max_d,
                                    ValueT derivative)
    {
        std::pair<WinSizeT, WinSizeT> adaptive_win;
        SizeT adaptive = this->adaptive_window(window, kz_d, max_d);

        // derivative = 0
        if (std::fabs(derivative) < tolerance) {
            adaptive_win.first = 
                this->normalize_left_win(adaptive, win_center);
            adaptive_win.second = 
                this->normalize_right_win(adaptive, win_center, data_size);
        } else
        if (derivative < 0) {
            adaptive_win.first = 
                this->normalize_left_win(adaptive, win_center);
            adaptive_win.second =
                this->normalize_right_win(window, win_center, data_size);
        } else {
            adaptive_win.second = 
                this->normalize_right_win(adaptive, win_center, data_size);
            adaptive_win.first = 
                this->normalize_left_win(window, win_center);
        }

        return adaptive_win;
    }

    /** 
     * \brief Calculates |Z(x+q2,y) - Z(x-q2,y)| and |Z(x,y+q1) - Z(x,y-q1)|
        \param kz 2D array of KZ2D filter result.
    */
    void calc_kz_difference(const std::vector<std::vector<ValueT>> &kz)
    {
        WinSizeT q1 = half_win_rows;
        WinSizeT q2 = half_win_cols;

        max_dx = 0;
        max_dy = 0;

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

    /** 
     * \brief Calculates discrete derivative of KZ2D filter result 
                in row and columns spaces  
    */
    void calc_kz_dirivative()
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
        x_derivative[rows-1][cols-1] = x_derivative[rows-1][cols-2];
        y_derivative[rows-1][cols-1] = y_derivative[rows-2][cols-1];
    }
};

#endif
