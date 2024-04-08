#ifndef KZ_GENERIC_H_SENTRY
#define KZ_GENERIC_H_SENTRY

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <vector>
#include <barrier>
#include <thread>
#include <functional>
#include <concepts>

typedef std::barrier<std::function<void()>> barrier_t;

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
class KZGeneric {
public:
    KZGeneric(WinSizeT window_size, const ValueT *data, SizeT data_size);
    virtual ~KZGeneric() = default;

    void perform_iterations(SizeT iterations);
    ValueT *get_ans();

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

    void update_prefix_sum(const std::vector<ValueT> &data);
  #endif

    void worker(ThreadData &t_data,
                SizeT iterations,
                barrier_t &sync_iteration);

    void start_threads(std::vector<ThreadData> &th,
                       SizeT iterations,
                       barrier_t &sync_iteration);

    ValueT average(SizeT start_idx, SizeT end_idx);

    virtual void perform_single_iteration(SizeT start_idx, SizeT end_idx) {}
};

  #ifdef PREFIX_SUM

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
KZGeneric<ValueT, SizeT, WinSizeT>::KZGeneric(
        WinSizeT window_size, const ValueT *data, SizeT data_size)
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

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
KZGeneric<ValueT, SizeT, WinSizeT>::KZGeneric(
        WinSizeT window_size, const ValueT *data, SizeT data_size)
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


template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
void KZGeneric<ValueT, SizeT, WinSizeT>::perform_iterations(SizeT iterations)
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

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
ValueT *KZGeneric<ValueT, SizeT, WinSizeT>::get_ans()
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

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
void KZGeneric<ValueT, SizeT, WinSizeT>::worker(
                                            KZGeneric::ThreadData &t_data, 
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

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
void KZGeneric<ValueT, SizeT, WinSizeT>::start_threads(
                              std::vector<KZGeneric::ThreadData> &th,
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

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
ValueT KZGeneric<ValueT, SizeT, WinSizeT>::average(SizeT start_idx, 
                                                   SizeT end_idx)
{
    /* (window sum) = (sum containig window) - (sum before window) */
    SizeT z = (pref_finite_cnt[end_idx + 1] - pref_finite_cnt[start_idx]);

    if (z == 0)
        return std::numeric_limits<ValueT>::quiet_NaN();

    return (pref_sum[end_idx+1] - pref_sum[start_idx]) / z;
}

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
void KZGeneric<ValueT, SizeT, WinSizeT>::update_prefix_sum(
                                        const std::vector<ValueT> &data)
{
    pref_sum[0] = 0;
    pref_finite_cnt[0] = 0;

    for (SizeT i = 1; size_t(i) <= data.size(); ++i) {
        bool is_finite_flag = std::isfinite(data[i-1]);
        pref_sum[i] = pref_sum[i-1] + is_finite_flag * data[i-1];
        pref_finite_cnt[i] = pref_finite_cnt[i-1] + is_finite_flag;
    }
}

  #else

template<class ValueT, std::unsigned_integral SizeT,
         std::unsigned_integral WinSizeT>
ValueT KZGeneric<ValueT, SizeT, WinSizeT>::average(SizeT start_idx, 
                                                   SizeT end_idx)
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

#endif
