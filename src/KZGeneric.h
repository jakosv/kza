#ifndef KZ_H_SENTRY
#define KZ_H_SENTRY

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <vector>
#include <barrier>
#include <thread>
#include <functional>

typedef std::barrier<std::function<void()>> barrier_t;

template<class ValueType, class SizeType>
class KZGeneric {
public:
    KZGeneric(SizeType window_size, const ValueType *data,
              SizeType data_size);
    virtual ~KZGeneric() = default;

    void perform_iterations(SizeType iterations);
    ValueType *get_ans();

protected:
    struct ThreadData {
        std::thread thread;
        SizeType start_idx, end_idx;
        SizeType iterations;
    };

    SizeType window_size;
    SizeType task_size;
    SizeType threads_cnt;
    std::vector<ValueType> ans;
    std::vector<ValueType> data;

#ifdef PREFIX_SUM
    std::vector<ValueType> pref_sum;
    std::vector<SizeType> pref_finite_cnt;
#endif

    void worker(ThreadData &t_data, barrier_t &sync_iteration);
    void start_threads(std::vector<ThreadData> &th,
                       barrier_t &sync_iteration, SizeType iterations);
    ValueType average(SizeType start_idx, SizeType end_idx);
    void update_prefix_sum(const std::vector<ValueType> &data);
    virtual void perform_single_iteration(SizeType start_idx,
                                          SizeType end_idx) 
        {}

};


  #ifdef PREFIX_SUM

template<class ValueType, class SizeType>
KZGeneric<ValueType, SizeType>::KZGeneric(SizeType window_size, 
                                          const ValueType *data, 
                                          SizeType data_size) : 
                                                window_size(window_size), 
                                                ans(data_size), 
                                                data(data, data+data_size), 
                                                pref_sum(data_size+1), 
                                                pref_finite_cnt(data_size+1)
{
    threads_cnt = std::thread::hardware_concurrency();
    #ifdef DEBUG
    printf("Number of cores: %d\n", threads_cnt);
    #endif

    task_size = (data_size + threads_cnt-1) / threads_cnt;

    update_prefix_sum(this->data); 
}

  #else

template<class ValueType, class SizeType>
KZGeneric<ValueType, SizeType>::KZGeneric(SizeType window_size, 
                                const ValueType *data, SizeType data_size) : 
                                        window_size(window_size),
                                        ans(data_size), 
                                        data(data, data+data_size)
{
    threads_cnt = std::thread::hardware_concurrency();
    #ifdef DEBUG
    printf("Number of cores: %d\n", threads_cnt);
    #endif

    task_size = (data_size + threads_cnt-1) / threads_cnt;
}
  #endif

template<class ValueType, class SizeType>
void KZGeneric<ValueType, SizeType>::perform_iterations(SizeType iterations)
{
    auto on_iteration_complete = [this]() {
  #ifdef PREFIX_SUM
        this->update_prefix_sum(ans);
  #else
        std::swap(data, ans);
  #endif
    };

    barrier_t sync_iteration(threads_cnt, on_iteration_complete);
    std::vector<ThreadData> th(threads_cnt);

    this->start_threads(th, sync_iteration, iterations);

    for (auto &thread_data : th)
        thread_data.thread.join();
}

template<class ValueType, class SizeType>
ValueType *KZGeneric<ValueType, SizeType>::get_ans()
{
    ValueType *answer;
    try {
        answer = new ValueType[data.size()];
    } 
    catch (const std::bad_alloc&) {
        return nullptr; 
    }

    std::copy(ans.begin(), ans.end(), answer);

    return answer;
}

template<class ValueType, class SizeType>
void KZGeneric<ValueType, SizeType>::worker(KZGeneric::ThreadData &t_data, 
                                            barrier_t &sync_iteration)
{
    for (SizeType i = 0; i < t_data.iterations-1; ++i) {
        this->perform_single_iteration(t_data.start_idx, t_data.end_idx);

        // waiting for other threads to complete iteration
        sync_iteration.arrive_and_wait();
    }
    // perform last iteration and finish
    this->perform_single_iteration(t_data.start_idx, t_data.end_idx);
} 

template<class ValueType, class SizeType>
void KZGeneric<ValueType, SizeType>::start_threads(
                              std::vector<KZGeneric::ThreadData> &th, 
                              barrier_t &sync_iteration, 
                              SizeType iterations)
{
    SizeType start_idx = 0;

    for (auto &t_data : th) {
        t_data.iterations = iterations;
        t_data.start_idx = start_idx;
        t_data.end_idx = (start_idx + task_size >= data.size()) ? 
                                                    data.size() - 1 : 
                                                    start_idx + task_size;
        start_idx = t_data.end_idx;

        t_data.thread = std::thread(&KZGeneric::worker, this, 
                                    std::ref(t_data),
                                    std::ref(sync_iteration));
    }
}


  #ifdef PREFIX_SUM

template<class ValueType, class SizeType>
ValueType KZGeneric<ValueType, SizeType>::average(SizeType start_idx, 
                                                  SizeType end_idx)
{
    /* (window sum) = (sum containig window) - (sum before window) */
    SizeType z = (pref_finite_cnt[end_idx+1] - pref_finite_cnt[start_idx]);

    if (z == 0)
        return std::numeric_limits<ValueType>::quiet_NaN();

    return (pref_sum[end_idx+1] - pref_sum[start_idx]) / z;
}

template<class ValueType, class SizeType>
void KZGeneric<ValueType, SizeType>::update_prefix_sum(
                                        const std::vector<ValueType> &data)
{
    pref_sum[0] = 0;
    pref_finite_cnt[0] = 0;

    for (SizeType i = 1; i <= data.size(); ++i) {
        SizeType is_finite_flag = std::isfinite(data[i-1]);
        pref_sum[i] = pref_sum[i-1] + is_finite_flag * data[i-1];
        pref_finite_cnt[i] = pref_finite_cnt[i-1] + is_finite_flag;
    }
}

  #else

template<class ValueType, class SizeType>
ValueType KZGeneric<ValueType, SizeType>::average(SizeType start_idx, 
                                                  SizeType end_idx)
{
    SizeType z = 0;
    ValueType s = 0;

    for (SizeType i = start_idx; i <= end_idx; ++i) {
        SizeType is_finite_flag = std::isfinite(data[i]);
        z += is_finite_flag;
        s += is_finite_flag * data[i];
    }

    if (z == 0)
        return std::numeric_limits<ValueType>::quiet_NaN();

    return s / z;
}

  #endif

#endif
