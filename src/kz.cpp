#include "kza.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <functional>

#ifdef TIMER
#include <chrono>
#endif

#include <thread>
#include <barrier>

typedef std::barrier<std::function<void()>> barrier_t;

struct ThreadData {
    std::thread t;
    long start_idx, end_idx;
};

class KZ {
public:

#ifdef PREFIX_SUM
    KZ(int window, int task_size, int iterations, const double *data, 
           int data_size) : 
                window(window), task_size(task_size), iterations(iterations),
                ans(data_size), data(data, data+data_size), 
                pref_sum(data_size+1), pref_finite_cnt(data_size+1)
    {
        update_prefix_sum(this->data); 
    }
#else
    KZ(int window, int task_size, int iterations, const double *data, 
           int data_size) : 
                window(window), task_size(task_size), iterations(iterations),
                ans(data_size), data(data, data+data_size)
    {}
#endif

    void perform_iterations(int threads_cnt)
    {
        auto on_iteration_complete = [this]() {
#ifdef PREFIX_SUM
            update_prefix_sum(ans);
#else
            std::swap(data, ans);
#endif
        };

        barrier_t sync_iteration(threads_cnt, on_iteration_complete);
        std::vector<ThreadData> th(threads_cnt);
        start_threads(th, sync_iteration);

        for (auto &thread : th)
            thread.t.join();
    }

    double *get_ans()
    {
        double *answer;
        try {
            answer = new double[data.size()];
        } 
        catch (const std::bad_alloc&) {
            return nullptr; 
        }

        std::copy(ans.begin(), ans.end(), answer);

        return answer;
    }

private:
    int window;
    int task_size;
    int iterations;
    std::vector<double> ans;
    std::vector<double> data;

#ifdef PREFIX_SUM
    std::vector<double> pref_sum;
    std::vector<int> pref_finite_cnt;
#endif

    void worker(ThreadData &thread, barrier_t &sync_iteration)
    {
        for (int k = 0; k < iterations-1; ++k) {
            perform_single_iteration(thread.start_idx, thread.end_idx);

            // waiting for other threads to complete iteration
            sync_iteration.arrive_and_wait();
        }
        // perform last iteration and finish
        perform_single_iteration(thread.start_idx, thread.end_idx);
    } 

    void start_threads(std::vector<ThreadData> &th, barrier_t &sync_iteration)
    {
        unsigned long start_idx = 0;

        for (auto &thread : th) {
            thread.start_idx = start_idx;
            thread.end_idx = (start_idx + task_size >= data.size()) ? 
                                                        data.size() : 
                                                        start_idx + task_size;
            start_idx = thread.end_idx;

            thread.t = std::thread(&KZ::worker, this, std::ref(thread),
                                   std::ref(sync_iteration));
        }
    }


    inline void calc_window_bounds(long *start_idx, long *end_idx, 
                                   int win_center)
    {
        // window length is 2*win_size+1
        *start_idx = win_center - window;
        *start_idx *= (*start_idx > 0);

        *end_idx = win_center + window + 1;
        *end_idx = (*end_idx >= data.size()) ? data.size() : *end_idx;
    }

#ifdef PREFIX_SUM
    double mavg1d(int win_center)
    {
        long start_idx, end_idx;
        calc_window_bounds(&start_idx, &end_idx, win_center);

        /* (window sum) = (sum containig window) - (sum before window) */
        int z = (pref_finite_cnt[end_idx] - pref_finite_cnt[start_idx]);

        if (z == 0)
            return std::numeric_limits<double>::quiet_NaN();

        return (pref_sum[end_idx] - pref_sum[start_idx]) / z;
    }

    void update_prefix_sum(const std::vector<double> &data)
    {
        pref_sum[0] = 0;
        pref_finite_cnt[0] = 0;

        for (int i = 1; i <= data.size(); ++i) {
            int is_finite_flag = std::isfinite(data[i-1]);
            pref_sum[i] = pref_sum[i-1] + is_finite_flag * data[i-1];
            pref_finite_cnt[i] = pref_finite_cnt[i-1] + is_finite_flag;
        }
    }

#else
    double mavg1d(int win_center)
    {
        long start_idx, end_idx;
        calc_window_bounds(&start_idx, &end_idx, win_center);

        long z = 0;
        double s = 0;
        for (int i = start_idx; i < end_idx; ++i) {
            int is_finite_flag = std::isfinite(data[i]);
            z += is_finite_flag;
            s += is_finite_flag * data[i];
        }

        if (z == 0)
            return std::numeric_limits<double>::quiet_NaN();

        return s / (double)z;
    }
#endif

    void perform_single_iteration(int start_idx, int end_idx)
    {
        for (int i = start_idx; i < end_idx; ++i) {
#ifdef PREFIX_SUM
            ans[i] = mavg1d(i);
#else
            ans[i] = mavg1d(i);
#endif
        }
    }

};

double *kz1d(const double *data, int data_size, int window, int iterations)
{
    int threads_cnt = std::thread::hardware_concurrency();
#ifdef DEBUG
    printf("Number of cores: %d\n", threads_cnt);
#endif

    int task_size = (data_size + threads_cnt-1) / threads_cnt;
    KZ kz_data(window, task_size, iterations, data, data_size);

    kz_data.perform_iterations(threads_cnt);

    return kz_data.get_ans();
}

double *kz(const double *data, int dimention, const int *data_size,
           const int *window, int iterations)
{
    double *ans = NULL;

#ifdef TIMER
    std::chrono::steady_clock::time_point begin, end;
    long int diff;
#endif

    switch (dimention) {
    case 1:
#ifdef TIMER
        begin = std::chrono::steady_clock::now();
#endif
        ans = kz1d(data, data_size[0], 0.5 * window[0], iterations); 
#ifdef TIMER
        end = std::chrono::steady_clock::now();
        diff = 
            std::chrono::duration_cast<std::chrono::milliseconds>(end-begin).count();
        std::cout << "kz(): " << diff << "ms" << std::endl;
#endif
        break;
    default:
        fprintf(stderr, "kza: Too many dimensions\n");
        break;
    }

    return ans;
}

void kz_free(double *data)
{
    delete[] data;
}
