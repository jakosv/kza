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

struct KzData {
    int window;
    int task_size;
    int iterations;
    std::vector<double> ans;
    std::vector<double> data;

#ifdef PREFIX_SUM
    std::vector<double> pref_sum;
    std::vector<int> pref_finite_cnt;
#endif

#ifdef PREFIX_SUM
    KzData(int window, int task_size, int iterations, const double *data, 
           int data_size) : 
                window(window), task_size(task_size), iterations(iterations),
                ans(data_size), data(data, data+data_size), 
                pref_sum(data_size+1), pref_finite_cnt(data_size+1)
    {}
#else
    KzData(int window, int task_size, int iterations, const double *data, 
           int data_size) : 
                window(window), task_size(task_size), iterations(iterations),
                ans(data_size), data(data, data+data_size)
    {}
#endif
};

struct ThreadData {
    std::thread t;
    long start_idx, end_idx;
};

static inline void calc_window_bounds(long *start_idx, long *end_idx, 
                                      const std::vector<double> &data, 
                                      int win_center, int win_size)
{
    // window length is 2*win_size+1
    *start_idx = win_center - win_size;
    *start_idx *= (*start_idx > 0);

    *end_idx = win_center + win_size + 1;
    if (*end_idx > data.size())
        *end_idx = data.size();
}

#ifdef PREFIX_SUM
static double mavg1d(const KzData &kz_data, int win_center, int win_size)
{
    long start_idx, end_idx;
    calc_window_bounds(&start_idx, &end_idx, data, win_center, win_size);

    /* (window sum) = (sum containig window) - (sum before window) */
    int s = (task.pref_sum[end_idx] - task.pref_sum[start_idx]);
    int z = (task.pref_finite_cnt[end_idx] - task.pref_finite_cnt[start_idx]);

    if (z == 0)
        return std::numeric_limits<double>::quiet_NaN();

    return s / (double)z;
}

static void calc_prefix_sum(const std::vector<double> &data, 
                            std::vector<double> &pref_sum,
                            std::vector<double> &pref_finite_cnt)
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
static double mavg1d(const std::vector<double> &data, int win_center, 
                     int win_size)
{
    long start_idx, end_idx;
    calc_window_bounds(&start_idx, &end_idx, data, win_center, win_size);

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

static void perform_task_iteration(KzData &task, int start_idx, int end_idx)
{
    for (int i = start_idx; i < end_idx; ++i) {
#ifdef PREFIX_SUM
        task.ans[i] = mavg1d(task, i, task.window);
#else
        task.ans[i] = mavg1d(task.data, i, task.window);
#endif
    }
}

static void worker(ThreadData &thread, KzData &kz_data,
                   barrier_t &sync_iteration)
{
    for (int k = 0; k < kz_data.iterations-1; ++k) {
        perform_task_iteration(kz_data, thread.start_idx, thread.end_idx);

        // waiting for other threads to complete iteration
        sync_iteration.arrive_and_wait();
    }
    // perform last iteration and finish
    perform_task_iteration(kz_data, thread.start_idx, thread.end_idx);
} 

static void start_threads(std::vector<ThreadData> &th, KzData &kz_data,
                          barrier_t &sync_iteration)
{
    unsigned long start_idx = 0;

    for (auto &thread : th) {
        thread.start_idx = start_idx;
        thread.end_idx = start_idx + kz_data.task_size;

        if (thread.end_idx > kz_data.data.size())
            thread.end_idx = kz_data.data.size();
        else
            start_idx = thread.end_idx;

        thread.t = std::thread(worker, std::ref(thread), std::ref(kz_data),
                               std::ref(sync_iteration));
    }
}

static void perform_iterations(int threads_cnt, KzData &kz_data)
{
    auto on_iteration_complete = [&kz_data]() {
#ifdef PREFIX_SUM
        calc_prefix_sum(kz_data.ans, kz_data.pref_sum, 
                        kz_data.pref_finite_cnt);
#else
        std::swap(kz_data.data, kz_data.ans);
#endif
    };

    barrier_t sync_iteration(threads_cnt, on_iteration_complete);
    std::vector<ThreadData> th(threads_cnt);
    start_threads(th, kz_data, sync_iteration);

    for (auto &thread : th)
        thread.t.join();
}

double *kz1d(const double *data, int data_size, int window, int iterations)
{
    int threads_cnt = std::thread::hardware_concurrency();
#ifdef DEBUG
    printf("Number of cores: %d\n", threads_cnt);
#endif

    int task_size = (data_size + threads_cnt-1) / threads_cnt;
    KzData kz_data(window, task_size, iterations, data, data_size);
#ifdef PREFIX_SUM
    calc_prefix_sum(kz_data.data, kz_data.pref_sum, kz_data.pref_finite_cnt);
#endif

    perform_iterations(threads_cnt, kz_data);

    double *ans = new double[data_size];
    std::copy(kz_data.ans.begin(), kz_data.ans.end(), ans);

    return ans;
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
        ans = kz1d(data, data_size[0], window[0] / 2.0, iterations); 
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
