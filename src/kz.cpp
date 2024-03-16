#include "kza.h"

#include <iostream>
#include <cmath>
#include <limits>

#ifdef TIMER
#include <chrono>
#endif

#include <thread>
#include <barrier>

typedef std::barrier<std::__empty_completion> barrier_t;

struct KzData {
    int window;
    int task_size;
    int iterations;
    int data_size;
    double *ans;
    double *data;

#ifdef PREFIX_SUM
    double *pref_sum;
    int *pref_finite_cnt;
#endif

    KzData(int window, int task_size, int iterations,
           const double *data, int data_size) : 
                window(window), task_size(task_size), iterations(iterations),
                data_size(data_size)
    {
        this->data = new double[data_size];
        std::copy(data, data + data_size, this->data);
        this->ans = new double[data_size];
#ifdef PREFIX_SUM
        this->pref_sum = new double[data_size+1];
        this->pref_finite_cnt = new double[data_size+1];
#endif
    }

    ~KzData()
    {
        delete[] data;
        delete[] ans;
#ifdef PREFIX_SUM
        delete[] pref_sum;
        delete[] pref_finite_cnt;
#endif
    }
    
    double *extract_answer()
    {
        double *ans = this->ans;
        this->ans = nullptr;
        return ans;
    }
};

struct ThreadData {
    std::thread t;
    int start_idx, end_idx;
};

#ifdef PREFIX_SUM
static double mavg1d(const KzData &task, int window_center, int w)
{
    // window length is 2*w+1
    int start_idx = (window_center+1) - w - 1; // index of value before window
    start_idx *= (start_idx > 0);

    end_idx = (window_center+1) + w; // the last window value index
    if (end_idx > task.data.size())
        end_idx = task.data.size();

    /* (window sum) = (sum containig window) - (sum before window) */
    int s = (task.pref_sum[end_idx] - task.pref_sum[start_idx]);
    int z = (task.pref_finite_cnt[end_idx] - task.pref_finite_cnt[start_idx]);

    if (z == 0)
        return std::numeric_limits<double>::quiet_NaN();

    return s / (double)z;
}

static void calc_prefix_sum(const double *data, int data_size,
                            double *pref_sum, double *pref_finite_cnt)
{
    pref_sum[0] = 0;
    pref_finite_cnt[0] = 0;
    for (int i = 1; i <= data_size; ++i) {
        int is_finite_flag = std::isfinite(data[i-1]);
        pref_sum[i] = pref_sum[i-1] + is_finite_flag * data[i-1];
        pref_finite_cnt[i] = pref_finite_cnt[i-1] + is_finite_flag;
    }
}

#else
static double mavg1d(const double *x, int data_size, int window_center, int w)
{
    // window length is 2*w+1
    int start_idx = window_center - w;
    start_idx *= (start_idx > 0);

    int end_idx = window_center + w + 1;
    if (end_idx > data_size)
        end_idx = data_size;

    long z = 0;
    double s = 0;
    for (int i = start_idx; i < end_idx; ++i) {
        int is_finite_flag = std::isfinite(x[i]);
        z += is_finite_flag;
        s += is_finite_flag * x[i];
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
        task.ans[i] = mavg1d(task.data, task.data_size, i, task.window);
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

        // waiting for data update on server
        sync_iteration.arrive_and_wait();
    }
    // perform last iteration and finish
    perform_task_iteration(kz_data, thread.start_idx, thread.end_idx);
} 

static void start_threads(ThreadData *th, int threads_cnt, KzData &kz_data,
                          barrier_t &sync_iteration)
{
    int start_idx = 0;

    for (int i = 0; i < threads_cnt; ++i) {
        th[i].start_idx = start_idx;
        th[i].end_idx = start_idx + kz_data.task_size;

        if (th[i].end_idx > kz_data.data_size)
            th[i].end_idx = kz_data.data_size;
        else
            start_idx = th[i].end_idx;

        th[i].t = std::thread(worker, std::ref(th[i]), std::ref(kz_data),
                              std::ref(sync_iteration));
    }
}

static void wait_threads(ThreadData *th, int threads_cnt)
{
    for (int i = 0; i < threads_cnt; ++i)
        th[i].t.join();
}

static void threads_server_loop(int threads_cnt, KzData &kz_data)
{
    ThreadData *th = new ThreadData[threads_cnt];
    std::barrier sync_iteration(threads_cnt + 1);

    start_threads(th, threads_cnt, kz_data, sync_iteration);

    for (int iter = 0; iter < kz_data.iterations-1; ++iter) {
        // waiting for workers to complete iteration
        sync_iteration.arrive_and_wait();
#ifdef DEBUG
        std::cout << "Server updates data" << std::endl;
#endif

#ifdef PREFIX_SUM
        calc_prefix_sum(kz_data.ans, kz_data.data_size, 
                        kz_data.pref_sum, kz_data.pref_finite_cnt);
#else
        std::swap(kz_data.data, kz_data.ans);
#endif
        // allow workers to perfrom next iteration 
        sync_iteration.arrive_and_wait();
    } 

    // waiting for the workers to complete the last iteration
    wait_threads(th, threads_cnt);
    delete[] th;
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

    threads_server_loop(threads_cnt, kz_data);

    return kz_data.extract_answer();
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
