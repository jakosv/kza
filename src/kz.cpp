#include "kza.h"

#include <vector>
#include <cstring>
#include <iostream>
#include <cmath>
#include <limits>
#include <algorithm>

#ifdef TIMER
#include "timer.h"
#endif

#include <thread>
#include <barrier>

struct TaskData {
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
    TaskData(int window, int task_size, int iterations,
             const double *data, int data_size) : 
                window(window), task_size(task_size), iterations(iterations),
                data(data, data + data_size), ans(data_size),
                pref_sum(data_size+1), pref_finite_cnt(data_size+1)
    {} 
#else
    TaskData(int window, int task_size, int iterations,
             const double *data, int data_size) : 
                window(window), task_size(task_size), iterations(iterations),
                data(data, data + data_size), ans(data_size)
    {}
#endif
};

struct ThreadData {
    std::thread t;
    int id;
    /*
    TaskData &task;
    */
    int start_idx, end_idx;

    /*
    ThreadData(TaskData &task) : task(task)
    {}
    */
};

#ifdef PREFIX_SUM
static double mavg1d(const TaskData &task, int window_center, int w)
{
    double s;
    int z;
    int start_idx, end_idx;

    /* window length is 2*w+1 */
    start_idx = (window_center+1) - w - 1; /* index of value before window */
    start_idx *= (start_idx > 0);

    end_idx = (window_center+1) + w; /* the last window value index */
    if (end_idx > task.data.size())
        end_idx = task.data.size();

    /* (window sum) = (sum containig window) - (sum before window) */
    s = (task.pref_sum[end_idx] - task.pref_sum[start_idx]);
    z = (task.pref_finite_cnt[end_idx] - task.pref_finite_cnt[start_idx]);

    if (z == 0)
        return std::numeric_limits<double>::quiet_NaN();

    return s / (double)z;
}

static void calc_prefix_sum(const std::vector<double> &data,
                            std::vector<double> &pref_sum, 
                            std::vector<int> &pref_finite_cnt)
{
    int i;

    pref_sum[0] = 0;
    pref_finite_cnt[0] = 0;
    for (i = 1; i <= data.size(); i++) {
        int is_finite_flag = std::isfinite(data[i-1]);
        pref_sum[i] = pref_sum[i-1] + is_finite_flag * data[i-1];
        pref_finite_cnt[i] = pref_finite_cnt[i-1] + is_finite_flag;
    }
}

#else
static double mavg1d(const std::vector<double> &x, int window_center, int w)
{
    double s = 0;
    long z;
    int i, start_idx, end_idx;

    /* window length is 2*w+1 */
    start_idx = window_center - w;
    start_idx *= (start_idx > 0);

    end_idx = window_center + w + 1;
    if (end_idx > x.size())
        end_idx = x.size();

    z = 0;
    for (i = start_idx; i < end_idx; i++) {
        int is_finite_flag = std::isfinite(x[i]);
        z += is_finite_flag;
        s += is_finite_flag * x[i];
    }

    if (z == 0)
        return std::numeric_limits<double>::quiet_NaN();

    return s / (double)z;
}
#endif

static void perform_task_iteration(TaskData &task, int start_idx, int end_idx)
{
    for (int i = start_idx; i < end_idx; i++) {
#ifdef PREFIX_SUM
        task.ans[i] = mavg1d(task, i, task.window);
#else
        task.ans[i] = mavg1d(task.data, i, task.window);
#endif
    }
}

static void worker(std::barrier<std::__empty_completion> &sync_iteration,
                   ThreadData &thread, TaskData &task)
{
    for (int k = 0; k < task.iterations-1; k++) {
        perform_task_iteration(task, thread.start_idx, thread.end_idx);

        // waiting for other threads to complete iteration
#ifdef DEBUG
        std::cout << "Worker " << thread.id << " iter " << k << std::endl;
#endif
        sync_iteration.arrive_and_wait();

        // waiting for data update on server
        sync_iteration.arrive_and_wait();
    }
    // perform last iteration and finish
    perform_task_iteration(task, thread.start_idx, thread.end_idx);
#ifdef DEBUG
    std::cout << "Worker " << thread.id << " last iter "<< std::endl;
#endif
} 

static void start_threads(std::vector<ThreadData> &th, TaskData &task,
                          std::barrier<std::__empty_completion> &sync_iteration)
{
    int start_idx = 0;

    for (int i = 0; i < th.size(); i++) {
        th[i].id = i;
        th[i].start_idx = start_idx;
        th[i].end_idx = start_idx + task.task_size;

        if (th[i].end_idx > task.data.size())
            th[i].end_idx = task.data.size();
        else
            start_idx = th[i].end_idx;

        th[i].t = std::thread(worker, std::ref(sync_iteration), 
                              std::ref(th[i]), std::ref(task));
    }
}

static void wait_threads(std::vector<ThreadData> &th)
{
    for (auto &thread : th) {
#ifdef DEBUG
        std::cout << thread.id << " was joined\n";
#endif
        thread.t.join();
    }
}

static void threads_server_loop(int threads_cnt, TaskData &task)
{
    std::vector<ThreadData> th(threads_cnt);
    std::barrier sync_iteration(th.size() + 1);

    start_threads(th, task, sync_iteration);

    for (int iter = 0; iter < task.iterations-1; iter++) {
        // waiting for workers to complete iteration
        sync_iteration.arrive_and_wait();
#ifdef DEBUG
        std::cout << "Server updates data" << std::endl;
#endif

#ifdef PREFIX_SUM
        calc_prefix_sum(task.ans, task.pref_sum, task.pref_finite_cnt);
#else
        std::swap(task.data, task.ans);
#endif
        // allow workers to perfrom next iteration 
        sync_iteration.arrive_and_wait();
    } 

    // waiting for the workers to complete the last iteration
    wait_threads(th);
}

double *kz1d(const double *data, int data_size, int window, int iterations)
{
    int threads_cnt = std::thread::hardware_concurrency();
#ifdef DEBUG
    printf("Number of cores: %d\n", threads_cnt);
#endif

    int task_size = (data_size + threads_cnt-1) / threads_cnt;
    TaskData task(window, task_size, iterations, data, data_size);
#ifdef PREFIX_SUM
    calc_prefix_sum(task.data, task.pref_sum, task.pref_finite_cnt);
#endif

    threads_server_loop(threads_cnt, task);

    double *ans = (double *)malloc(sizeof(double) * data_size);
    if (!ans)
        return NULL;

    std::copy(task.ans.begin(), task.ans.end(), ans);

    return ans;
}

double *kz(const double *data, int dimention, const int *data_size,
           const int *window, int iterations)
{
    double *ans = NULL;

#ifdef TIMER
    struct timeval t_start;
#endif

    switch (dimention) {
    case 1:
#ifdef TIMER
        timestamp(&t_start);
#endif
        ans = kz1d(data, data_size[0], window[0] / 2.0, iterations); 
#ifdef TIMER
        print_timer(&t_start, "kz():");
#endif

        break;
    default:
        fprintf(stderr, "kza: Too many dimensions\n");
        break;
    }

    return ans;
}
