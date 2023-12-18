#include "kza.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef TIMER
#include <sys/time.h>
#endif

#ifdef KZ_PARALLEL
#include <pthread.h>

#ifdef _WIN32
#include <windows.h>
#elif MACOS
#include <sys/param.h>
#include <sys/sysctl.h>
#else
#include <unistd.h>
#endif

int get_num_cores() {
#ifdef WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;
#elif MACOS
    int nm[2];
    size_t len = 4;
    uint32_t count;

    nm[0] = CTL_HW; nm[1] = HW_AVAILCPU;
    sysctl(nm, 2, &count, &len, NULL, 0);

    if (count < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &count, &len, NULL, 0);
        if (count < 1) {
            count = 1;
        }
    }
    return count;
#else
    return sysconf(_SC_NPROCESSORS_ONLN);
#endif
}

struct task_data {
    int start_idx, end_idx;
    int window;
    int task_size;
    int data_size;
    double *ans;
    double *pref_sum;
    int *pref_finite_cnt;
};

#endif

static double mavg1d(const double *pref_sum, const int *pref_finite_cnt, 
                     int length, int col, int w)
{
    double s;
    int z;
    int start_col, end_col;

    if (!pref_sum) {
        fprintf(stderr, "mavg1d(): Incorrect params\n");
        return 0;
    }

    start_col = col - w;
    if (start_col < 0)
        start_col = 0;

    end_col = col + w + 1;
    if (end_col > length)
        end_col = length;

    s = (pref_sum[end_col] - pref_sum[start_col]);
    z = (pref_finite_cnt[end_col] - pref_finite_cnt[start_col]);
    /*
    if (z == 0) 
        return nan("");
    */
    return s/z;
}

#ifdef KZ_PARALLEL

static void *thread_task(void *data)
{
    int i;
    struct task_data *args = data;

    for (i = args->start_idx; i < args->end_idx; i++) {
        args->ans[i] = mavg1d(args->pref_sum, args->pref_finite_cnt,
                              args->data_size, i, args->window);
    }

    return NULL;
} 

static void init_tasks(struct task_data **tasks, int tasks_cnt, int task_size,
                       int window, int data_size, double *ans, 
                       double *pref_sum, int *pref_finite_cnt)
{
    int i, start_idx;
    start_idx = 0;
    for (i = 0; i < tasks_cnt; i++) {
        struct task_data *task;

        task = malloc(sizeof(struct task_data));
        task->start_idx = start_idx;
        task->end_idx = start_idx + task_size;
        task->window = window;
        task->data_size = data_size;
        task->ans = ans;
        task->pref_sum = pref_sum;
        task->pref_finite_cnt = pref_finite_cnt;

        tasks[i] = task;

        start_idx = task->end_idx;
    }

}

static void start_threads(pthread_t *th, int threads_cnt,
                          struct task_data **tasks)
{
    int i;
    for (i = 0; i < threads_cnt; i++) {
        int res;

        res = pthread_create(&th[i], NULL, thread_task, tasks[i]);
        if (res != 0) {
            perror("thread_create");
            exit(1);
        }
    }
}

static void wait_threads(const pthread_t *th, int threads_cnt)
{
    int i;
    for (i = 0; i < threads_cnt; i++) {
        pthread_join(th[i], NULL);
    }
}

static void free_tasks(struct task_data **tasks, int tasks_cnt)
{
    int i;
    for (i = 0; i < tasks_cnt; i++) {
        free(tasks[i]);
    }
}

#endif

static void calc_prefix_sum(const double *data, int size, double *pref_sum, 
                            int *pref_finite_cnt)
{
    int i;

    pref_sum[0] = 0;
    pref_finite_cnt[0] = 0;
    for (i = 1; i <= size; i++) {
        pref_sum[i] = pref_sum[i-1];
        pref_finite_cnt[i] = pref_finite_cnt[i-1];
        if (isfinite(data[i-1])) {
            pref_sum[i] += data[i-1];
            ++pref_finite_cnt[i];
        }
    }
}

static double *kz1d(const double *x, int length, int window, int iterations)
{
    int i, k;
    int mem_size;
    double *data = NULL, *ans = NULL, *pref_sum = NULL;
    int *pref_finite_cnt = NULL;

#ifdef KZ_PARALLEL
    pthread_t *th;
    struct task_data **tasks;
    int thread_task_size, max_threads_cnt, tasks_cnt;

    max_threads_cnt = get_num_cores();

#ifdef DEBUG
    printf("Number of cores: %d\n", max_threads_cnt);
#endif

    tasks_cnt = max_threads_cnt - 1;
    th = malloc(tasks_cnt * sizeof(pthread_t));
    tasks = malloc(tasks_cnt * sizeof(struct task_data*));

#endif

    mem_size = length * sizeof(double);
    ans = malloc(mem_size);
    data = malloc(mem_size);
    pref_sum = malloc((length+1) * sizeof(double));
    pref_finite_cnt = malloc((length+1) * sizeof(int));
    if (!ans || !data || !pref_sum || !pref_finite_cnt)
        goto quit;

    memcpy(data, x, mem_size);
    calc_prefix_sum(data, length, pref_sum, pref_finite_cnt);

#ifdef KZ_PARALLEL
    thread_task_size = length / tasks_cnt;
    
    init_tasks(tasks, tasks_cnt, thread_task_size, window, length,
               ans, pref_sum, pref_finite_cnt);
#endif

    for (k = 0; k < iterations; k++) {
#ifdef KZ_PARALLEL
        start_threads(th, tasks_cnt, tasks);

        for (i = (tasks_cnt)*thread_task_size; i < length; i++) {
            ans[i] = mavg1d(pref_sum, pref_finite_cnt, length, i, window);
        }
#else
        for (i = 0; i < length; i++) {
            ans[i] = mavg1d(pref_sum, pref_finite_cnt, length, i, window);
        }
#endif

#ifdef KZ_PARALLEL
        wait_threads(th, tasks_cnt);
#endif

        memcpy(data, ans, mem_size); 
        calc_prefix_sum(data, length, pref_sum, pref_finite_cnt);
    }

#ifdef KZ_PARALLEL
    free_tasks(tasks, tasks_cnt);
#endif

quit:
    free(data);
    free(pref_sum);
    free(pref_finite_cnt);

#ifdef KZ_PARALLEL
    free(th);
    free(tasks);
#endif

    return ans;
}

double *kz(const double *x, int dim, const int *size, const int *window,
           int iterations)
{
    int i;
    int *win_size;
    double *ans = NULL;

#ifdef TIMER
    enum { usecs_in_sec = 1000000, usecs_in_msec = 1000 };
    struct timeval t_start, t_end;
    long elapsed_time, secs, msecs;
#endif

    if (!x || !size || !window) {
        fprintf(stderr, "kz(): Incorrect params\n");
        return NULL;
    }

    win_size = malloc(dim * sizeof(int));
    if (!win_size)
        goto quit;

    for (i = 0; i < dim; i++)
        win_size[i] = floor(window[i] / 2.0);

    switch (dim) {
    case 1:
        ans = malloc(size[0] * sizeof(double));
        if (!ans)
            goto quit;

#ifdef TIMER
        gettimeofday(&t_start, NULL);
#endif
        ans = kz1d(x, size[0], win_size[0], iterations); 
#ifdef TIMER
        gettimeofday(&t_end, NULL);
        elapsed_time = (t_end.tv_sec - t_start.tv_sec) * 1e6 + 
            t_end.tv_usec - t_start.tv_usec;
        secs = elapsed_time / usecs_in_sec;
        msecs = elapsed_time % usecs_in_sec;
        msecs /= usecs_in_msec;
        printf("elapsed time %lds, %ldms (%ldusec)\n", secs, msecs, elapsed_time);
#endif
        break;
    default:
        fprintf(stderr, "kza: Too many dimensions\n");
    }

quit:
    free(win_size);

    return ans;
}
