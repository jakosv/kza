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

enum { max_threads_cnt = 4, min_thread_task_size = 10 };

struct task_data {
    int start_idx, end_idx;
    int window;
    int task_size;
    int data_size;
    double *data;
    double *ans;
};

struct thread_data {
    pthread_t id;
    struct task_args* args;
};

#endif

static double mavg1d(const double *v, int length, int col, int w)
{
    double s;
    int i, z;
    int start_col, end_col;

    if (!v) {
        fprintf(stderr, "mavg1d(): Incorrect params\n");
        return 0;
    }

    start_col = col - w;
#ifdef SPEED_HACKS
    start_col *= (start_col > 0);
#else
    if (start_col < 0)
        start_col = 0;
#endif

    end_col = col + w + 1;
    if (end_col > length)
        end_col = length;

    s = 0;
    z = 0;
    for (i = start_col; i < end_col; i++) {
        if (isfinite(v[i])) {
            z++;
            s += v[i];
        }
    }

    if (z == 0) 
        return nan("");
    return s/z;
}

#ifdef KZ_PARALLEL

static void *thread_task(void *data)
{
    int i;
    struct task_data *args = data;

    for (i = args->start_idx; i < args->end_idx; i++) {
        args->ans[i] = mavg1d(args->data, args->data_size, i, args->window);
    }

    return NULL;
} 

static void init_tasks(struct task_data **tasks, int tasks_cnt,
                       int task_size, int window, int data_size)
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
        task->data = NULL;
        task->ans = NULL;

        tasks[i] = task;

        start_idx = task->end_idx;
    }

}

static void start_threads(pthread_t *th, int threads_cnt,
                          struct task_data **tasks, double *data, double *ans)
{
    int i;
    for (i = 0; i < threads_cnt; i++) {
        int res;

        tasks[i]->data = data;
        tasks[i]->ans = ans;

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

static double *kz1d(const double *x, int length, int window, int iterations)
{
    int i, k;
    int mem_size;
    double *tmp = NULL, *ans = NULL;

#ifdef KZ_PARALLEL
    pthread_t th[max_threads_cnt-1];
    struct task_data *tasks[max_threads_cnt-1];
    int threads_cnt, thread_task_size;
    threads_cnt = max_threads_cnt;
#endif

    mem_size = length * sizeof(double);
    ans = malloc(mem_size);
    tmp = malloc(mem_size);
    if (!ans || !tmp)
        goto quit;

    memcpy(tmp, x, mem_size);

#ifdef KZ_PARALLEL
    thread_task_size = length / threads_cnt;
    
    if (thread_task_size < min_thread_task_size) {
        threads_cnt = 1;
    } else {
        init_tasks(tasks, threads_cnt-1, thread_task_size, window, length);
    }
#endif

    for (k = 0; k < iterations; k++) {
#ifdef KZ_PARALLEL
        if (threads_cnt > 1) {
            start_threads(th, threads_cnt-1, tasks, tmp, ans);
        } 

        for (i = (threads_cnt-1)*thread_task_size; i < length; i++) {
            ans[i] = mavg1d(tmp, length, i, window);
        }
#else
        for (i = 0; i < length; i++) {
            ans[i] = mavg1d(tmp, length, i, window);
        }
#endif

#ifdef KZ_PARALLEL
        if (threads_cnt > 1)
            wait_threads(th, threads_cnt-1);
#endif

        memcpy(tmp, ans, mem_size); 
    }

#ifdef KZ_PARALLEL
        if (threads_cnt > 1)
            free_tasks(tasks, threads_cnt-1);
#endif

quit:
    if (tmp)
        free(tmp);

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
    if (win_size)
        free(win_size);

    return ans;
}
