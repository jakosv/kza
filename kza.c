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


#include "kza.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX(A, B) ((A) > (B) ? (A) : (B))

#ifdef TIMER
#include <sys/time.h>
#endif

#ifdef KZA_THREADS_STRATEGY_1

#include <pthread.h>
#include <semaphore.h>
#include <sys/sysinfo.h>

struct task_data {
    int window;
    int iterations;
    int task_size;
    int data_size;
    int min_window_len;
    double max_d;
    double tolerance;
    double *ans;
    double *data;
    double *d;
    double *dprime;
};

struct thread_data {
    pthread_t id;
    struct task_data *task;
    int start_idx, end_idx;
    sem_t *done_work_sem;
    sem_t can_work_sem;
};

#endif

static double mavg1d(const double *x, int a, int b)
{
    double s = 0;
    long i, z;

    for (i = a, z = 0; i < b; i++) {
        if (isfinite(x[i])) {
            z++;
            s += x[i];
        }
    }

    if (z == 0) 
        return nan("");
    return s/z;
}

static double adaptive(double d, double m)
{
    return 1 - d/m;
}

static void normalize_window(int *left_win, int *right_win, int n, 
                             int window_center, int min_window_len)
{
    int max_right_bound;

    if (*left_win < min_window_len)
        *left_win = min_window_len;
    if (*right_win < min_window_len)
        *right_win = min_window_len;

    /* check bounds */
    max_right_bound = n - window_center - 1;
    if (*right_win > max_right_bound)
        *right_win = max_right_bound;
    if (*left_win > window_center)
        *left_win = window_center;
}

void calc_adaptive_windows(int *left_win, int *right_win,
                           const double *d, const double *dprime, int t,
                           int window, double eps, double max_d)  
{
    int window_adaptive;

    window_adaptive = (int)floor(window * adaptive(d[t], max_d));

    if (fabs(dprime[t]) < eps) { /* dprime[t] = 0 */
        *left_win = window_adaptive;
        *right_win = window_adaptive;
    } else if (dprime[t] < 0) {
        *right_win = window;
        *left_win = window_adaptive;
    } else {
        *right_win = window_adaptive;
        *left_win = window;
    }
}

void get_window_bounds(int *left_bound, int *right_bound,
                       int left_win, int right_win, int t, 
                       int data_size, int min_window_len)
{
    normalize_window(&left_win, &right_win, data_size, t, min_window_len);
    *left_bound = t - left_win;
    *right_bound = t + right_win;
}


#ifdef KZA_THREADS_STRATEGY_1

static void perform_task_iteration(struct task_data *task, int start_idx,
                                   int end_idx)
{
    int t;
    for (t = start_idx; t < end_idx; t++) {
        int left_win, right_win, left_bound, right_bound;

        calc_adaptive_windows(&left_win, &right_win, 
                              task->d, task->dprime, t,
                              task->window, task->tolerance, task->max_d);
        get_window_bounds(&left_bound, &right_bound, left_win, right_win, 
                          t, task->data_size, task->min_window_len);

        task->ans[t] = mavg1d(task->data, left_bound, right_bound+1);
    }
}

static void *worker(void *data)
{
    int k;

    struct thread_data *thread_data = data;
    struct task_data *task = thread_data->task;

    for (k = 0; k < task->iterations; k++) {
        sem_wait(&thread_data->can_work_sem);

        perform_task_iteration(task, thread_data->start_idx,
                               thread_data->end_idx);

        sem_post(thread_data->done_work_sem);
    }

    return NULL;
} 

static void start_threads(struct thread_data *th, int threads_cnt,
                          struct task_data *task, sem_t *done_work_sem)
{
    int i, start_idx;

    start_idx = 0;
    for (i = 0; i < threads_cnt; i++) {
        int res;

        th[i].task = task;
        th[i].done_work_sem = done_work_sem;
        th[i].start_idx = start_idx;
        th[i].end_idx = start_idx + task->task_size;

        if (th[i].end_idx > task->data_size)
            th[i].end_idx = task->data_size;
        else
            start_idx = th[i].end_idx;

        sem_init(&th[i].can_work_sem, 0, 1);

        res = pthread_create(&th[i].id, NULL, worker, (void*)&th[i]);
        if (res != 0) {
            perror("thread_create");
            exit(1);
        }
    }
}

static void wait_threads(struct thread_data *th, int threads_cnt)
{
    int i;
    for (i = 0; i < threads_cnt; i++) {
        pthread_join(th[i].id, NULL);
        sem_destroy(&th[i].can_work_sem);
    }
}

static void update_iteration_data(double *new, const double *old, int size)
{
    memcpy(new, old, size * sizeof(double)); 
}

static void threads_server_loop(struct thread_data *th, int threads_cnt,
                                struct task_data *task)
{
    int iter, idle_workers_count;
    sem_t finished_workers_sem;

    sem_init(&finished_workers_sem, 0, 0);

    start_threads(th, threads_cnt, task, &finished_workers_sem);

    idle_workers_count = 0;
    iter = 0;
    while (iter < task->iterations) {
        sem_wait(&finished_workers_sem);
        idle_workers_count++; 
        if (idle_workers_count == threads_cnt) {
            int i;
            update_iteration_data(task->data, task->ans, task->data_size);
            idle_workers_count = 0;
            for (i = 0; i < threads_cnt; i++)
                sem_post(&th[i].can_work_sem);
            iter++;
        }
    } 

    wait_threads(th, threads_cnt);
    sem_destroy(&finished_workers_sem);
}

#endif

static double maximum(const double *v, int length) 
{
    double m;
    int i;

    for (i = 0, m = v[0]; i < length; i++)
        if (isfinite(v[i])) 
            m = MAX(v[i], m);

    return m;
}

static void differenced(const double *y, double *d, double *dprime,
                        int length, int q)
{
    int i;
    long n;

    n = length;    

    /* calculate d = |Z(i+q) - Z(i-q)| */
    for (i = 0; i < q; i++)
        d[i] = fabs(y[i+q] - y[0]);
    for (i = q; i < n-q; i++)
        d[i] = fabs(y[i+q] - y[i-q]);
    for (i = n-q; i < n; i++)
        d[i] = fabs(y[n-1] - y[i-q]);

    /* d'(t) = d(i+1)-d(i) */
    for (i = 0; i < n-1; i++)
        dprime[i] = d[i+1] - d[i];
    dprime[n-1] = dprime[n-2];
}

static double *kza1d(const double *v, int n, const double *y, int window,
                     int iterations, int min_window_len, double tolerance)
{
    int mem_size;
    double m;
    double *d = NULL;
    double *dprime = NULL;
    double *tmp = NULL;
    double *ans = NULL;

#ifdef KZA_THREADS_STRATEGY_1
    struct thread_data *th = NULL;
    struct task_data task;
#else
    int i;
#endif

#ifdef KZA_THREADS_STRATEGY_1
    int task_size, threads_cnt;

    threads_cnt = get_nprocs();

#  ifdef DEBUG
    printf("kza(): Number of cores: %d\n", threads_cnt);
#  endif

    th = malloc(threads_cnt * sizeof(struct thread_data));
    if (!th)
        goto quit;

#endif

    mem_size = n * sizeof(double);

    d = malloc(mem_size);
    dprime = malloc(mem_size);
    if (!d || !dprime)
        goto quit;

    memset(dprime, 0, mem_size);
    differenced(y, d, dprime, n, window);

    m = maximum(d, n);

    tmp = malloc(mem_size);
    if (!tmp)
        goto quit;

    memcpy(tmp, v, mem_size);

    ans = malloc(mem_size);
    if (!ans)
        goto quit;

#ifdef KZA_THREADS_STRATEGY_1
    task_size = (n + threads_cnt-1) / threads_cnt;

    task.task_size = task_size;
    task.iterations = iterations;
    task.window = window;
    task.min_window_len = min_window_len;
    task.data_size = n;
    task.tolerance = tolerance;
    task.data = tmp;
    task.ans = ans;
    task.d = d;
    task.dprime = dprime;
    task.max_d = m;

    threads_server_loop(th, threads_cnt, &task);

#else

    for (i = 0; i < iterations; i++) {
        int t;
        for (t = 0; t < n; t++) {
            int left_win, right_win, left_bound, right_bound;

            calc_adaptive_windows(&left_win, &right_win, d, dprime, t,
                                  window, tolerance, m);
            get_window_bounds(&left_bound, &right_bound,
                              left_win, right_win, t, n, min_window_len);

            ans[t] = mavg1d(tmp, left_bound, right_bound+1); 
        }
        memcpy(tmp, ans, mem_size);
    }

#endif

quit:
#ifdef KZA_THREADS_STRATEGY_1
    free(th);
#endif
    free(tmp);
    free(d);
    free(dprime);

    return ans;
}

double *kza(const double *x, int dim, const int *size, const double *y, 
            const int *window, int iterations, double min_window, 
            double tolerance)
{
    double *kz_ans = NULL;
    double *kza_ans = NULL;
    int mem_size;

#ifdef TIMER
    enum { usecs_in_sec = 1000000, usecs_in_msec = 1000 };
    struct timeval t_start, t_end;
    long elapsed_time, secs, msecs;
#endif

    if (!x || !size || !window) {
        fprintf(stderr, "kza: Incorrect params\n");
        return NULL;
    }

    if (dim > 1) {
        fprintf(stderr, "kza: Not yet implemented\n");
        return NULL;
    }

    mem_size = size[0] * sizeof(double);
    kz_ans = malloc(mem_size);
    if (!kz_ans)
        return NULL;

    if (!y) {
        memcpy(kz_ans, x, mem_size);
        kz(kz_ans, dim, size, window, iterations);
    } else {
        memcpy(kz_ans, y, mem_size);
    }

#ifdef TIMER
    gettimeofday(&t_start, NULL);
#endif
    kza_ans = kza1d(x, size[0], kz_ans, window[0], iterations, min_window,
                    tolerance);
#ifdef TIMER
        gettimeofday(&t_end, NULL);
        elapsed_time = (t_end.tv_sec - t_start.tv_sec) * 1e6 + 
            t_end.tv_usec - t_start.tv_usec;
        secs = elapsed_time / usecs_in_sec;
        msecs = elapsed_time % usecs_in_sec;
        msecs /= usecs_in_msec;
        printf("kza(): elapsed time %lds, %ldms (%ld usec)\n", 
               secs, msecs, elapsed_time);
#endif

    free(kz_ans);

    return kza_ans;
}

void kza_free(double *x)
{
    free(x);
}
