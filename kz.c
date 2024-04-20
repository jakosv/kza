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

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef TIMER
#include "timer.h"
#endif

#include <pthread.h>
#include <sys/sysinfo.h>

struct task_data {
    int window;
    int task_size;
    int data_size;
    int start_idx, end_idx;
    double *ans;

#ifdef PREFIX_SUM
    double *pref_sum;
    int *pref_finite_cnt;
#else
    double *data;
#endif
};

#ifdef PREFIX_SUM
static double mavg1d(const double *pref_sum, const int *pref_finite_cnt,
                     int data_size, int window_center, int w)
{
    double s;
    int z;
    int start_idx, end_idx;

    /* window length is 2*w+1 */
    start_idx = (window_center+1) - w - 1; /* index of value before window */
    start_idx *= (start_idx > 0);

    end_idx = (window_center+1) + w; /* the last window value index */
    if (end_idx > data_size)
        end_idx = data_size;

    /* (window sum) = (sum containig window) - (sum before window) */
    s = (pref_sum[end_idx] - pref_sum[start_idx]);
    z = (pref_finite_cnt[end_idx] - pref_finite_cnt[start_idx]);

    return s / (double)z;
}

static void calc_prefix_sum(const double *data, int size, double *pref_sum, 
                            int *pref_finite_cnt)
{
    int i;

    pref_sum[0] = 0;
    pref_finite_cnt[0] = 0;
    for (i = 1; i <= size; i++) {
        int is_finite_flag = isfinite(data[i-1]);
        pref_sum[i] = pref_sum[i-1] + is_finite_flag * data[i-1];
        pref_finite_cnt[i] = pref_finite_cnt[i-1] + is_finite_flag;
    }
}

#else
static double mavg1d(const double *x, int data_size, int window_center, int w)
{
    double s = 0;
    long z;
    int i, start_idx, end_idx;

    /* window length is 2*w+1 */
    start_idx = window_center - w;
    start_idx *= (start_idx > 0);

    end_idx = window_center + w + 1;
    if (end_idx > data_size)
        end_idx = data_size;

    z = 0;
    for (i = start_idx; i < end_idx; i++) {
        int is_finite_flag = isfinite(x[i]);
        z += is_finite_flag;
        s += is_finite_flag * x[i];
    }

    return s / (double)z;
}
#endif


static void perform_task_iteration(struct task_data *task, int start_idx,
                                   int end_idx)
{
    int i;
    for (i = start_idx; i < end_idx; i++) {
#ifdef PREFIX_SUM
        task->ans[i] = mavg1d(task->pref_sum, task->pref_finite_cnt,
                              task->data_size, i, task->window);
#else
        task->ans[i] = mavg1d(task->data, task->data_size, i, task->window);
#endif
    }
}

static void *worker(void *data)
{
    struct task_data *task = (struct task_data *)data;
    perform_task_iteration(task, task->start_idx, task->end_idx);

    pthread_exit(NULL);
} 

static void init_tasks(struct task_data **tasks, int tasks_cnt,
                       const struct task_data *task)
{
    int i, start_idx;
    start_idx = 0;
    for (i = 0; i < tasks_cnt; i++) {
        tasks[i] = malloc(sizeof(struct task_data));
        memcpy(tasks[i], task, sizeof(struct task_data));
        tasks[i]->start_idx = start_idx;
        tasks[i]->end_idx = start_idx + task->task_size;
        start_idx = tasks[i]->end_idx;
    }
}


static void start_threads(pthread_t *th, int threads_cnt,
                          struct task_data **tasks)
{
    int i;
    for (i = 0; i < threads_cnt; i++) {
        int res;

        res = pthread_create(&th[i], NULL, worker, tasks[i]);
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


static double *kz1d(const double *x, int length, int window, int iterations)
{
    int mem_size;
    double *ans = NULL;
#ifdef PREFIX_SUM
    double *pref_sum = NULL;
    int *pref_finite_cnt = NULL;
#else
    double *data = NULL;
#endif
    int i, k, tasks_cnt, task_size, threads_cnt;
    pthread_t *th = NULL;
    struct task_data **tasks = NULL;
    struct task_data task;

    threads_cnt = get_nprocs();

#  ifdef DEBUG
    printf("Number of cores: %d\n", threads_cnt);
#  endif

    tasks_cnt = threads_cnt - 1;
    th = malloc(tasks_cnt * sizeof(pthread_t));
    tasks = malloc(tasks_cnt * sizeof(struct task_data*));
    if (!th || !tasks)
        goto quit;

    mem_size = length * sizeof(double);
    ans = malloc(mem_size);

#ifdef PREFIX_SUM
    pref_sum = malloc((length+1) * sizeof(double));
    pref_finite_cnt = malloc((length+1) * sizeof(int));
    if (!ans || !pref_sum || !pref_finite_cnt)
        goto quit;

    calc_prefix_sum(x, length, pref_sum, pref_finite_cnt);
#else
    data = malloc(mem_size);
    if (!ans || !data)
        goto quit;

    memcpy(data, x, mem_size);
#endif

    task.window = window;
    task.data_size = length;
    task.ans = ans;
#ifdef PREFIX_SUM
    task.pref_sum = pref_sum;
    task.pref_finite_cnt = pref_finite_cnt;
#else
    task.data = data;
#endif

    task_size = length / tasks_cnt;
    task.task_size = task_size;
    task.start_idx = 0;
    
    init_tasks(tasks, tasks_cnt, &task);

    for (k = 0; k < iterations; k++) {
        start_threads(th, tasks_cnt, tasks);

        for (i = (tasks_cnt)*task_size; i < length; i++) {
#ifdef PREFIX_SUM
            ans[i] = mavg1d(pref_sum, pref_finite_cnt, length, i, window);
#else
            ans[i] = mavg1d(data, length, i, window); 
#endif
        }

        wait_threads(th, tasks_cnt);

#ifdef PREFIX_SUM
        calc_prefix_sum(ans, length, pref_sum, pref_finite_cnt);
#else
        memcpy(data, ans, mem_size);
#endif
    }


quit:
    free(th);
    if (tasks != NULL) {
        free_tasks(tasks, tasks_cnt);
        free(tasks);
    }

#ifdef PREFIX_SUM
    free(pref_sum);
    free(pref_finite_cnt);
#else
    free(data);
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
    struct timeval t_start;
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
        timestamp(&t_start);
#endif
        ans = kz1d(x, size[0], win_size[0], iterations); 
#ifdef TIMER
        print_timer(&t_start, "kz():");
#endif
        break;
    default:
        fprintf(stderr, "kza: Too many dimensions\n");
    }

quit:
    free(win_size);

    return ans;
}
