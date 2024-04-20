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
#include <sys/time.h>
#endif


#if defined(KZ_THREADS_STRATEGY_1) || defined(KZ_THREADS_STRATEGY_2)
#include <pthread.h>

#  ifdef KZ_THREADS_STRATEGY_1
#include <semaphore.h>
#  endif

#  ifdef _WIN32
#include <windows.h>

#  elifdef MACOS
#include <sys/param.h>
#include <sys/sysctl.h>

#  else
#include <unistd.h>
#  endif

int get_num_cores() {
#  ifdef WIN32
    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    return sysinfo.dwNumberOfProcessors;

#  elifdef MACOS
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

#  else
    return sysconf(_SC_NPROCESSORS_ONLN);
#  endif
}

struct task_data {
    int window;
    int task_size;
    int data_size;
    double *ans;
    double *pref_sum;
    int *pref_finite_cnt;

#  ifdef KZ_THREADS_STRATEGY_1
    int iterations;
    double *data;

#  elifdef KZ_THREADS_STRATEGY_2
    int start_idx, end_idx;
#  endif
};

#endif

static double mavg1d(const double *pref_sum, const int *pref_finite_cnt,
                     int data_size, int window_center, int w)
{
    double s;
    int z;
    int start_idx, end_idx;

    /* window length is 2*w+1 */
    start_idx = (window_center+1) - w; /* the first window value index */
    if (start_idx < 0)
        start_idx = 0;

    end_idx = (window_center+1) + w; /* the last window value index */
    if (end_idx > data_size)
        end_idx = data_size;

    /* (window sum) = (sum containig window) - (sum before window) */
    s = (pref_sum[end_idx] - pref_sum[start_idx-1]);
    z = (pref_finite_cnt[end_idx] - pref_finite_cnt[start_idx-1]);
    /*
    if (z == 0) 
        return nan("");
    */
    return s/z;
}

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


#if defined(KZ_THREADS_STRATEGY_1) || defined(KZ_THREADS_STRATEGY_2)

static void perform_task_iteration(struct task_data *task, int start_idx,
                                   int end_idx)
{
    int i;
    for (i = start_idx; i < end_idx; i++) {
        task->ans[i] = mavg1d(task->pref_sum, task->pref_finite_cnt,
                              task->data_size, i, task->window);
    }
}

#endif


#ifdef KZ_THREADS_STRATEGY_1

struct thread_data {
    pthread_t id;
    struct task_data *task;
    int start_idx, end_idx;
    sem_t *done_work;
    sem_t can_work;
};

static void *worker(void *data)
{
    int k;

    struct thread_data *thread_data = data;
    struct task_data *task = thread_data->task;

    for (k = 0; k < task->iterations; k++) {
        
        sem_wait(&thread_data->can_work);

        perform_task_iteration(task, thread_data->start_idx,
                               thread_data->end_idx);

        sem_post(thread_data->done_work);
    }

    return NULL;
} 

static void start_threads(struct thread_data *th, int threads_cnt,
                          struct task_data *task, sem_t *done_work)
{
    int i, start_idx;

    start_idx = 0;
    for (i = 0; i < threads_cnt; i++) {
        int res;

        th[i].task = task;
        th[i].done_work = done_work;
        th[i].start_idx = start_idx;
        th[i].end_idx = start_idx + task->task_size;
        if (th[i].end_idx > task->data_size)
            th[i].end_idx = task->data_size;
        else
            start_idx = th[i].end_idx;
        sem_init(&th[i].can_work, 0, 1);
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
        sem_destroy(&th[i].can_work);
    }
}

static void update_iteration_data(double *new, const double *old, int size,
                                  double *pref_sum, int *pref_finite_cnt)
{
    memcpy(new, old, size * sizeof(double)); 
    calc_prefix_sum(new, size, pref_sum, pref_finite_cnt);
}

static void threads_server_loop(struct thread_data *th, int threads_cnt,
                                struct task_data *task)
{
    int k, done_workers;
    sem_t done_work_sem;

    sem_init(&done_work_sem, 0, 0);

    start_threads(th, threads_cnt, task, &done_work_sem);

    done_workers = 0;
    k = 0;
    while (k < task->iterations) {
        sem_wait(&done_work_sem);
        done_workers++; 
        if (done_workers == threads_cnt) {
            int i;
            update_iteration_data(task->data, task->ans, task->data_size,
                                  task->pref_sum, task->pref_finite_cnt);
            done_workers = 0;
            for (i = 0; i < threads_cnt; i++)
                sem_post(&th[i].can_work);
            k++;
        }
    } 

    wait_threads(th, threads_cnt);
    sem_destroy(&done_work_sem);
}

#elifdef KZ_THREADS_STRATEGY_2

static void *worker(void *data)
{
    struct task_data *task = (struct task_data *)data;
    perform_task_iteration(task, task->start_idx, task->end_idx);

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

#endif


static double *kz1d(const double *x, int length, int window, int iterations)
{
    int mem_size;
    double *ans = NULL, *pref_sum = NULL;
    int *pref_finite_cnt = NULL;

#ifdef KZ_THREADS_STRATEGY_1
    struct thread_data *th;
    struct task_data task;
    double *data;
#elifdef KZ_THREADS_STRATEGY_2
    int i, k, tasks_cnt;
    pthread_t *th;
    struct task_data **tasks;
#else
    int i, k;
#endif

#if defined(KZ_THREADS_STRATEGY_1) || defined(KZ_THREADS_STRATEGY_2)
    int task_size, threads_cnt;

    threads_cnt = get_num_cores();

#  ifdef DEBUG
    printf("Number of cores: %d\n", threads_cnt);
#  endif

#  ifdef KZ_THREADS_STRATEGY_1
    th = malloc(threads_cnt * sizeof(struct thread_data));
    data = malloc(length * sizeof(double));
    if (!data)
        goto quit;

#  elifdef KZ_THREADS_STRATEGY_2
    tasks_cnt = threads_cnt - 1;
    th = malloc(tasks_cnt * sizeof(pthread_t));
    tasks = malloc(tasks_cnt * sizeof(struct task_data*));
#  endif

#endif

    mem_size = length * sizeof(double);
    ans = malloc(mem_size);
    pref_sum = malloc((length+1) * sizeof(double));
    pref_finite_cnt = malloc((length+1) * sizeof(int));
    if (!ans || !pref_sum || !pref_finite_cnt)
        goto quit;

    calc_prefix_sum(x, length, pref_sum, pref_finite_cnt);

#ifdef KZ_THREADS_STRATEGY_1
    task_size = (length + threads_cnt-1) / threads_cnt;

    task.task_size = task_size;
    task.iterations = iterations;
    task.window = window;
    task.data_size = length;
    task.data = data;
    task.ans = ans;
    task.pref_sum = pref_sum;
    task.pref_finite_cnt = pref_finite_cnt;

    threads_server_loop(th, threads_cnt, &task);

#elifdef KZ_THREADS_STRATEGY_2
    task_size = length / tasks_cnt;
    
    init_tasks(tasks, tasks_cnt, task_size, window, length,
               ans, pref_sum, pref_finite_cnt);

    for (k = 0; k < iterations; k++) {
        start_threads(th, tasks_cnt, tasks);

        for (i = (tasks_cnt)*task_size; i < length; i++)
            ans[i] = mavg1d(pref_sum, pref_finite_cnt, length, i, window);

        wait_threads(th, tasks_cnt);

        calc_prefix_sum(ans, length, pref_sum, pref_finite_cnt);
    }

#else
    for (k = 0; k < iterations; k++) {
        for (i = 0; i < length; i++)
            ans[i] = mavg1d(pref_sum, pref_finite_cnt, length, i, window);

        calc_prefix_sum(ans, length, pref_sum, pref_finite_cnt);
    }
#endif

quit:
    free(pref_sum);
    free(pref_finite_cnt);

#if defined(KZ_THREADS_STRATEGY_1) || defined(KZ_THREADS_STRATEGY_2)
    free(th);
#  ifdef KZ_THREADS_STRATEGY_1
    free(data);
#  elifdef KZ_THREADS_STRATEGY_2
    if (tasks != NULL) {
        free_tasks(tasks, tasks_cnt);
        free(tasks);
    }
#  endif
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
        printf("kz(): elapsed time %lds, %ldms (%ld usec)\n", 
               secs, msecs, elapsed_time);
#endif
        break;
    default:
        fprintf(stderr, "kza: Too many dimensions\n");
    }

quit:
    free(win_size);

    return ans;
}
