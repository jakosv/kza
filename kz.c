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
#include <semaphore.h>

enum { max_threads_cnt = 2, min_thread_task_size = 10 };

struct task_data {
    int thread_task_size;
    int iterations;
    int window;
    int data_size;
    double *ans;
    double *data;
    double *pref_sum;
    int *pref_finite_cnt;
};

struct thread_data {
    pthread_t id;
    struct task_data *task;
    int start_idx, end_idx;
    sem_t *done_work;
    sem_t can_work;
};

#endif

static double mavg1d(const double *pref_sum, const int *pref_finite_cnt,
                     int length, int col, int w)
{
    double s;
    int z;
    int start_col, end_col;

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


static void update_iteration_data(double *new, const double *old, int size,
                                  double *pref_sum, int *pref_finite_cnt)
{
    memcpy(new, old, size * sizeof(double)); 
    calc_prefix_sum(new, size, pref_sum, pref_finite_cnt);
}

#ifdef KZ_PARALLEL

static void *worker(void *data)
{
    int k;
    struct thread_data *thread_data = data;
    struct task_data *task = thread_data->task;

    for (k = 0; k < task->iterations; k++) {
        int i;
        
        sem_wait(&thread_data->can_work);

        for (i = thread_data->start_idx; i < thread_data->end_idx; i++) {
            task->ans[i] = mavg1d(task->pref_sum, task->pref_finite_cnt,
                                  task->data_size, i, task->window);
        }

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
        th[i].end_idx = start_idx + task->thread_task_size;
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

static void free_threads(struct thread_data *th, int threads_cnt)
{
    int i;
    for (i = 0; i < threads_cnt; i++) {
        pthread_join(th[i].id, NULL);
        sem_destroy(&th[i].can_work);
    }
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

    free_threads(th, threads_cnt);
    sem_destroy(&done_work_sem);
}

#endif

static double *kz1d(const double *x, int length, int window, int iterations)
{
    int mem_size;
    double *data = NULL, *ans = NULL, *pref_sum = NULL;
    int *pref_finite_cnt = NULL;

#ifdef KZ_PARALLEL
    struct thread_data th[max_threads_cnt];
    int threads_cnt, thread_task_size;
#else
    int i, k;
#endif

    mem_size = length * sizeof(double);
    ans = malloc(mem_size);
    data = malloc(mem_size);
    pref_sum = malloc((length+1) * sizeof(double));
    pref_finite_cnt = malloc((length+1) * sizeof(int));
    if (!ans || !data || !pref_sum || !pref_finite_cnt)
        goto quit;

    update_iteration_data(data, x, length, pref_sum, pref_finite_cnt);

#ifdef KZ_PARALLEL
    threads_cnt = max_threads_cnt;
    thread_task_size = (length + threads_cnt-1) / threads_cnt;

    struct task_data task;
    task.thread_task_size = thread_task_size;
    task.iterations = iterations;
    task.window = window;
    task.data_size = length;
    task.data = data;
    task.ans = ans;
    task.pref_sum = pref_sum;
    task.pref_finite_cnt = pref_finite_cnt;

    threads_server_loop(th, threads_cnt, &task);

#else
    for (k = 0; k < iterations; k++) {
        for (i = 0; i < length; i++) {
            ans[i] = mavg1d(pref_sum, pref_finite_cnt, length, i, window);
        }
        update_iteration_data(data, ans, length, pref_sum, pref_finite_cnt);
    }
#endif


quit:
    free(data);
    free(pref_sum);
    free(pref_finite_cnt);

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
