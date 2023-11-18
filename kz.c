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

enum { max_threads_cnt = 4, min_thread_task_size = 10 };

struct task_data {
    int thread_task_size;
    int iterations;
    int window;
    int data_size;
    double *data;
    double *ans;
};

struct thread_data {
    pthread_t id;
    struct task_data *task;
    int start_idx, end_idx;
    sem_t *done_work;
    sem_t can_work;
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
    if (start_col < 0)
        start_col = 0;

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

    /*
    if (z == 0) 
        return nan("");
    */
    return s/z;
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
            task->ans[i] = mavg1d(task->data, task->data_size, i, 
                                  task->window);
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
    int k, done_workers, data_mem_size;
    sem_t done_work_sem;

    sem_init(&done_work_sem, 0, 0);

    start_threads(th, threads_cnt, task, &done_work_sem);

    data_mem_size = task->data_size * sizeof(double);
    done_workers = 0;
    k = 0;
    while (k < task->iterations) {
        sem_wait(&done_work_sem);
        done_workers++; 
        if (done_workers == threads_cnt) {
            int i;
            memcpy(task->data, task->ans, data_mem_size); 
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
    int i, k;
    int mem_size;
    double *tmp = NULL, *ans = NULL;

#ifdef KZ_PARALLEL
    struct thread_data th[max_threads_cnt];
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
    thread_task_size = (length + threads_cnt-1) / threads_cnt;

    if (thread_task_size < min_thread_task_size) {
        threads_cnt = 1;
        for (k = 0; k < iterations; k++) {
            for (i = 0; i < length; i++) {
                ans[i] = mavg1d(tmp, length, i, window);
            }
            memcpy(tmp, ans, mem_size); 
        }
    } else {
        struct task_data task;
        task.thread_task_size = thread_task_size;
        task.iterations = iterations;
        task.window = window;
        task.data_size = length;
        task.data = tmp;
        task.ans = ans;
        threads_server_loop(th, threads_cnt, &task);
    }
#else
    for (k = 0; k < iterations; k++) {
        for (i = 0; i < length; i++) {
            ans[i] = mavg1d(tmp, length, i, window);
        }
        memcpy(tmp, ans, mem_size); 
    }
#endif


quit:
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
    free(win_size);

    return ans;
}
