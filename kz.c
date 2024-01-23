#include "kza.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef TIMER
#include <sys/time.h>
#endif


#if defined(KZ_THREADS_CLIENT_SERVER) || defined(KZ_THREADS_LOOP)
#include <pthread.h>
#include <sys/sysinfo.h>

#  ifdef KZ_THREADS_CLIENT_SERVER
#include <semaphore.h>
#  endif

struct task_data {
    int window;
    int task_size;
    int data_size;
    double *ans;

#  ifdef PREFIX_SUM
    double *pref_sum;
    int *pref_finite_cnt;
#  else
    double *data;
#  endif

#  ifdef KZ_THREADS_CLIENT_SERVER
    int iterations;

#  elif KZ_THREADS_LOOP
    int start_idx, end_idx;
#  endif
};

#endif

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


#if defined(KZ_THREADS_CLIENT_SERVER) || defined(KZ_THREADS_LOOP)

static void perform_task_iteration(struct task_data *task, int start_idx,
                                   int end_idx)
{
    int i;
    for (i = start_idx; i < end_idx; i++) {
#  ifdef PREFIX_SUM
        task->ans[i] = mavg1d(task->pref_sum, task->pref_finite_cnt,
                              task->data_size, i, task->window);
#  else
        task->ans[i] = mavg1d(task->data, task->data_size, i, task->window);
#  endif
    }
}

#endif


#ifdef KZ_THREADS_CLIENT_SERVER

struct thread_data {
    pthread_t id;
    struct task_data *task;
    int start_idx, end_idx;
    sem_t *finished_workers_sem;
    sem_t can_work_sem;
};

static void *worker(void *data)
{
    int k;

    struct thread_data *thread_data = data;
    struct task_data *task = thread_data->task;

    for (k = 0; k < task->iterations; k++) {
        
        sem_wait(&thread_data->can_work_sem);

        perform_task_iteration(task, thread_data->start_idx,
                               thread_data->end_idx);

        sem_post(thread_data->finished_workers_sem);
    }

    pthread_exit(NULL);
} 

static void start_threads(struct thread_data *th, int threads_cnt,
                          struct task_data *task, sem_t *finished_workers_sem)
{
    int i, start_idx;

    start_idx = 0;
    for (i = 0; i < threads_cnt; i++) {
        int res;

        th[i].task = task;
        th[i].finished_workers_sem = finished_workers_sem;
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
#  ifdef PREFIX_SUM
            calc_prefix_sum(task->ans, task->data_size, task->pref_sum, 
                            task->pref_finite_cnt);
#  else
            memcpy(task->data, task->ans, task->data_size * sizeof(double));
#  endif
            idle_workers_count = 0;
            for (i = 0; i < threads_cnt; i++)
                sem_post(&th[i].can_work_sem);
            iter++;
        }
    } 

    wait_threads(th, threads_cnt);
    sem_destroy(&finished_workers_sem);
}

#elif KZ_THREADS_LOOP

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

#endif


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

#ifdef KZ_THREADS_CLIENT_SERVER
    struct thread_data *th = NULL;
#elif KZ_THREADS_LOOP
    int i, k, tasks_cnt;
    pthread_t *th = NULL;
    struct task_data **tasks = NULL;
#else
    int i, k;
#endif

#if defined(KZ_THREADS_CLIENT_SERVER) || defined(KZ_THREADS_LOOP)
    struct task_data task;
    int task_size, threads_cnt;

    threads_cnt = get_nprocs();

#  ifdef DEBUG
    printf("Number of cores: %d\n", threads_cnt);
#  endif

#  ifdef KZ_THREADS_CLIENT_SERVER
    th = malloc(threads_cnt * sizeof(struct thread_data));
    if (!th)
        goto quit;

#  elif KZ_THREADS_LOOP
    tasks_cnt = threads_cnt - 1;
    th = malloc(tasks_cnt * sizeof(pthread_t));
    tasks = malloc(tasks_cnt * sizeof(struct task_data*));
    if (!th || !tasks)
        goto quit;
#  endif

#endif

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


#if defined(KZ_THREADS_CLIENT_SERVER) || defined(KZ_THREADS_LOOP)
    task.window = window;
    task.data_size = length;
    task.ans = ans;
#  ifdef PREFIX_SUM
    task.pref_sum = pref_sum;
    task.pref_finite_cnt = pref_finite_cnt;
#  else
    task.data = data;
#  endif
#endif

#ifdef KZ_THREADS_CLIENT_SERVER
    task_size = (length + threads_cnt-1) / threads_cnt;
    task.task_size = task_size;
    task.iterations = iterations;

    threads_server_loop(th, threads_cnt, &task);

#elif KZ_THREADS_LOOP
    task_size = length / tasks_cnt;
    task.task_size = task_size;
    task.start_idx = 0;
    
    init_tasks(tasks, tasks_cnt, &task);

    for (k = 0; k < iterations; k++) {
        start_threads(th, tasks_cnt, tasks);

        for (i = (tasks_cnt)*task_size; i < length; i++) {
#  ifdef PREFIX_SUM
            ans[i] = mavg1d(pref_sum, pref_finite_cnt, length, i, window);
#  else
            ans[i] = mavg1d(data, length, i, window); 
#  endif
        }

        wait_threads(th, tasks_cnt);

#  ifdef PREFIX_SUM
        calc_prefix_sum(ans, length, pref_sum, pref_finite_cnt);
#  else
        memcpy(data, ans, mem_size);
#  endif
    }

#else
    for (k = 0; k < iterations; k++) {
        for (i = 0; i < length; i++) {
#  ifdef PREFIX_SUM
            ans[i] = mavg1d(pref_sum, pref_finite_cnt, length, i, window);
#  else
            ans[i] = mavg1d(data, length, i, window); 
#  endif
        }

#  ifdef PREFIX_SUM
        calc_prefix_sum(ans, length, pref_sum, pref_finite_cnt);
#  else
        memcpy(data, ans, mem_size);
#  endif
    }
#endif

quit:
#if defined(KZ_THREADS_CLIENT_SERVER) || defined(KZ_THREADS_LOOP)
    free(th);
#  ifdef KZ_THREADS_LOOP
    if (tasks != NULL) {
        free_tasks(tasks, tasks_cnt);
        free(tasks);
    }
#  endif
#endif

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
