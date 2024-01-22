#include "kza.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX(A, B) ((A) > (B) ? (A) : (B))

#ifdef TIMER
#include <sys/time.h>
#endif

#if defined(KZA_THREADS_CLIENT_SERVER) || defined(KZA_THREADS_LOOP)
#include <pthread.h>
#include <sys/sysinfo.h>

#  ifdef KZA_THREADS_CLIENT_SERVER
#include <semaphore.h>
#  endif

struct task_data {
    int window;
    int task_size;
    int data_size;
    int min_window_len;
    double max_d;
    double tolerance;
    double *data;
    double *ans;
    double *d;
    double *dprime;

#  ifdef KZA_THREADS_CLIENT_SERVER
    int iterations;

#  elif KZA_THREADS_LOOP
    int start_idx, end_idx;
#  endif
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


#if defined(KZA_THREADS_CLIENT_SERVER) || defined(KZA_THREADS_LOOP)

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

#endif

#ifdef KZA_THREADS_CLIENT_SERVER

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

    return NULL;
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

#elif KZA_THREADS_LOOP

static void *worker(void *data)
{
    struct task_data *task = (struct task_data *)data;
    perform_task_iteration(task, task->start_idx, task->end_idx);

    return NULL;
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

#ifdef KZA_THREADS_CLIENT_SERVER
    struct thread_data *th = NULL;
    struct task_data task;
#elif KZA_THREADS_LOOP
    int i, tasks_cnt;
    pthread_t *th = NULL;
    struct task_data task;
    struct task_data **tasks = NULL;
#else
    int i;
#endif

#if defined(KZA_THREADS_CLIENT_SERVER) || defined(KZA_THREADS_LOOP)
    int task_size, threads_cnt;

    threads_cnt = get_nprocs();

#  ifdef DEBUG
    printf("Number of cores: %d\n", threads_cnt);
#  endif

#  ifdef KZA_THREADS_CLIENT_SERVER
    th = malloc(threads_cnt * sizeof(struct thread_data));
    if (!th)
        goto quit;

#  elif KZA_THREADS_LOOP
    tasks_cnt = threads_cnt - 1;
    th = malloc(tasks_cnt * sizeof(pthread_t));
    tasks = malloc(tasks_cnt * sizeof(struct task_data*));
    if (!th || !tasks)
        goto quit;
#  endif

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

#if defined(KZA_THREADS_CLIENT_SERVER) || defined(KZA_THREADS_LOOP)
    task.window = window;
    task.min_window_len = min_window_len;
    task.tolerance = tolerance;
    task.ans = ans;
    task.data = tmp;
    task.data_size = n;
    task.d = d;
    task.dprime = dprime;
    task.max_d = m;
#endif

#ifdef KZA_THREADS_CLIENT_SERVER
    task_size = (n + threads_cnt-1) / threads_cnt;

    task.task_size = task_size;
    task.iterations = iterations;
    task.data = tmp;

    threads_server_loop(th, threads_cnt, &task);

#elif KZA_THREADS_LOOP
    task_size = n / tasks_cnt;
    
    task.task_size = task_size;
    task.start_idx = 0;

    init_tasks(tasks, tasks_cnt, &task);

    for (i = 0; i < iterations; i++) {
        int t;

        start_threads(th, tasks_cnt, tasks);
        for (t = (tasks_cnt)*task_size; t < n; t++) {
            int left_win, right_win, left_bound, right_bound;

            calc_adaptive_windows(&left_win, &right_win, d, dprime, t,
                                  window, tolerance, m);
            get_window_bounds(&left_bound, &right_bound,
                              left_win, right_win, t, n, min_window_len);

            ans[t] = mavg1d(tmp, left_bound, right_bound+1); 
        }

        wait_threads(th, tasks_cnt);
        memcpy(tmp, ans, mem_size);
    }

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
#if defined(KZA_THREADS_CLIENT_SERVER) || defined(KZA_THREADS_LOOP)
    free(th);
#  ifdef KZA_THREADS_LOOP
    if (tasks != NULL) {
        free_tasks(tasks, tasks_cnt);
        free(tasks);
    }
#  endif
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
