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

struct thread_args {
    int start_col, end_col;
    const double *v;
};

struct thread_res {
    double s;
    int z;
};

static void *thread_task(void *data)
{
    int i;
    struct thread_args *args = data;
    struct thread_res *res;

    res = malloc(sizeof(*res));
    res->s = 0;
    res->z = 0;
    for (i = args->start_col; i < args->end_col; i++)
        if (isfinite(args->v[i])) {
            res->s += args->v[i];
            (res->z)++;
        }
    free(args);
    return res;
} 

static void start_threads(pthread_t *th, int threads_cnt, int start_col,
                          int task_size, const double *v)
{
    int i;
    for (i = 0; i < threads_cnt; i++) {
        int res;
        struct thread_args *args; 

        args = malloc(sizeof(struct thread_args));
        args->start_col = start_col;
        args->end_col = start_col + task_size;
        args->v = v;

        res = pthread_create(&th[i], NULL, thread_task, args);
        if (res != 0) {
            perror("thread_create");
            exit(1);
        }
        start_col += task_size;
    }
}

static void get_threads_result(const pthread_t *th, int threads_cnt, 
                               double *s, int *z)
{
    int i;
    for (i = 0; i < threads_cnt; i++) {
        struct thread_res *res;
        pthread_join(th[i], (void*)&res);
        *s += res->s;
        *z += res->z;
        free(res);
    }
}

#endif

static double mavg1d(const double *v, int length, int col, int w)
{
    double s;
    int i, z;
    int start_col, end_col;

#ifdef KZ_PARALLEL
    pthread_t th[max_threads_cnt];
    int threads_cnt, thread_task_size;
    threads_cnt = max_threads_cnt;
#endif

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

#ifdef KZ_PARALLEL
    thread_task_size = (end_col - start_col) / threads_cnt;
    
    if (thread_task_size < min_thread_task_size) {
        threads_cnt = 1;
    } else {
        start_threads(th, threads_cnt-1, start_col, thread_task_size, v);
        start_col += (threads_cnt-1) * thread_task_size;
    }
#endif

    s = 0;
    z = 0;
    for (i = start_col; i < end_col; i++) {
        if (isfinite(v[i])) {
            z++;
            s += v[i];
        }
    }

#ifdef KZ_PARALLEL
    if (threads_cnt > 1)
        get_threads_result(th, threads_cnt-1, &s, &z);
#endif

    if (z == 0) 
        return nan("");
    return s/z;
}

static double *kz1d(const double *x, int length, int window, int iterations)
{
    int i, k;
    int mem_size;
    double *tmp = NULL, *ans = NULL;

    mem_size = length * sizeof(double);
    ans = malloc(mem_size);
    tmp = malloc(mem_size);
    if (!ans || !tmp)
        goto quit;

    memcpy(tmp, x, mem_size);

    for (k = 0; k < iterations; k++) {
        for (i = 0; i < length; i++) {
            ans[i] = mavg1d(tmp, length, i, window);
        }
        memcpy(tmp, ans, mem_size); 
    }

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
        printf("elapsed time %lds, %ldms\n", secs, msecs);
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
