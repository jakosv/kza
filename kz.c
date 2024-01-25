#include "kza.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef TIMER
#include "timer.h"
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


static double *kz1d(const double *x, int length, int window, int iterations)
{
    int i, k, mem_size;
    double *ans = NULL;
#ifdef PREFIX_SUM
    double *pref_sum = NULL;
    int *pref_finite_cnt = NULL;
#else
    double *data = NULL;
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

    for (k = 0; k < iterations; k++) {
        for (i = 0; i < length; i++) {
#ifdef PREFIX_SUM
            ans[i] = mavg1d(pref_sum, pref_finite_cnt, length, i, window);
#else
            ans[i] = mavg1d(data, length, i, window); 
#endif
        }

#ifdef PREFIX_SUM
        calc_prefix_sum(ans, length, pref_sum, pref_finite_cnt);
#else
        memcpy(data, ans, mem_size);
#endif
    }

quit:
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
