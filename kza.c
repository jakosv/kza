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
#include "timer.h"
#endif

#ifdef PREFIX_SUM
static double mavg1d(const double *pref_sum, const int *pref_finite_cnt,
                     int left_bound, int right_bound)
{
    double s;
    int z;
    int start_idx, end_idx;

    /* pref[start_idx] = data[0] + data[1] + ... + data[left_bound] */
    start_idx = left_bound;
    /* pref[end_idx] = data[0] + data[1] + ... + data[right_bound] */
    end_idx = right_bound;

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
static double mavg1d(const double *x, int a, int b)
{
    double s = 0;
    long i, z;

    for (i = a, z = 0; i < b; i++) {
        int is_finite_flag = isfinite(x[i]);
        z += is_finite_flag;
        s += is_finite_flag * x[i];
    }

    return s / (double)z;
}
#endif

static void normalize_window(int *left_win, int *right_win, int n, 
                             int window_center, int min_window_len)
{
    int max_right_bound, new_left, new_right;

    new_left = (*left_win < min_window_len) ? min_window_len : *left_win;
    new_right = (*right_win < min_window_len) ? min_window_len : *right_win;

    /* check bounds */
    max_right_bound = n - window_center - 1;
    new_right = (new_right > max_right_bound) ? max_right_bound : new_right;
    new_left = (new_left > window_center) ? window_center : new_left;

    *left_win = new_left;
    *right_win = new_right;
}

void calc_adaptive_windows(int *left_win, int *right_win,
                           const double *d, const double *dprime, int t,
                           int window, double eps, double max_d)  
{
    int window_adaptive;

    window_adaptive = (int)floor(window * (1 - d[t]/max_d));

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
    int i, mem_size;
    double m;
    double *d = NULL;
    double *dprime = NULL;
    double *ans = NULL;
#ifdef PREFIX_SUM
    double *pref_sum = NULL;
    int *pref_finite_cnt = NULL;
#else
    double *tmp = NULL;
#endif

    mem_size = n * sizeof(double);

    d = malloc(mem_size);
    dprime = malloc(mem_size);
    if (!d || !dprime)
        goto quit;

    memset(dprime, 0, mem_size);
    differenced(y, d, dprime, n, window);

    m = maximum(d, n);

#ifdef PREFIX_SUM
    pref_sum = malloc((n+1) * sizeof(double));
    pref_finite_cnt = malloc((n+1) * sizeof(int));
    if (!pref_sum || !pref_finite_cnt)
        goto quit;

    calc_prefix_sum(v, n, pref_sum, pref_finite_cnt);
#else
    tmp = malloc(mem_size);
    if (!tmp)
        goto quit;

    memcpy(tmp, v, mem_size);
#endif

    ans = malloc(mem_size);
    if (!ans)
        goto quit;

    for (i = 0; i < iterations; i++) {
        int t;
        for (t = 0; t < n; t++) {
            int left_win, right_win, left_bound, right_bound;

            calc_adaptive_windows(&left_win, &right_win, d, dprime, t,
                                  window, tolerance, m);
            get_window_bounds(&left_bound, &right_bound,
                              left_win, right_win, t, n, min_window_len);

#ifdef PREFIX_SUM
            ans[t] = mavg1d(pref_sum, pref_finite_cnt, left_bound, 
                            right_bound+1);
#else
            ans[t] = mavg1d(tmp, left_bound, right_bound+1); 
#endif
        }
#ifdef PREFIX_SUM
        calc_prefix_sum(ans, n, pref_sum, pref_finite_cnt);
#else
        memcpy(tmp, ans, mem_size);
#endif
    }

quit:
#ifdef PREFIX_SUM
    free(pref_sum);
    free(pref_finite_cnt);
#else
    free(tmp);
#endif

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
    struct timeval t_start;
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
    timestamp(&t_start);
#endif
    kza_ans = kza1d(x, size[0], kz_ans, window[0], iterations, min_window,
                    tolerance);
#ifdef TIMER
    print_timer(&t_start, "kza():");
#endif

    free(kz_ans);

    return kza_ans;
}

/* it is just a wrapper */
void kza_free(double *x)
{
    free(x);
}
