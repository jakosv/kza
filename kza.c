#include "kza.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MAX(A, B) ((A) > (B) ? (A) : (B))

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

static double adaptive(double d, double m)
{
    return 1 - d/m;
}

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

static void normalize_window(int *left_window, int *right_window, int n,
                             int time, int min_window_len)
{
    int right_bound;

    if (*left_window < min_window_len)
        *left_window = min_window_len;
    if (*right_window < min_window_len)
        *right_window = min_window_len;

    /* check bounds */
    right_bound = n - time - 1;
    if (*right_window > right_bound)
        *right_window = right_bound; /* head past end of series */
    if (*left_window > time)
        *left_window = time;
}

static double *kza1d(const double *v, int n, const double *y, int window,
                     int iterations, int min_window_len, double tolerance)
{
    int i, mem_size;
    double m, eps;
    double *d = NULL;
    double *dprime = NULL;
    double *tmp = NULL;
    double *ans = NULL;

    eps = tolerance;

    mem_size = n * sizeof(double);

    d = malloc(mem_size);
    dprime = malloc(mem_size);
    if (!d || !dprime)
        goto quit;

    differenced(y, d, dprime, n, window);

    m = maximum(d, n);

    tmp = malloc(mem_size);
    if (!tmp)
        goto quit;

    memcpy(tmp, v, mem_size);

    ans = malloc(mem_size);
    if (!ans)
        goto quit;

    for (i = 0; i < iterations; i++) {
        int t;
        for (t = 0; t < n; t++) {
            int window_adaptive, left_window, right_window;
            int left_bound, right_bound;
            window_adaptive = (int)floor(window * adaptive(d[t], m));
            if (fabs(dprime[t]) < eps) { /* dprime[t] = 0 */
                right_window = window_adaptive;
                left_window = window_adaptive;
            } else if (dprime[t] < 0) {
                right_window = window;
                left_window = window_adaptive;
            } else {
                right_window = window_adaptive;
                left_window = window;
            }

            normalize_window(&left_window, &right_window, n, t, 
                             min_window_len);
            left_bound = t - left_window;
            right_bound = t + right_window + 1;

            ans[t] = mavg1d(tmp, left_bound, right_bound); 
        }
        memcpy(tmp, ans, mem_size);
    }

quit:
    free(tmp);
    free(d);
    free(dprime);

    return ans;
}

double *kza(const double *x, int dim, const int *size, const double *y, 
            const int *window, int iterations, double min_window, double tol)
{
    double *kz_ans = NULL;
    double *kza_ans = NULL;
    int mem_size;

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

    kza_ans = kza1d(x, size[0], kz_ans, window[0], iterations, min_window,
                    tol);

    free(kz_ans);

    return kza_ans;
}

void kza_free(double *x)
{
    free(x);
}
