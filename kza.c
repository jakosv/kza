#include "kz.h"

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

static double adaptive(double d, double m)
{
    return 1 - d/m;
}

static double mavg1d(const double *x, int a, int b)
{
    double s=0.00;
    long i, z;

    for(i = a, z = 0; i < b; i++)
    {
        if (isfinite(x[i])) {
            z++;
            s += x[i];
        }
    }
    if (z == 0) 
        return nan("");
    return s/z;
}

void kza1d(double *v, int n, const double *y, int window, int iterations,
           int min_windown_len, double tolerance)
{
    int i, mem_size;
    double *d, *dprime;
    double m;
    long t, qh, qt;
    double *tmp, *ans;
    double eps;

    eps = tolerance;

    mem_size = n * sizeof(double);
    d = malloc(mem_size);
    dprime = malloc(mem_size);

    differenced(y, d, dprime, n, window);

    m = maximum(d, n);

    tmp = malloc(mem_size);
    memcpy(tmp, v, mem_size);

    ans = malloc(mem_size);

    for(i = 0; i < iterations; i++) {
        for (t = 0; t < n; t++) {
            if (fabs(dprime[t]) < eps) { /* dprime[t] = 0 */
                qh = (int)floor(window * adaptive(d[t], m));
                qt = (int)floor(window * adaptive(d[t], m));
            } else if (dprime[t] < 0) {
                qh = window;
                qt = (int)floor(window * adaptive(d[t], m));
            } else {
                qh = (int)floor(window * adaptive(d[t], m));
                qt = window;
            }
            qt = ((qt) < min_windown_len) ? min_windown_len : qt;
            qh = ((qh) < min_windown_len) ? min_windown_len : qh;

            /* check bounds */
            qh = (qh > n-t-1) ? n-t-1 : qh; /* head past end of series */
            qt = (qt > t) ? t : qt;  	        		
            ans[t] = mavg1d(tmp, t-qt, t+qh+1); 
        }
        memcpy(tmp, ans, mem_size);
    }

    memcpy(v, ans, mem_size);

    free(tmp);
    free(d);
    free(dprime);
}

void kza(double *x, int dim, const int *size, const double *y, 
         const int *window, int iterations, double min_window, double tol)
{
    double *kz_ans;
    int mem_size;

    mem_size = size[0] * sizeof(double);
    kz_ans = malloc(mem_size);
    if (!y) {
        memcpy(kz_ans, x, mem_size);
        kz(kz_ans, dim, size, window, iterations);
    } else {
        memcpy(kz_ans, y, mem_size);
    }

    kza1d(x, size[0], kz_ans, window[0], iterations, min_window, tol);

    free(kz_ans);
}

