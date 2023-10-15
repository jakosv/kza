#include "kza.h"

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

    for(i = start_col, z = 0, s = 0.0; i < end_col; i++) {
        if (isfinite(v[i])) {
            z++;
            s += v[i];
        }
    }

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

    for(k = 0; k < iterations; k++) {
        for(i = 0; i < length; i++) {
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
    int *m;
    double *ans = NULL;

    if (!x || !size || !window) {
        fprintf(stderr, "kz(): Incorrect params\n");
        return NULL;
    }

    m = malloc(dim * sizeof(int));
    if (!m)
        goto quit;

    for (i = 0; i < dim; i++)
        m[i] = floor(window[i] / 2.0);

    switch (dim) {
    case 1:
        ans = malloc(size[0] * sizeof(double));
        if (!ans)
            goto quit;
        ans = kz1d(x, size[0], m[0], iterations); 
        break;
    default:
        fprintf(stderr, "kza: Too many dimensions\n");
    }

quit:
    if (m)
        free(m);

    return ans;
}
