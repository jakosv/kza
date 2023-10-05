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
    double *tmp, *ans;

    mem_size = length * sizeof(double);
    ans = malloc(mem_size);
    tmp = malloc(mem_size);
    memcpy(tmp, x, mem_size);

    for(k = 0; k < iterations; k++) {
        for(i = 0; i < length; i++) {
            ans[i] = mavg1d(tmp, length, i, window);
        }
        memcpy(tmp, ans, mem_size); 
    }

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
    for (i = 0; i < dim; i++)
        m[i] = floor(window[i] / 2.0);

    if (dim > 3) {
        fprintf(stderr, "kza: Too many dimensions\n");
    } else if (dim == 3) {
        /*
        kz3d(x, size, m, iterations); 
        */
        fprintf(stderr, "kza: Not yet implemented\n");
    } else if (dim == 2) {
        /*
        kz2d(x, size, m, iterations); 
        */
        fprintf(stderr, "kza: Not yet implemented\n");
    } else if (dim == 1) {
        ans = malloc(size[0] * sizeof(double));
        ans = kz1d(x, size[0], m[0], iterations); 
    }

    free(m);

    return ans;
}
