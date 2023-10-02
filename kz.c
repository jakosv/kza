#include "kz.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

static double mavg1d(double *v, int length, int col, int w)
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

void kz1d(double *x, int length, int window, int iterations)
{
    int i, k;
    int arr_size;
    double *tmp;

    arr_size = length * sizeof(double);
    tmp = malloc(arr_size);
    memcpy(tmp, x, arr_size);

    for(k = 0; k < iterations; k++) {
        for(i = 0; i < length; i++) {
            x[i] = mavg1d(tmp, length, i, window);
        }
        memcpy(tmp, x, arr_size); 
    }

    free(tmp);
}

void kz(double *x, int dim, int *size, int *window, int iterations)
{
    int i;

    if (!x || !size || !window) {
        fprintf(stderr, "kz(): Incorrect params\n");
        return;
    }

    for (i = 0; i < dim; i++)
        window[i] = floor(window[i] / 2.0);

    if (dim > 3) {
        fprintf(stderr, "Too many dimensions\n");
        return;
    } else if (dim == 3) {
        /*
        kz3d(x, size, window, iterations); 
        */
    } else if (dim == 2) {
        /*
        kz2d(x, size, window, iterations); 
        */
    } else if (dim == 1) {
        kz1d(x, size[0], window[0], iterations); 
    }
}
