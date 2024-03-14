#ifndef KZA_H_SENTRY
#define KZA_H_SENTRY

double *kz(const double *x, int dim, const int *size, const int *window,
        int iterations);
double *kza(const double *x, int dim, const int *size, const double *y, 
         const int *window, int iterations, double min_window, double tol);
void kza_free(double *x);

#endif
