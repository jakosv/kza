#ifndef KZA_H_SENTRY
#define KZA_H_SENTRY

void kz(double *x, int dim, const int *size, const int *window,
        int iterations);
void kza(double *x, int dim, const int *size, const double *y, 
         const int *window, int iterations, double min_window, double tol);
double maximum(const double *v, int length);
void differenced(const double *y, double *d, double *dprime, int length, int q);

#endif
