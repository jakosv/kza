#ifndef KZA_H_SENTRY
#define KZA_H_SENTRY

#ifdef __cplusplus
extern "C" {
#endif

double *kz(const double *data, int dimention, const int *data_size,
           const int *window, int iterations);

#ifdef __cplusplus
}
#endif

double *kza(const double *x, int dim, const int *size, const double *y, 
         const int *window, int iterations, double min_window, double tol);
void kza_free(double *x);


#endif
