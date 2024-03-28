#ifndef KZA_H_SENTRY
#define KZA_H_SENTRY

#ifdef __cplusplus
extern "C" {
#endif

double *kz(const double *data, int dimention, const int *data_size,
           const int *window, int iterations);
void kz_free(double *data);

double *kza(const double *x, int dim, const int *size, const double *y, 
         const int *window, int iterations, double min_window, double tol);
void kza_free(double *x);

#ifdef __cplusplus
}
#endif



#endif
