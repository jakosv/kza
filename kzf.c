#include "kz.h"

#include <math.h>

#define MAX(A, B) ((A) > (B) ? (A) : (B))

double maximum(const double *v, int length) 
{
	double m;
	int i;
	
	for(i = 0, m = v[0]; i < length; i++) {
        if (isfinite(v[i])) 
            m = MAX(v[i], m);
    }
	return m;
}

void differenced(const double *y, double *d, double *dprime, int length, int q)
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
	for(i = 0; i < n-1; i++)
        dprime[i] = d[i+1] - d[i];
	dprime[n-1] = dprime[n-2];
}
