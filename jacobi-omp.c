/* Solving linear system use Jacobi iteration, with criterion 10^-4 to initial residual */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "util.h"
double residual(int m, double h2, double **x, double rhs);

int main (int argc, char **argv)
{
  long  i,j;
  double d, l, r, rhs, tmp;// idiagonal, right, left, right-hand-side
  double  **x, **x0; //  variable
  double prod,h,h2,a,b,res,rres; // grid size, initial, end, resdiual
  long m, m2,m1, maxiter, n = 0;// size, max iteration, iteration index
  double crit; // critera 

        a = 0;
        b = 1;
//	m = atol(argv[1]);
        m = 100;
	m2 = m+2 ;
	m1 = m+1 ;
        maxiter = 100;
//	maxiter = atol(argv[2]);

	x = malloc(m2 * sizeof(double *));
        x0 = malloc(m2 * sizeof(double *));
	for(i = 0; i < m2; i++)
	{
		x[i] = malloc(m2 * sizeof(double));
		x0[i] = malloc(m2 * sizeof(double));		
	}

	h = (b-a)/((double)m+1.0) ;h2 = h*h; rhs = 1;  
        // initial matrix
  	for (j = 0; j < m2; ++j) {
		for (i = 0; i < m2; ++i) {
			x[j][i] = 0.0;
			x0[j][i] = 0.0;
		}
	}
  	timestamp_type time1, time2;
  	get_timestamp(&time1);
	int myid = 0;	
	#pragma omp parallel shared(x,x0,m,m1,m2,h,h2,rhs,crit,maxiter) private(myid,n,i,j,tmp,rres)
	{
#ifdef _OPENMP
  	myid = omp_get_thread_num();
#endif
        res = 0.0; // boundary condition on (0)
	#pragma omp barrier
	#pragma omp for schedule(dynamic,10) reduction(+:res)
		for (j = 1; j < m1; ++j) {
        		for (i = 1; i < m1; ++i) {
                	        tmp = rhs + (x[j-1][i] + x[j][i-1] - 4*x[j][i] + x[j][i+1] + x[j+1][i])/h2;
        		        res = res + tmp*tmp;
        		}
		}
	rres = sqrt(res); crit = 1.e-4*rres; n = 0;
        printf("myid = %li, residul = %10e \n", myid, rres);
	while (n < maxiter && rres > crit ) {
	n++ ; 
	#pragma omp barrier
	res = 0.0;

        #pragma omp for schedule(dynamic,10)
		for (j = 1; j < m1; ++j) {
			for (i = 1; i < m1; ++i) {
				x[j][i] = .25*( rhs*h*h + x0[j-1][i]+ x0[j][i-1] + x0[j][i+1]+ x0[j+1][i]);
			}	
		}

        #pragma omp for schedule(dynamic,10) reduction(+:res)
	       for (j = 1; j < m1; ++j) {
                        for (i = 1; i < m1; ++i) {
                                tmp = rhs + (x[j-1][i] + x[j][i-1] - 4*x[j][i] + x[j][i+1] + x[j+1][i])/h2;
                                res = res + tmp*tmp;
                        }
                }
	rres = sqrt(res);
        #pragma omp for schedule(dynamic,10) 
                for (j = 1; j < m1; ++j) {
                        for (i = 1; i < m1; ++i) {
				x0[j][i] = x[j][i];// also copy the array
                        }
                }
		if (myid == 0 ) printf("niter = %li, residul = %10e \n", n, rres);		
	}

	}

	get_timestamp(&time2);
  	double elapsed = timestamp_diff_in_seconds(time1,time2);
	
	// print final output
	printf("Numer of iteration %li, residual = %10e, Time elapsed is %f secs. \n", n, sqrt(res), elapsed);

        for(i = 0; i < m2; i++)
        {
                free(x[i]);
                free(x0[i]);
        }

	free(x);
	free(x0);

  return 0;
}

double residual (int m, double h2, double **x, double rhs) {

int i,j;
double tmp,res;

res = 0.0; // boundary condition on (0)

for (j = 1; j < m; ++j) {
	for (i = 1; i < m; ++i) {
			tmp = rhs + (x[j-1][i] + x[j][i-1] - 4*x[j][i] + x[j][i+1] + x[j+1][i])/h2;
		res = res + tmp*tmp;
	}
}

res = sqrt(res);

return res;
}
