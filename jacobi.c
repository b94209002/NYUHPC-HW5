/* Solving linear system use Jacobi iteration, with criterion 10^-4 to initial residual */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void restriction(double **fine, double **crse, int N);
void prolongation(double **crse, double **fine, int N);
void  vcycle(double ***u, double ***rhs, int lv, int nlv, int *m, double* hsq, double *invhsq );
void jacobi(double **u, double **rhs, int N, double hsq, int ssteps);
double residual (double **x, double **rhs, int N, double invhsq);

int main (int argc, char **argv)
{
  int nlevel;
  int i,j,n;
  double d, l, r, tmp;// idiagonal, right, left, right-hand-side
  double ***x, ***rhs; //  variable
  double prod,h,a,b,res,rres; // grid size, initial, end, resdiual
  int m0,m1,m2,maxiter;// size, max iteration, iteration index
  double crit; // critera 

  int *m; 
  double *hsq, *invhsq;

        a = 0.; 
        b = 1.;
	maxiter = 10; 
        nlevel = 10; m0 =2;
	m = malloc(nlevel * sizeof(int));
 	hsq = malloc(nlevel * sizeof(double));
	invhsq = malloc(nlevel * sizeof(double));
	x = malloc(nlevel * sizeof(double **));
	rhs = malloc(nlevel * sizeof(double **));
        for (i = 0; i<nlevel; i++){
		m[i] = m0*pow(2,nlevel-i) - 1;
                m2 = m[i]+2;
		//printf("number of pt =  %i, at level %i \n", m2,i);
		h = (b-a)/((double)m2 - 1.0) ;
		hsq[i] = h*h; 
		invhsq[i] = 1.0/h/h;
        	x[i] = malloc(m2 * sizeof(double *));
		rhs[i] = malloc(m2 * sizeof(double *));
		for (j =0; j<m2; j++ ){
			x[i][j] = malloc(m2 * sizeof(double));
			rhs[i][j] = malloc(m2 * sizeof(double));
		}

	} 
        // initial matrix (include bc)
        for (n = 0; n < nlevel; n++) {
  		for (j = 0; j < m[n]+2; ++j) {
			for (i = 0; i < m[n]+2; ++i) {
				x[n][j][i] = 0.0;
				rhs[n][j][i] = 1.0;
			}	
		}
	}
//  	timestamp_type time1, time2;
//  	get_timestamp(&time1);
	int myid = 0;	
 	res = residual(x[0],rhs[0],m[0],invhsq[0]); crit = 1.e-4*res; n = 0;
        printf("myid = %li, residul = %10e \n", myid, res);
      
  	vcycle(x,rhs,0,nlevel, m, hsq,invhsq);

//	get_timestamp(&time2);
//  	double elapsed = timestamp_diff_in_seconds(time1,time2);
	
	// print final output
//	printf("Numer of iteration %li, residual = %10e, Time elapsed is %f secs. \n", n, sqrt(res), elapsed);

        for (n = 0; n < nlevel ; n++) {
		for(i = 0; i < m2; i++) { 
			free(x[n][i]);free(rhs[n][i]);
		}
		free(x[n]);free(rhs[n]);
	}
	free(x);free(rhs);
	free(m);free(hsq);free(invhsq);


  return 0;
}

void restriction(double **fine, double **crse, int N) { // coarse N
  int i,j,it,jt;
  double a = 4.0/16.0, b = 2.0 / 16.0, c = 1.0/16.0; 

  for (j = 1; j < N+1 ; ++j ){
	jt = 2*j;
	for (i = 1; i < N + 1 ; ++i){
		it = 2*i;
		crse[j][i] = a*fine[jt][it]+b*(fine[jt][it-1]+fine[jt][it+1]+fine[jt-1][it]+ fine[jt+1][it]);
		crse[j][i] = crse[j][i] + c*(fine[jt-1][it-1]+fine[jt-1][it+1]+fine[jt+1][it-1]+ fine[jt+1][it+1]);
	}
  }

}

void prolongation(double **crse, double **fine, int N) { // fine N
  int i,j,it,jt;
  double c1 = .5 ,c2 = .25;
  for (j = 1; j < N +1; ++j ){
        jt = j/2;
        for (i = 1; i < N +1; ++i){                
                it = i/2;
                    if (i%2==0 && j%2==0) // even x, even y
                        fine[j][i] = crse[jt][it];
                    else if (i%2==0 && j%2==1) // even x, odd y
                        fine[j][i] = c1 * (crse[jt][it] + crse[jt+1][it]) ;
                    else if (i%2==1 && j%2==0) // odd x, even y
                        fine[j][i] = c1 * (crse[jt][it] + crse[jt][it+1]) ;
                    else if (i%2==1 && j%2==1) // odd x, odd y
                        fine[j][i] = c2 * (crse[jt][it] + crse[jt][it+1] + crse[jt+1][it] + crse[jt+1][it+1]) ;
        }
  }

}
void  vcycle(double ***u, double ***rhs, int lv, int nlv, int *m, double* hsq, double *invhsq )
{ double res;

if ( lv  == nlv - 1) { 
// maximum level
        res = residual(u[lv],rhs[lv],m[lv],invhsq[lv]);
        printf("Multigrid level %i, sec residual = %10e. \n", lv+1, res);
        jacobi(u[lv],rhs[lv],m[lv],hsq[lv],100);
        res = residual(u[lv],rhs[lv],m[lv],invhsq[lv]);
        printf("Multigrid level %i, residual = %10e. \n", lv+1, res);
} else {

	res = residual(u[lv],rhs[lv],m[lv],invhsq[lv]);
        printf("Multigrid level %i, sec residual = %10e. \n", lv+1, res);
	jacobi(u[lv],rhs[lv],m[lv],hsq[lv],3);
        res = residual(u[lv],rhs[lv],m[lv],invhsq[lv]);
        printf("Multigrid level %i, before residual = %10e. \n", lv+1, res);
	restriction(u[lv], u[lv+1], m[lv+1]);
	vcycle(u,rhs,lv+1,nlv,m, hsq,invhsq);
	prolongation(u[lv+1], u[lv], m[lv]);
        res = residual(u[lv],rhs[lv],m[lv],invhsq[lv]);
        printf("Multigrid level %i, sec residual = %10e. \n", lv+1, res);
	jacobi(u[lv],rhs[lv],m[lv],hsq[lv],10*(nlv-lv));
	res = residual(u[lv],rhs[lv],m[lv],invhsq[lv]);
        printf("Multigrid level %i, residual = %10e. \n", lv+1, res);

}
}


void jacobi(double **u, double **rhs, int N, double hsq, int ssteps)
{
  int i, j,n;
  /* Jacobi damping parameter -- plays an important role in MG */
  double omega = 2./3.;
  double **unew = calloc(sizeof(double*), N+2);
  for (i = 1; i < N+1; i++) unew[i] = calloc(sizeof(double), N+2);
  		
  for (n = 0; n < ssteps; ++n) {
	for (j = 1; j < N+1; j++){
		for (i = 1; i < N+1; i++){
			unew[j][i]  = .25*( rhs[j][i]*hsq + u[j-1][i]+ u[j][i-1] + u[j][i+1]+ u[j+1][i]);
		}	
	}
	for (j = 1; j < N+1; j++ ) {
     		memcpy(u[j], unew[j], (N+2)*sizeof(double));
	}
  }
  for (i = 1; i < N+1; i++) {
	free(unew[i]); 
  }
  free (unew);
}

double residual (double **x, double **rhs, int N, double invhsq) {

int i,j;
double tmp,res;

res = 0.0; // boundary condition on (0)

for (j = 1; j < N + 1; ++j) {
	for (i = 1; i < N + 1; ++i) {
		tmp = rhs[j][i] + (x[j-1][i] + x[j][i-1] - 4*x[j][i] + x[j][i+1] + x[j+1][i])*invhsq;
		res = res + tmp*tmp;
	}
}

return sqrt(res);
}
