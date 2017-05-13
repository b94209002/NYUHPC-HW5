/* Solving linear system use Jacobi iteration, with criterion 10^-4 to initial residual */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mpi.h"
#define COMM MPI_COMM_WORLD

struct node get_info(int num);
struct node
       {
        int rank, n, s, e, w;
        MPI_Comm comm;
       };
void restriction(double **fine, double **crse, int N);
void prolongation(double **crse, double **fine, int N);
void vcycle(double ***u, double ***rhs, int lv, int nlv, int *m, double* hsq, double *invhsq, struct node info);
void jacobi(double **u, double **rhs, int N, double hsq, int maxit, double crit, struct node info);
double residual (double **x, double **rhs, int N, double invhsq);
void non_blocking_communication(double **x, int m,struct node info, double** recvbuf, double** sendbuf);
void blocking_communication(double **x, int m,struct node info, double** recvbuf, double** sendbuf);

int main (int argc, char **argv)
{
  int nlevel;
  int i,j,n,num;
  double d, l, r, tmp;// idiagonal, right, left, right-hand-side
  double ***x, ***rhs; //  variable
  double prod,h,a,b,res,rres; // grid size, initial, end, resdiual
  int m0,m1,m2,maxiter;// size, max iteration, iteration index
  double crit; // critera 
  struct node info;

  int *m; 
  double *hsq, *invhsq;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(COMM, &num);
  info = get_info(sqrt(num));

        a = 0.; 
        b = 1.;
	maxiter = 10; 
        nlevel = 2; m0 = 10;
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
        double  time1,time2;
        time1 = MPI_Wtime();

 	res = residual(x[0],rhs[0],m[0],invhsq[0]); crit = 1.e-4*res; n = 0;
        printf("myid = %li, residul = %10e \n", info.rank, res);
   
 	vcycle(x,rhs,0,nlevel, m, hsq,invhsq,info);

        time2 = MPI_Wtime();
        double elapsed = time2 - time1;
        // print final output
        printf("Numer of iteration %li, residual = %10e, Time elapsed is %f secs. \n", n, res, elapsed);

	

        for (n = 0; n < nlevel ; n++) {
		for(i = 0; i < m2; i++) { 
			free(x[n][i]);free(rhs[n][i]);
		}
		free(x[n]);free(rhs[n]);
	}
	free(x);free(rhs);
	free(m);free(hsq);free(invhsq);

  MPI_Finalize();
  return 0;
}

struct node get_info(int num){

struct node info;
MPI_Comm_rank(COMM, &info.rank);

int dim[2], bc[2];
dim[0] = num; dim[1] = num; bc[0] = 0; bc[1] =0;
// define communication information
   MPI_Cart_create(COMM, 2, dim, bc, 0 , &info.comm);
   MPI_Cart_shift(info.comm, 1, -1, &info.e, &info.w );
   MPI_Cart_shift(info.comm, 1, 1, &info.w, &info.e );
   MPI_Cart_shift(info.comm, 0, -1, &info.n, &info.s );
   MPI_Cart_shift(info.comm, 0, 1, &info.s, &info.n );
//

   return info;
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
void vcycle(double ***u, double ***rhs, int lv, int nlv, int *m, double* hsq, double *invhsq, struct node info)
{ double res, crit = 1.e-5;

	res = residual(u[lv],rhs[lv],m[lv],invhsq[lv]); 
        printf("Multigrid level %i, residual = %10e. \n", lv+1, res);

	crit = crit * res;

if ( lv  == nlv - 1) { 
// maximum level
        jacobi(u[lv],rhs[lv],m[lv],hsq[lv],100,crit*1e-5,info);
        res = residual(u[lv],rhs[lv],m[lv],invhsq[lv]);
        printf("Multigrid level %i, residual = %10e. \n", lv+1, res);
} else {

	jacobi(u[lv],rhs[lv],m[lv],hsq[lv],5,crit,info);
	restriction(u[lv], u[lv+1], m[lv+1]);
	vcycle(u,rhs,lv+1,nlv,m, hsq,invhsq,info);
	prolongation(u[lv+1], u[lv], m[lv]);
	jacobi(u[lv],rhs[lv],m[lv],hsq[lv],10000*(lv+1),crit,info);
	res = residual(u[lv],rhs[lv],m[lv],invhsq[lv]);
        printf("Multigrid level %i, residual = %10e. \n", lv+1, res);

}
}

void jacobi(double **u, double **rhs, int N, double hsq, int maxit, double crit, struct node info)
{
  int i, j,n = 0;
  /* Jacobi damping parameter -- plays an important role in MG */
  double omega = 2./3.;
  double res;
  double **unew = calloc(sizeof(double*), N+2);
  double **sbuf , **rbuf;
  sbuf = calloc(sizeof(double*), 4);
  rbuf = calloc(sizeof(double*), 4); 

  for (i = 0; i < N+2; i++) unew[i] = calloc(sizeof(double), N+2);
  for (i = 0; i < 4; i++) sbuf[i] = calloc(sizeof(double), N);
  for (i = 0; i < 4; i++) rbuf[i] = calloc(sizeof(double), N);

  res = residual(u,rhs,N,1.0/hsq);

  while (n < maxit && res > crit) { 
	n++;
	for (j = 1; j < N+1; j++){
		for (i = 1; i < N+1; i++){
			unew[j][i]  = u[j][i] +omega*.25*( rhs[j][i]*hsq + u[j-1][i]+ u[j][i-1] + u[j][i+1]+ u[j+1][i]- 4*u[j][i]);
		}	
	}

        non_blocking_communication(unew, N, info, rbuf, sbuf);
	for (j = 0; j < N+2; j++ ) {
     		memcpy(u[j], unew[j], (N+2)*sizeof(double));
	}

	if ( n%10 == 0) res = residual(u,rhs,N,1.0/hsq);
  }

  if (res < crit){ 
	printf("Jaocbi iteration converges to criterion = %10e. \n", res);
  } else {
	printf("Jaocbi iteration reaches max itereation = %i \n", maxit);
  }

  for (i = 0; i < N+2; i++) free(unew[i]); 
  for (i = 0; i < 4; i++) free(sbuf[i]);
  for (i = 0; i < 4; i++) free(rbuf[i]);

  free(unew);free(sbuf);free(rbuf);
}

double residual (double **x, double **rhs, int N, double invhsq) {

int i,j,ierr;
double tmp,res;

res = 0.0; // boundary condition on (0)

for (j = 1; j < N + 1; ++j) {
	for (i = 1; i < N + 1; ++i) {
		tmp = rhs[j][i] + (x[j-1][i] + x[j][i-1] - 4*x[j][i] + x[j][i+1] + x[j+1][i])*invhsq;
		res = res + tmp*tmp;
	}
}

ierr = MPI_Allreduce(&res , &tmp, 1, MPI_DOUBLE, MPI_SUM, COMM);
res = sqrt(tmp);

return sqrt(res);
}

void non_blocking_communication(double **x, int m,struct node info, double** recvbuf, double** sendbuf) {

MPI_Request reqs[8];
MPI_Status stats[8];
int i,tag = 1;


MPI_Irecv(&recvbuf[0][0], m, MPI_DOUBLE, info.w,tag, info.comm, &reqs[0]);
MPI_Irecv(&recvbuf[1][0], m, MPI_DOUBLE, info.e,tag, info.comm, &reqs[1]);
MPI_Irecv(&recvbuf[2][0], m, MPI_DOUBLE, info.s,tag, info.comm, &reqs[2]);
MPI_Irecv(&recvbuf[3][0], m, MPI_DOUBLE, info.n,tag, info.comm, &reqs[3]);

if (info.w != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        sendbuf[0][i] = x[i+1][1];
}
if (info.e != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        sendbuf[1][i] = x[i+1][m];
}

if (info.s != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        sendbuf[2][i] = x[1][i+1];
}

if (info.n != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        sendbuf[3][i] = x[m][i+1];
}
MPI_Isend(&sendbuf[0][0], m, MPI_DOUBLE, info.w, tag, info.comm, &reqs[4]);
MPI_Isend(&sendbuf[1][0], m, MPI_DOUBLE, info.e, tag, info.comm, &reqs[5]);
MPI_Isend(&sendbuf[2][0], m, MPI_DOUBLE, info.s, tag, info.comm, &reqs[6]);
MPI_Isend(&sendbuf[3][0], m, MPI_DOUBLE, info.n, tag, info.comm, &reqs[7]);

MPI_Waitall(8, reqs, stats);

if (info.w != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        x[i+1][0] = recvbuf[0][i];
}
if (info.e != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        x[i+1][m+1] = recvbuf[1][i];
}

if (info.s != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        x[0][i+1] = recvbuf[2][i];
}

if (info.n != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        x[m+1][i+1] = recvbuf[3][i];
}


}
void blocking_communication(double **x, int m,struct node info, double** recvbuf, double** sendbuf) {
MPI_Status stats[4];
int i,tag = 1;

if (info.w != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        sendbuf[0][i] = x[i+1][1];
}
if (info.e != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        sendbuf[1][i] = x[i+1][m];
}

if (info.s != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        sendbuf[2][i] = x[1][i+1];
}

if (info.n != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        sendbuf[3][i] = x[m][i+1];
}

MPI_Sendrecv(&sendbuf[0][0], m, MPI_DOUBLE, info.w, tag, &recvbuf[1][0], m, MPI_DOUBLE, info.e,tag, info.comm, &stats[0]);
MPI_Sendrecv(&sendbuf[1][0], m, MPI_DOUBLE, info.e, tag, &recvbuf[0][0], m, MPI_DOUBLE, info.w,tag, info.comm, &stats[1]);
MPI_Sendrecv(&sendbuf[2][0], m, MPI_DOUBLE, info.s, tag, &recvbuf[3][0], m, MPI_DOUBLE, info.n,tag, info.comm, &stats[2]);
MPI_Sendrecv(&sendbuf[3][0], m, MPI_DOUBLE, info.n, tag, &recvbuf[2][0], m, MPI_DOUBLE, info.s,tag, info.comm, &stats[3]);

if (info.w != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        x[i+1][0] = recvbuf[0][i];
}
if (info.e != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        x[i+1][m+1] = recvbuf[1][i];
}

if (info.s != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        x[0][i+1] = recvbuf[2][i];
}

if (info.n != MPI_PROC_NULL) {
for (i=0; i<m; i++)
        x[m+1][i+1] = recvbuf[3][i];
}

}

