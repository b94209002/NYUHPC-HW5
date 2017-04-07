/* Solving linear system use Jacobi iteration, with criterion 10^-4 to initial residual */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mpi.h"
#define COMM MPI_COMM_WORLD
double residual(int m, double h2, double **x, double rhs);
//void non_blocking_communication(double **x, int m,struct node info, double** recvbuf, double** sendbuf);
struct node get_info(int num);
struct node
       {
        int rank, n, s, e, w;
        MPI_Comm comm;
       };
void non_blocking_communication(double **x, int m,struct node info, double** recvbuf, double** sendbuf);
void blocking_communication(double **x, int m,struct node info, double** recvbuf, double** sendbuf);

int main (int argc, char **argv)
{
  long  i,j,k;
  int i1,i2,num;
  double d, l, r, rhs, tmp;// idiagonal, right, left, right-hand-side
  double ***x; //  variable
  double ***buf; // send and recv buf
  double prod,h,h2,a,b,res; // grid size, initial, end, resdiual
  long m, m2,m1, maxiter, n = 0;// size, max iteration, iteration index
  double crit; // critera 
  struct node info;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(COMM, &num);
  info = get_info(sqrt(num));

        a = 0;
        b = 1;
 	m = atol(argv[1]);
	m2 = m+2 ;
	m1 = m+1 ;
 	maxiter = atol(argv[2]);

	x = malloc(2 * sizeof(double **));
	x[0] = malloc(m2 * sizeof(double *));
        x[1] = malloc(m2 * sizeof(double *));
	for(k = 0; k < 2 ; k++){
		for(i = 0; i < m2; i++){
			x[k][i] = malloc(m2 * sizeof(double));
		}
	}

	// initialize send and recv buffer. 
        buf = malloc(2 * sizeof(double**));
        buf[0] = malloc(4 * sizeof(double *));
        buf[1] = malloc(4 * sizeof(double *));
	for(k = 0; k < 2; k++){
        	for(i = 0; i < 4; i++)
        		buf[k][i] = malloc(m * sizeof(double));
	}

	h = (b-a)/((double)m+1.0) ;h2 = h*h; rhs = 1;  
        // initial matrix
  	for (j = 0; j < m2; ++j) {
		for (i = 0; i < m2; ++i) {
			x[0][j][i] = 0.0;
		}
	}
	i2= 1; i1= 0;

        double  time1,time2;
        time1 = MPI_Wtime();

	res = residual(m,h*h,x[i1],rhs);

	crit = 1.e-4*res; n = 0;
        if (info.rank == 0) printf("myid = %i, residul = %10e \n", info.rank, res);

	while (n < maxiter && res > crit ) {
	n++ ; i = i2; i2 = i1; i1 = i; // use the different data

	for (j = 1; j < m1; ++j) {
		for (i = 1; i < m1; ++i) {
			x[i1][j][i] = .25*( rhs*h*h + x[i2][j-1][i]+ x[i2][j][i-1] + x[i2][j][i+1]+ x[i2][j+1][i]);
		}	
	}

	non_blocking_communication(x[i1], m, info, buf[0], buf[1]); 

	res = residual(m,h*h,x[i1],rhs);
	
	if (info.rank == 0 ) printf("niter = %li, residul = %10e \n", n, res);		
	}

        time2 = MPI_Wtime();
        double elapsed = time2 - time1;
	// print final output
	printf("Numer of iteration %li, residual = %10e, Time elapsed is %f secs. \n", n, res, elapsed);




        for(k = 0; k < 2 ; k++){
                for(i = 0; i < 4; i++)
                free(buf[k][i]);
        }
	free(buf[0]);free(buf[1]);
        free(buf);
	for(k = 0; k < 2 ; k++){
        	for(i = 0; i < m2; i++)
                free(x[k][i]);
	}
	free(x[0]);free(x[1]);
	free(x);

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

return info;
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

double residual (int m, double h2, double **x, double rhs) {

int i,j,ierr;
double tmp,res;

res = 0.0; // boundary condition on (0)

for (j = 1; j < m+1; ++j) {
	for (i = 1; i < m+1; ++i) {
			tmp = rhs + (x[j-1][i] + x[j][i-1] - 4*x[j][i] + x[j][i+1] + x[j+1][i])/h2;
		res = res + tmp*tmp;
	}
}

ierr = MPI_Allreduce(&res , &tmp, 1, MPI_DOUBLE, MPI_SUM, COMM);
res = sqrt(tmp);

return res;
}



