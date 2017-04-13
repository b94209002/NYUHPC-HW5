/* Parallel sample sort
 */
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>
#include <stdlib.h>


static int compare(const void *a, const void *b)
{
  int *da = (int *)a;
  int *db = (int *)b;

  if (*da > *db)
    return 1;
  else if (*da < *db)
    return -1;
  else
    return 0;
}

int main( int argc, char *argv[])
{
  int rank,num;
  int i, N, id, m ,n,n0;
  int *vec, *vec3, *vec2, *sendcount, *sdispls, *recvcount, *rdispls;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Number of random numbers per processor (this should be increased
   * for actual tests or could be passed in through the command line */

  for (n = 0; n < 5;n++){
  N = 1000;
  for (i = 0; i<n;i++) N = N*10;
  n0=N;
  m = N/10;
  MPI_Barrier( MPI_COMM_WORLD );
  double  time1,time2;
  time1 = MPI_Wtime();

  vec = calloc(N, sizeof(int));
  /* seed random number generator differently on every core */
  srand((unsigned int) (rank + 393919));

  /* fill vector with random integers */
  for (i = 0; i < N; ++i) {
    vec[i] = rand();
  }
//  printf("rank: %d, first entry: %d\n", rank, vec[0]);

  /* sort locally */
  qsort(vec, N, sizeof(int), compare);

  /* randomly sample s entries from vector or select local splitters,
   * i.e., every N/P-th entry of the sorted vector */
    vec2 = calloc(m, sizeof(int));	

  for (i = 0; i < m; ++i) {
    id = rand()%N ;
    vec2[i] = vec[id];
  }

  /* every processor communicates the selected entries
   * to the root processor; use for instance an MPI_Gather */

  if (0 == rank) vec3 = calloc(m*num, sizeof(int));

  MPI_Gather( &vec2[0],m, MPI_INT, &vec3[0], m, MPI_INT, 0, MPI_COMM_WORLD);
  
  free(vec2); vec2 = calloc(num, sizeof(int));

  /* root processor does a sort, determinates splitters that
   * split the data into P buckets of approximately the same size */

  if (rank ==0 ){
  qsort( vec3, m*num,sizeof(int), compare);

  for (i = 0; i < num; ++i)
  	vec2[i] = vec3[m*i];
  free(vec3);
  }
  vec2[num-1] = 2147483647; //INT_MAX in case to over count the sdispls.

  /* root process broadcasts splitters */
  MPI_Bcast( vec2, num, MPI_INT, 0, MPI_COMM_WORLD );

  /* every processor uses the obtained splitters to decide
   * which integers need to be sent to which other processor (local bins) */

  sendcount = calloc(num, sizeof(int));
  sdispls = calloc(num, sizeof(int));
  id = 0; sdispls[0] = 0;

  for (i=0; i < N; i++){
	if (vec[i] >= vec2[id] ){
                id++;
	 	sdispls[id] = i; 
	}
  }

  for (i=0; i < num; i++)
      sendcount[i] = sdispls[i+1] - sdispls[i] ;   
  sendcount[num-1] = N - sdispls[num-1] ;
  
  if (63 == rank) {
    for (i=0;i<num;i++)
    printf("test for sendcount %i, ndispls = %i  \n",sendcount[i], sdispls[i]);
  }



  /* send and receive: either you use MPI_AlltoallV, or
   * (and that might be easier), use an MPI_Alltoall to share
   * with every processor how many integers it should expect,
   * and then use MPI_Send and MPI_Recv to exchange the data */
  recvcount = calloc(num, sizeof(int));

  MPI_Alltoall(&sendcount[0], 1, MPI_INT, &recvcount[0], 1,MPI_INT, MPI_COMM_WORLD);

  N = 0;
  for (i = 0; i < num; i++ ){
	N += recvcount[i];
  }

  rdispls = calloc(num, sizeof(int));
  rdispls[0] = 0;
  for (i = 0; i < num-1; i++ ){
        rdispls[i+1] = rdispls[i] + recvcount[i];
  }

  vec3 = calloc(N, sizeof(int));

  MPI_Alltoallv(&vec[0], sendcount, sdispls, MPI_INT, &vec3[0], recvcount, rdispls, MPI_INT, MPI_COMM_WORLD);

  /* do a local sort */

  qsort(vec3, N, sizeof(int), compare);

  /* every processor writes its result to a file */
  {
    FILE* fd = NULL;
    char filename[256];
    snprintf(filename, 256, "output%02d.txt", rank);
    fd = fopen(filename,"w+");

    if(NULL == fd)
    {
      printf("Error opening file \n");
      return 1;
    }

    for(i = 0; i < N; ++i)
      fprintf(fd, "  %i \n", vec3[i]);

    fclose(fd);
  }
        time2 = MPI_Wtime();
        double elapsed = time2 - time1;
        // print final output
        if (0 == rank) printf("For N = %i, Time elapsed is %f secs. \n", n0, elapsed);

  

  free(vec); free(vec2); free(vec3) ; 
  free(sendcount);free(recvcount); 
  free(sdispls); free(rdispls);

  }

  MPI_Finalize();
  return 0;
}
