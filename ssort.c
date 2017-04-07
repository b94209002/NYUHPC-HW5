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
  int i, N, id, m  ;
  int *vec, *vec2, *vec3, *vec4;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &num);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  /* Number of random numbers per processor (this should be increased
   * for actual tests or could be passed in through the command line */
  N = 100; m = N/10;

  vec = calloc(N, sizeof(int));
  /* seed random number generator differently on every core */
  srand((unsigned int) (rank + 393919));

  /* fill vector with random integers */
  for (i = 0; i < N; ++i) {
    vec[i] = rand();
  }
  printf("rank: %d, first entry: %d\n", rank, vec[0]);

  /* sort locally */
  qsort(vec, N, sizeof(int), compare);

  /* randomly sample s entries from vector or select local splitters,
   * i.e., every N/P-th entry of the sorted vector */

  vec2 = calloc(num*m, sizeof(int));
  for (i = 0; i < m; ++i) {
    id = rand()%num ;
    vec2[i] = vec[id];
  }
 



  /* every processor communicates the selected entries
   * to the root processor; use for instance an MPI_Gather */

  MPI_Gather( &vec2[0],m, MPI_INT, vec2, m, MPI_INT, 0, MPI_COMM_WORLD);

  /* root processor does a sort, determinates splitters that
   * split the data into P buckets of approximately the same size */

  if (rank ==0 )
  qsort( vec2, m*num,sizeof(int), compare);



  /* root process broadcasts splitters */
  MPI_Bcast( vec2, m*num, MPI_INT, 0, MPI_COMM_WORLD );

  /* every processor uses the obtained splitters to decide
   * which integers need to be sent to which other processor (local bins) */

  vec3 = calloc(num, sizeof(int));
  id = 0;
  for (i=0; i < N; i++){
	if (vec[i] >= vec2[id] ){
	 	vec3[id] = i; 
		id++ ;
	}
  }
  vec2[0] = vec3[0];
  for (i=1; i < num; i++)
  	vec2[i] = vec3[i] - vec3[i-1];


  /* send and receive: either you use MPI_AlltoallV, or
   * (and that might be easier), use an MPI_Alltoall to share
   * with every processor how many integers it should expect,
   * and then use MPI_Send and MPI_Recv to exchange the data */

  MPI_Alltoall(&vec3, 1, MPI_INT, vec3,1,MPI_INT, MPI_COMM_WORLD);
  N = 0;
  for (i = 0; i < N; i++ )
	N += vec[2];
  

  MPI_Alltoallv(&vec, vec3, vec2, MPI_INT, vec4, vec3, vec2, MPI_INT, MPI_COMM_WORLD);

  /* do a local sort */

  qsort(vec4, N, sizeof(int), compare);

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
      fprintf(fd, "  %f\n", vec4[i]);

    fclose(fd);
  }

  free(vec); free(vec2); free(vec3) ; free(vec4);
  MPI_Finalize();
  return 0;
}
