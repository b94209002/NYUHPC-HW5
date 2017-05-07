mpc=mpicc
CC=gcc
FLAGS=-O3
EXECS= MultiGrid_v_cycle MultiGrid_v_cycle_omp jacobi-mpi2D 

all: ${EXECS}

jacobi-mpi2D: jacobi-mpi2D.c
	${mpc} ${FLAGS} $^ -o $@

MultiGrid_v_cycle: MultiGrid_v_cycle.c
	$(CC) $^ $(CFLAGS) -lm -o $@

MultiGrid_v_cycle_omp: MultiGrid_v_cycle_omp.c 
	$(CC) $^ $(CFLAGS) -lm -o $@ -fopenmp

clean:
	rm -f ${EXECS}
