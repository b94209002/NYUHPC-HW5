mpc=mpicc
CC=gcc
FLAGS=-O3
EXECS= MultiGrid_v_cycle MultiGrid_v_cycle_omp MultiGrid_v_cycle_mpi

all: ${EXECS}

MultiGrid_v_cycle: MultiGrid_v_cycle.c
	$(CC) $^ $(CFLAGS) -lm -o $@

MultiGrid_v_cycle_omp: MultiGrid_v_cycle_omp.c 
	$(CC) $^ $(CFLAGS) -lm -o $@ -fopenmp


MultiGrid_v_cycle_mpi: MultiGrid_v_cycle_mpi.c
	$(mpc) $^ $(CFLAGS) -lm -o $@



clean:
	rm -f ${EXECS}
