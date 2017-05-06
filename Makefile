mpc=mpicc
CC=gcc
FLAGS=-O3
EXECS= MultiGrid_v_cycle jacobi-mpi2D 

all: ${EXECS}

jacobi-mpi2D: jacobi-mpi2D.c
	${mpc} ${FLAGS} $^ -o $@

MultiGrid_v_cycle: MultiGrid_v_cycle.c
	$(CC) $^ $(CFLAGS) -lm -o $@

clean:
	rm -f ${EXECS}
