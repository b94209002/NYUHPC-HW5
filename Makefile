mpc=mpicc
CC=gcc
FLAGS=-O3
EXECS= jacobi jacobi-mpi2D 

all: ${EXECS}

jacobi-mpi2D: jacobi-mpi2D.c
	${mpc} ${FLAGS} $^ -o $@

jacobi: jacobi.c
	$(CC) $^ $(CFLAGS) -lm -o $@

clean:
	rm -f ${EXECS}
