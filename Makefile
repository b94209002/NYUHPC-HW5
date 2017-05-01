CC=mpicc
FLAGS=-O3
EXECS= jacobi-mpi2D 

all: ${EXECS}

jacobi-mpi2D: jacobi-mpi2D.c
	${CC} ${FLAGS} $^ -o $@

clean:
	rm -f ${EXECS}
