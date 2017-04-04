CC=mpicc
FLAGS=-O3
EXECS= mpi_solved1

all: ${EXECS}

mpi_solved1: mpi_bug1.c
	${CC} ${FLAGS} $^ -o mpi_solved1

clean:
	rm -f ${EXECS}
