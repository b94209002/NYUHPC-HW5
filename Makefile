CC=mpicc
FLAGS=-O3
EXECS= mpi_solved1 mpi_solved2 mpi_solved3

all: ${EXECS}

mpi_solved1: mpi_bug1.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved2: mpi_bug2.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved3: mpi_bug3.c
	${CC} ${FLAGS} $^ -o $@


clean:
	rm -f ${EXECS}
