CC=mpicc
FLAGS=-O3
EXECS= mpi_solved1 mpi_solved2 mpi_solved3 mpi_solved4 mpi_solved5 mpi_solved6 mpi_solved7

all: ${EXECS}

mpi_solved1: mpi_bug1.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved2: mpi_bug2.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved3: mpi_bug3.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved4: mpi_bug4.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved5: mpi_bug5.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved6: mpi_bug6.c
	${CC} ${FLAGS} $^ -o $@

mpi_solved7: mpi_bug7.c
	${CC} ${FLAGS} $^ -o $@

clean:
	rm -f ${EXECS}
