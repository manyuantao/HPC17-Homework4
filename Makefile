CC = mpicc
EXECS = mpi_solved1 mpi_solved2 mpi_solved3 mpi_solved4 mpi_solved5 mpi_solved6 mpi_solved7 jacobi2D-mpi jacobi2D-mpi-nonblocking ssort

all: ${EXECS}

mpi_solved1: mpi_solved1.c
	${CC} mpi_solved1.c -o mpi_solved1

mpi_solved2: mpi_solved2.c
	${CC} mpi_solved2.c -o mpi_solved2

mpi_solved3: mpi_solved3.c
	${CC} mpi_solved3.c -o mpi_solved3

mpi_solved4: mpi_solved4.c
	${CC} mpi_solved4.c -o mpi_solved4

mpi_solved5: mpi_solved5.c
	${CC} mpi_solved5.c -o mpi_solved5

mpi_solved6: mpi_solved6.c
	${CC} mpi_solved6.c -o mpi_solved6

mpi_solved7: mpi_solved7.c
	${CC} mpi_solved7.c -o mpi_solved7

jacobi2D-mpi: jacobi2D-mpi.c
	${CC} jacobi2D-mpi.c -o jacobi-mpi2D -lm

jacobi2D-mpi-nonblocking: jacobi2D-mpi-nonblocking.c
	${CC} jacobi2D-mpi-nonblocking.c -o nonblocking -lm

ssort: ssort.c
	${CC} ssort.c -o ssort

clean:
	rm -f ${EXECS}