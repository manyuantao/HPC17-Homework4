# HPC17-Homework4
Homework4, MPI blocking and non-blocking, parallel sample sort.

## Compile
Type `make` in a terminal.

## Find MPI bugs

Enter `mpirun -np 4 ./mpi_solved{1,2,...,7}`.

## Run MPI-parallel 2D Jacobi smoother
Enter `mpirun -np 16 ./jacobi-mpi2D N max_iters` and `mpirun -np 16 ./nonblocking N max_iters`,

where `N` is the total number of discretization points, `max_iters` is the maximum number of iterations.

## Run parallel sample sort
Enter `mpirun -np 10 ./ssort N`,

where `N` is the number of random numbers to be sorted per processor.

`N` is set by default to be `100` if not input.
