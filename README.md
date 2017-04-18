# HPC17-Homework4
Homework4, MPI blocking and non-blocking, parallel sample sort.

## Compile
Type `make` in a terminal.

## Find MPI bugs

Enter `mpirun -np 4 ./mpi_solved{1,2,...,7}`.

## Run MPI-parallel 2D Jacobi smoother
Enter `mpirun -np 16 ./jacobi-mpi2D N max_iters` and `mpirun -np 16 ./nonblocking N max_iters`,

where `N` is the total number of discretization points, `max_iters` is the maximum number of iterations.
