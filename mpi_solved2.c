/******************************************************************************
 * FILE: mpi_solved2.c
 * DESCRIPTION:
 *   This program has a bug that causes wrong answers and/or termination - depends
 *   upon the MPI library and platform. Fixed
 * AUTHOR: Manyuan Tao
 * LAST REVISED: 04/12/17
 ******************************************************************************/

/*
 BUG: datatypes in MPI_Isend and MPI_Irecv do not match each other.
 SOLUTION: change datatype of 'beta' from float to int.
*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
    int numtasks, rank, tag=1, alpha, i;
    int beta;   //!BUG: datatype of 'beta' didn't match that of 'alpha'
    MPI_Request reqs[10];
    MPI_Status stats[10];
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0) {
        if (numtasks > 2)
            printf("Numtasks=%d. Only 2 needed. Ignoring extra...\n",numtasks);
        for (i=0; i<10; i++) {
            alpha = i*10;
            MPI_Isend(&alpha, 1, MPI_INT, 1, tag, MPI_COMM_WORLD, &reqs[i]);
            MPI_Wait(&reqs[i], &stats[i]);
            printf("Task %d sent = %d\n",rank,alpha);
        }
    }
    
    if (rank == 1) {
        for (i=0; i<10; i++) {
            MPI_Irecv(&beta, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &reqs[i]);   //!BUG: datatype should be 'MPI_INT'
            MPI_Wait(&reqs[i], &stats[i]);
            printf("Task %d received = %d\n",rank,beta);    //!BUG: should be '%d' instead of '%f'
        }
    }
    
    MPI_Finalize();
}
