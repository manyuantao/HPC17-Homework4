/******************************************************************************
 * FILE: mpi_solved5.c
 * DESCRIPTION:
 *   This is an "unsafe" program. It's behavior varies depending upon the
 *   platform and MPI library. Fixed.
 * AUTHOR: Manyuan Tao
 * LAST REVISED: 04/12/17
 * Hint: If possible, try to run the program on two different machines,
 * which are connected through a network. You should see uneven timings;
 * try to understand/explain them.
 ******************************************************************************/

/*
 BUG: rank 0 uses a blocking send(MPI_Send) that returns to the program once the buffer stores the message, so it continues to send data until the buffer becomes full; while rank 1 does more work(slower) than the send task, thus when the buffer storing sending data and waiting for being received is full, the MPI_Send call blocks and has to wait until there is space in the buffer. This causes uneven timings.
 SOLUTION: use a synchronous blocking send(MPI_Ssend) where the sender will not return until matching receive is posted.
*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define MSGSIZE 2000

int main (int argc, char *argv[])
{
    int        numtasks, rank, i, tag=111, dest=1, source=0, count=0;
    char       data[MSGSIZE];
    double     start, end, result;
    MPI_Status status;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    if (rank == 0) {
        printf ("mpi_bug5 has started...\n");
        if (numtasks > 2)
            printf("INFO: Number of tasks= %d. Only using 2 tasks.\n", numtasks);
    }
    
    /******************************* Send task **********************************/
    if (rank == 0) {
        
        /* Initialize send data */
        for(i=0; i<MSGSIZE; i++)
            data[i] =  'x';
        
        start = MPI_Wtime();
        while (1) {
            MPI_Ssend(data, MSGSIZE, MPI_BYTE, dest, tag, MPI_COMM_WORLD);   //!BUG: use MPI_Ssend() for synchronization
            count++;
            if (count % 10 == 0) {
                end = MPI_Wtime();
                printf("Count= %d  Time= %f sec.\n", count, end-start);
                start = MPI_Wtime();
            }
        }
    }
    
    /****************************** Receive task ********************************/
    
    if (rank == 1) {
        while (1) {
            MPI_Recv(data, MSGSIZE, MPI_BYTE, source, tag, MPI_COMM_WORLD, &status);
            /* Do some work  - at least more than the send task */
            result = 0.0;
            for (i=0; i < 1000000; i++) 
                result = result + (double)random();
        }
    }
    
    MPI_Finalize();
}
