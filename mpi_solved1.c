/******************************************************************************
 * FILE: mpi_solved1.c
 * DESCRIPTION:
 *   This program has a bug that causes it to hang. Fixed
 * AUTHOR: Manyuan Tao
 * LAST REVISED: 04/12/17
 ******************************************************************************/

/*
 BUG: tags of MPI_Send and MPI_Recv do not match each other, while message sent with a tag will be received by a process looking for a message with the same tag.
 SOLUTION: set 'tag' to be the same value for each process, here I set 'tag = 99'.
*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

int main (int argc, char *argv[])
{
    int numtasks, rank, dest, tag, source, rc, count;
    char inmsg, outmsg='x';
    MPI_Status Stat;
    
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    printf("Task %d starting...\n",rank);
    
    tag = 99;   //!BUG: set 'tag' to be the same value for each process
    
    if (rank == 0) {
        if (numtasks > 2)
            printf("Numtasks=%d. Only 2 needed. Ignoring extra...\n",numtasks);
        dest = rank + 1;
        source = dest;
        //tag = rank;   //!BUG: tags of send/receive do not match
        rc = MPI_Send(&outmsg, 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
        printf("Sent to task %d...\n",dest);
        rc = MPI_Recv(&inmsg, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, &Stat);
        printf("Received from task %d...\n",source);
    }
    
    else if (rank == 1) {
        dest = rank - 1;
        source = dest;
        //tag = rank;   //!BUG: tags of send/receive do not match
        rc = MPI_Recv(&inmsg, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, &Stat);
        printf("Received from task %d...\n",source);
        rc = MPI_Send(&outmsg, 1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
        printf("Sent to task %d...\n",dest);
    }
    
    if (rank < 2) {
        rc = MPI_Get_count(&Stat, MPI_CHAR, &count);
        printf("Task %d: Received %d char(s) from task %d with tag %d \n",
               rank, count, Stat.MPI_SOURCE, Stat.MPI_TAG);
    }
    
    MPI_Finalize();
}
