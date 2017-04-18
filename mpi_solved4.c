/******************************************************************************
 * FILE: mpi_solved4.c
 * DESCRIPTION:
 *   This program gives the wrong result for Final sum - compare to mpi_array. Fixed.
 *   The number of MPI tasks must be divisible by 4.
 *   Hint: Every process needs to be included in a collective communication call.
 * AUTHOR: Manyuan Tao
 * LAST REVISED: 04/12/17
 ******************************************************************************/

/*
 BUG1: missing MPI_Reduce() in the Master process(rank=0).
 SOLUTION: call MPI_Reduce() in the Master process.
 BUG2: the precision of 'float' is not enough for the large summation.
 SOLUTION: change 'float' to 'double', for higher precision.
*/

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#define  ARRAYSIZE	16000000
#define  MASTER		0

double  data[ARRAYSIZE];    //!BUG2: change datatype from float to double

int main (int argc, char *argv[])
{
    int   numtasks, taskid, rc, dest, offset, i, j, tag1, tag2, source, chunksize;
    double mysum, sum;
    double update(int myoffset, int chunk, int myid);
    MPI_Status status;
    
    /***** Initializations *****/
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    if (numtasks % 4 != 0) {
        printf("Quitting. Number of MPI tasks must be divisible by 4.\n");
        MPI_Abort(MPI_COMM_WORLD, rc);
        exit(0);
    }
    MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
    printf ("MPI task %d has started...\n", taskid);
    chunksize = (ARRAYSIZE / numtasks);
    tag2 = 1;
    tag1 = 2;
    
    /***** Master task only ******/
    if (taskid == MASTER){
        
        /* Initialize the array */
        sum = 0;
        for(i=0; i<ARRAYSIZE; i++) {
            data[i] =  i * 1.0;
            sum = sum + data[i];
        }
        printf("Initialized array sum = %e\n",sum);
        
        /* Send each task its portion of the array - master keeps 1st part */
        offset = chunksize;
        for (dest=1; dest<numtasks; dest++) {
            MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
            MPI_Send(&data[offset], chunksize, MPI_DOUBLE, dest, tag2, MPI_COMM_WORLD);
            printf("Sent %d elements to task %d offset= %d\n",chunksize,dest,offset);
            offset = offset + chunksize;
        }
        
        /* Master does its part of the work */
        offset = 0;
        mysum = update(offset, chunksize, taskid);
        
        /* Wait to receive results from each task */
        for (i=1; i<numtasks; i++) {
            source = i;
            MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
            MPI_Recv(&data[offset], chunksize, MPI_DOUBLE, source, tag2,
                     MPI_COMM_WORLD, &status);
        }
        
        /* Get final sum and print sample results */
        
        MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);   //!BUG1: MASTER needs to call MPI_Reduce()
        
        printf("Sample results: \n");
        offset = 0;
        for (i=0; i<numtasks; i++) {
            for (j=0; j<5; j++)
                printf("  %e",data[offset+j]);
            printf("\n");
            offset = offset + chunksize;
        }
        printf("*** Final sum= %e ***\n",sum);
        
    }  /* end of master section */
    
    
    
    /***** Non-master tasks only *****/
    
    if (taskid > MASTER) {
        
        /* Receive my portion of array from the master task */
        source = MASTER;
        MPI_Recv(&offset, 1, MPI_INT, source, tag1, MPI_COMM_WORLD, &status);
        MPI_Recv(&data[offset], chunksize, MPI_DOUBLE, source, tag2,
                 MPI_COMM_WORLD, &status);
        
        mysum = update(offset, chunksize, taskid);
        
        /* Send my results back to the master task */
        dest = MASTER;
        MPI_Send(&offset, 1, MPI_INT, dest, tag1, MPI_COMM_WORLD);
        MPI_Send(&data[offset], chunksize, MPI_DOUBLE, MASTER, tag2, MPI_COMM_WORLD);
        
        MPI_Reduce(&mysum, &sum, 1, MPI_DOUBLE, MPI_SUM, MASTER, MPI_COMM_WORLD);
        
    } /* end of non-master */
    
    MPI_Finalize();
    
}   /* end of main */


double update(int myoffset, int chunk, int myid) {
    int i; 
    double mysum;
    /* Perform addition to each of my array elements and keep my sum */
    mysum = 0;
    for(i=myoffset; i < myoffset + chunk; i++) {
        data[i] = data[i] + i * 1.0;
        mysum = mysum + data[i];
    }
    printf("Task %d mysum = %e\n",myid,mysum);
    return(mysum);
}
