// The MPI-parallel nonblocking version of Jacobi method

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "util.h"
#include <mpi.h>
#include <assert.h>

#define index(i,j) (i)*(lN+2)+j

double residual(int lN, double h, double *lu);

int main(int argc, char *argv[])
{
    int mpirank, p, sp, i, j, N, lN, iter, max_iters;
    double h, gres, gres_init;
    MPI_Status status;
    MPI_Request request_out1, request_out2, request_out3, request_out4;
    MPI_Request request_in1, request_in2, request_in3, request_in4;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    /* get name of host running MPI process */
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);
    printf("Rank %d/%d running on %s.\n", mpirank, p, processor_name);
    
    sscanf(argv[1], "%d", &N);
    sscanf(argv[2], "%d", &max_iters);
    h = 1./(N + 1);
    
    sp = floor(sqrt(p + 0.5));
    assert(p == sp * sp);
    
    /* compute number of unknowns handled by each process */
    lN = N / sp;
    if ((N % sp != 0) && mpirank == 0) {
        printf("N: %d, local N: %d\n", N, lN);
        printf("Exiting. N must be a multiple of sqrt(p)\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }
    
    /* allocation of vectors, including boundary ghost points */
    double * lu_pre = calloc((lN+2)*(lN+2), sizeof(double));
    double * lu_next = calloc((lN+2)*(lN+2), sizeof(double));
    double * ghost_left = calloc(lN, sizeof(double));
    double * ghost_right = calloc(lN, sizeof(double));
    double * ghost_up = calloc(lN, sizeof(double));
    double * ghost_down = calloc(lN, sizeof(double));
    
    /* compute initial residual */
    gres_init = residual(lN, h, lu_pre);
    gres = gres_init;
    
    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    timestamp_type time1, time2;
    get_timestamp(&time1);
    
    for (iter = 1; iter <= max_iters && gres_init/gres < 1e4; iter++) {
        /* interleaf computation and communication: compute the left, right,
         * up and down values, which are communicated with non-blocking
         * send/recv. During that communication, do all the local work */
        
        /* Jacobi step for the left, right, up and down points */
        for (i = 1; i <= lN; i++) {
            lu_next[index(i,1)] = 0.25 * (h * h + lu_pre[index(i-1,1)] + lu_pre[index(i,0)] + lu_pre[index(i+1,1)] + lu_pre[index(i,2)]);
            lu_next[index(i,lN)] = 0.25 * (h * h + lu_pre[index(i-1,lN)] + lu_pre[index(i,lN-1)] + lu_pre[index(i+1,lN)] + lu_pre[index(i,lN+1)]);
        }
        for (j = 1; j <= lN; j++) {
            lu_next[index(1,j)] = 0.25 * (h * h + lu_pre[index(0,j)] + lu_pre[index(1,j-1)] + lu_pre[index(2,j)] + lu_pre[index(1,j+1)]);
            lu_next[index(lN,j)] = 0.25 * (h * h + lu_pre[index(lN-1,j)] + lu_pre[index(lN,j-1)] + lu_pre[index(lN+1,j)] + lu_pre[index(lN,j+1)]);
        }
        
        /* communicate ghost values */
        if (mpirank % sp < sp - 1) {
            for (i = 1; i <= lN; i++)
                ghost_left[i-1] = lu_next[index(i,lN)];
            MPI_Isend(ghost_left, lN, MPI_DOUBLE, mpirank+1, 124, MPI_COMM_WORLD, &request_out1);
            MPI_Irecv(ghost_right, lN, MPI_DOUBLE, mpirank+1, 123, MPI_COMM_WORLD, &request_in1);
        }
        if (mpirank % sp > 0) {
            for (i = 1; i <= lN; i++)
                ghost_right[i-1] = lu_next[index(i,1)];
            MPI_Isend(ghost_right, lN, MPI_DOUBLE, mpirank-1, 123, MPI_COMM_WORLD, &request_out2);
            MPI_Irecv(ghost_left, lN, MPI_DOUBLE, mpirank-1, 124, MPI_COMM_WORLD, &request_in2);
        }
        if (mpirank < p - sp) {
            for (j = 1; j <= lN; j++)
                ghost_down[j-1] = lu_next[index(lN,j)];
            MPI_Isend(ghost_down, lN, MPI_DOUBLE, mpirank+sp, 126, MPI_COMM_WORLD, &request_out3);
            MPI_Irecv(ghost_up, lN, MPI_DOUBLE, mpirank+sp, 125, MPI_COMM_WORLD, &request_in3);
        }
        if (mpirank > sp - 1) {
            for (j = 1; j <= lN; j++)
                ghost_up[j-1] = lu_next[index(1,j)];
            MPI_Isend(ghost_up, lN, MPI_DOUBLE, mpirank-sp, 125, MPI_COMM_WORLD, &request_out4);
            MPI_Irecv(ghost_down, lN, MPI_DOUBLE, mpirank-sp, 126, MPI_COMM_WORLD, &request_in4);
        }
        
        /* Jacobi step for all the inner points */
        for (i = 2; i < lN; i++)
            for (j = 2; j < lN; j++)
                lu_next[index(i,j)]= 0.25 * (h * h + lu_pre[index(i-1,j)] + lu_pre[index(i,j-1)] + lu_pre[index(i+1,j)] + lu_pre[index(i,j+1)]);
        
        /* check if Isend/Irecv are done */
        if (mpirank % sp < sp - 1) {
            MPI_Wait(&request_out1, &status);
            MPI_Wait(&request_in1, &status);
            for (i = 1; i <= lN; i++)
                lu_next[index(i,lN+1)] = ghost_right[i-1];
        }
        if (mpirank % sp > 0) {
            MPI_Wait(&request_out2, &status);
            MPI_Wait(&request_in2, &status);
            for (i = 1; i <= lN; i++)
                lu_next[index(i,0)] = ghost_left[i-1];
        }
        if (mpirank < p - sp) {
            MPI_Wait(&request_out3, &status);
            MPI_Wait(&request_in3, &status);
            for (j = 1; j <= lN; j++)
                lu_next[index(lN+1,j)] = ghost_up[j-1];
        }
        if (mpirank > sp - 1) {
            MPI_Wait(&request_out4, &status);
            MPI_Wait(&request_in4, &status);
            for (j = 1; j <= lN; j++)
                lu_next[index(0,j)] = ghost_down[j-1];
        }
        
        /* compute residual for each iteration */
        gres = residual(lN, h, lu_next);
        if (mpirank == 0)
            printf("Iteration: %d;\t Residual: %f.\n", iter, gres);
        
        /* copy lu_next onto lu_pre */
        memcpy(lu_pre, lu_next, (lN+2)*(lN+2)*sizeof(double));
    }
    
    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    get_timestamp(&time2);
    double elapsed = timestamp_diff_in_seconds(time1, time2);
    if (mpirank == 0)
        printf("Time elapsed is %f seconds.\n", elapsed);
    
    /* clean up */
    free(lu_pre);
    free(lu_next);
    free(ghost_left);
    free(ghost_right);
    free(ghost_up);
    free(ghost_down);
        
    MPI_Finalize();
    return 0;
}

/* define funtion to compute residual */
double residual(int lN, double h, double *lu)
{
    int i, j;
    double temp, gres = 0.0, lres = 0.0;
    
    for (i = 1; i <= lN; i++)
        for (j = 1; j <= lN; j++){
            temp = (4 * lu[index(i,j)] - lu[index(i-1,j)] - lu[index(i,j-1)] - lu[index(i+1,j)] - lu[index(i,j+1)]) / (h * h) - 1;
            lres += temp * temp;
        }
    
    MPI_Allreduce(&lres, &gres, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return sqrt(gres);
}
