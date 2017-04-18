// Parallel sample sort

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>

static int compare(const void *a, const void *b)
{
    int *da = (int *)a;
    int *db = (int *)b;
    
    if (*da > *db)
        return 1;
    else if (*da < *db)
        return -1;
    else
        return 0;
}

int main( int argc, char *argv[])
{
    int rank, P, root = 0;
    int i, n, N, S, newN;
    int *vec, *newvec, *sample, *gather, *splitter;
    int *SendCounts, *RecvCounts, *sdispls, *rdispls;
    double time1, time2;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);
    
    /* Number of random numbers per processor (this should be increased
     * for actual tests or could be passed in through the command line */
    if (argc > 1)
        sscanf(argv[1], "%d", &N);
    else
        N = 100;
    if (rank == 0){
        printf("Number of integers to be sorted per processor is %d.\n", N);
    }
    
    /* seed random number generator differently on every core */
    srand((unsigned int) (rank + 393919));
    
    /* fill vector with random integers */
    vec = calloc(N, sizeof(int));
    for (i = 0; i < N; i++)
        vec[i] = rand();
    //printf("rank: %d, first entry: %d\n", rank, vec[0]);
    
    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    time1 = MPI_Wtime();
    
    /* sort locally */
    qsort(vec, N, sizeof(int), compare);
    
    /* randomly sample s entries from vector or select local splitters,
     * i.e., every N/P-th entry of the sorted vector */
    S = P - 1;
    sample = calloc(S, sizeof(int));
    for(i = 1; i <= S; i++)
        sample[i-1] = vec[i*N/P];
    
    /* every processor communicates the selected entries
     * to the root processor; use for instance an MPI_Gather */
    if (rank == root)
        gather = calloc(S * P, sizeof(int));
    MPI_Gather(sample, S, MPI_INT, gather, S, MPI_INT, root, MPI_COMM_WORLD);
    
    /* root processor does a sort, determinates splitters that
     * split the data into P buckets of approximately the same size */
    splitter = calloc(P - 1, sizeof(int));
    if (rank == root) {
        qsort(gather, S * P, sizeof(int), compare);
        for (i = 1; i <= P - 1; i++)
            splitter[i-1] = gather[i*S-1];
        free(gather);
    }
    
    /* root process broadcasts splitters */
    MPI_Bcast(splitter, P - 1, MPI_INT, root, MPI_COMM_WORLD);
    
    /* every processor uses the obtained splitters to decide
     * which integers need to be sent to which other processor (local bins) */
    SendCounts = calloc(P, sizeof(int));
    RecvCounts = calloc(P, sizeof(int));
    i = 0;
    for (n = 0; n < N; n++) {
        if (vec[n] <= splitter[i])
            SendCounts[i]++;
        else if (i < P - 2) {
            i++;
            SendCounts[i]++;
        }
        else
            SendCounts[P - 1]++;
    }
    
    /* send and receive: either you use MPI_AlltoallV, or
     * (and that might be easier), use an MPI_Alltoall to share
     * with every processor how many integers it should expect,
     * and then use MPI_Send and MPI_Recv to exchange the data */
    MPI_Alltoall(SendCounts, 1, MPI_INT, RecvCounts, 1, MPI_INT, MPI_COMM_WORLD);
    
    newN = 0;
    for (i = 0; i < P; i++)
        newN += RecvCounts[i];
    newvec = calloc(newN, sizeof(int));
    
    sdispls = calloc(P, sizeof(int));
    rdispls = calloc(P, sizeof(int));
    for (i = 1; i < P; i++) {
        sdispls[i] = sdispls[i-1] + SendCounts[i-1];
        rdispls[i] = rdispls[i-1] + RecvCounts[i-1];
    }
    
    MPI_Alltoallv(vec, SendCounts, sdispls, MPI_INT, newvec, RecvCounts, rdispls, MPI_INT, MPI_COMM_WORLD);
    
    /* do a local sort */
    qsort(newvec, newN, sizeof(int), compare);
    
    /* timing */
    MPI_Barrier(MPI_COMM_WORLD);
    time2 = MPI_Wtime();
    
    /* every processor writes its result to a file */
    FILE *fp1 = NULL, *fp2 = NULL;
    char filename1[256];
    char filename2[256];
    snprintf(filename1, 256, "output%02d.txt", rank);
    snprintf(filename2, 256, "timing%d.txt", N);
    fp1 = fopen(filename1, "w+");
    fp2 = fopen(filename2, "a+");
    
    if(NULL == fp1 || NULL == fp2) {
        printf("Error opening file\n");
        return 1;
    }
    
    fprintf(fp1, "rank %d sorted %d numbers:\n\n", rank, newN);
    for(i = 0; i < newN; i++)
        fprintf(fp1, "%d\n", newvec[i]);
    fprintf(fp2, "rank %02d takes %f seconds.\n", rank, time2 - time1);
    
    fclose(fp1);
    fclose(fp2);
    
    /* Use barrier for clean output */
    for (i = 0; i < P; i++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (rank == i)
            printf("rank %02d sorted %d numbers.\tTime elapsed is %f seconds.\tmin = %d, max = %d.\n", rank, newN, time2 - time1, newvec[0], newvec[newN-1]);
    }
    
    free(vec);
    free(newvec);
    free(sample);
    free(splitter);
    free(SendCounts);
    free(RecvCounts);
    free(sdispls);
    free(rdispls);
    MPI_Finalize();
    return 0;
}
