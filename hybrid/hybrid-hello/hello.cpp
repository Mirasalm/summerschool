#include <cstdio>
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    int my_id, omp_rank;
    int provided, required=MPI_THREAD_FUNNELED;

    /* TODO: Initialize MPI with thread support. */
    MPI_Init_thread(&argc, &argv, required, &provided);
    
    #pragma omp parallel private(omp_rank, my_id)
    {

    /* TODO: Find out the MPI rank and thread ID of each thread and print
     *       out the results. */
    MPI_Comm_rank(MPI_COMM_WORLD, &omp_rank);
    my_id = omp_get_thread_num();
    printf("Hello from thread %i with MPI rank %i\n", my_id, omp_rank);

    /* TODO: Investigate the provided thread support level. */

    if (my_id == 0) {
    printf("\nProvided thread support level: %d\n", provided);
    printf("  %d - MPI_THREAD_SINGLE\n", MPI_THREAD_SINGLE);
    printf("  %d - MPI_THREAD_FUNNELED\n", MPI_THREAD_FUNNELED);
    printf("  %d - MPI_THREAD_SERIALIZED\n", MPI_THREAD_SERIALIZED);
    printf("  %d - MPI_THREAD_MULTIPLE\n", MPI_THREAD_MULTIPLE);

    }
    }
    MPI_Finalize();
    return 0;
}
