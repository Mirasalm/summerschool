#include <cstdio>
#include <mpi.h>
#include <omp.h>

int main(int argc, char *argv[])
{
    int my_id, omp_rank, ntasks;
    int provided, required=MPI_THREAD_MULTIPLE;

    MPI_Init_thread(&argc, &argv, required, &provided);
    MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
    
    #pragma omp parallel shared(omp_rank) private(my_id)
    {

    int recv_id;

    #pragma omp single
    MPI_Comm_rank(MPI_COMM_WORLD, &omp_rank);
    my_id = omp_get_thread_num();

    if (omp_rank == 0) {
        for (int i = 1; i < ntasks; i++) {
            MPI_Send(&my_id, 1, MPI_INT, i, my_id, MPI_COMM_WORLD);
            printf("Thread %i on rank %i sent to thread %i on rank %i\n", my_id, omp_rank, my_id, i);
        }

    } else {
        MPI_Recv(&recv_id, 1, MPI_INT, 0, my_id, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("Thread %i on rank %i received from thread %i on rank %i\n", my_id, omp_rank, recv_id, 0);

    }
    }
    MPI_Finalize();
    return 0;
}
