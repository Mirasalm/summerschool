#include <iostream>
# include <mpi.h>

int main(int argc, char *argv[])
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    char name[MPI_MAX_PROCESSOR_NAME];
    int namelen;

    MPI_Get_processor_name(name, &namelen);

    std::cout << "Hello from rank " << rank << " running on " << name  << "!";
    
    if (rank == 0) {
        std::cout << " Number of MPI processes is " << size;
    }

    std::cout << std::endl;
    MPI_Finalize();
}
