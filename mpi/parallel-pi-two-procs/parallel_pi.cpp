#include <cstdio>
#include <cmath>
#include <mpi.h>

constexpr int n = 840;

int main(int argc, char** argv)
{
  int rank, ntasks;
  
  //int arraysize = 100000;
  //int msgsize = 100000;
  //int *message;
  //int *receiveBuffer;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  printf("Computing approximation to pi with N=%d\n", n);
    
  int istart = 1 + (rank * n);
  int istop = n + (rank * n);
  int N = n * ntasks;

  double pi = 0.0;

  for (int i=istart; i <= istop; i++) {
    double x = (i - 0.5) / N;
    pi += 1.0 / (1.0 + x*x);
  }

  if (rank == 1) {
    MPI_Send(&pi, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  else if (rank == 0) {
    double pi_recv;
    MPI_Recv(&pi_recv, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, &status);
    pi += pi_recv;
    pi *= 4.0 / N;
    printf("Approximate pi=%18.16f (exact pi=%10.8f)\n", pi, M_PI);
  }

  MPI_Finalize();

}
