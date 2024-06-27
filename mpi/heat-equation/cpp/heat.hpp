#pragma once
#include <string>
#include "matrix.hpp"
#include <mpi.h>

// Class for basic parallelization information
struct ParallelData {
    int size;            // Number of MPI tasks
    int rank;
    int nup, ndown;      // Ranks of neighbouring MPI tasks
    MPI_Request requests[4];
    MPI_Comm comm;

    ParallelData() {      // Constructor

      // TODO start: query number of MPI tasks and store it in
      // the size attribute of the class
      MPI_Comm_size(MPI_COMM_WORLD, &size);
      int dims[1] = {size};
      int periods[1] = {0};

      MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 1, &comm);
      MPI_Comm_rank(comm, &rank);
      MPI_Cart_shift(comm, 0, 1, &nup, &ndown);
      // TODO end

    };

};

// Class for temperature field
struct Field {
    // nx and ny are the true dimensions of the field. The temperature matrix
    // contains also ghost layers, so it will have dimensions nx+2 x ny+2
    int nx;                     // Local dimensions of the field
    int ny;
    int nx_full;                // Global dimensions of the field
    int ny_full;                // Global dimensions of the field
    double dx = 0.01;           // Grid spacing
    double dy = 0.01;

    Matrix<double> temperature;

    void setup(int nx_in, int ny_in, ParallelData parallel);

    void generate(ParallelData parallel);

    // standard (i,j) syntax for setting elements
    double& operator()(int i, int j) {return temperature(i, j);}

    // standard (i,j) syntax for getting elements
    const double& operator()(int i, int j) const {return temperature(i, j);}

};

// Function declarations
void initialize(int argc, char *argv[], Field& current,
                Field& previous, int& nsteps, ParallelData parallel);

void exchange_init(Field& field, ParallelData& parallel);

void exchange_finalize(ParallelData& parallel);

void evolve_interior(Field& curr, const Field& prev, const double a, const double dt);

void evolve_edges(Field& curr, const Field& prev, const double a, const double dt);

void write_field(const Field& field, const int iter, const ParallelData parallel);

void read_field(Field& field, std::string filename,
                const ParallelData parallel);

double average(const Field& field, const ParallelData parallel);
