// Main solver routines for heat equation solver

#include <mpi.h>

#include "heat.hpp"

// Exchange the boundary values
void exchange(Field& field, const ParallelData& parallel)
{

    // TODO start: implement halo exchange

    // You can utilize the data() method of the Matrix class to obtain pointer
    // to element, e.g. field.temperature.data(i, j)
    double* sbufup = field.temperature.data(1,0);
    double* rbufdown = field.temperature.data(field.nx+1, 0);
    double* sbufdown = field.temperature.data(field.nx,0);
    double* rbufup = field.temperature.data(0,0);

    // Send to up, receive from down
    MPI_Isend(sbufup, field.ny+2, MPI_DOUBLE, parallel.nup, 0, MPI_COMM_WORLD, parallel.requests[0]);
    MPI_Irecv(rbufdown, field.ny+2, MPI_DOUBLE, parallel.ndown, 0, MPI_COMM_WORLD, parallel.requests[1]);
    // Send to down, receive from up
    MPI_Isend(sbufdown, field.ny+2, MPI_DOUBLE, parallel.ndown, 0, MPI_COMM_WORLD, parallel.requests[2]);
    MPI_Irecv(rbufup, field.ny+2, MPI_DOUBLE, parallel.nup , 0, MPI_COMM_WORLD, parallel.requests[3]);

    // TODO end
}

// Update the temperature values using five-point stencil */
void evolve(Field& curr, const Field& prev, const double a, const double dt, ParallelData& parallel)
{

  // Compilers do not necessarily optimize division to multiplication, so make it explicit
  auto inv_dx2 = 1.0 / (prev.dx * prev.dx);
  auto inv_dy2 = 1.0 / (prev.dy * prev.dy);

  // Determine the temperature field at next time step
  // As we have fixed boundary conditions, the outermost gridpoints
  // are not updated.
  if (curr.nx > 3) {
    for (int i = 3; i < curr.nx-1; i++) {
      for (int j = 1; j < curr.ny+1; j++) {
              curr(i, j) = prev(i, j) + a * dt * (
	         ( prev(i + 1, j) - 2.0 * prev(i, j) + prev(i - 1, j) ) * inv_dx2 +
	         ( prev(i, j + 1) - 2.0 * prev(i, j) + prev(i, j - 1) ) * inv_dy2
                 );
      }
    }
  
    // Wait for ghost layer exchange to finish before calculating outer values.
    MPI_Waitall(4, *parallel.requests, MPI_STATUSES_IGNORE);

    for (int i = 1; i < 3; i++) {
      for (int j = 1; j < curr.ny+1; j++) {
              curr(i, j) = prev(i, j) + a * dt * (
	         ( prev(i + 1, j) - 2.0 * prev(i, j) + prev(i - 1, j) ) * inv_dx2 +
	         ( prev(i, j + 1) - 2.0 * prev(i, j) + prev(i, j - 1) ) * inv_dy2
                 );
      }
    }
    for (int i = curr.nx-2; i < curr.nx+1; i++) {
      for (int j = 1; j < curr.ny+1; j++) {
              curr(i, j) = prev(i, j) + a * dt * (
	         ( prev(i + 1, j) - 2.0 * prev(i, j) + prev(i - 1, j) ) * inv_dx2 +
	         ( prev(i, j + 1) - 2.0 * prev(i, j) + prev(i, j - 1) ) * inv_dy2
                 );
      }
    }
  } else {
    // No regions that don't neighbour unsafe ones.
    MPI_Waitall(4, *parallel.requests, MPI_STATUSES_IGNORE);

    for (int i = 1; i < curr.nx+1; i++) {
      for (int j = 1; j < curr.ny+1; j++) {
              curr(i, j) = prev(i, j) + a * dt * (
	         ( prev(i + 1, j) - 2.0 * prev(i, j) + prev(i - 1, j) ) * inv_dx2 +
	         ( prev(i, j + 1) - 2.0 * prev(i, j) + prev(i, j - 1) ) * inv_dy2
                 );
      }
    }
  }
}
