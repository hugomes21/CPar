#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include <omp.h>
#include <mpi.h>
#include <cmath>
#include <cstdint>

void vel_step(int M, int N, int O, float *u, float *v, float *w, float *u0, float *v0, float *w0, float visc, float dt, int rank, int size, int start_i, int end_i, int start_j, int end_j, int start_k, int end_k, const int dims[3], MPI_Comm cart_comm);

void dens_step(int M, int N, int O, float *x, float *x0, float *u, float *v, float *w, float diff, float dt, int rank, int size, int start_i, int end_i, int start_j, int end_j, int start_k, int end_k, const int dims[3], MPI_Comm cart_comm);
void distribute_workload(int M, int N, int O, int size, int rank, int &start_i,
                         int &end_i, int &start_j, int &end_j, int &start_k,
                         int &end_k, int dims[3], MPI_Comm &cart_comm);

#endif // FLUID_SOLVER_H