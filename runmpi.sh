#!/bin/sh

module load gcc/11.2.0

# Get the number of processors from the command line argument
#np=${1:-10}
make runmpi

mpirun -np 1 perf stat -e instructions,cycles,cache-references,cache-misses ./fluid_sim_mpi > "outputmpi_ntasks${np}_rank${OMPI_COMM_WORLD_RANK}.txt"