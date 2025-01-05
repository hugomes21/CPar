#!/bin/bash
#SBATCH --output=fluid_sim_%j.out
#SBATCH --error=fluid_sim_%j.err

# Load necessary modules
module load openmpi

# Get the number of processors from the command line argument
np=$1

# Run the MPI program with perf stat and redirect output to a file
mpirun -np $np perf stat -e instructions,cycles,cache-references,cache-misses ./fluid_sim_mpi > "outputmpi_ntasks${np}.txt"