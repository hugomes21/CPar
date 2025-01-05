#!/bin/bash

# Load necessary modules
module load gcc/11.2.0

# Array of processor counts to test
processor_counts=(20 25 30 35)

for ntasks in "${processor_counts[@]}"; do
    # Submit the job to the scheduler with the desired number of tasks
    sbatch --job-name="mpi_test_ntasks${ntasks}" --partition=day --constraint=c24 --ntasks=$ntasks --time=5:00 runmpi_test.sh $ntasks
done