#!/bin/bash

module load gcc/11.2.0

processor_counts=(20 25 30 35)

for ntasks in "${processor_counts[@]}"; do
    sbatch --job-name="mpi_test_ntasks${ntasks}" \
           --partition=cpar \
           --ntasks=$ntasks \
           --time=5:00 \
           runmpi.sh $ntasks
done
