#!/bin/sh
#
#SBATCH --exclusive      # exclusive node for the job
#SBATCH --time=02:00     # allocation for 2 minutes

# Run the sequential version first
echo "Running fluid_sim_seq (1 thread)..."
(time ./fluid_sim_seq)

# Run fluid_sim with different thread counts
for threads in 2 4 6 8 12 16 20 24 30 40
do
    echo "Running fluid_sim with $threads threads..." 
    export OMP_NUM_THREADS=$threads
    (time ./fluid_sim)
done
