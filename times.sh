#!/bin/sh
#
#SBATCH --exclusive 	  # exclusive node for the job
#SBATCH --time=02:00      # allocation for 2 minutes
 
time ./fluid_sim_seq
export OMP_NUM_THREADS=1
time ./fluid_sim
export OMP_NUM_THREADS=2
time ./fluid_sim
export OMP_NUM_THREADS=4
time ./fluid_sim
export OMP_NUM_THREADS=6
time ./fluid_sim
export OMP_NUM_THREADS=8
time ./fluid_sim
export OMP_NUM_THREADS=12
time ./fluid_sim
export OMP_NUM_THREADS=16
time ./fluid_sim
export OMP_NUM_THREADS=20
time ./fluid_sim