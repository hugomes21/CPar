CPP = g++ -Wall -Ofast -msse4.1 -funroll-loops -ftree-vectorize -ffast-math -fno-math-errno -mfpmath=sse
SRCS_SEQ = main.cpp fluid_solver.cpp EventManager.cpp
SRCS_PAR = main.cpp fluid_solver.cpp EventManager.cpp
SRCS_MPI = main_mpi.cpp fluid_solver_mpi.cpp EventManager.cpp
MPI_CPP = mpic++
NUM_THREADS = 20
MPI_INCLUDE = -I/usr/lib/x86_64-linux-gnu/openmpi/include

all: fluid_sim_seq fluid_sim fluid_sim_mpi

fluid_sim_seq: $(SRCS_SEQ)
    @echo Compiling with command: $(CPP) $(SRCS_SEQ) -o fluid_sim_seq
    $(CPP) $(SRCS_SEQ) -o fluid_sim_seq

fluid_sim: $(SRCS_PAR)
    @echo Compiling with OpenMP support
    $(CPP) -fopenmp $(SRCS_PAR) -o fluid_sim

fluid_sim_mpi: $(SRCS_MPI)
    @echo Compiling with MPI support
    $(MPI_CPP) $(MPI_INCLUDE) $(SRCS_MPI) -o fluid_sim_mpi

clean:
    @echo Cleaning up...
    @rm -f fluid_sim fluid_sim_seq fluid_sim_mpi
    @echo Done.

run: all
    @echo Submitting job to the cluster with sbatch
    @sbatch runmpi.sh

.PHONY: all clean run