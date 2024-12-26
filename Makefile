CPP = g++ -Wall -Ofast -msse4.1 -funroll-loops -ftree-vectorize -ffast-math -fno-math-errno -mfpmath=sse
SRCS_SEQ = main.cpp fluid_solver_par.cpp EventManager.cpp
SRCS_PAR = main.cpp fluid_solver_par.cpp EventManager.cpp
SRCS_MPI = main.cpp fluid_solver.cpp EventManager.cpp
NUM_THREADS = 20

all: runseq runpar runmpi

clean:
	@echo Cleaning up...
	@rm -f fluid_sim fluid_sim_seq fluid_sim_mpi
	@echo Done.

runseq:
	module load gcc/11.2.0; $(CPP) $(SRCS_SEQ) -o fluid_sim_seq

runpar:
	module load gcc/11.2.0; $(CPP) -fopenmp $(SRCS_PAR) -o fluid_sim
	@echo Executing with $(NUM_THREADS) threads
	OMP_NUM_THREADS=$(NUM_THREADS) ./fluid_sim

runmpi:
	module load gcc/11.2.0; mpic++ $(SRCS_MPI) -o fluid_sim_mpi
	@echo Executing with MPI
	mpirun -np $(NUM_THREADS) ./fluid_sim_mpi