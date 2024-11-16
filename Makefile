CPP = g++ -Wall -Ofast -msse4.1 -funroll-loops -ftree-vectorize -ffast-math -fno-math-errno -mfpmath=sse
SRCS = main.cpp fluid_solver.cpp EventManager.cpp

all: runseq runpar

clean:
	@echo Cleaning up...
	@rm -f fluid_sim fluid_sim_seq
	@echo Done.

runseq:
	$(CPP) $(SRCS) -o fluid_sim_seq

runpar:
	$(CPP) -fopenmp $(SRCS) -o fluid_sim  
