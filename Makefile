CPP = g++ -Wall -fopenmp -Ofast -msse4.1 -funroll-loops -ftree-vectorize -ffast-math -fno-math-errno -mfpmath=sse
SRCS = main.cpp fluid_solver.cpp EventManager.cpp

all: phase2

phase2:
	$(CPP) $(SRCS) -o fluid_sim

clean:
	@echo Cleaning up...
	@rm fluid_sim
	@echo Done.

runseq:
	./fluid_sim_seq

runpar:
	./fluid_sim
