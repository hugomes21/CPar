CPP = g++ -Wall -pg -Ofast -march=native -msse4.1 -funroll-loops -ftree-vectorize -ffast-math -flto -fno-math-errno -mfpmath=sse
SRCS = main.cpp fluid_solver.cpp EventManager.cpp

all:
	$(CPP) $(SRCS) -o fluid_sim

clean:
	@echo Cleaning up...
	@rm fluid_sim
	@echo Done.
