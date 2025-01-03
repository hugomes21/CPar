#include "EventManager.h"
#include "fluid_solver_mpi.h"
#include <iostream>
#include <vector>
#include <mpi.h>

#define SIZE 168

#define IX(i, j, k) ((i) + (M + 2) * (j) + (M + 2) * (N + 2) * (k))

// Globals for the grid size
static int M = SIZE;
static int N = SIZE;
static int O = SIZE;
static float dt = 0.1f;      // Time delta
static float diff = 0.0001f; // Diffusion constant
static float visc = 0.0001f; // Viscosity constant

// Fluid simulation arrays
static float *u, *v, *w, *u_prev, *v_prev, *w_prev;
static float *dens, *dens_prev;

// Function to allocate simulation data
int allocate_data() {
  int size = (M + 2) * (N + 2) * (O + 2);
  u = new float[size];
  v = new float[size];
  w = new float[size];
  u_prev = new float[size];
  v_prev = new float[size];
  w_prev = new float[size];
  dens = new float[size];
  dens_prev = new float[size];

  if (!u || !v || !w || !u_prev || !v_prev || !w_prev || !dens || !dens_prev) {
    std::cerr << "Cannot allocate memory" << std::endl;
    return 0;
  }
  return 1;
}

// Function to clear the data (set all to zero)
void clear_data() {
  int size = (M + 2) * (N + 2) * (O + 2);
  for (int i = 0; i < size; i++) {
    u[i] = v[i] = w[i] = u_prev[i] = v_prev[i] = w_prev[i] = dens[i] =
        dens_prev[i] = 0.0f;
  }
}

// Free allocated memory
void free_data() {
  delete[] u;
  delete[] v;
  delete[] w;
  delete[] u_prev;
  delete[] v_prev;
  delete[] w_prev;
  delete[] dens;
  delete[] dens_prev;
}

// Apply events (source or force) for the current timestep
void apply_events(const std::vector<Event> &events) {
  for (const auto &event : events) {
    if (event.type == ADD_SOURCE) {
      // Apply density source at the center of the grid
      int i = M / 2, j = N / 2, k = O / 2;
      dens[IX(i, j, k)] = event.density;
    } else if (event.type == APPLY_FORCE) {
      // Apply forces based on the event's vector (fx, fy, fz)
      int i = M / 2, j = N / 2, k = O / 2;
      u[IX(i, j, k)] = event.force.x;
      v[IX(i, j, k)] = event.force.y;
      w[IX(i, j, k)] = event.force.z;
    }
  }
}

// Function to sum the total density
float sum_density(int M, int N, int O, float *x) {
    float total_density = 0.0f;
    for (int i = 0; i < (M + 2) * (N + 2) * (O + 2); i++) {
        total_density += x[i];
    }
    //printf("Local total density: %f\n", total_density); // Debug message
    return total_density;
}

// Simulation loop
void simulate(EventManager &eventManager, int timesteps, int rank, int size) {
  for (int t = 0; t < timesteps; t++) {
    // Get the events for the current timestep
    std::vector<Event> events = eventManager.get_events_at_timestamp(t);

    // Apply events to the simulation
    apply_events(events);

    // Perform the simulation steps
    vel_step(M, N, O, u, v, w, u_prev, v_prev, w_prev, visc, dt, rank, size);
    dens_step(M, N, O, dens, dens_prev, u, v, w, diff, dt, rank, size);
  }
}

int main() {
    // Inicialização do MPI
    int mpi = MPI_Init(NULL, NULL);
    if (mpi != MPI_SUCCESS) {
        printf("Error starting MPI program. Terminating.\n");
        MPI_Abort(MPI_COMM_WORLD, mpi);
    }

    // Obter rank e tamanho do MPI
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    printf("Rank %d: MPI initialized with size %d\n", rank, size);

    // Initialize EventManager
    EventManager eventManager;
    eventManager.read_events("events.txt");

    // Get the total number of timesteps from the event file
    int timesteps = eventManager.get_total_timesteps();

    // Allocate and clear data
    if (!allocate_data())
        return -1;
    clear_data();
    
    // Run simulation with events
    simulate(eventManager, timesteps, rank, size);

    // Sincronize all processes
    MPI_Barrier(MPI_COMM_WORLD);

    // Print total density at the end of simulation
    float local_density = sum_density(M, N, O, dens);
    printf("Rank %d: Local density = %f\n", rank, local_density);
    float global_density;
    MPI_Reduce(&local_density, &global_density, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) { // Apenas o processo com rank 0 imprime o resultado
        std::cout << "BOAAAAS Total density after " << timesteps
                  << " timesteps: " << global_density << std::endl;
    }

    // Free memory
    free_data();

    // Finalização do MPI
    MPI_Finalize();

    printf("Rank %d: MPI finalized\n", rank);

    return 0;
}