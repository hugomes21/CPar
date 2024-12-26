#include "fluid_solver.h"
#include <mpi.h>
#include <cmath>
#include <cstdint>

#define IX(i, j, k) ((i) + (M + 2) * (j) + (M + 2) * (N + 2) * (k))
void SWAP(float *&x0, float *&x) { float *tmp = x0; x0 = x; x = tmp;}
float MAX(float a, float b) { return (a > b) ? a : b; }
#define LINEARSOLVERTIMES 20

// Função para inicializar o MPI
void initialize_mpi(int argc, char **argv, int &rank, int &size) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
}

// Função para finalizar o MPI
void finalize_mpi() {
    MPI_Finalize();
}

// Add sources (density or velocity) usando MPI
void add_source(int M, int N, int O, float *x, float *s, float dt, int rank, int size) {
  int local_M = M / size; // Dividir o domínio em partes iguais
  int start_i = rank * local_M + 1;
  int end_i = (rank + 1) * local_M;
  
  // Adicionar fontes localmente
  for (int k = 0; k <= O + 1; k++) {
    for (int j = 0; j <= N + 1; j++) {
      for (int i = start_i; i <= end_i; i++) {
          int idx = IX(i, j, k);
          x[idx] += dt * s[idx];
      }
    }
  }
}

// Set boundary conditions
void set_bnd(int M, int N, int O, int b, float *x, int rank, int size) {
    int i, j;
    int local_M = M / size;
    int start_i = rank * local_M + 1;
    int end_i = (rank + 1) * local_M;

    // Atualizar os limites das faces Z (k = 0 e k = O + 1)
    for (j = 1; j <= N; j++) {
        for (i = start_i; i <= end_i; i++) {
            x[IX(i, j, 0)] = (b == 3) ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
            x[IX(i, j, O + 1)] = (b == 3) ? -x[IX(i, j, O)] : x[IX(i, j, O)];
        }
    }

    // Atualizar os limites das faces X (i = 0 e i = M + 1)
    if (rank == 0) {
        for (j = 1; j <= O; j++) {
            for (i = 1; i <= N; i++) {
                x[IX(0, i, j)] = (b == 1) ? -x[IX(1, i, j)] : x[IX(1, i, j)];
                x[IX(M + 1, i, j)] = (b == 1) ? -x[IX(M, i, j)] : x[IX(M, i, j)];
            }
        }
    }

    // Atualizar os limites das faces Y (j = 0 e j = N + 1)
    for (i = start_i; i <= end_i; i++) {
        for (j = 1; j <= O; j++) {
            x[IX(i, 0, j)] = (b == 2) ? -x[IX(i, 1, j)] : x[IX(i, 1, j)];
            x[IX(i, N + 1, j)] = (b == 2) ? -x[IX(i, N, j)] : x[IX(i, N, j)];
        }
    }

    // Troca de dados entre processos
    if (rank > 0) {
        MPI_Send(x + IX(start_i, 0, 0), (N + 2) * (O + 2), MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD);
        MPI_Recv(x + IX(start_i - 1, 0, 0), (N + 2) * (O + 2), MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
    if (rank < size - 1) {
        MPI_Recv(x + IX(end_i + 1, 0, 0), (N + 2) * (O + 2), MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Send(x + IX(end_i, 0, 0), (N + 2) * (O + 2), MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD);
    }

    // Atualizar os cantos
    if (rank == 0) {
        x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
        x[IX(M + 1, 0, 0)] = 0.33f * (x[IX(M, 0, 0)] + x[IX(M + 1, 1, 0)] + x[IX(M + 1, 0, 1)]);
        x[IX(0, N + 1, 0)] = 0.33f * (x[IX(1, N + 1, 0)] + x[IX(0, N, 0)] + x[IX(0, N + 1, 1)]);
        x[IX(M + 1, N + 1, 0)] = 0.33f * (x[IX(M, N + 1, 0)] + x[IX(M + 1, N, 0)] + x[IX(M + 1, N + 1, 1)]);
    }
}

// Aux Function of lin_solve
float calculate_new_value(int i, int j, int k, float *x, float *x0, float a, float c, int M, int N, int O) {
  int M2 = M + 2;
  int N2 = N + 2;
  int MN2 = M2 * N2;

  int index = i + j * M2 + k * MN2;
  int idx_left = index - 1;
  int idx_right = index + 1;
  int idx_below = index - M2;
  int idx_above = index + M2;
  int idx_back = index - MN2;
  int idx_front = index + MN2;

  float sum_neighbors = x[idx_left] + x[idx_right] + x[idx_below] + x[idx_above] + x[idx_back] + x[idx_front];

  return (x0[index] + a * sum_neighbors) / c;
}

// Função lin_solve com MPI
void lin_solve(int M, int N, int O, int b, float *x, float *x0, float a, float c, int rank, int size) {
    float tol = 1e-7, max_c, old_x, change;
    int l = 0;

    int local_M = M / size; // Dividir o domínio em partes iguais

    do {
        max_c = 0.0f;

        // Red phase
        for (int kk = 1; kk <= O; kk++) {
            for (int jj = 1; jj <= N; jj++) {
                for (int k = kk; k <= kk; k++) {
                    for (int j = jj; j <= jj; j++) {
                        int start_i = 1 + (j + k) % 2;
                        for (int i = start_i + rank * local_M; i < (rank + 1) * local_M; i += 2) {
                            int idx = IX(i, j, k);
                            old_x = x[idx];

                            x[idx] = calculate_new_value(i, j, k, x, x0, a, c, M, N, O);
                            change = fabs(x[idx] - old_x);
                            if (change > max_c) max_c = change;
                        }
                    }
                }
            }
        }

        // Troca de dados entre processos
        MPI_Allreduce(MPI_IN_PLACE, &max_c, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

        // Black phase
        for (int kk = 1; kk <= O; kk++) {
            for (int jj = 1; jj <= N; jj++) {
                for (int k = kk; k <= kk; k++) {
                    for (int j = jj; j <= jj; j++) {
                        int start_i = 2 - (j + k) % 2;
                        for (int i = start_i + rank * local_M; i < (rank + 1) * local_M; i += 2) {
                            int idx = IX(i, j, k);
                            old_x = x[idx];

                            x[idx] = calculate_new_value(i, j, k, x, x0, a, c, M, N, O);
                            change = fabs(x[idx] - old_x);
                            if (change > max_c) max_c = change;
                        }
                    }
                }
            }
        }

        // Troca de dados entre processos
        MPI_Allreduce(MPI_IN_PLACE, &max_c, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);

        set_bnd(M, N, O, b, x, rank, size);
    } while (max_c > tol && ++l < LINEARSOLVERTIMES);
}

// Diffusion step (uses implicit method)
void diffuse(int M, int N, int O, int b, float *x, float *x0, float diff, float dt, int rank, int size) {
  int max = MAX(MAX(M, N), O);
  float a = dt * diff * max * max;
  lin_solve(M, N, O, b, x, x0, a, 1 + 6 * a, rank, size);
}

// Advection step (uses velocity field to move quantities)
void advect(int M, int N, int O, int b, float *d, float *d0, float *u, float *v, float *w, float dt, int rank, int size) {
    float dtX = dt * M, dtY = dt * N, dtZ = dt * O;

    int local_M = M / size;
    int start_i = rank * local_M + 1;
    int end_i = (rank + 1) * local_M;

    for (int k = 1; k <= O; k++) {
        for (int j = 1; j <= N; j++) {
            for (int i = start_i; i <= end_i; i++) {
                int idx = IX(i, j, k);
                float x = i - dtX * u[idx];
                float y = j - dtY * v[idx];
                float z = k - dtZ * w[idx];

                // Clamping nas fronteiras da grade
                if (x < 0.5f) x = 0.5f;
                if (x > M + 0.5f) x = M + 0.5f;
                if (y < 0.5f) y = 0.5f;
                if (y > N + 0.5f) y = N + 0.5f;
                if (z < 0.5f) z = 0.5f;
                if (z > O + 0.5f) z = O + 0.5f;

                int i0 = (int)x, i1 = i0 + 1;
                int j0 = (int)y, j1 = j0 + 1;
                int k0 = (int)z, k1 = k0 + 1;

                float s1 = x - i0, s0 = 1 - s1;
                float t1 = y - j0, t0 = 1 - t1;
                float u1 = z - k0, u0 = 1 - u1;

                d[idx] = s0 * (t0 * (u0 * d0[IX(i0, j0, k0)] + u1 * d0[IX(i0, j0, k1)]) +
                              t1 * (u0 * d0[IX(i0, j1, k0)] + u1 * d0[IX(i0, j1, k1)])) +
                        s1 * (t0 * (u0 * d0[IX(i1, j0, k0)] + u1 * d0[IX(i1, j0, k1)]) +
                              t1 * (u0 * d0[IX(i1, j1, k0)] + u1 * d0[IX(i1, j1, k1)]));
            }
        }
    }

    set_bnd(M, N, O, b, d, rank, size);
}

// Projection step to ensure incompressibility (make the velocity field divergence-free)
void project(int M, int N, int O, float *u, float *v, float *w, float *p, float *div, int rank, int size) {
  int M2 = M + 2;
  int N2 = N + 2;
  int MN2 = M2 * N2;

  // Primeira parte: cálculo de div e inicialização de p
  for (int k = 1; k <= O; k++) {
    for (int j = 1; j <= N; j++) {
      for (int i = 1 + rank * (M / size); i <= (rank + 1) * (M / size); i++) {
        int idx = i + j * M2 + k * MN2;
        div[idx] = -0.5f * (u[IX(i + 1, j, k)] - u[IX(i - 1, j, k)] +
                            v[IX(i, j + 1, k)] - v[IX(i, j - 1, k)] +
                            w[IX(i, j, k + 1)] - w[IX(i, j, k - 1)]) /
                  MAX(M, MAX(N, O));
        p[idx] = 0;
      }
    }
  }

  set_bnd(M, N, O, 0, div, rank, size);
  set_bnd(M, N, O, 0, p, rank, size);
  lin_solve(M, N, O, 0, p, div, 1, 6, rank, size);

  // Segunda parte: atualização de u, v e w
  for (int k = 1; k <= O; k++) {
    for (int j = 1; j <= N; j++) {
      for (int i = 1 + rank * (M / size); i <= (rank + 1) * (M / size); i++) {
        int idx = i + j * M2 + k * MN2;
        u[idx] -= 0.5f * (p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]);
        v[idx] -= 0.5f * (p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]);
        w[idx] -= 0.5f * (p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]);
      }
    }
  }

  set_bnd(M, N, O, 1, u, rank, size);
  set_bnd(M, N, O, 2, v, rank, size);
  set_bnd(M, N, O, 3, w, rank, size);
}



// Step function for density
void dens_step(int M, int N, int O, float *x, float *x0, float *u, float *v, float *w, float diff, float dt) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    add_source(M, N, O, x, x0, dt, rank, size);
    SWAP(x0, x);
    diffuse(M, N, O, 0, x, x0, diff, dt, rank, size);
    SWAP(x0, x);
    advect(M, N, O, 0, x, x0, u, v, w, dt, rank, size);
}

// Step function for velocity
void vel_step(int M, int N, int O, float *u, float *v, float *w, float *u0, float *v0, float *w0, float visc, float dt) {
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    add_source(M, N, O, u, u0, dt, rank, size);
    add_source(M, N, O, v, v0, dt, rank, size);
    add_source(M, N, O, w, w0, dt, rank, size);
    SWAP(u0, u);
    diffuse(M, N, O, 1, u, u0, visc, dt, rank, size);
    SWAP(v0, v);
    diffuse(M, N, O, 2, v, v0, visc, dt, rank, size);
    SWAP(w0, w);
    diffuse(M, N, O, 3, w, w0, visc, dt, rank, size);
    project(M, N, O, u, v, w, u0, v0, rank, size);
    SWAP(u0, u);
    SWAP(v0, v);
    SWAP(w0, w);
    advect(M, N, O, 1, u, u0, u0, v0, w0, dt, rank, size);
    advect(M, N, O, 2, v, v0, u0, v0, w0, dt, rank, size);
    advect(M, N, O, 3, w, w0, u0, v0, w0, dt, rank, size);
    project(M, N, O, u, v, w, u0, v0, rank, size);
}