#include "fluid_solver_mpi.h"
#include <omp.h>
#include <mpi.h>
#include <cmath>
#include <cstdint>
#include <cstdio>

#define IX(i, j, k) ((i) + (M + 2) * (j) + (M + 2) * (N + 2) * (k))
void SWAP(float *&x0, float *&x) { float *tmp = x0; x0 = x; x = tmp;}
float MAX(float a, float b) { return (a > b) ? a : b; }
#define LINEARSOLVERTIMES 20

// Add sources (density or velocity)
void add_source(int M, int N, int O, float *x, float *s, float dt) {
  int size = (M + 2) * (N + 2) * (O + 2);

  // Unroll the loop by a factor of 16
  for (int i = 0; i <= size - 16; i += 16) {
    x[i] += dt * s[i];
    x[i + 1] += dt * s[i + 1];
    x[i + 2] += dt * s[i + 2];
    x[i + 3] += dt * s[i + 3];
    x[i + 4] += dt * s[i + 4];
    x[i + 5] += dt * s[i + 5];
    x[i + 6] += dt * s[i + 6];
    x[i + 7] += dt * s[i + 7];
    x[i + 8] += dt * s[i + 8];
    x[i + 9] += dt * s[i + 9];
    x[i + 10] += dt * s[i + 10];
    x[i + 11] += dt * s[i + 11];
    x[i + 12] += dt * s[i + 12];
    x[i + 13] += dt * s[i + 13];
    x[i + 14] += dt * s[i + 14];
    x[i + 15] += dt * s[i + 15];
  }

  // Handle the remaining elements
  for (int i = size - (size % 16); i < size; i++) {
    x[i] += dt * s[i];
  }
}

// Set boundary conditions
void set_bnd(int M, int N, int O, int b, float *x) {
  int i, j;

  // Atualizar os limites das faces Z (k = 0 e k = O + 1)
  for (j = 1; j <= N; j++) {
    for (i = 1; i <= M; i += 4) { // Loop Unrolling com passo de 4
      x[IX(i, j, 0)] = (b == 3) ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
      x[IX(i + 1, j, 0)] = (b == 3) ? -x[IX(i + 1, j, 1)] : x[IX(i + 1, j, 1)];
      x[IX(i + 2, j, 0)] = (b == 3) ? -x[IX(i + 2, j, 1)] : x[IX(i + 2, j, 1)];
      x[IX(i + 3, j, 0)] = (b == 3) ? -x[IX(i + 3, j, 1)] : x[IX(i + 3, j, 1)];

      x[IX(i, j, O + 1)] = (b == 3) ? -x[IX(i, j, O)] : x[IX(i, j, O)];
      x[IX(i + 1, j, O + 1)] = (b == 3) ? -x[IX(i + 1, j, O)] : x[IX(i + 1, j, O)];
      x[IX(i + 2, j, O + 1)] = (b == 3) ? -x[IX(i + 2, j, O)] : x[IX(i + 2, j, O)];
      x[IX(i + 3, j, O + 1)] = (b == 3) ? -x[IX(i + 3, j, O)] : x[IX(i + 3, j, O)];
    }
  }

  // Atualizar os limites das faces X (i = 0 e i = M + 1)
  for (j = 1; j <= O; j++) {
    for (i = 1; i <= N; i += 2) {
      x[IX(0, i, j)] = (b == 1) ? -x[IX(1, i, j)] : x[IX(1, i, j)];
      x[IX(0, i + 1, j)] = (b == 1) ? -x[IX(1, i + 1, j)] : x[IX(1, i + 1, j)];

      x[IX(M + 1, i, j)] = (b == 1) ? -x[IX(M, i, j)] : x[IX(M, i, j)];
      x[IX(M + 1, i + 1, j)] = (b == 1) ? -x[IX(M, i + 1, j)] : x[IX(M, i + 1, j)];
    }
  }


  // Atualizar os limites das faces Y (j = 0 e j = N + 1)
  for (i = 1; i <= M; i++) {
    for (j = 1; j <= O; j += 2) {
      x[IX(i, 0, j)] = (b == 2) ? -x[IX(i, 1, j)] : x[IX(i, 1, j)];
      x[IX(i, 0, j + 1)] = (b == 2) ? -x[IX(i, 1, j + 1)] : x[IX(i, 1, j + 1)];

      x[IX(i, N + 1, j)] = (b == 2) ? -x[IX(i, N, j)] : x[IX(i, N, j)];
      x[IX(i, N + 1, j + 1)] = (b == 2) ? -x[IX(i, N, j + 1)] : x[IX(i, N, j + 1)];
    }
  }



  // Atualizar os cantos
  x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
  x[IX(M + 1, 0, 0)] = 0.33f * (x[IX(M, 0, 0)] + x[IX(M + 1, 1, 0)] + x[IX(M + 1, 0, 1)]);
  x[IX(0, N + 1, 0)] = 0.33f * (x[IX(1, N + 1, 0)] + x[IX(0, N, 0)] + x[IX(0, N + 1, 1)]);
  x[IX(M + 1, N + 1, 0)] = 0.33f * (x[IX(M, N + 1, 0)] + x[IX(M + 1, N, 0)] + x[IX(M + 1, N + 1, 1)]);
}

// Aux Function of lin_solve
float calculate_new_value(int i, int j, int k, float *x, float *x0, float a, float c, int M, int N, int O) {
  int M2 = M + 2;
  int N2 = N + 2;
  int MN2 = M2 * N2;

  int index = i + j * M2 + k * MN2;

  if (fabs(c) < 1e-6) {
    printf("Error: Division by near-zero value at (%d, %d, %d)\n", i, j, k);
    return x0[index]; // Or handle appropriately
  }

  int idx_left = (i > 1) ? index - 1 : index;
  int idx_right = (i < M) ? index + 1 : index;
  int idx_below = (j > 1) ? index - M2 : index;
  int idx_above = (j < N) ? index + M2 : index;
  int idx_back = (k > 1) ? index - MN2 : index;
  int idx_front = (k < O) ? index + MN2 : index;

  float sum_neighbors = x[idx_left] + x[idx_right] + x[idx_below] + x[idx_above] + x[idx_back] + x[idx_front];
  float new_value = (x0[index] + a * sum_neighbors) / (c + 1e-7); 

  if (fabs(new_value - x[index]) < 1e-8) return x[index];

  return new_value;
}

int distribute_workload(int N, int size, int rank, int &start_index, int &end_index) {
    // Dividir o intervalo do eixo j entre os ranks
    int rows_per_process = N / size;
    int remaining_rows = N % size;

    // Calcular start_index e end_index para o rank atual
    start_index = rank * rows_per_process + std::min(rank, remaining_rows);
    end_index = start_index + rows_per_process - 1 + (rank < remaining_rows ? 1 : 0);

    // Garantir que os índices estão dentro dos limites
    if (start_index > end_index || start_index < 0 || end_index >= N) {
        printf("Erro: Índices inválidos! Rank %d tem start_index=%d e end_index=%d\n",
               rank, start_index, end_index);
        start_index = -1;
        end_index = -1;
    } else {
        printf("Rank %d: start_index=%d, end_index=%d\n", rank, start_index, end_index);
    }

    // Verificação e retorno da quantidade de trabalho atribuída ao rank
    return (start_index != -1 && end_index != -1) ? end_index - start_index + 1 : 0;
}

//void distribute_workload(int N, int size, int rank, int &start_index, int &end_index) {
//    if (size == 0) {
//        printf("Error: Division by zero in distribute_workload. Variable size is zero.\n");
//        MPI_Abort(MPI_COMM_WORLD, 1);
//    }
//
//    int base_chunk = N / size; // Verifique se size não é zero
//    int remainder = N % size;
//
//    start_index = rank * base_chunk + std::min(rank, remainder) + 1;
//    end_index = start_index + base_chunk - 1 + (rank < remainder ? 1 : 0);
//}

// Função para trocar fronteiras entre processos
void exchange_boundaries(float *x, int M, int N, int O, int start_index, int end_index, int rank, int size) {
    MPI_Request requests[4];
    MPI_Status statuses[4];

    int slice_size = (M + 2) * (O + 2); // Dimensão de uma fatia (XY plano)

    // Inicializar requests como nulos
    requests[0] = requests[1] = requests[2] = requests[3] = MPI_REQUEST_NULL;

    // Top boundary
    if (rank > 0) {
        MPI_Irecv(&x[IX(1, start_index - 1, 1)], slice_size, MPI_FLOAT, rank - 1, 1, MPI_COMM_WORLD, &requests[1]);
        MPI_Isend(&x[IX(1, start_index, 1)], slice_size, MPI_FLOAT, rank - 1, 0, MPI_COMM_WORLD, &requests[0]);
    }

    // Bottom boundary
    if (rank < size - 1) {
        MPI_Irecv(&x[IX(1, end_index + 1, 1)], slice_size, MPI_FLOAT, rank + 1, 0, MPI_COMM_WORLD, &requests[3]);
        MPI_Isend(&x[IX(1, end_index, 1)], slice_size, MPI_FLOAT, rank + 1, 1, MPI_COMM_WORLD, &requests[2]);
    }

    // Esperar todas as comunicações serem finalizadas
    MPI_Waitall(4, requests, statuses);
    MPI_Barrier(MPI_COMM_WORLD);
}



void lin_solve(int M, int N, int O, int b, float *x, float *x0, float a, float c, int rank, int size) {
    float tol = 1e-7, max_c, global_max_c, old_x, change;
    int l = 0, blockSize = 8;

    // Divisão do domínio por processo
    int start_index, end_index;
    distribute_workload(N, size, rank, start_index, end_index);

    do {
        max_c = 0.0f;

        // Red phase
        for (int kk = 1; kk <= O; kk += blockSize) {
            for (int jj = start_index; jj < end_index; jj += blockSize) {
                int k_max = std::min(kk + blockSize, O + 1);
                int j_max = std::min(jj + blockSize, N + 1);

                for (int k = kk; k < k_max; k++) {
                    for (int j = jj; j < j_max; j++) {
                        int start_i_red = 1 + (j + k) % 2;
                        for (int i = start_i_red; i <= M; i += 2) {
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

        exchange_boundaries(x, M, N, O, start_index, end_index, rank, size);

        // Black phase
        for (int kk = 1; kk <= O; kk += blockSize) {
            for (int jj = start_index; jj < end_index; jj += blockSize) {
                int k_max = std::min(kk + blockSize, O + 1);
                int j_max = std::min(jj + blockSize, N + 1);

                for (int k = kk; k < k_max; k++) {
                    for (int j = jj; j < j_max; j++) {
                        int start_i_black = 2 - (j + k) % 2;
                        for (int i = start_i_black; i <= M; i += 2) {
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

        exchange_boundaries(x, M, N, O, start_index, end_index, rank, size);

        // Redução global do max_c
        MPI_Allreduce(&max_c, &global_max_c, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
        set_bnd(M, N, O, b, x);

        if (rank == 0) {
            printf("Iteration %d, max_c = %f\n", l, global_max_c);
        } else {
            printf("Rank %d: Iteration %d, max_c = %f\n", rank, l, global_max_c);
        }
    } while (global_max_c > tol && ++l < 20);
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

    // Troca de fronteiras antes da advecção
    exchange_boundaries(d0, M, N, O, 1, N, rank, size);

    int block_size_i = 8;
    int block_size_j = 8;
    int block_size_k = 8;

    for (int kk = 1; kk <= O; kk += block_size_k) {
        for (int jj = 1; jj <= N; jj += block_size_j) {
            for (int ii = 1; ii <= M; ii += block_size_i) {
                // Iterar sobre os elementos dentro de cada bloco

                for (int k = kk; k < kk + block_size_k && k <= O; k++) {
                    for (int j = jj; j < jj + block_size_j && j <= N; j++) {
                        for (int i = ii; i < ii + block_size_i && i <= M; i++) {
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
            }
        }
    }

    // Troca de fronteiras após a advecção
    exchange_boundaries(d, M, N, O, 1, N, rank, size);
    set_bnd(M, N, O, b, d);
}


// Projection step to ensure incompressibility (make the velocity field divergence-free)
void project(int M, int N, int O, float *u, float *v, float *w, float *p, float *div, int rank, int size) {
  int M2 = M + 2;
  int N2 = N + 2;
  int MN2 = M2 * N2;

  int block_size_i = 8;
  int block_size_j = 8;
  int block_size_k = 8;

  // Troca de fronteiras antes da projeção
  exchange_boundaries(u, M, N, O, 1, N, rank, size);
  exchange_boundaries(v, M, N, O, 1, N, rank, size);
  exchange_boundaries(w, M, N, O, 1, N, rank, size);

  // Primeira parte: cálculo de div e inicialização de p com loop unrolling
  for (int kk = 1; kk <= O; kk += block_size_k) {
    for (int jj = 1; jj <= N; jj += block_size_j) {
      for (int ii = 1; ii <= M; ii += block_size_i) {

        for (int k = kk; k < kk + block_size_k && k <= O; k++) {
          for (int j = jj; j < jj + block_size_j && j <= N; j++) {
            int i;
            // Loop unrolling por um fator de 8
            for (i = ii; i <= M - 7 && i < ii + block_size_i; i += 8) {
              for (int offset = 0; offset < 8; ++offset) {
                int idx = (i + offset) + j * M2 + k * MN2;
                div[idx] = -0.5f * (u[(i + offset + 1) + j * M2 + k * MN2] - u[(i + offset - 1) + j * M2 + k * MN2] +
                                    v[(i + offset) + (j + 1) * M2 + k * MN2] - v[(i + offset) + (j - 1) * M2 + k * MN2] +
                                    w[(i + offset) + j * M2 + (k + 1) * MN2] - w[(i + offset) + j * M2 + (k - 1) * MN2]) 
                          / MAX(M, MAX(N, O));
                p[idx] = 0;
              }
            }
            // Elementos restantes
            for (; i <= M && i < ii + block_size_i; i++) {
              int idx = i + j * M2 + k * MN2;
              div[idx] = -0.5f * (u[(i + 1) + j * M2 + k * MN2] - u[(i - 1) + j * M2 + k * MN2] +
                                  v[i + (j + 1) * M2 + k * MN2] - v[i + (j - 1) * M2 + k * MN2] +
                                  w[i + j * M2 + (k + 1) * MN2] - w[i + j * M2 + (k - 1) * MN2]) 
                        / MAX(M, MAX(N, O));
              p[idx] = 0;
            }
          }
        }
      }
    }
  }

  set_bnd(M, N, O, 0, div);
  set_bnd(M, N, O, 0, p);
  lin_solve(M, N, O, 0, p, div, 1, 6, rank, size);

  // Segunda parte: atualização de u, v e w com loop unrolling
  for (int kk = 1; kk <= O; kk += block_size_k) {
    for (int jj = 1; jj <= N; jj += block_size_j) {
      for (int ii = 1; ii <= M; ii += block_size_i) {
        for (int k = kk; k < kk + block_size_k && k <= O; k++) {
          for (int j = jj; j < jj + block_size_j && j <= N; j++) {
            int i;
            // Loop unrolling por um fator de 8
            for (i = ii; i <= M - 7 && i < ii + block_size_i; i += 8) {
              for (int offset = 0; offset < 8; ++offset) {
                int idx = (i + offset) + j * M2 + k * MN2;
                u[idx] -= 0.5f * (p[(i + offset + 1) + j * M2 + k * MN2] - p[(i + offset - 1) + j * M2 + k * MN2]);
                v[idx] -= 0.5f * (p[(i + offset) + (j + 1) * M2 + k * MN2] - p[(i + offset) + (j - 1) * M2 + k * MN2]);
                w[idx] -= 0.5f * (p[(i + offset) + j * M2 + (k + 1) * MN2] - p[(i + offset) + j * M2 + (k - 1) * MN2]);
              }
            }
            // Elementos restantes
            for (; i <= M && i < ii + block_size_i; i++) {
              int idx = i + j * M2 + k * MN2;
              u[idx] -= 0.5f * (p[(i + 1) + j * M2 + k * MN2] - p[(i - 1) + j * M2 + k * MN2]);
              v[idx] -= 0.5f * (p[i + (j + 1) * M2 + k * MN2] - p[i + (j - 1) * M2 + k * MN2]);
              w[idx] -= 0.5f * (p[i + j * M2 + (k + 1) * MN2] - p[i + j * M2 + (k - 1) * MN2]);
            }
          }
        }
      }
    }
  }

  set_bnd(M, N, O, 1, u);
  set_bnd(M, N, O, 2, v);
  set_bnd(M, N, O, 3, w);

  // Troca de fronteiras após a projeção
  exchange_boundaries(u, M, N, O, 1, N, rank, size);
  exchange_boundaries(v, M, N, O, 1, N, rank, size);
  exchange_boundaries(w, M, N, O, 1, N, rank, size);
}


// Step function for density
void dens_step(int M, int N, int O, float *x, float *x0, float *u, float *v, float *w, float diff, float dt, int rank, int size) {
    add_source(M, N, O, x, x0, dt);
    SWAP(x0, x);
    diffuse(M, N, O, 0, x, x0, diff, dt, rank, size);
    SWAP(x0, x);
    advect(M, N, O, 0, x, x0, u, v, w, dt, rank, size);
}


// Step function for velocity
void vel_step(int M, int N, int O, float *u, float *v, float *w, float *u0, float *v0, float *w0, float visc, float dt, int rank, int size) {
  add_source(M, N, O, u, u0, dt);
  add_source(M, N, O, v, v0, dt);
  add_source(M, N, O, w, w0, dt);
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