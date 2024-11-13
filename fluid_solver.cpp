#include "fluid_solver.h"
#include <omp.h>
#include <cmath>

#define IX(i, j, k) ((i) + (M + 2) * (j) + (M + 2) * (N + 2) * (k))
inline void SWAP(float *&x0, float *&x)                                        
  {                                                                           
    float *tmp = x0;                                                           
    x0 = x;                                                                    
    x = tmp;                                                                   
  }
// Usar uma função inline em vez de uma macro para o MAX
inline float MAX(float a, float b) {
  return (a > b) ? a : b;
}
#define LINEARSOLVERTIMES 20

// Add sources (density or velocity)
void add_source(int M, int N, int O, float *x, float *s, float dt) {
  int size = (M + 2) * (N + 2) * (O + 2);
  int i = 0;

  // Unroll the loop by a factor of 16
  for (; i <= size - 16; i += 16) {
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
  for (; i < size; i++) {
    x[i] += dt * s[i];
  }
}

// Set boundary conditions
void set_bnd(int M, int N, int O, int b, float *x) {
  int i, j;

  // Atualizar os limites das faces Z (k = 0 e k = O + 1)
  for (j = 1; j <= N; j++) {
    for (i = 1; i <= M; i++) {
      x[IX(i, j, 0)] = (b == 3) ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
      x[IX(i, j, O + 1)] = (b == 3) ? -x[IX(i, j, O)] : x[IX(i, j, O)];
    }
  }

  // Atualizar os limites das faces X (i = 0 e i = M + 1)
  for (j = 1; j <= O; j++) {
    for (i = 1; i <= N; i++) {
      x[IX(0, i, j)] = (b == 1) ? -x[IX(1, i, j)] : x[IX(1, i, j)];
      x[IX(M + 1, i, j)] = (b == 1) ? -x[IX(M, i, j)] : x[IX(M, i, j)];
    }
  }

  // Atualizar os limites das faces Y (j = 0 e j = N + 1)
  for (i = 1; i <= M; i++) {
    for (j = 1; j <= O; j++) {
      x[IX(i, 0, j)] = (b == 2) ? -x[IX(i, 1, j)] : x[IX(i, 1, j)];
      x[IX(i, N + 1, j)] = (b == 2) ? -x[IX(i, N, j)] : x[IX(i, N, j)];
    }
  }

  // Atualizar os cantos
  x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
  x[IX(M + 1, 0, 0)] = 0.33f * (x[IX(M, 0, 0)] + x[IX(M + 1, 1, 0)] + x[IX(M + 1, 0, 1)]);
  x[IX(0, N + 1, 0)] = 0.33f * (x[IX(1, N + 1, 0)] + x[IX(0, N, 0)] + x[IX(0, N + 1, 1)]);
  x[IX(M + 1, N + 1, 0)] = 0.33f * (x[IX(M, N + 1, 0)] + x[IX(M + 1, N, 0)] + x[IX(M + 1, N + 1, 1)]);
}


/*
inline float calculate_new_value(int i, int j, int k, float *x, float *x0, float a, float c, int M, int N, int O) {
  int M2 = M + 2;
  int N2 = N + 2;
  int MN2 = M2 * N2;

  // Calcular o índice da posição atual
  int index = i + j * M2 + k * MN2;

  // Acessar os vizinhos diretamente
  float left = x[index - 1];
  float right = x[index + 1];
  float below = x[index - M2];
  float above = x[index + M2];
  float back = x[index - MN2];
  float front = x[index + MN2];

  // Calcular a soma dos vizinhos
  float sum_neighbors = left + right + below + above + back + front;

  // Calcular e retornar o novo valor
  return (x0[index] + a * sum_neighbors) / c;
}

void lin_solve(int M, int N, int O, int b, float *x, float *x0, float a, float c) {
  int M2 = M + 2;
  int N2 = N + 2;
  int MN2 = M2 * N2;

  int block_size_i = 8;  
  int block_size_j = 8;  
  int block_size_k = 8;  

  for (int l = 0; l < LINEARSOLVERTIMES; l++) {

    // Loop over blocks
    for (int kk = 1; kk <= O; kk += block_size_k) {
      for (int jj = 1; jj <= N; jj += block_size_j) {
        for (int ii = 1; ii <= M; ii += block_size_i) {

          // Dentro de cada bloco, iterar sobre os elementos
          for (int k = kk; k < kk + block_size_k && k <= O; k++) {
            for (int j = jj; j < jj + block_size_j && j <= N; j++) {
              for (int i = ii; i < ii + block_size_i && i <= M; i+=4) {
                x[i + j * M2 + k * MN2] = calculate_new_value(i, j, k, x, x0, a, c, M, N, O);
                x[i + 1 + j * M2 + k * MN2] = calculate_new_value(i + 1, j, k, x, x0, a, c, M, N, O);
                x[i + 2 + j * M2 + k * MN2] = calculate_new_value(i + 2, j, k, x, x0, a, c, M, N, O);
                x[i + 3 + j * M2 + k * MN2] = calculate_new_value(i + 3, j, k, x, x0, a, c, M, N, O);
              }
            }
          }
        }
      }
    }
    set_bnd(M, N, O, b, x);
  }
}
*/


// red-black solver with convergence check
void lin_solve(int M, int N, int O, int b, float *x, float *x0, float a, float c) {
    float tol = 1e-7, max_c, old_x, change;
    int l = 0, blockSize = 8;

    do {
        max_c = 0.0f;

        // Primeira fase: Células vermelhas
        #pragma omp parallel for collapse(3) reduction(max:max_c) private(old_x, change)
        for (int kk = 1; kk <= O; kk += blockSize) {
            for (int jj = 1; jj <= N; jj += blockSize) {
                for (int ii = 1; ii <= M; ii += blockSize) {
                    for (int k = kk; k < kk + blockSize && k <= O; k++) {
                        for (int j = jj; j < jj + blockSize && j <= N; j++) {
                            for (int i = ii + (j + k) % 2; i < ii + blockSize && i <= M; i += 2) {
                                old_x = x[IX(i, j, k)];
                                x[IX(i, j, k)] = (x0[IX(i, j, k)] +
                                                  a * (x[IX(i - 1, j, k)] + x[IX(i + 1, j, k)] +
                                                       x[IX(i, j - 1, k)] + x[IX(i, j + 1, k)] +
                                                       x[IX(i, j, k - 1)] + x[IX(i, j, k + 1)])) / c;
                                change = fabs(x[IX(i, j, k)] - old_x);
                                if (change > max_c) max_c = change;
                            }
                        }
                    }
                }
            }
        }

        // Segunda fase: Células pretas
        #pragma omp parallel for collapse(3) reduction(max:max_c) private(old_x, change)
        for (int kk = 1; kk <= O; kk += blockSize) {
            for (int jj = 1; jj <= N; jj += blockSize) {
                for (int ii = 1; ii <= M; ii += blockSize) {
                    for (int k = kk; k < kk + blockSize && k <= O; k++) {
                        for (int j = jj; j < jj + blockSize && j <= N; j++) {
                            for (int i = ii + (j + k + 1) % 2; i < ii + blockSize && i <= M; i += 2) {
                                old_x = x[IX(i, j, k)];
                                x[IX(i, j, k)] = (x0[IX(i, j, k)] +
                                                  a * (x[IX(i - 1, j, k)] + x[IX(i + 1, j, k)] +
                                                       x[IX(i, j - 1, k)] + x[IX(i, j + 1, k)] +
                                                       x[IX(i, j, k - 1)] + x[IX(i, j, k + 1)])) / c;
                                change = fabs(x[IX(i, j, k)] - old_x);
                                if (change > max_c) max_c = change;
                            }
                        }
                    }
                }
            }
        }

        set_bnd(M, N, O, b, x);
    } while (max_c > tol && ++l < 20);
}




// Diffusion step (uses implicit method)
void diffuse(int M, int N, int O, int b, float *x, float *x0, float diff, float dt) {
  int max = MAX(MAX(M, N), O);
  float a = dt * diff * max * max;
  lin_solve(M, N, O, b, x, x0, a, 1 + 6 * a);
}

// Advection step (uses velocity field to move quantities)
void advect(int M, int N, int O, int b, float *d, float *d0, float *u, float *v, float *w, float dt) {
  float dtX = dt * M, dtY = dt * N, dtZ = dt * O;

  // Definir tamanhos dos blocos
  int block_size_i = 8;
  int block_size_j = 8;
  int block_size_k = 8;

  for (int kk = 1; kk <= O; kk += block_size_k) {
    for (int jj = 1; jj <= N; jj += block_size_j) {
      for (int ii = 1; ii <= M; ii += block_size_i) {
        // Dentro de cada bloco, iterar sobre os elementos
        for (int k = kk; k < kk + block_size_k && k <= O; k++) {
          for (int j = jj; j < jj + block_size_j && j <= N; j++) {
            for (int i = ii; i < ii + block_size_i && i <= M; i++) {
              int idx = IX(i, j, k);
              float x = i - dtX * u[idx];
              float y = j - dtY * v[idx];
              float z = k - dtZ * w[idx];

              // Clamping nas fronteiras da grade
              if (x < 0.5f)
                x = 0.5f;
              if (x > M + 0.5f)
                x = M + 0.5f;
              if (y < 0.5f)
                y = 0.5f;
              if (y > N + 0.5f)
                y = N + 0.5f;
              if (z < 0.5f)
                z = 0.5f;
              if (z > O + 0.5f)
                z = O + 0.5f;

              int i0 = (int)x, i1 = i0 + 1;
              int j0 = (int)y, j1 = j0 + 1;
              int k0 = (int)z, k1 = k0 + 1;

              float s1 = x - i0, s0 = 1 - s1;
              float t1 = y - j0, t0 = 1 - t1;
              float u1 = z - k0, u0 = 1 - u1;

              d[idx] =
                  s0 * (t0 * (u0 * d0[IX(i0, j0, k0)] + u1 * d0[IX(i0, j0, k1)]) +
                        t1 * (u0 * d0[IX(i0, j1, k0)] + u1 * d0[IX(i0, j1, k1)])) +
                  s1 * (t0 * (u0 * d0[IX(i1, j0, k0)] + u1 * d0[IX(i1, j0, k1)]) +
                        t1 * (u0 * d0[IX(i1, j1, k0)] + u1 * d0[IX(i1, j1, k1)]));
            }
          }
        }
      }
    }
  }
  set_bnd(M, N, O, b, d);
}


// Projection step to ensure incompressibility (make the velocity field
// divergence-free)
void project(int M, int N, int O, float *u, float *v, float *w, float *p, float *div) {
  int M2 = M + 2;
  int N2 = N + 2;
  int MN2 = M2 * N2;

  int block_size_i = 8;
  int block_size_j = 8;
  int block_size_k = 8;

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
  lin_solve(M, N, O, 0, p, div, 1, 6);

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
}




// Step function for density
void dens_step(int M, int N, int O, float *x, float *x0, float *u, float *v,
               float *w, float diff, float dt) {
  add_source(M, N, O, x, x0, dt);
  SWAP(x0, x);
  diffuse(M, N, O, 0, x, x0, diff, dt);
  SWAP(x0, x);
  advect(M, N, O, 0, x, x0, u, v, w, dt);
}

// Step function for velocity
void vel_step(int M, int N, int O, float *u, float *v, float *w, float *u0,
              float *v0, float *w0, float visc, float dt) {
  add_source(M, N, O, u, u0, dt);
  add_source(M, N, O, v, v0, dt);
  add_source(M, N, O, w, w0, dt);
  SWAP(u0, u);
  diffuse(M, N, O, 1, u, u0, visc, dt);
  SWAP(v0, v);
  diffuse(M, N, O, 2, v, v0, visc, dt);
  SWAP(w0, w);
  diffuse(M, N, O, 3, w, w0, visc, dt);
  project(M, N, O, u, v, w, u0, v0);
  SWAP(u0, u);
  SWAP(v0, v);
  SWAP(w0, w);
  advect(M, N, O, 1, u, u0, u0, v0, w0, dt);
  advect(M, N, O, 2, v, v0, u0, v0, w0, dt);
  advect(M, N, O, 3, w, w0, u0, v0, w0, dt);
  project(M, N, O, u, v, w, u0, v0);
}
