#include "fluid_solver.h"
#include <cmath>

#define IX(i, j, k) ((i) + (M + 2) * (j) + (M + 2) * (N + 2) * (k))
#define SWAP(x0, x)                                                            \
  {                                                                            \
    float *tmp = x0;                                                           \
    x0 = x;                                                                    \
    x = tmp;                                                                   \
  }
#define MAX(a, b) (((a) > (b)) ? (a) : (b))
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
  int i = 0, j = 0;

  // Set boundary on faces com loop unrolling fator 16
  for (i = 1; i <= M; i++) {
    for (j = 1; j <= N - 15; j += 16) {
      // Processar 16 elementos consecutivos
      x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
      x[IX(i, j + 1, 0)] = b == 3 ? -x[IX(i, j + 1, 1)] : x[IX(i, j + 1, 1)];
      x[IX(i, j + 2, 0)] = b == 3 ? -x[IX(i, j + 2, 1)] : x[IX(i, j + 2, 1)];
      x[IX(i, j + 3, 0)] = b == 3 ? -x[IX(i, j + 3, 1)] : x[IX(i, j + 3, 1)];
      x[IX(i, j + 4, 0)] = b == 3 ? -x[IX(i, j + 4, 1)] : x[IX(i, j + 4, 1)];
      x[IX(i, j + 5, 0)] = b == 3 ? -x[IX(i, j + 5, 1)] : x[IX(i, j + 5, 1)];
      x[IX(i, j + 6, 0)] = b == 3 ? -x[IX(i, j + 6, 1)] : x[IX(i, j + 6, 1)];
      x[IX(i, j + 7, 0)] = b == 3 ? -x[IX(i, j + 7, 1)] : x[IX(i, j + 7, 1)];
      x[IX(i, j + 8, 0)] = b == 3 ? -x[IX(i, j + 8, 1)] : x[IX(i, j + 8, 1)];
      x[IX(i, j + 9, 0)] = b == 3 ? -x[IX(i, j + 9, 1)] : x[IX(i, j + 9, 1)];
      x[IX(i, j + 10, 0)] = b == 3 ? -x[IX(i, j + 10, 1)] : x[IX(i, j + 10, 1)];
      x[IX(i, j + 11, 0)] = b == 3 ? -x[IX(i, j + 11, 1)] : x[IX(i, j + 11, 1)];
      x[IX(i, j + 12, 0)] = b == 3 ? -x[IX(i, j + 12, 1)] : x[IX(i, j + 12, 1)];
      x[IX(i, j + 13, 0)] = b == 3 ? -x[IX(i, j + 13, 1)] : x[IX(i, j + 13, 1)];
      x[IX(i, j + 14, 0)] = b == 3 ? -x[IX(i, j + 14, 1)] : x[IX(i, j + 14, 1)];
      x[IX(i, j + 15, 0)] = b == 3 ? -x[IX(i, j + 15, 1)] : x[IX(i, j + 15, 1)];
    }
    // Processar elementos restantes
    for (; j <= N; j++) {
      x[IX(i, j, 0)] = b == 3 ? -x[IX(i, j, 1)] : x[IX(i, j, 1)];
    }
  }

  // Aplicar estratégia semelhante aos outros loops
  for (i = 1; i <= N; i++) {
    for (j = 1; j <= O - 15; j += 16) {
      x[IX(0, i, j)] = b == 1 ? -x[IX(1, i, j)] : x[IX(1, i, j)];
      x[IX(0, i, j + 1)] = b == 1 ? -x[IX(1, i, j + 1)] : x[IX(1, i, j + 1)];
      x[IX(0, i, j + 2)] = b == 1 ? -x[IX(1, i, j + 2)] : x[IX(1, i, j + 2)];
      x[IX(0, i, j + 3)] = b == 1 ? -x[IX(1, i, j + 3)] : x[IX(1, i, j + 3)];
      x[IX(0, i, j + 4)] = b == 1 ? -x[IX(1, i, j + 4)] : x[IX(1, i, j + 4)];
      x[IX(0, i, j + 5)] = b == 1 ? -x[IX(1, i, j + 5)] : x[IX(1, i, j + 5)];
      x[IX(0, i, j + 6)] = b == 1 ? -x[IX(1, i, j + 6)] : x[IX(1, i, j + 6)];
      x[IX(0, i, j + 7)] = b == 1 ? -x[IX(1, i, j + 7)] : x[IX(1, i, j + 7)];
      x[IX(0, i, j + 8)] = b == 1 ? -x[IX(1, i, j + 8)] : x[IX(1, i, j + 8)];
      x[IX(0, i, j + 9)] = b == 1 ? -x[IX(1, i, j + 9)] : x[IX(1, i, j + 9)];
      x[IX(0, i, j + 10)] = b == 1 ? -x[IX(1, i, j + 10)] : x[IX(1, i, j + 10)];
      x[IX(0, i, j + 11)] = b == 1 ? -x[IX(1, i, j + 11)] : x[IX(1, i, j + 11)];
      x[IX(0, i, j + 12)] = b == 1 ? -x[IX(1, i, j + 12)] : x[IX(1, i, j + 12)];
      x[IX(0, i, j + 13)] = b == 1 ? -x[IX(1, i, j + 13)] : x[IX(1, i, j + 13)];
      x[IX(0, i, j + 14)] = b == 1 ? -x[IX(1, i, j + 14)] : x[IX(1, i, j + 14)];
      x[IX(0, i, j + 15)] = b == 1 ? -x[IX(1, i, j + 15)] : x[IX(1, i, j + 15)];
    }
  }

  // Faça o mesmo para o terceiro loop
  for (i = 1; i <= M; i++) {
    for (j = 1; j <= O - 15; j += 16) {
      x[IX(i, 0, j)] = b == 2 ? -x[IX(i, 1, j)] : x[IX(i, 1, j)];
      x[IX(i, 0, j + 1)] = b == 2 ? -x[IX(i, 1, j + 1)] : x[IX(i, 1, j + 1)];
      x[IX(i, 0, j + 2)] = b == 2 ? -x[IX(i, 1, j + 2)] : x[IX(i, 1, j + 2)];
      x[IX(i, 0, j + 3)] = b == 2 ? -x[IX(i, 1, j + 3)] : x[IX(i, 1, j + 3)];
      x[IX(i, 0, j + 4)] = b == 2 ? -x[IX(i, 1, j + 4)] : x[IX(i, 1, j + 4)];
      x[IX(i, 0, j + 5)] = b == 2 ? -x[IX(i, 1, j + 5)] : x[IX(i, 1, j + 5)];
      x[IX(i, 0, j + 6)] = b == 2 ? -x[IX(i, 1, j + 6)] : x[IX(i, 1, j + 6)];
      x[IX(i, 0, j + 7)] = b == 2 ? -x[IX(i, 1, j + 7)] : x[IX(i, 1, j + 7)];
      x[IX(i, 0, j + 8)] = b == 2 ? -x[IX(i, 1, j + 8)] : x[IX(i, 1, j + 8)];
      x[IX(i, 0, j + 9)] = b == 2 ? -x[IX(i, 1, j + 9)] : x[IX(i,1, j + 9)];
      x[IX(i, 0, j + 10)] = b == 2 ? -x[IX(i, 1, j + 10)] : x[IX(i, 1, j + 10)];
      x[IX(i, 0, j + 11)] = b == 2 ? -x[IX(i, 1, j + 11)] : x[IX(i, 1, j + 11)];
      x[IX(i, 0, j + 12)] = b == 2 ? -x[IX(i, 1, j + 12)] : x[IX(i, 1, j + 12)];
      x[IX(i, 0, j + 13)] = b == 2 ? -x[IX(i, 1, j + 13)] : x[IX(i, 1, j + 13)];
      x[IX(i, 0, j + 14)] = b == 2 ? -x[IX(i, 1, j + 14)] : x[IX(i, 1, j + 14)];
      x[IX(i, 0, j + 15)] = b == 2 ? -x[IX(i, 1, j + 15)] : x[IX(i, 1, j + 15)];
    }
  }


  // Set corners
  x[IX(0, 0, 0)] = 0.33f * (x[IX(1, 0, 0)] + x[IX(0, 1, 0)] + x[IX(0, 0, 1)]);
  x[IX(M + 1, 0, 0)] =
      0.33f * (x[IX(M, 0, 0)] + x[IX(M + 1, 1, 0)] + x[IX(M + 1, 0, 1)]);
  x[IX(0, N + 1, 0)] =
      0.33f * (x[IX(1, N + 1, 0)] + x[IX(0, N, 0)] + x[IX(0, N + 1, 1)]);
  x[IX(M + 1, N + 1, 0)] = 0.33f * (x[IX(M, N + 1, 0)] + x[IX(M + 1, N, 0)] +
                                    x[IX(M + 1, N + 1, 1)]);
}

// Linear solve for implicit methods (diffusion)
void lin_solve(int M, int N, int O, int b, float *x, float *x0, float a, float c) {
  for (int l = 0; l < LINEARSOLVERTIMES; l++) {
    for (int k = 1; k <= O; k++) {
      for (int j = 1; j <= N; j++) {
        int i;
        // Loop unrolling com fator 8
        for (i = 1; i <= M - 7; i += 8) {
          // Índices para os 8 elementos consecutivos
          int idx1 = IX(i, j, k);
          int idx2 = IX(i + 1, j, k);
          int idx3 = IX(i + 2, j, k);
          int idx4 = IX(i + 3, j, k);
          int idx5 = IX(i + 4, j, k);
          int idx6 = IX(i + 5, j, k);
          int idx7 = IX(i + 6, j, k);
          int idx8 = IX(i + 7, j, k);

          // Atualização dos 8 elementos com os cálculos dos vizinhos
          x[idx1] = (x0[idx1] + a * (x[IX(i - 1, j, k)] + x[IX(i + 1, j, k)] +
                                     x[IX(i, j - 1, k)] + x[IX(i, j + 1, k)] +
                                     x[IX(i, j, k - 1)] + x[IX(i, j, k + 1)])) / c;

          x[idx2] = (x0[idx2] + a * (x[IX(i, j, k)] + x[IX(i + 2, j, k)] +
                                     x[IX(i + 1, j - 1, k)] + x[IX(i + 1, j + 1, k)] +
                                     x[IX(i + 1, j, k - 1)] + x[IX(i + 1, j, k + 1)])) / c;

          x[idx3] = (x0[idx3] + a * (x[IX(i + 1, j, k)] + x[IX(i + 3, j, k)] +
                                     x[IX(i + 2, j - 1, k)] + x[IX(i + 2, j + 1, k)] +
                                     x[IX(i + 2, j, k - 1)] + x[IX(i + 2, j, k + 1)])) / c;

          x[idx4] = (x0[idx4] + a * (x[IX(i + 2, j, k)] + x[IX(i + 4, j, k)] +
                                     x[IX(i + 3, j - 1, k)] + x[IX(i + 3, j + 1, k)] +
                                     x[IX(i + 3, j, k - 1)] + x[IX(i + 3, j, k + 1)])) / c;

          x[idx5] = (x0[idx5] + a * (x[IX(i + 3, j, k)] + x[IX(i + 5, j, k)] +
                                     x[IX(i + 4, j - 1, k)] + x[IX(i + 4, j + 1, k)] +
                                     x[IX(i + 4, j, k - 1)] + x[IX(i + 4, j, k + 1)])) / c;

          x[idx6] = (x0[idx6] + a * (x[IX(i + 4, j, k)] + x[IX(i + 6, j, k)] +
                                     x[IX(i + 5, j - 1, k)] + x[IX(i + 5, j + 1, k)] +
                                     x[IX(i + 5, j, k - 1)] + x[IX(i + 5, j, k + 1)])) / c;

          x[idx7] = (x0[idx7] + a * (x[IX(i + 5, j, k)] + x[IX(i + 7, j, k)] +
                                     x[IX(i + 6, j - 1, k)] + x[IX(i + 6, j + 1, k)] +
                                     x[IX(i + 6, j, k - 1)] + x[IX(i + 6, j, k + 1)])) / c;

          x[idx8] = (x0[idx8] + a * (x[IX(i + 6, j, k)] + x[IX(i + 8, j, k)] +
                                     x[IX(i + 7, j - 1, k)] + x[IX(i + 7, j + 1, k)] +
                                     x[IX(i + 7, j, k - 1)] + x[IX(i + 7, j, k + 1)])) / c;
        }

        // Processar qualquer elemento restante (se M não for múltiplo de 8)
        for (; i <= M; i++) {
          int idx = IX(i, j, k);
          int idx_i1 = IX(i - 1, j, k);
          int idx_i2 = IX(i + 1, j, k);
          int idx_j1 = IX(i, j - 1, k);
          int idx_j2 = IX(i, j + 1, k);
          int idx_k1 = IX(i, j, k - 1);
          int idx_k2 = IX(i, j, k + 1);
          x[idx] = (x0[idx] + a * (x[idx_i1] + x[idx_i2] +
                                   x[idx_j1] + x[idx_j2] +
                                   x[idx_k1] + x[idx_k2])) / c;
        }
      }
    }
    set_bnd(M, N, O, b, x);
  }
}



// Diffusion step (uses implicit method)
void diffuse(int M, int N, int O, int b, float *x, float *x0, float diff,
             float dt) {
  int max = MAX(MAX(M, N), O);
  float a = dt * diff * max * max;
  lin_solve(M, N, O, b, x, x0, a, 1 + 6 * a);
}

// Advection step (uses velocity field to move quantities)
void advect(int M, int N, int O, int b, float *d, float *d0, float *u, float *v,
            float *w, float dt) {
  float dtX = dt * M, dtY = dt * N, dtZ = dt * O;

  for (int i = 1; i <= M; i++) {
    for (int j = 1; j <= N; j++) {
      for (int k = 1; k <= O; k++) {
        float x = i - dtX * u[IX(i, j, k)];
        float y = j - dtY * v[IX(i, j, k)];
        float z = k - dtZ * w[IX(i, j, k)];

        // Clamp to grid boundaries
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

        d[IX(i, j, k)] =
            s0 * (t0 * (u0 * d0[IX(i0, j0, k0)] + u1 * d0[IX(i0, j0, k1)]) +
                  t1 * (u0 * d0[IX(i0, j1, k0)] + u1 * d0[IX(i0, j1, k1)])) +
            s1 * (t0 * (u0 * d0[IX(i1, j0, k0)] + u1 * d0[IX(i1, j0, k1)]) +
                  t1 * (u0 * d0[IX(i1, j1, k0)] + u1 * d0[IX(i1, j1, k1)]));
      }
    }
  }
  set_bnd(M, N, O, b, d);
}

// Projection step to ensure incompressibility (make the velocity field
// divergence-free)
void project(int M, int N, int O, float *u, float *v, float *w, float *p, float *div) {
  // Primeira parte: cálculo de div e inicialização de p
  for (int i = 1; i <= M; i++) {
    for (int j = 1; j <= N; j++) {
      int k;
      // Unroll o loop de k com fator 16
      for (k = 1; k <= O - 15; k += 16) {
        // Atualizar 16 elementos de div e p
        for (int offset = 0; offset < 16; offset++) {
          int idx = IX(i, j, k + offset);
          div[idx] = -0.5f * (u[IX(i + 1, j, k + offset)] - u[IX(i - 1, j, k + offset)] +
                              v[IX(i, j + 1, k + offset)] - v[IX(i, j - 1, k + offset)] +
                              w[IX(i, j, k + offset + 1)] - w[IX(i, j, k + offset - 1)]) / MAX(M, MAX(N, O));
          p[idx] = 0;
        }
      }
      // Processar elementos restantes (se houver)
      for (; k <= O; k++) {
        int idx = IX(i, j, k);
        div[idx] = -0.5f * (u[IX(i + 1, j, k)] - u[IX(i - 1, j, k)] +
                            v[IX(i, j + 1, k)] - v[IX(i, j - 1, k)] +
                            w[IX(i, j, k + 1)] - w[IX(i, j, k - 1)]) / MAX(M, MAX(N, O));
        p[idx] = 0;
      }
    }
  }

  set_bnd(M, N, O, 0, div);
  set_bnd(M, N, O, 0, p);
  lin_solve(M, N, O, 0, p, div, 1, 6);

  // Segunda parte: atualização de u, v e w
  for (int i = 1; i <= M; i++) {
    for (int j = 1; j <= N; j++) {
      int k;
      // Unroll o loop de k com fator 16
      for (k = 1; k <= O - 15; k += 16) {
        // Atualizar 16 elementos de u, v e w
        for (int offset = 0; offset < 16; offset++) {
          int idx = IX(i, j, k + offset);
          u[idx] -= 0.5f * (p[IX(i + 1, j, k + offset)] - p[IX(i - 1, j, k + offset)]);
          v[idx] -= 0.5f * (p[IX(i, j + 1, k + offset)] - p[IX(i, j - 1, k + offset)]);
          w[idx] -= 0.5f * (p[IX(i, j, k + offset + 1)] - p[IX(i, j, k + offset - 1)]);
        }
      }
      // Processar elementos restantes (se houver)
      for (; k <= O; k++) {
        int idx = IX(i, j, k);
        u[idx] -= 0.5f * (p[IX(i + 1, j, k)] - p[IX(i - 1, j, k)]);
        v[idx] -= 0.5f * (p[IX(i, j + 1, k)] - p[IX(i, j - 1, k)]);
        w[idx] -= 0.5f * (p[IX(i, j, k + 1)] - p[IX(i, j, k - 1)]);
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
