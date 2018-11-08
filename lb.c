#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define _USE_MATH_DEFINES
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef N
#define N 128
#endif

#ifndef M
#define M 48
#endif

static const int cx = N / 4;
static const int cy = M / 2;
static const int r = M / 9;

static const double w[9] = {
   4./9. ,
   1./9. ,
   1./9. ,
   1./9. ,
   1./36.,
   1./36.,
   1./9. ,
   1./36.,
   1./36.
};

static const int c[9][2] = {
   { 0,  0},
   { 0, -1},
   { 0,  1},
   {-1,  0},
   {-1, -1},
   {-1,  1},
   { 1,  0},
   { 1, -1},
   { 1,  1}
};

static const int noslip[9] = {0, 2, 1, 6, 8, 7, 3, 5, 4};

struct vec {
   double x, y;
};

static inline int idx(int i, int j) {
   int id = M * i + j;
   assert(id >= 0);
   assert(id < N*M);
   return id;
}

void update_rho(double* restrict rho, double* restrict f) {
#pragma omp parallel for
   for (int i = 0; i < N*M; i++) {
      double sum = 0;
      for (int q = 0; q < 9; q++)
         sum += f[9*i + q];

      rho[i] = sum;
   }
}

void update_u(struct vec * restrict u, double* restrict rho, double* restrict f) {
#pragma omp parallel for
   for (int i = 0; i < N*M; i++) {
      double sum1 = 0, sum2 = 0;
      for (int q = 0; q < 9; q++) {
         sum1 += f[9*i + q] * c[q][0];
         sum2 += f[9*i + q] * c[q][1];
      }

      u[i].x = sum1 / rho[i];
      u[i].y = sum2 / rho[i];
   }
}

void update_feq(double* restrict feq, double* restrict rho, struct vec * restrict u, double cs) {
   double oocs2 = 1. / (cs*cs);
   double oocs4 = oocs2 * oocs2;
#pragma omp parallel for collapse(2)
   for (int i = 0; i < N*M; i++) {
      for (int q = 0; q < 9; q++) {
         double cu = c[q][0] * u[i].x + c[q][1] * u[i].y;
         double cu2 = cu * cu;
         double u2 = u[i].x * u[i].x + u[i].y * u[i].y;

         feq[9*i + q] = w[q] * rho[i]
            * (1 + oocs2 * cu + 0.5 * oocs4 * cu2 - 0.5 * oocs2 * u2);
      }
   }
}

void collide(double* restrict ft, double* restrict f, double* restrict feq, int obstacle[], double omega) {
#pragma omp parallel for collapse(2)
   for (int i = 0; i < N*M; i++) {
      for (int q = 0; q < 9; q++) {
         int id = 9*i + q;
         if (!obstacle[i])
            ft[id] = f[id] - omega * (f[id] - feq[id]);
         else
            ft[id] = f[9*i + noslip[q]];
      }
   }
}

void stream(double* restrict f, double* restrict ft) {
#pragma omp parallel for collapse(3)
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         for (int q = 0; q < 9; q++) {
            f[9*idx((N + (i + c[q][0]) % N) % N, (M + (j + c[q][1]) % M) % M) + q] = ft[9*idx(i, j) + q];
         }
      }
   }
}

static inline double inlet_vel(int k, int i, double ulb) {
   double x = (double)i / (M-1);
   return (1. - k) * ulb * (1. + 1.e-4 * sin(2.*M_PI * x));
}

void update(double* restrict f, double* restrict ft, double* restrict feq, double* restrict rho, struct vec * restrict u, int obstacle[], double cs, double omega, double ulb) {
   // outflow
#pragma omp parallel for collapse(2)
   for (int j = 0; j < M; j++)
      for (int q = 3; q < 6; q++)
         f[9*idx(N-1, j) + q] = f[9*idx(N-2, j) + q];

   update_rho(rho, f);
   update_u(u, rho, f);

   // inflow
#pragma omp parallel for
   for (int j = 0; j < M; j++) {
      double sum1 = 0, sum2 = 0;
      for (int q = 0; q < 6; q++) {
         if (q < 3)
            sum2 += f[9*j + q];
         else
            sum1 += f[9*j + q];
      }
      u[j].x = inlet_vel(0, j, ulb);
      u[j].y = inlet_vel(1, j, ulb);

      rho[j] = 1. / (1. - u[j].x) * (sum2 + 2. * sum1);
   }
   update_feq(feq, rho, u, cs);
#pragma omp parallel for collapse(2)
   for (int j = 0; j < M; j++)
      for (int q = 6; q < 9; q++)
         f[9*j + q] = feq[9*j + q];

   collide(ft, f, feq, obstacle, omega);
   stream(f, ft);
}

void print_u(struct vec * restrict u, int t) {
   printf("%d", t);
   for (int i = 0; i < N*M; i++) {
      printf(" %E", sqrt(u[i].x * u[i].x + u[i].y * u[i].y));
   }
   putchar('\n');
}

void print_rho(double* restrict rho, int t) {
   printf("%d", t);
   for (int i = 0; i < N*M; i++) {
      printf(" %E", rho[i]);
   }
   putchar('\n');
}

int main() {
   double cs = 1./sqrt(3.);
   double Re = 160.;
   double ulb = 0.04;
   double nu = ulb * r / Re;
   double omega = 1. / (3. * nu + 0.5);

   double* restrict f = malloc(9 * N * M * sizeof(double));
   double* restrict ft = malloc(9 * N * M * sizeof(double));
   double* restrict feq = malloc(9 * N * M * sizeof(double));
   double* restrict rho = malloc(N * M * sizeof(double));
   struct vec* restrict u = malloc(N * M * sizeof(struct vec));

   int* restrict obstacle = malloc(N*M * sizeof(int));
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         rho[M*i + j] = 1.;
         if ((i-cx)*(i-cx) + (j-cy)*(j-cy) < r*r)
            obstacle[idx(i, j)] = 1;
         else
            obstacle[idx(i, j)] = 0;
      }
   }

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         u[idx(i, j)].x = inlet_vel(0, j, ulb);
         u[idx(i, j)].y = inlet_vel(1, j, ulb);
      }
   }

   update_feq(feq, rho, u, cs);
   for (int i = 0; i < N*M; i++)
      for (int q = 0; q < 9; q++)
         f[9*i + q] = feq[9*i + q];

   int T = 20000;
   //int barsteps = (int)ceil(T / 100.);
   //int barwidth = 40;
   for (int t = 0; t < T; t++) {
      /*
      if (!(t % 10)) {
         float progress = (float)t / (T - 2);

         putchar('[');
         int pos = barwidth * progress;
         for (int i = 0; i < barwidth; ++i) {
            if (i < pos)
               putchar('=');
            else if (i == pos)
               putchar('>');
            else
               putchar(' ');
         }
         printf("] %d \r", (int)(progress * 100.0));
         fflush(stdout);
      }
      */
      if (!(t % 100)) {
         print_u(u, t);
         //print_rho(rho, t);
      }
      update(f, ft, feq, rho, u, obstacle, cs, omega, ulb);
   }

   free(f);
   free(ft);
   free(feq);
   free(rho);
   free(u);
   free(obstacle);

   return 0;
}
