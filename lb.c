#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#define _USE_MATH_DEFINES
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#ifndef N
#define N 130
#endif

#ifndef M
#define M 45
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

static double f[9][N * M];
static double ft[9][N * M];
static double feq[9][N * M];
static double u[2][N * M];
static double rho[N * M];

static inline int idx(int i, int j) {
   int idx = M * i + j;
   assert(idx >= 0);
   assert(idx < N*M);
   return idx;
}

void update_rho(double rho[], double f[][N*M]) {
#pragma omp parallel for
   for (int i = 0; i < N*M; i++) {
      double sum = 0;
      for (int q = 0; q < 9; q++)
         sum += f[q][i];

      rho[i] = sum;
   }
}

void update_u(double u[][N*M], double rho[], double f[][N*M]) {
#pragma omp parallel for collapse(2)
   for (int j = 0; j < 2; j++) {
      for (int i = 0; i < N*M; i++) {
         double sum = 0;
         for (int q = 0; q < 9; q++)
            sum += f[q][i] * c[q][j];

         u[j][i] = sum / rho[i];
      }
   }
}

void update_feq(double feq[][N*M], double rho[], double u[][N*M], double cs) {
   double oocs2 = 1. / (cs*cs);
   double oocs4 = oocs2 * oocs2;
#pragma omp parallel for collapse(2)
   for (int q = 0; q < 9; q++) {
      for (int i = 0; i < N*M; i++) {
         double cu = 0, cu2 = 0, u2 = 0;
         for (int j = 0; j < 2; j++) {
            cu += c[q][j] * u[j][i];
            u2 += u[j][i] * u[j][i];
         }
         cu2 = cu * cu;

         feq[q][i] = w[q] * rho[i]
            * (1 + oocs2 * cu + 0.5 * oocs4 * cu2 - 0.5 * oocs2 * u2);
      }
   }
}

void collide(double ft[][N*M], double f[][N*M], double feq[][N*M], int obstacle[], double omega) {
#pragma omp parallel for collapse(2)
   for (int q = 0; q < 9; q++) {
      for (int i = 0; i < N*M; i++) {
         ft[q][i] = f[q][i] - omega * (f[q][i] - feq[q][i]);
      }
   }

#pragma omp parallel for collapse(2)
   for (int q = 0; q < 9; q++) {
      for (int i = 0; i < N*M; i++) {
         if (obstacle[i]) {
            ft[q][i] = f[noslip[q]][i];
         }
      }
   }
}

void stream(double f[][N*M], double ft[][N*M]) {
#pragma omp parallel for collapse(3)
   for (int q = 0; q < 9; q++) {
      for (int i = 0; i < N; i++) {
         for (int j = 0; j < M; j++) {
            f[q][idx((N + (i + c[q][0]) % N) % N, (M + (j + c[q][1]) % M) % M)] = ft[q][idx(i, j)];
         }
      }
   }
}

static inline double vel(int k, int i, double ulb) {
   double x = (double)i / (M-1);
   return (1. - k) * ulb * (1. + 1.e-4 * sin(2.*M_PI * x));
}

void update(double f[][N*M], double ft[][N*M], double feq[][N*M], double rho[], double u[][N*M], int obstacle[], double cs, double omega, double ulb) {
   // outflow
#pragma omp parallel for collapse(2)
   for (int q = 3; q < 6; q++)
      for (int j = 0; j < M; j++)
         f[q][idx(N-1, j)] = f[q][idx(N-2, j)];

   update_rho(rho, f);
   update_u(u, rho, f);

   // inflow
#pragma omp parallel for
   for (int j = 0; j < M; j++) {
      double sum1 = 0, sum2 = 0;
      for (int q = 0; q < 6; q++) {
         if (q < 3)
            sum2 += f[q][j];
         else
            sum1 += f[q][j];
      }
      for (int k = 0; k < 2; k++)
         u[k][j] = vel(k, j, ulb);

      rho[j] = 1. / (1. - u[0][j]) * (sum2 + 2. * sum1);
   }
   update_feq(feq, rho, u, cs);
#pragma omp parallel for collapse(2)
   for (int q = 6; q < 9; q++)
      for (int j = 0; j < M; j++)
         f[q][j] = feq[q][j];

   collide(ft, f, feq, obstacle, omega);
   stream(f, ft);
}

void print_u(double u[][N*M], int t) {
   printf("%d", t);
   for (int i = 0; i < N*M; i++) {
      printf(" %E", sqrt(u[0][i] * u[0][i] + u[1][i] * u[1][i]));
   }
   putchar('\n');
}

void print_rho(double rho[], int t) {
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

   int *obstacle = malloc(N*M * sizeof(int));
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         rho[M*i + j] = 1.;
         if ((i-cx)*(i-cx) + (j-cy)*(j-cy) < r*r)
            obstacle[idx(i, j)] = 1;
         else
            obstacle[idx(i, j)] = 0;
      }
   }

   for (int k = 0; k < 2; k++)
      for (int i = 0; i < N; i++)
         for (int j = 0; j < M; j++)
            u[k][idx(i, j)] = vel(k, j, ulb);

   update_feq(feq, rho, u, cs);
   for (int q = 0; q < 9; q++)
      for (int i = 0; i < N*M; i++)
         f[q][i] = feq[q][i];

   int T = 30000;
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

   free(obstacle);

   return 0;
}
