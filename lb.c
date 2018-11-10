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

static inline int idx(int i, int j) {
   int id = M * i + j;
   assert(id >= 0);
   assert(id < N*M);
   return id;
}

static inline int mod(int a, int b) {
   return (b + a%b) % b;
}

static inline double inlet_vel(int k, int i, double ulb) {
   double x = (double)i / (M-1);
   return (1. - k) * ulb * (1. + 1.e-4 * sin(2.*M_PI * x));
}

double rho(int i, int j, double* restrict f) {
   double sum = 0;
   for (int q = 0; q < 9; q++)
      sum += f[9*idx(i,j) + q];

   return sum;
}

double u(int k, int i, int j, double* restrict f) {
   double sum = 0;
   for (int q = 0; q < 9; q++)
      sum += f[9*idx(i,j) + q] * c[q][k];

   return sum / rho(i,j, f);
}

double feq(int i, int j, int q, double* restrict f) {
   double u_x = u(0,i,j, f);
   double u_y = u(1,i,j, f);

   double cu = c[q][0] * u_x + c[q][1] * u_y;
   double u2 = u_x * u_x + u_y * u_y;

   return w[q] * rho(i,j, f) * (1 + 3 * cu + 0.5 * 9 * cu*cu - 0.5 * 3 * u2);
}

void boundary(double* restrict f, double ulb) {
   // outflow
   for (int j = 0; j < M; j++)
      for (int q = 3; q < 6; q++)
         f[9*idx(N-1, j) + q] = f[9*idx(N-2, j) + q];

   // inflow
   for (int j = 0; j < M; j++) {
      double u_x = inlet_vel(0, j, ulb);
      double u_y = inlet_vel(1, j, ulb);
      double u2 = u_x*u_x + u_y*u_y;
      double sum1 = 0, sum2 = 0;
      for (int q = 0; q < 6; q++) {
         if (q < 3)
            sum2 += f[9*j + q];
         else
            sum1 += f[9*j + q];
      }

      double rho_j = 1. / (1 - u_x) * (sum2 + 2*sum1);

      for (int q = 6; q < 9; q++) {
         double cu = c[q][0] * u_x + c[q][1] * u_y;
         f[9*j + q] = w[q] * rho_j * (1 + 3 * cu + 0.5 * 9 * cu*cu - 0.5 * 3 * u2);
      }
   }
}

void collstream(double* restrict fnew, double* restrict fold, int* restrict obstacle, double omega) {
#pragma omp parallel for
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         for (int q = 0; q < 9; q++) {
            int id = 9*idx(i,j) + q;
            int idnew = 9*idx(mod(i + c[q][0], N), mod(j + c[q][1], M)) + q;
            if (!obstacle[idx(i,j)])
               fnew[idnew] = (1 - omega) * fold[id] + omega*feq(i,j,q, fold);
            else
               fnew[idnew] = fold[9*idx(i,j) + noslip[q]];
         }
      }
   }
}

void update(double* restrict fnew, double* restrict fold, int* restrict obstacle, double omega, double ulb) {
   collstream(fnew, fold, obstacle, omega);
   boundary(fnew, ulb);
}

static inline void print_u(double* restrict f, int t) {
   if (t % 100)
      return;

   printf("%d", t);
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         double ux = u(0, i, j, f);
         double uy = u(1, i, j, f);
         printf(" %E", sqrt(ux*ux + uy*uy));
      }
   }
   putchar('\n');
}

int main() {
   double Re = 160.;
   double ulb = 0.04;
   double nu = ulb * r / Re;
   double omega = 1. / (3. * nu + 0.5);

   double* restrict fnew = malloc(9 * N * M * sizeof(double));
   double* restrict fold = malloc(9 * N * M * sizeof(double));

   int* restrict obstacle = malloc(N*M * sizeof(int));
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         if ((i-cx)*(i-cx) + (j-cy)*(j-cy) < r*r)
            obstacle[idx(i, j)] = 1;
         else
            obstacle[idx(i, j)] = 0;
      }
   }

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         double ux = inlet_vel(0, j, ulb);
         double uy = inlet_vel(1, j, ulb);
         double u2 = ux*ux + uy*uy;
         for (int q = 0; q < 9; q++) {
            double cu = c[q][0] * ux + c[q][1] * uy;
            fold[9*idx(i,j) + q] = w[q] * 1.0 * (1 + 3 * cu + 0.5 * 9 * cu*cu - 0.5 * 3 * u2);
         }
      }
   }

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
      if (t % 2) {
         print_u(fnew, t);
         update(fold, fnew, obstacle, omega, ulb);
      } else {
         print_u(fold, t);
         update(fnew, fold, obstacle, omega, ulb);
      }
   }

   free(fnew);
   free(fold);
   free(obstacle);

   return 0;
}
