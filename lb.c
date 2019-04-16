#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "param.h"

#define _USE_MATH_DEFINES
#ifndef M_PI
#define M_PI 3.14159265358979323846
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
   return id;
}

static inline int mod(int a, int b) {
   return (b + a%b) % b;
}

static inline double inlet_vel(int i) {
   double x = (double)i / (M-1);
   return ULB * (1. + 1.e-4 * sin(2.*M_PI * x));
}

static inline double rho(int i, int j, double* restrict f) {
   double sum = 0;
   for (int q = 0; q < 9; q++)
      sum += f[9*idx(i,j) + q];

   return sum;
}

static inline double u(int k, int i, int j, double* restrict f) {
   double sum = 0;
   for (int q = 0; q < 9; q++)
      sum += f[9*idx(i,j) + q] * c[q][k];

   return sum / rho(i,j, f);
}

void boundary(double* restrict f) {
   // outflow
   for (int j = 0; j < M; j++)
      for (int q = 3; q < 6; q++)
         f[9*idx(N-1, j) + q] = f[9*idx(N-2, j) + q];

   // inflow
   for (int j = 0; j < M; j++) {
      double ux = inlet_vel(j);
      double u2 = ux*ux;
      double sum1 = 0, sum2 = 0;
      for (int q = 0; q < 6; q++) {
         if (q < 3)
            sum2 += f[9*j + q];
         else
            sum1 += f[9*j + q];
      }

      double rho_j = 1. / (1 - ux) * (sum2 + 2*sum1);

      for (int q = 6; q < 9; q++) {
         double cu = c[q][0] * ux;
         f[9*j + q] = w[q] * rho_j * (1 + 3 * cu + 0.5 * 9 * cu*cu - 0.5 * 3 * u2);
      }
   }
}

void collstream(double* restrict fnew, double* restrict fold, int* restrict obstacle, double omega) {
#pragma omp parallel for collapse(2)
   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         int id = idx(i,j);
         double ux = 0;
         double uy = 0;
         double rho_ij = 0;
         for (int q = 0; q < 9; q++) {
            rho_ij += fold[9*id + q];
            ux += fold[9*id + q] * c[q][0];
            uy += fold[9*id + q] * c[q][1];
         }
         ux /= rho_ij;
         uy /= rho_ij;
         double u2 = ux*ux + uy*uy;
         for (int q = 0; q < 9; q++) {
            int idnew = 9*idx(mod(i + c[q][0], N), mod(j + c[q][1], M)) + q;
            double cu = c[q][0] * ux + c[q][1] * uy;
            if (!obstacle[id]) {
               double feq_ijq = w[q] * rho_ij * (1 + 3 * cu + 0.5 * 9 * cu*cu - 0.5 * 3 * u2);
               fnew[idnew] = (1 - omega) * fold[9*id + q] + omega*feq_ijq;
            } else
               fnew[idnew] = fold[9*id + noslip[q]];
         }
      }
   }
}

void update(double* restrict fnew, double* restrict fold, int* restrict obstacle, double omega) {
   collstream(fnew, fold, obstacle, omega);
   boundary(fnew);
}

static void print_u(double* restrict f, int t) {
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
   double nu = ULB * r / RE;
   double omega = 1. / (3. * nu + 0.5);

   double* restrict fnew = malloc(9 * N * M * sizeof(double));
   double* restrict fold = malloc(9 * N * M * sizeof(double));

   int* restrict obstacle = malloc(N*M * sizeof(int));
   for (int i = 0; i < N; i++)
      for (int j = 0; j < M; j++)
         obstacle[idx(i,j)] = ((i-cx)*(i-cx) + (j-cy)*(j-cy) < r*r);

   for (int i = 0; i < N; i++) {
      for (int j = 0; j < M; j++) {
         double ux = inlet_vel(j);
         double u2 = ux*ux;
         for (int q = 0; q < 9; q++) {
            double cu = c[q][0] * ux;
            fold[9*idx(i,j) + q] = w[q] * 1.0 * (1 + 3 * cu + 0.5 * 9 * cu*cu - 0.5 * 3 * u2);
         }
      }
   }

   for (int t = 0; t < T; t++) {
      if (t % 2) {
         update(fold, fnew, obstacle, omega);
      } else {
         if(!(t % 100))
            print_u(fold, t);
         update(fnew, fold, obstacle, omega);
      }
   }

   free(fnew);
   free(fold);
   free(obstacle);

   return 0;
}
