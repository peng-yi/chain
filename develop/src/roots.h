/*
    program:	roots.h
    author:	Pieter J. in 't Veld
    date:	May 10, 2001
    purpose:	roots of upto fourth order polynomials
*/
#ifndef __ROOTS_HEADER
#define __ROOTS_HEADER

#include <math.h>
#include <complex.h>


#ifdef __ROOTS_MODULE

#else


//extern long R_ComplexCount(long n, double complex *z);
//extern long R_SolveFirstOrder(double *c, double complex *z);
//extern long R_SolveSecondOrder(double *c, double complex *z);
//extern long R_SolveThirdOrder(double *c, double complex *z);
//extern long R_SolveFourthOrder(double *c, double complex *z);
//extern long R_SolveAnalytical(long n, double *c, double complex *z);
extern long R_NewtonRaphson(
  double y(double), double x_low, double x_high, double *root);
extern long R_Brent(
  double y(double), double x_low, double x_high, double *root);
extern long R_SolveNumerical(
  double y(double), double x_low, double x_high, long n_mesh, long depth,
  double *roots, long n_max);

#endif

#endif

