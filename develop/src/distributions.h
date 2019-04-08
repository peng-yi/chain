/*
    program:	distributions.h
    author:	Pieter J. in 't Veld for UT at Austin
    date:	June 12, 1998
    purpose:	header file for distributions.c
*/
#ifndef __DISTRIBUTIONS_HEADER
#define __DISTRIBUTIONS_HEADER

#include "header.h"

#define D_NAVERAGE	6		// do average to x, x^2, x^3 up to x^6

#ifndef __DISTRIBUTIONS_MODULE

extern void D_Allocate(diststruct *d, long n);
extern void D_Reset(diststruct *d);
extern void D_Copy(diststruct *dest, diststruct *src);
extern void D_Add(diststruct *dest, diststruct *src);
extern void D_Subtr(diststruct *dest, diststruct *src);
extern void D_Submit(diststruct *d, double *x, double *y, double *weight);
extern void D_Print(diststruct *d, int dist_type, double *L);
extern void D_PrintMath(diststruct *d, int dist_type, double *L);
extern void D_Init(
  diststruct **dist, char *header, long nlevels, double *binsize);

extern void PutInDistribution(diststruct *d, double x, double y, double weight);
extern void PrintDistribution(diststruct *d);
extern void InitDist(long flag, diststruct **dist, char *s, long nlevels, double *binsize);

#endif

#endif

