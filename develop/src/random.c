/*
    Name:       random.h
    Author:     Peng Yi at MIT
    Date:       October 23, 2006
    Purpose:    Using code from 'Numerical Recipes', chapter 7.
*/
#define __RANDOM_MODULE
#include "random.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) 
    {
      if (-(*idum)< 1) *idum=1;
      else *idum = -(*idum);
      for (j=NTAB+7; j>=0; j--)
	{
	  k=(*idum)/IQ;
	  *idum=IA*(*idum-k*IQ) - IR*k;
	  if(*idum < 0) *idum += IM;
	  if (j < NTAB) iv[j] = *idum;
	}
      iy=iv[0];
    }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=*idum;
  if ((temp=AM*iy)> RNMX) return RNMX;
  else return temp;
}


double gauss(double std, double mean)		// Frenkel and Smit, Algorithm 44
{
   double	r, x, v1, v2;
   
   r	=	2.0;
   while (r >= 1.0) {
      v1	=	2.0 * ran1(seed) - 1.0;
      v2	=	2.0 * ran1(seed) - 1.0;
      r		=	v1*v1 + v2*v2;
   }
   x	=	v1 * sqrt(-2.0*log(r)/r);
   x	=	mean + std * x;
   return	x;
}




