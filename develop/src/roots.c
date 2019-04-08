/*
    program:	roots.c
    author:	Pieter J. in 't Veld
    date:	May 10, 2001
    purpose:	roots of upto fourth order polynomials
*/
#define __ROOTS_MODULE
#include "roots.h"
#include <stdio.h>
#define R_VERSION	"3.2"

#define R_ZERO		1e-14
#define R_CZERO		1e-8

long R_Round(long n, double complex *z)
{
  long			i;

  for (i=0; i<n; ++i)
  {
    if ((creal(z[i])<R_CZERO)&&(creal(z[i])>-R_CZERO)) 
       z[i]	= 0.0 + cimag(z[i]) * I;	// z[i].re = 0
    if ((cimag(z[i])<R_CZERO)&&(cimag(z[i])>-R_CZERO)) 
       z[i]	= creal(z[i]) + 0.0 * I;	// z[i].im = 0
  }
  return n;
}


long R_ComplexCount(long n, double complex *z)
{
  long			i, count = 0;

  for (i=0; i<n; ++i)
    if (cimag(z[i])) ++count;
  return count;
}


long R_Sort(long n, double complex *z)
{
  long			i, j;
  double complex		zt;
  
  for (i=0; i<n-1; ++i)
    for (j=i+1; j<n; ++j)
      if (creal(z[j]) < creal(z[i]))
      {
        zt		= z[i];
	z[i]		= z[j];
	z[j]		= zt;
      }
      else
        if ((creal(z[i])==creal(z[j]))&&(cimag(z[j])<cimag(z[i])))
	{
          zt		= z[i];
	  z[i]		= z[j];
	  z[j]		= zt;
	}
  return n;
}


long R_SolveFirstOrder(double *c, double complex *z)	// solve c[0] + c[1]*z = 0
{
  if (fabs(c[1])<R_ZERO)
    return 0;

  *z	=	-c[0]/c[1] + 0 * I;
  return 1;
}


long R_SolveSecondOrder(double *c, double complex *z)	// c[0] + c[1] * z + c[2] * z^2 = 0
{
  double complex	zR = c[1]*c[1]-4.0*c[0]*c[2] + 0.0 * I;

  if (fabs(c[2])<R_ZERO)
    return R_SolveFirstOrder(c, z);
  
  if (fabs(creal(zR))<R_ZERO)
    zR	=	0.0;	// zR.re = 0.0;

  zR	=	cpow(zR, 0.5);

  zR	/=	2.0 * c[2];
  z[0]	=	-c[1]/c[2]/2.0 + 0.0 * I;
  z[1]	=	z[0];

  z[0]	-=	zR;
  z[1]	+=	zR;

  return R_Round(2, z);
}


long R_SolveThirdOrder(double *c, double complex *z)	// c[0]+c[1]z+c[2]z^2+c[3]z^3=0
{
  double		f1, f2, f3;
  double complex	A, A1, B, z1, z2;

  if (fabs(c[3])<R_ZERO)
    return R_SolveSecondOrder(c, z);
  
  // real constants
  
  f1			= 3.0*c[1]*c[3]-c[2]*c[2];
  f2			= (c[2]*(9.0*c[1]*c[3]-2.0*c[2]*c[2])
  			  -27.0*c[0]*c[3]*c[3])/2.0;
  f3			= -c[2]/(3.0*c[3]);
 
  // complex constants
  A1		=	0.0 + 0.0 * I; 
  if (fabs(f1)<1e-3)					// Limit f1->0
  {
    A		=	(f2<0.0 ? pow(-2.0*f2, 1.0/3.0)/(3.0*c[3]) : 0.0) + 0.0 * I;
  } 
  else if (fabs(f2)<R_ZERO)				// Limit f2->0
  {
    A		=	sqrt(fabs(f1))/(3.0*c[3]) + 0.0 * I;
    if (f1 < 0.0) 
       A	=	creal(A)*sqrt(3.0)/2.0 - creal(A)/2.0 * I;
    else
       A	=	creal(A) + 0.0 * I;
  }
  else
  {
    A1	=	f1*f1*f1 + f2*f2 + 0.0 * I;
    A1	=	cpow(A1, 0.5);
    A1	+=	f2;
    A1	=	cpow(A1, 1.0/3.0);
    *z	=	f1/(3.0*c[3]) + 0.0 * I;
    A	=	(*z)/A1;
  }
  B	=	A/(3.0*c[3]);
  z1	=	1.0/2.0 + sqrt(3.0)/2.0 * I;
  z2	=	conj(z1);
  
  // roots
  z[0]	=	f3+creal(B)-creal(A) + (cimag(B)-cimag(A))*I;
  z[1]	=	z1 * A;
  A1	=	z2 * B;
  z[1]	+=	f3 - A1;
  z[2]	=	z2 * A;
  A1	=	z1 * B;
  z[2]	+=	f3 - A1;
  
  return R_Round(3, z);
}


long R_SolveFourthOrder(double *c, double complex *z)
{
  double		f1, f2, f3, f4 ;
  double complex	A, B, 
  			A1, A11, A111, B1, B2;

  if (fabs(c[4])<R_ZERO)
    return R_SolveThirdOrder(c, z);
  
  // real constant
  
  f1			= c[2]*c[2]-3.0*c[1]*c[3]+12.0*c[0]*c[4];
  f2			= (c[2]*(2.0*c[2]*c[2]-9.0*c[1]*c[3]-72.0*c[0]*c[4])+
  			  27.0*(c[0]*c[3]*c[3]+c[1]*c[1]*c[4]))/2.0;
  f3			= c[3]*c[3]/(4.0*c[4]*c[4])-2.0*c[2]/(3.0*c[4]);
  f4			= -c[3]*c[3]*c[3]/(c[4]*c[4]*c[4])
  		          +4.0*c[2]*c[3]/(c[4]*c[4])-8.0*c[1]/c[4];

  // complex constants
 
  A11	=	0.0 + 0.0 * I; 
  if (fabs(f1)<R_ZERO)				// limit f1->0
     A11	=	(f2>0.0 ? 2.0*f2 : 0.0) + 0.0 * I;
  else if (fabs(f2)<R_ZERO)			// limit f2->0
     A11	=	(f1>0.0 ? sqrt(f1/3.0) : 0.0) + 0.0 * I;
  else {
     A111	=	f2*f2 - f1*f1*f1 + 0.0 * I;
     A111	=	cpow(A111, 0.5);
     A111	+=	f2;
     A111	=	cpow(A111, 1.0/3.0);
     A		=	A111 * A111;
     A		+=	f1;
     A111	*=	3.0*c[4];
     A11	=	A/A111;
  }
  A1	=	A11+f3;
  A1	=	cpow(A1, 0.5);
  B1	=	2.0*f3-A11;

  B2	=	0.0 * I;
  if ((fabs(creal(A1))<5e-6)&&(fabs(cimag(A1))<5e-6)) {		// Limit A1->0
    A	=	c[0]/c[4] + cimag(A)*I;
    B2	=	creal(A)>=0.0 ? 0.0 : -4.0*sqrt(-creal(A));
  }
  else {
    A	=	f4/4.0 + 0.0 * I;
    B2	=	A/A1;
  }
  
  // roots

  A	=	-c[3]/(4.0*c[4])+A1/2.0;
  B	=	B1 + B2;
  B	=	cpow(B, 0.5);
  B	/=	2.0;
  z[0]	=	A+B;
  z[1]	=	A-B;
  
  A	=	-c[3]/(4.0*c[4])-A1/2.0;
  B	=	B1 - B2;
  B	=	cpow(B, 0.5);
  B	/=	2.0;
  z[2]	=	A+B;
  z[3]	=	A-B;
  return R_Round(4, z);
}


long R_SolveAnalytical(long n, double *c, double complex *z)
{
  switch (n)
  {
    case 1: return R_SolveFirstOrder(c, z);
    case 2: return R_SolveSecondOrder(c, z);
    case 3: return R_SolveThirdOrder(c, z);
    case 4: return R_SolveFourthOrder(c, z);
    default: return 0;
  }
}


#define R_MAXNROOTS	20
#define R_D2		1e-5
#define R_MAXITER	60
#define R_EPS		2.5e-13
#define sign(a)		((a<0) ? -1 : 1)

long R_NewtonRaphson(
  double y(double), double x_low, double x_high, double *root)
{
  static long		iter;
  static double		x0, x1, x2,
  			y0, y1, y2;

  x0			= x_low;
  y0			= y(x_low);
  x1			= x_high;
  y1			= y(x_high);
  iter			= 0;
  do
  {
    x2			= (x1*y0-x0*y1)/(y0-y1);
    y2			= y(x2);
    if (sign(y2)==sign(y1))
    {
      x0		= x1;
      y0		= y1;
      x1		= x2;
      y1		= y2;
    }
    else
    {
      x1		= x0;
      y1		= y0;
      x0		= x2;
      y0		= y2;
    }
  } while((fabs(x0-x1)>R_EPS)&&(++iter<R_MAXITER));
  *root			= x2;
  return (iter<R_MAXITER);
}

// Adapted from 
//   W.H. Press, Numerical Recipes in C, Chapter 9.3, Cambridge University
//   Press, 1997 (p. 361)

//long			R_NYCALL = 0;

long R_BrentInt(
  double y(double), double x_low, double x_high, 
  double y_low, double y_high, double tol, double *root)
{
  static long		iter;
  static double		a, b, c, d, e, min1,
  			min2, fa, fb, fc, p, q, r, s, tol1,
			xm;

  a			= x_low;
  b			= x_high;
  c			= x_high;
  fa			= y_low; 
  fb			= y_high;
  //if (sign(fa)!=sign(fb))
  {
    fc			= fb;
    for (iter=0; iter<R_MAXITER; ++iter)
    {
      if (((fb>0.0)&&(fc>0.0))||((fb<0.0)&&(fc<0.0)))
      {
        c		= a;
	fc		= fa;
	e		= d = b-a;
      }
      if (fabs(fc)<fabs(fb))
      {
        a		= b;
	b		= c;
	c		= a;
	fa		= fb;
	fb		= fc;
	fc		= fa;
      }
      tol1		= 2.0*R_EPS*fabs(b)+0.5*tol;
      xm		= 0.5*(c-b);
      if ((fabs(xm)<=tol1)||(fb==0.0))
      {
	//if (fabs(fb)<1e-10)
	{
	  *root		= b;
          return 1;
	}
	//else 
	//  return 0;
      }
      if ((fabs(e)>=tol1)&&(fabs(fa)>fabs(fb)))
      {
        s		= fb/fa;
	if (a==c)
	{
	  p		= 2.0*xm*s;
	  q		= 1.0-s;
	}
	else
	{
	  q		= fa/fc;
	  r		= fb/fc;
	  p		= s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	  q		= (q-1.0)*(r-1.0)*(s-1.0);
	}  
	if (p>0.0)
	  q		= -q;
	p		= fabs(p);
	min1		= 3.0*xm*q-fabs(tol1*q);
	min2		= fabs(e*q);
	if (2.0*p<(min1<min2 ? min1 : min2))
	{
	  e		= d;
	  d		= p/q;
	}
	else
	{
	  d		= xm;
	  e		= d;
	}
      }
      else
      {
        d		= xm;
	e		= d;
      }
      a			= b;
      fa		= fb;
      if (fabs(d)>tol1)
        b		+= d;
      else
        b		+= sign(xm)*tol1;
      fb		= y(b);
      //++R_NYCALL;
    }
  }
  return 0;
}


long R_Brent(
  double y(double), double x_low, double x_high, double tol, double *root)
{
  return R_BrentInt(y, x_low, x_high, y(x_low), y(x_high), tol, root);
}


#define R_DX		1e-6
#define R_TOL		1e-12

long R_SolveNumericalInt(
  long recursive,
  double y(double), 
  double x_low, double x_high, 
  double y_low, double y_high, double dy_low, double dy_high,
  long n_mesh, long depth,
  double *roots, long n_max)
{
  long			i, nr, k, abort = 0;
  double	 	x, x_old, x_step, 
			yx, yx_old, 
			dx, dy, dy_old;
  
  if (depth<=0) 
    return 0;
  
  x_step		= (x_high-x_low)/n_mesh;
  dx			= R_DX*x_step;
  x			= x_low;
  yx			= y_low;
  dy			= dy_low;
  if (fabs(yx)<R_TOL)
  {
    *roots		= x;
    nr			= 1;
    x			+= x_step;
    yx			= y(x);
    dy			= y(x+dx)-yx;
    //R_NYCALL		+= 2;
  }
  else
    nr			= 0;
  for (i=nr ? 2 : 1; (i<=n_mesh)&&(nr<R_MAXNROOTS)&&(!abort); ++i)
  {
    // Function scan
    
    x_old		= x;
    yx_old		= yx;
    dy_old		= dy;
    if (i<n_mesh)
    {
      x			= x_low+i*x_step;
      yx		= y(x);
      dy		= y(x+dx)-yx;
      //R_NYCALL		+= 2;
    }
    else
    {
      x			= x_high;
      yx		= y_high;
      dy		= dy_high;
    }
    if (fabs(yx)>R_TOL)
    { 
      // Combined bisection and Newton-Raphson scheme

      if (sign(yx_old) != sign(yx))
      {
        //fprintf(stderr, "d = %d, x E [%g, %g], y E [%g, %g]",
	//  depth, x_old, x, yx_old, yx);
        if (R_BrentInt(y, x_old, x, yx_old, yx, R_TOL, roots+nr))
	{
	  //fprintf(stderr, ", x[%d] = %g", nr, roots[nr]);
	  ++nr;
	}
	//fprintf(stderr, "\n");
      }
      else
        if (sign(dy_old) != sign(dy))
	{
	  nr		+= k = R_SolveNumericalInt(
	  			1, y, x_old, x, yx_old, yx, dy_old, dy, 
				4, depth-1, roots+nr, n_max);
	  abort		=  k==0 ? recursive : 0;
	}
    }
    else
    {
      if (nr<n_max)
        roots[nr++]	= x;
      else
	return nr;
    }
  }
  return nr;
}


long R_SolveNumerical(
  double y(double), double x_low, double x_high, long n_mesh, long depth,
  double *roots, long n_max)
{
  long			nr;
  double		dx, y_low, y_high;

  //R_NYCALL		= 0;
  dx			= R_DX*(x_high-x_low)/n_mesh;
  x_low			+= dx;
  x_high		-= dx;
  y_low			= y(x_low);
  y_high		= y(x_high);
  dx			= R_DX*(x_high-x_low)/n_mesh;
  nr 			= R_SolveNumericalInt(0, y, x_low, x_high, 
    y_low, y_high, y(x_low+dx)-y_low, y_high-y(x_high-dx), 
    n_mesh, depth, roots, n_max);
  //fprintf(stderr, "R_NYCALL = %d\n", R_NYCALL);
  return nr;
}

