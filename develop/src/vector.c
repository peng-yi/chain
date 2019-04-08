/*
    program:	vector.c
    author:	Pieter J. in 't Veld for UT at Austin
    date:	April 1, 1999
    purpose:	3D vector algebra
*/
#define __VECTOR_MODULE
#include "vector.h"

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>

void V_Null(vector *a)
{
  a->x			= 0.0;
  a->y			= 0.0;
  a->z			= 0.0;
}


void M_Null(matrix *a)
{
  V_Null(&(a->x));
  V_Null(&(a->y));
  V_Null(&(a->z));
}


void V_Unit(vector *a)
{
  a->x			= 1.0;
  a->y			= 1.0;
  a->z			= 1.0;
}


void V_Negate(vector *a)
{
  a->x			= -a->x;
  a->y			= -a->y;
  a->z			= -a->z;
}


void M_Unit(matrix *a)
{
  M_Null(a);
  a->x.x		= 1.0;
  a->y.y		= 1.0;
  a->z.z		= 1.0;
}


void M_Negate(matrix *a)
{
  a->x.x		= -a->x.x;
  a->x.y		= -a->x.y;
  a->x.z		= -a->x.z;
  a->y.x		= -a->y.x;
  a->y.y		= -a->y.y;
  a->y.z		= -a->y.z;
  a->z.x		= -a->z.x;
  a->z.y		= -a->z.y;
  a->z.z		= -a->z.z;
}


void V_PV(vector a)
{
  printf("{%g, %g, %g}", a.x, a.y, a.z);
//  fprintf(stderr, "{%g, %g, %g}", a.x, a.y, a.z);
}


void V_Print(vector a)
{
  V_PV(a);
  printf("\n");
//  fprintf(stderr, "\n");
}


void M_Print(matrix a)
{
  printf("{"); V_PV(a.x);
  printf(", "); V_PV(a.y);
  printf(", "); V_PV(a.z);
  printf("}\n");
/*
  fprintf(stderr, "{"); V_PV(a.x);
  fprintf(stderr, ", "); V_PV(a.y);
  fprintf(stderr, ", "); V_PV(a.z);
  fprintf(stderr, "}\n");
*/
}


double V_Dot(vector *a, vector *b)
{
  return a->x*b->x+a->y*b->y+a->z*b->z;
}

vector V_Add(vector *a, vector *b)
{
  vector		c;

  c.x			= a->x+b->x;
  c.y			= a->y+b->y;
  c.z			= a->z+b->z;
  return c;
}

vector V_Subtr(vector *a, vector *b)
{
  vector		c;

  c.x			= a->x-b->x;
  c.y			= a->y-b->y;
  c.z			= a->z-b->z;
  return c;
}

vector V_Cross(vector *a, vector *b)
{
  vector		c;

  c.x			= a->y*b->z - a->z*b->y;
  c.y			= a->z*b->x - a->x*b->z;
  c.z			= a->x*b->y - a->y*b->x;
  return c;
}


vector V_Mult(double f, vector *a)
{
  vector		b;

  b.x			= f*a->x;
  b.y			= f*a->y;
  b.z			= f*a->z;
  return b;
}


void V_Swap(vector *a, vector *b)
{
  vector		c = *a;

  *a			= *b;
  *b			= c;
}


vector MV_Dot(matrix *a, vector *b)
{
  vector		c;

  c.x			= a->x.x*b->x + a->y.x*b->y + a->z.x*b->z;
  c.y			= a->x.y*b->x + a->y.y*b->y + a->z.y*b->z;
  c.z			= a->x.z*b->x + a->y.z*b->y + a->z.z*b->z;
  return c;
}


matrix M_Rotation(double alpha, double beta)
{
  double		cosa = cos(alpha), cosb = cos(beta),
  			sina = sin(alpha), sinb = sin(beta);
  matrix		R;

  R.x.x			= cosa;
  R.x.y			= cosb*sina;
  R.x.z			= sinb*sina;
  R.y.x			= -sina;
  R.y.y			= cosb*cosa;
  R.y.z			= sinb*cosa;
  R.z.x			= 0.0;
  R.z.y			= -sinb;
  R.z.z			= cosb;
  return R;
}


vector V_Rotate(vector *a, double alpha, double beta)
{
  matrix		R = M_Rotation(alpha, beta);

  return MV_Dot(&R, a);
}


matrix M_Transpose(matrix *x)
{
  matrix		m;

  m			= *x;
  m.x.y			= x->y.x;
  m.x.z			= x->z.x;
  m.y.x			= x->x.y;
  m.y.z			= x->z.y;
  m.z.x			= x->x.z;
  m.z.y			= x->y.z;
  return m;
}


matrix M_Dot(matrix *x, matrix *y)
{
  matrix		m;

  m.x			= MV_Dot(x, &(y->x));
  m.y			= MV_Dot(x, &(y->y));
  m.z			= MV_Dot(x, &(y->z));
  return m;
}


double M_Det(matrix *x)
{
  vector		y = V_Cross(&(x->x), &(x->y));
  
  return V_Dot(&y, &(x->z));
}


int M_Inverse(matrix *a, matrix *b)
{
  double		f = M_Det(a);

  if (fabs(f)<1e-14) return V_E_SINGULAR;
  b->x			= V_Cross(&(a->y), &(a->z));
  b->y			= V_Cross(&(a->z), &(a->x));
  b->z			= V_Cross(&(a->x), &(a->y));
  b->x			= V_Mult(f = 1.0/f, &(b->x));
  b->y			= V_Mult(f, &(b->y));
  b->z			= V_Mult(f, &(b->z));
  *b			= M_Transpose(b);
  return 0;
}


matrix M_Mult(double f, matrix *x)
{
  matrix		y;

  y.x			= V_Mult(f, &(x->x));
  y.y			= V_Mult(f, &(x->y));
  y.z			= V_Mult(f, &(x->z));
  return y;
}


matrix M_Add(matrix *x, matrix *y)
{
  matrix		z;

  z.x			= V_Add(&(x->x), &(y->x));
  z.y			= V_Add(&(x->y), &(y->y));
  z.z			= V_Add(&(x->z), &(y->z));
  return z;
}


matrix M_Subtr(matrix *x, matrix *y)
{
  matrix		z;
  
  z.x			= V_Subtr(&(x->x), &(y->x));
  z.y			= V_Subtr(&(x->y), &(y->y));
  z.z			= V_Subtr(&(x->z), &(y->z));
  return z;
}


matrix M_XNormalRotate(vector *x)
{
  matrix		m;
  vector		r;
  double		f = V_Dot(x, x), 
  			sina, cosa, sinb, cosb;

  M_Unit(&m);
  if (f)
  {
    r			= V_Mult(1.0/sqrt(f), x);
    f			= (r.y<0) ? -sqrt(r.y*r.y+r.z*r.z) :
    				sqrt(r.y*r.y+r.z*r.z);
    if (f)
    {
      sina		= f;
      cosa		= r.x;
      sinb		= r.z/f;
      cosb		= -r.y/f;
      m.x.x		= cosa;
      m.x.y		= -sina*cosb;
      m.x.z		= sina*sinb;
      m.y.x		= sina;
      m.y.y		= cosa*cosb;
      m.y.z		= -cosa*sinb;
      m.z.x		= 0;
      m.z.y		= sinb;
      m.z.z		= cosb;
    }
  }
  return m;
}


// Determines the rotation matrix for rotating x to y

matrix M_Rotate(vector *x, vector *y)
{
  matrix		m1, m2;

  m1			= M_XNormalRotate(x);
  m2			= M_XNormalRotate(y);
  m1			= M_Transpose(&m1);
  m1			= M_Dot(&m2, &m1);
  return m1;
}


// Determines the orientation of a plane made up out of x and y

matrix M_Orientation(vector *x, vector *y)
{
  matrix		m;

  m.x			= V_Mult(1.0/sqrt(V_Dot(y, y)), y);
  m.z			= V_Cross(x, y);
  m.z			= V_Mult(1.0/sqrt(V_Dot(&(m.z),&(m.z))), &(m.z));
  m.y			= V_Cross(&(m.z), &(m.x));
  return m;
}


vector M_eig(matrix M)	// only for 3x3 REAL matrix, from Numerical Recipes
{			// note: real matrix can still have complex eigenvalues
   double	a, b, c, temp;
   double	Q, R, theta;
   double complex	A, B;
   vector	eigenvalue;

   a	=	-(M.x.x + M.y.y + M.z.z);
   b	=	M.x.x * M.y.y + M.x.x * M.z.z + M.y.y * M.z.z - M.z.y * M.y.z - M.x.y * M.y.x - M.x.z * M.z.x;
   c	=	-1.0* M.x.x * (M.y.y * M.z.z - M.y.z * M.z.y) 
		+ M.x.y * (M.y.x * M.z.z - M.y.z * M.z.x) 
		- M.x.z * (M.y.x * M.z.y - M.z.x * M.y.y);

   Q	=	(a*a - 3*b) /9;
   R	=	(2*a*a*a - 9*a*b + 27*c) / 54;

   if (R*R < Q*Q*Q) {
      theta	=	acos(R/sqrt(Q*Q*Q));
      eigenvalue.x	=	-2.0 * sqrt(Q) * cos(theta/3) - a/3;
      eigenvalue.y	=	-2.0 * sqrt(Q) * cos((theta+2*M_PI)/3) - a/3;	
      eigenvalue.z	=	-2.0 * sqrt(Q) * cos((theta-2*M_PI)/3) - a/3;

      //printf("%f\t%f\t%f\n", eigenvalue.x, eigenvalue.y, eigenvalue.z);
   }
   else {

      if (R<0)
         A	=	-cpow(R-sqrt(R*R-Q*Q*Q), 1.0/3);
      else
	 A	=	-cpow(R+sqrt(R*R-Q*Q*Q), 1.0/3);

      if (fabs(A)<1e-14)
         B	=	0;
      else
	 B	=	Q/A;

      //printf("%f\t%f\n", creal(A+B-a/3), cimag(A+B-a/3));
      //printf("%f\t%f\n",creal( -0.5*(A+B) - a/3 + 1.732050807569 * 0.5 * (A-B) * I ), cimag( -0.5*(A+B) - a/3 + 1.732050807569 * 0.5 * (A-B) * I ));
      //printf("%f\t%f\n",creal( -0.5*(A+B) - a/3 - 1.732050807569 * 0.5 * (A-B) * I ), cimag( -0.5*(A+B) - a/3 - 1.732050807569 * 0.5 * (A-B) * I ));

      eigenvalue.x	=	creal( A+B-a/3 );
      eigenvalue.y	=	creal( -0.5*(A+B) - a/3 + 1.732050807569 * 0.5 * (A-B) * I );
      eigenvalue.z	=	creal( -0.5*(A+B) - a/3 - 1.732050807569 * 0.5 * (A-B) * I );
   }
   return	eigenvalue;
}


vector V_eig(matrix M, double eigvalue) // the eigenvector with eigvalue
{
   vector       eigvector;
   double       xx, xy, xz, yx, yy, yz, zx, zy, zz;
   
   float	zero=1e-14;		// in case off-diagonal element very small
   
   xx   =       M.x.x;          xy      =       M.x.y + zero;   xz      =       M.x.z + zero;
   yx   =       M.y.x + zero;   yy      =       M.y.y;          yz      =       M.y.z + zero;
   zx   =       M.z.x + zero;   zy      =       M.z.y + zero;   zz      =       M.z.z;

   eigvector.x  =       1.0;
   eigvector.y  =       -(zx * xz - (zz-eigvalue)*(xx-eigvalue))/(zy*xz-(zz-eigvalue)*xy);
   eigvector.z  =       -((xx-eigvalue)/xz + xy/xz*eigvector.y);

   return       eigvector;
}

//======================================================//
//	eigensol(): solve eigenvalues and eigenvectors	//
//		    of a 3x3, real, symmetric matrix	//
//======================================================//
void eigensol(matrix *M, vector *e_val, matrix *e_vec)	// coded on 4/11/2013, only for 3x3 matrix
{							// output eval in asending abs value
   int		i, j, k;				// evec corresponds to eval

   gsl_matrix *m	=	gsl_matrix_alloc (3, 3);
   gsl_matrix_set (m, 0, 0, M->x.x);
   gsl_matrix_set (m, 0, 1, M->x.y);
   gsl_matrix_set (m, 0, 2, M->x.z);
   gsl_matrix_set (m, 1, 0, M->y.x);
   gsl_matrix_set (m, 1, 1, M->y.y);
   gsl_matrix_set (m, 1, 2, M->y.z);
   gsl_matrix_set (m, 2, 0, M->z.x);
   gsl_matrix_set (m, 2, 1, M->z.y);
   gsl_matrix_set (m, 2, 2, M->z.z);
	 
   gsl_vector *eval		=	gsl_vector_alloc (3);
   gsl_matrix *evec		=	gsl_matrix_alloc (3, 3);
   gsl_eigen_symmv_workspace *w	=	gsl_eigen_symmv_alloc(3);
   gsl_eigen_symmv(m, eval, evec, w);
   gsl_eigen_symmv_free(w);
   gsl_eigen_symmv_sort(eval, evec, GSL_EIGEN_SORT_ABS_ASC);	// asending in abs values of eval

   e_val->x	=	gsl_vector_get(eval, 0);
   e_val->y	=	gsl_vector_get(eval, 1);
   e_val->z	=	gsl_vector_get(eval, 2);

   e_vec->x.x	=	gsl_matrix_get(evec, 0, 0); 	// column vector is e_vector
   e_vec->x.y	=	gsl_matrix_get(evec, 1, 0);    
   e_vec->x.z	=	gsl_matrix_get(evec, 2, 0);    
   e_vec->y.x	=	gsl_matrix_get(evec, 0, 1);     
   e_vec->y.y	=	gsl_matrix_get(evec, 1, 1);     
   e_vec->y.z	=	gsl_matrix_get(evec, 2, 1);     
   e_vec->z.x	=	gsl_matrix_get(evec, 0, 2);     
   e_vec->z.y	=	gsl_matrix_get(evec, 1, 2);     
   e_vec->z.z	=	gsl_matrix_get(evec, 2, 2);     

   return;
} 

//======================================================================//
//	Rot_matrix(): rotation matrix of angle theta about axis v.	//
//		      4/15/2013						//
//======================================================================//
matrix matrix_rot(vector *v, float theta)
{
   matrix	R;
   float	costheta, sintheta;
   float	r, r2;
   vector	u;

   r2		=	v->x * v->x + v->y * v->y + v->z * v->z;

   if ( fabsf(r2-1) > FLT_EPSILON) {		// not unit vector
      r		=	sqrtf(r2);
      u.x	=	v->x/r;
      u.y	=	v->y/r;
      u.z	=	v->z/r;
   }
   else {
      u.x	=	v->x;
      u.y	=	v->y;
      u.z	=	v->z;
   }

   costheta	=	cos(theta);
   sintheta	=	sin(theta);

   R.x.x	=	costheta + u.x * u.x * (1-costheta);
   R.y.y	=	costheta + u.y * u.y * (1-costheta);
   R.z.z	=	costheta + u.z * u.z * (1-costheta);

   R.x.y	=	u.x * u.y * (1-costheta) - u.z * sintheta;
   R.x.z	=	u.x * u.z * (1-costheta) + u.y * sintheta;
   R.y.z	=	u.y * u.z * (1-costheta) - u.x * sintheta;

   R.y.x	=	u.y * u.x * (1-costheta) + u.z * sintheta;
   R.z.x	=	u.z * u.x * (1-costheta) - u.y * sintheta;
   R.z.y	=	u.z * u.y * (1-costheta) + u.x * sintheta;
   return	R;
}

//======================================================================//
//	angle2vector(): calculate x,y,z component of unit vector 	//
//			from angles theta and phi			//
//======================================================================//
vector angle2vector(float theta, float phi)
{
   vector	r;
   r.x	=	sin(theta)*cos(phi);
   r.y	=	sin(theta)*sin(phi);
   r.z	=	cos(theta);
   return	r;
}
