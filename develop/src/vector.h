/*
    program:	vector.h
    author:	Pieter J. in 't Veld for UT at Austin
    date:	April 1, 1999
    purpose:	3D vector algebra
*/

#ifndef __VECTOR_HEADER
#define __VECTOR_HEADER

#define V_E_SINGULAR	1
/*
typedef struct {
  double	x, y, z; } vector;
typedef struct {
  vector	x, y, z; } matrix;
*/
#ifdef __VECTOR_MODULE

#include "header.h"
/*
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "types.h"
*/
#else

// Vector operations

extern void V_Null(vector *);
extern void V_Unit(vector *);
extern void V_Negate(vector *);
extern double V_Dot(vector *, vector *);
extern vector V_Add(vector *, vector *);
extern vector V_Subtr(vector *, vector *);
extern vector V_Cross(vector *, vector *);
extern vector V_Mult(double, vector *);
extern vector V_Rotate(vector *, double, double);
extern void V_Print(vector);

// Matrix operations

extern void M_Null(matrix *);
extern void M_Unit(matrix *);
extern void M_Negate(matrix *);
extern matrix M_Dot(matrix *, matrix *);
extern matrix M_Add(matrix *, matrix *);
extern matrix M_Subtr(matrix *, matrix *);
extern matrix M_Mult(double, matrix *);
extern matrix M_Rotation(double, double);
extern matrix M_XNormalRotate(vector *);
extern matrix M_Rotate(vector *, vector *);
extern matrix M_Orientation(vector *, vector *);
extern matrix M_Transpose(matrix *);
extern int M_Inverse(matrix *, matrix *);
extern double M_Det(matrix *);
extern void M_Print(matrix);
extern vector M_eig(matrix);
extern vector V_eig(matrix, double);

extern void eigensol(matrix *, vector *, matrix *);	// eigen-solver
extern matrix matrix_rot(vector *v, float theta);	// rotation matrix
extern vector angle2vector(float theta, float phi);

// Hybrid operarions

extern vector MV_Dot(matrix *, vector *);


#endif

#endif

