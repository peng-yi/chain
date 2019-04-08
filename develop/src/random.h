/*
    Name:       random.h
    Author:     Peng Yi at MIT
    Date:       October 23, 2006
    Purpose:    Using code from 'Numerical Recipes', chapter 7.
*/
#ifndef __RANDOM_HEADER
#define __RANDOM_HEADER

#ifdef __RANDOM_MODULE
/*
#include <stdio.h>
#include <math.h>
#include <time.h>
*/

#include "header.h"

#else

extern float 	ran1(long *idum);
extern double 	gauss(double, double);

#endif

#endif

