/*
    program:	rebridge.c
    author:	Pieter J. in 't Veld for MIT Boston
    date:	May 4, 2001
    purpose:	Trimer rebridging module, using Mavrantzas et al.,
    		Macromolecules 32, 5072 (1999).
*/
#ifndef __REBRIDGE_HEADER
#define __REBRIDGE_HEADER

#include "header.h"
/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"
#include "roots.h"
#include "random.h"
#include "types.h"
#include "globals.h"
#include "move.h"
#ifdef MPI
#include "mpicomm.h"
#endif
*/

#define B_MAXNROOTS	17

typedef struct {
  vector		r;
  sphere		s;
} position;

#endif

#ifndef __REBRIDGE_MODULE

// enter 7 positions (both cartesian and spheric)

extern long Frame(vector *, vector *, vector *, vector *, vector *, vector *);
extern long Feasible(vector *p, sphere *s);
extern long Rebridge(vector *p, sphere *s);
extern double Jacobian(vector *p, sphere *s);
extern long RebridgeSetup(
	molstruct *molm, long end, long reverse, vector *p, sphere *s);

#endif

