/*
    program:	units.h
    author:	Peng Yi at MIT
    date:	Nov. 7, 2007
    purpose:	Convert general SI units to program units
*/

#ifndef __UNITS_HEADER
#define __UNITS_HEADER

#define	N_NAV		6.0221367e23	/* Avogadro constant 1/mol */
#define	N_KB		1.380658e-23	/* Boltzmann constant J/K */
#define N_ATM		101325.024	/* pascals/atm */

typedef struct {
   double	LENGTH,			// SI length = system length * LENGTH
		ENERGY,
		PRESSURE,
		TEMP,
		ANGLE,
		MASS,
		DENSITY;
} unitstruct;

#ifdef __UNITS_MODULE

#include "header.h"

unitstruct		unit;
long			ConvertUnits;

#else

extern unitstruct	unit;
extern long		ConvertUnits;

extern double		CalcMass();
extern void		CalcUnits(long flag);
extern void		ConvertSIToSystem();
extern void		ConvertSystemToSI();
extern void		CoorSystem2SI();
extern void		CoorSI2System();
extern void		InitUnits();
extern void		PrintUnits();

#endif

#endif

