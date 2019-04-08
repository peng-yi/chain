/*
    program:	units.c
    author:	Peng Yi at MIT
    date:	Nov. 7, 2007
    purpose:	Units conversion from SI to program units
*/
#define __UNITS_MODULE
#include "units.h"

/*
			SI units used in input		Program units

	energy:		kJ/mol				epsilon[0]
	length:		Angstrom			sigma[0]
	temperature:	Kelvin				epsilon[0]
	pressure:	atm				epsilon[0]/sigma[0]^3
	angle:		degree				radian
	mass:		g/mol				1
	density:	g/cm^3				1/sigma[0]^3
*/

double CalcMass()		// calculate average mass of each site
{
   molstruct	*moli;
   long		i;
   double	avemass = 0.0;
/*
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         avemass	+=	type[moli->type[i]].M;

   avemass	/=	NSITES;
   return	avemass;
*/

   // We don't use above calc. because we first need units to read in coord.
   // This is for now, because I don't want to assume that I know how many types are there

   if (1==NTYPES)
      return	type[0].M;

   if (NTYPES>=2) {
      avemass	=	2 * NMOLS * type[0].M + (NSITES - 2* NMOLS) * type[1].M;
      avemass	/=	NSITES;
      return	avemass;
   }
   return	0;
}


void CalcUnits(long flag)	// calc. the link b/w system units and SI units
{
   switch(flag){
      case 1: 
      	unit.MASS	=	CalcMass();	// ave mass of each site in g/mol
      	unit.LENGTH	=	type[0].SIGMA;
      	unit.ENERGY	=	type[0].EPSILON;
      	unit.TEMP	=	type[0].EPSILON / (N_KB * N_NAV * 1e-3);
      	unit.PRESSURE	=	type[0].EPSILON / pow(type[0].SIGMA, 3) / (N_ATM * N_NAV * 1e-33);
      	unit.ANGLE	=	180.0 / M_PI;
      	unit.DENSITY	=	unit.MASS / pow(unit.LENGTH, 3) / (N_NAV * 1e-24);
        break;

      case 0:
        unit.MASS	=	1.0;
  	unit.LENGTH	=	1.0;
      	unit.ENERGY	=	1.0;
      	unit.TEMP	=	1.0;
      	unit.PRESSURE	=	1.0;
      	unit.ANGLE	=	1.0;
      	unit.DENSITY	=	1.0;
        break;

      default:
        break;
   }
}


void PrintUnits()
{
   fprintf(foutput, "\nUnits:\n\n");
   fprintf(foutput, "\tMASS (g/mol):		\t%f\n", unit.MASS);
   fprintf(foutput, "\tLENGTH (Angstrom):	\t%f\n", unit.LENGTH);
   fprintf(foutput, "\tENERGY (kJ/mol): 	\t%f\n", unit.ENERGY);
   fprintf(foutput, "\tTEMP (K):		\t%f\n", unit.TEMP);
   fprintf(foutput, "\tPRESSURE (atm):		\t%f\n", unit.PRESSURE);
   fprintf(foutput, "\tANGLE (degree):		\t%f\n", unit.ANGLE);
   fprintf(foutput, "\tDENSITY (g/cm^3):	\t%f\n", unit.DENSITY);
}


void ConvertSIToSystem()		// convert read-in data to system units
{
   long		i, j;
  
   kT		/=	unit.TEMP;
   P		/=	unit.PRESSURE;
   Rho		/=	unit.DENSITY;

   for (i=0; i<NTYPES; i++) {
      type[i].M		/=	unit.MASS;
      type[i].SIGMA	/=	unit.LENGTH;
      type[i].EPSILON	/=	unit.ENERGY;
      type[i].LSTRETCH	/=	unit.LENGTH;
      type[i].KSTRETCH	/=	unit.ENERGY;	// need to be fixed/improved ...
      type[i].THETA	/=	unit.ANGLE;
      type[i].KBENDING	/=	unit.ENERGY;
      type[i].HS	/=	unit.LENGTH;

      for (j=0; j<6; j++)
         type[i].TORSION[j]	/=	unit.ENERGY;
   }
}


void ConvertSystemToSI()		// convert output data to SI units
{
   molstruct	*moli;
   long		i, system;

   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->p[i].x		*=	unit.LENGTH;
         moli->p[i].y		*=	unit.LENGTH;
         moli->p[i].z		*=	unit.LENGTH;
         moli->s[i].d		*=	unit.LENGTH;
         moli->s[i].alpha	*=	unit.ANGLE;
         moli->s[i].beta	*=	unit.ANGLE;
      }
   }
}


void CoorSI2System()		// box size and particle coordinates unit conversion
{
   molstruct	*moli;
   long		i;

   for (i=0; i<NSYSTEMS; i++) {
      BOX[i].lx		/=	unit.LENGTH;
      BOX[i].ly		/=	unit.LENGTH;
      BOX[i].lz		/=	unit.LENGTH;
      BOX[i].lbox	/=	unit.LENGTH;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->p[i].x	/=	unit.LENGTH;
         moli->p[i].y	/=	unit.LENGTH;
         moli->p[i].z	/=	unit.LENGTH;
         moli->s[i].d	/=	unit.LENGTH;
         moli->s[i].alpha	/=	unit.ANGLE;
         moli->s[i].beta	/=	unit.ANGLE;
      }
   }
   return;
}


void CoorSystem2SI()		// box size and particle coordinates unit conversion
{
   molstruct	*moli;
   long		i;

   for (i=0; i<NSYSTEMS; i++) {
      BOX[i].lx		*=	unit.LENGTH;
      BOX[i].ly		*=	unit.LENGTH;
      BOX[i].lz		*=	unit.LENGTH;
      BOX[i].lbox	*=	unit.LENGTH;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->p[i].x	*=	unit.LENGTH;
         moli->p[i].y	*=	unit.LENGTH;
         moli->p[i].z	*=	unit.LENGTH;
         moli->s[i].d	*=	unit.LENGTH;
         moli->s[i].alpha	*=	unit.ANGLE;
         moli->s[i].beta	*=	unit.ANGLE;
      }
   }
   return;
}  


void InitUnits()
{
   CalcUnits(ConvertUnits);
   ConvertSIToSystem();
}
