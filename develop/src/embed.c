/*
	program:	embed.c
	author:		Peng Yi at MIT
	date:		May 18, 2008
	purpose:	embed a smaller system to the center of a larger system, do NO analysis
	note:		require setup file
*/

#define __MAIN_PROGRAM
#include "header.h"


int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   char		s[255], ff[255], dummy[80];
   molstruct	*moli, *mol_in;
   FILE		*fPtr;

   vector	com,
		Linner[MAXNSYSTEMS];		// inner BOX dimension
   long		i, j, site_in, n, system=0,
		nsystems,			// for one system only, 5/18/2008
		nmols1, nsites1,		// NMOLS and NSITES for inner box
		nmolstot, nsitestot,		// total NMOLS and NSITES
		id;
   double	lambda;
   long		overlap;

   if (argc<4) {
      printf("embed (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tembed file1 file2 lambda\n\n");
      printf("Notes:\n");
      printf("\t* file1 is the configuration file for inner box, file2 for outer box\n\n");
      printf("\t* lambda <= 1.0 is the overlap tuning parameter\n\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }
   lambda	=	atof(argv[3]);

   InitMols(MAXNMOLS, 1);	// allocate memory for molecules, MAXNMOLS must > the sum of NMOLS
   GetSetup(argv);		// mainly PBC I think

   // Read in the conf. of the inner system

   Read_Conf(argv[1]);
   for (moli=mol; moli<mol+NMOLS; moli++)	// arrange the inner system
      MolInBox2(moli);

   nsystems	=	NSYSTEMS;		// save inner system information
   nmols1	=	NMOLS;
   nsites1	=	NSITES;
   nmolstot	=	NMOLS;
   nsitestot	=	NSITES;

   for (i=0; i<NSYSTEMS; i++) {
      Linner[i].x	=	BOX[i].lx;
      Linner[i].y	=	BOX[i].ly;
      Linner[i].z	=	BOX[i].lz;
   }

   // Read in the conf. of the outer system

   if (!(fPtr=fopen(argv[2], "r")))		// open the conf. file of the outer system
      Exit("embed", "main", "big system conf. file failed to open.");

   fscanf(fPtr, "%s%ld", dummy, &TIMESTEP);			// read in timestep
   fscanf(fPtr, "%ld%ld%ld", &NSYSTEMS, &NMOLS, &NSITES);	// update NSYSTEMS, NMOLS, NSITES
   if (nsystems != NSYSTEMS) {}					// might do something in the future
   for (i=0; i<NSYSTEMS; i++)
      fscanf(fPtr, "%lf%lf%lf", &(BOX[i].lx), &(BOX[i].ly), &(BOX[i].lz));	// update BOX dimension

   moli	=	mol + nmols1;			// clear statement, although not necessary
   for (i=0; i<NMOLS; i++) {			// start reading in the molecules in the outer box
      fscanf(fPtr, "%ld%ld%ld", &id, &(moli->box), &(moli->nsites));
      for (j=0; j<moli->nsites; j++)
         fscanf(fPtr, "%ld%lf%lf%lf", &(moli->type[j]), &(moli->p[j].x), &(moli->p[j].y), &(moli->p[j].z));

      MolInBox2(moli);				// map back in outer box

      overlap	=	0;
      for (j=0; j<moli->nsites; j++) {		// check each bead in the outer box
         for (mol_in=mol; mol_in<mol+nmols1; mol_in++) {
            for (site_in=0; site_in<mol_in->nsites; site_in++) {
	       if (DistSQ(moli->p[j], mol_in->p[site_in], system) < type[0].SIGMA * type[0].SIGMA * lambda) {
		  overlap	=	1;
		  break;
	       }
	    }
	    if (overlap)
		break;
	 }
	 if (overlap)
	    break;
      }
      if (!overlap) {
         nmolstot	++;
	 nsitestot	+=	moli->nsites;
	 moli	++;
      }

/*
      overlap	=	0;			// is this molecule overlap with inner box?
      for (j=0; j<moli->nsites; j++) {		// check every bead
         if (fabs(moli->p[j].x)*lambda < 0.5*Linner[moli->box].x  && 
		fabs(moli->p[j].y)*lambda < 0.5*Linner[moli->box].y && 
		fabs(moli->p[j].z)*lambda < 0.5*Linner[moli->box].z ) {
            overlap	=	1;
            break;
	 }
      }
      if (!overlap) {
         nmolstot	++;
	 nsitestot	+=	moli->nsites;
	 moli	++;
      }
*/
/*
      com	=	CenterofMass(moli);	// center of mass of one chain
      com.x	*=	lambda;		// give some gap
      com.y	*=	lambda;
      com.z	*=	lambda;

      if (!(com.x < 0.5*Linner[moli->box].x && com.x >= -0.5*Linner[moli->box].x &&
	   com.y < 0.5*Linner[moli->box].y && com.y >= -0.5*Linner[moli->box].y &&
	   com.z < 0.5*Linner[moli->box].z && com.z >= -0.5*Linner[moli->box].z) ) {	// com NOT in the inner box
         nmolstot	++;
	 nsitestot	+=	moli->nsites;	
	 moli		++;
      }
*/
   }
   NMOLS	=	nmolstot;
   NSITES	=	nsitestot;

   // PRINT OUT CONFIGURATION

   system	=	0;			// only one system for now
   unit.LENGTH	=	1.0;			// not doing any analysis
   Write_Conf(-1);				// output combined configuration

   return	0;
}
