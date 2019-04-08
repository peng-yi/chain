/*
    program:    init.c
    author:     Peng Yi at MIT
    date:       October 23, 2006
    purpose:    Suggested initialization sequence
*/

#define __INIT_MODULE
#include "init.h"


void InitMols(long nmols, long noldmols)	// allocate memory for molecules
{
   if (nmols > 0)
      mol		=	(molstruct *) calloc(nmols, sizeof(molstruct));
   if (noldmols > 0)
      oldmol	=	(molstruct *) calloc(noldmols, sizeof(molstruct));
   if ( (nmols>0&&mol==NULL) || (noldmols>0&&oldmol==NULL) )
      Exit("init", "InitMols", "Out of memory.");
}


void GetCoordinates(char * filename)	// read coordinates from initial file
{
   long		i, j, id, system;
   long		nsystems, nmols, nsites;
   molstruct	*moli;
  
   Read_Conf(filename);			// read in WITHOUT unit conversion

   for (i=0; i<NMOLS; i++)		// convert coordinates to system units
      for (j=0; j<mol[i].nsites; j++)
         mol[i].p[j]		=	V_Mult(1.0/unit.LENGTH, mol[i].p+j);

   for (i=0; i<NSYSTEMS; i++) {		// convert box size to system units
      BOX[i].lx	/=	unit.LENGTH;
      BOX[i].ly	/=	unit.LENGTH;
      BOX[i].lz	/=	unit.LENGTH; 
   }

#ifdef MPI						// if parallelized
   system		=	CommRank -1 ;		// system id named by processor id

   nsites		=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {		// move mols in this processor
      if (moli->box == system) {			// up front in the mol array
         nsites		+=	moli->nsites;		// and change NMOLS and NSITES
         *oldmol	=	mol[nmols];		// such that in the following
         mol[nmols++]	=	*moli;			// this processor can only "see"
         *moli		=	*oldmol;		// these mols and sites
      }
   }
   NMOLS		=	nmols;			// # of mols in this processor
   NSITES		=	nsites;			// # of sites in this processor
   NMols[system]	=	nmols;
   NSites[system]	=	nsites;
#else							// not parallelized, but still multiple systems
   for (system=0; system<NSYSTEMS; system++) {
      NMols[system]	=	0;
      NSites[system]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if ( (system=moli->box) >= 0) {
         NMols[system]	++;				// total # of mols in certain system
         NSites[system]	+=	moli->nsites;		// total # of sites in certain system
      }
   }
#endif /* MPI */

   for (i=0; i<NSYSTEMS; i++) {
      BOX[i].lbox	=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
      BOX[i].vol	=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
      BOX[i].pres 	=	P;
      BOX[i].temp 	=	kT;
      BOX[i].drmax 	=	DRMAX;
      BOX[i].dlmax 	=	DLMAX;
      BOX[i].damax	=	DAMAX;
      BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
      BOX[i].rb		=	Rb;
      BOX[i].rv		=	Rv;
   } 

   for (i=0; i<NMOLS; i++) { 				// note NMOLS might change due to MPI
      for (j=0; j<mol[i].nsites; j++)  {
         mol[i].flags[j]	=	1;		// activate all the sites on this processor
         mol[i].parent[j]	=	j-1;		// initialize parent sites
      }
      mol[i].flip		=	0;		// flip to the original direction
      mol[i].origin		=	CenterofMass(mol+i);
   }
}


void InitSystem()
{
   if (!strcmp(INITCONF, "initconf"))			// read in coordinates from initconf
      GetCoordinates("initconf");
   else						
      Exit("init", "InitSystem", "initconf not found");

   AllSpherical();					// calculate spherical coordinates
}


void CheckSetup()					// check setup parameters
{
   if ( NTYPES>MAXNTYPES )	
      Exit("init", "CheckSetup", "NTYPES invalid!");
   if (NMOLS >MAXNMOLS)			Exit("init", "CheckSetup", "NMOLS invalid!");
   if (NSYSTEMS >MAXNSYSTEMS)		Exit("init", "CheckSetup", "NSYSTEMS invalid!");
}


void PrintInit()					// print out initial system info.
{
   long		i, j;

   PrintUnits();

   fprintf(foutput, "\nRandom number test:\n\n");
   for (i=0; i<5; i++)
      fprintf(foutput, "\trandom number #%ld:\t%f\n", i+1, ran1(seed));

   fprintf(foutput, "\nInitial box setup (system units):\n");
   for (i=0; i<NSYSTEMS; i++) {
      fprintf(foutput, "\n\tBOX[%ld]\n\n", i);
      fprintf(foutput, "\tNMols:\t%ld\n", NMols[i]);
      fprintf(foutput, "\tNSites:\t%ld\n", NSites[i]);
      fprintf(foutput, "\tLBOX:\t%f\n", BOX[i].lbox);
      fprintf(foutput, "\tLBOXX:\t%f\n", BOX[i].lx);
      fprintf(foutput, "\tLBOXY:\t%f\n", BOX[i].ly);
      fprintf(foutput, "\tLBOXZ:\t%f\n", BOX[i].lz);
      fprintf(foutput, "\tRc:\t%f\n", BOX[i].rc);
      fprintf(foutput, "\tkT:\t%f\n", BOX[i].temp);
#ifdef CELL_LIST
      fprintf(foutput, "\tCell[%ld]:\t%3d %3d %3d\n", i, M[i].x, M[i].y, M[i].z);
#endif
   }
   fprintf(foutput, "\nInitial potential energy and virial (system units):\n");
   for (i=0; i<NSYSTEMS; i++) {
      fprintf(foutput, "\n\tBOX[%ld]\n\n", i);
      fprintf(foutput, "\tlj energy (without long range):\t%f\n", v[i].lj);
      fprintf(foutput, "\tlj energy (long range part):\t%f\n", v[i].ljcorr);
      fprintf(foutput, "\tstretch energy:\t%f\n", v[i].stretch);
      fprintf(foutput, "\tbending energy:\t%f\n", v[i].bending);
      fprintf(foutput, "\ttorsion energy:\t%f\n", v[i].torsion);
      fprintf(foutput, "\ttotal energy:\t%f\n", v[i].tot);
      fprintf(foutput, "\tlj virial:\t%f\n", vir[i].lj);
      fprintf(foutput, "\tstretch virial:\t%f\n", vir[i].stretch);
      fprintf(foutput, "\ttotal virial:\t%f\n", vir[i].tot);
   }
   fflush(foutput);
}


void InitAll(char * argv[])
{
   InitMols(MAXNMOLS, MAXNMOLS);	// allocate max. memory for molecules
   InitFile();				// initialize the i/o files
   GetSetup(argv);			// read setup parameters
   CheckSetup();			// check setup parameters
 
   InitUnits();				// convert units, read in coord. will need units
   PrintSetup();			// print out setup file
   InitSystem();			// initialize systems and molecule coordinates
   					// and convert the units to system units
  
#ifdef CELL_LIST
   CL_Init();				// Cell list should be ready before Verlet list
   CL_Build();
#endif	/* CELL_LIST */

   InitForcefield();			// initialize mixing rule
   CalcV();				// calculate total energy and virial
   if (UMBRELLA) {
      New_Qlm(l_of_Ylm);		// calculate initial qlm and Qlm
      New_Clist();			// set up connected neighbor list
      Find_Nuclei(dynvar);
      Calc_Qlm(6);			// calc. Q6
   }
   InitSample(argv);			// initialize distributions
   PrintInit();				// print initial system info.
}


// below is some old stuff...

void InitConfig()   //initialize the whole system, old stuff
{
  long		i, j, k;
  FILE 		*fconf;
  molstruct	*moli;

  S_STOP	=	0;  

  if (V_SCALECUTOFF) {
     Rc	=	MIN(0.5 * LBOX, Rc / LBOX0 * LBOX);	// initialize cutoff radii, with Rc0, Rb0 and Rv0 ...
     Rb	=	Rb / LBOX0 * LBOX;			// ... the desired value at equilibrium, and LBOX0
     Rv	=	Rv / LBOX0 * LBOX;			// ... is the estimated LBOX at equilibrium
  }
#ifdef CELL_LIST
  CL_Init();
  CL_Build();
#endif

#ifdef VERLET_LIST
  vlistplus		=	calloc(2*MAXVERLETNEIGH+1, sizeof(long));	// allocate for move-affected Verlet list
  New_Vlist();				// get Verlet list for each particle, using array
  maxnverlet = 0;
#endif	/* VERLET_LIST */
}
