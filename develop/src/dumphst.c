/*
    program:    dumphst.c
    author:     Peng Yi at MIT
    date:       January 12, 2008
    purpose:    get info from dump files and do analysis
    note:	previously I have a program hst.c to deal with the binary history file
		produced by Pieter's code.  However later I realize that I didn't need
		that much information to be stored in the history file.  Even it is 
		in binary format, it still take a huge amount of disk space.
		I also noticed that lammps dumps in human-readable format.  Of course 
		as an MD code, lammps can dump coordinates as velocity and force and 
		etc.  So I decide to only dump coordinates in my MC code and in human-
		readable format.  The analysis program is thus change to dumphst.c
*/

#define __MAIN_PROGRAM
#include "header.h"

#include "correlation.h"

long	timestep = 0,
	min=MAXNMOLS, 
	max=0,
	nnmax[MAXNMOLS],		// nnmax[i] is the counts of conf. with nmax=i
	n[1000][MAXNMOLS];		// n[i][j] is the counts of nuclei of size j when nmax=i
double	pn[1000][MAXNMOLS],		// pn[i][j] is the prob. of nuclei of size j when nmax=i
	Gn[MAXNMOLS],			// free energy of one nucleus
	G[MAXNMOLS];			// free energy of the system

// some functions declared here
void printout(FILE *);
void hst2conf();
void hst2lammpsdump(FILE *fPtr);

long Dump_GetNext(FILE *fPtr)		// get next configuration from dump file, and primary processing
{					// very much like GetCoordinate() function in io.c
   char			name[80], dummy[80];
   long			i, j, id, system;
   molstruct		*moli;

   // read in data

   fscanf(fPtr, "%s", name);
   fscanf(fPtr, "%ld", &timestep);			// read in timestep
   fscanf(fPtr, "%ld%ld%ld", &NSYSTEMS, &NMOLS, &NSITES);
   NBOX	=	NSYSTEMS;		// temporary
   for (i=0; i<NSYSTEMS; i++)
      fscanf(fPtr, "%lf%lf%lf", &(BOX[i].lx), &(BOX[i].ly), &(BOX[i].lz));
   for (i=0; i<NMOLS; i++) {
      fscanf(fPtr, "%ld%ld%ld", &id, &(mol[i].box), &(mol[i].nsites));
      //fscanf(fconf, "%ld%ld%ld%ld%ld", &id, &(mol[i].box), &(mol[i].nsites), &(mol[i].fix), &(mol[i].flip));
      for (j=0; j<mol[i].nsites; j++)
         fscanf(fPtr, "%ld%lf%lf%lf", &(mol[i].type[j]), &(mol[i].p[j].x), &(mol[i].p[j].y), &(mol[i].p[j].z));
   }
   fgets(dummy, sizeof(dummy), fPtr);

   // some primary data processing

   CoorSI2System();				// the dump file is in real units (8/7/08)
   for (i=0; i<NSYSTEMS; i++) {
      NMols[i]	=	0;
      NSites[i]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if ( (system=moli->box) >= 0) {
         NMols[system]	++;				// total # of mols in certain system
         NSites[system]	+=	moli->nsites;		// total # of sites in certain system
      }
   }
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

   // return 

   return feof(fPtr)==0;
}


int main(int argc, char * argv[])
{
   FILE		*fhist, *fp, *fp1, *fsize, *fconf, *frdf, *flammps;
   long		Nsplit, 
		Nblock = -1;		// # of records in each block, could be set by command line
   long		version, system;
   char		file_hst1[256],
		comments[1024];
   double	dummy;
   long		i, j, k; 
   molstruct	*moli;
   double	var, bias; 
   static long	init = 1;

   if (argc<2) {
      printf("dumphst (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tdumphst dumpfile [N]\n\n");
      printf("Notes:\n");
      printf("\t* require setup file to get e.g. Temp for analysis\n\n");
      printf("\t\t if N is given, then do analysis in blocks each with N records, otherwise as one big block\n\n");
      exit(-1);
   }

   if (!(fp=fopen(argv[1], "r"))) {			// open dump file
      printf("dumphst: dump file cannot open!\n");
      exit(-1);
   }
   if (argc == 3)
      Nblock		=	atol(argv[2]);		// otherwise = -1

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate max. memory for molecules
   GetSetup(argv);			// most of the variables read from setup 
					// will be updated by the reading from history file
   InitUnits();				// calculate the units needed for system<->SI conversion
					// and convert setup variable to system units 
   InitForcefield();
 
   if (!(fsize = fopen("Psize","w"))) {			// open output file for size distribution information
      printf("dumphst: Size prob. distribution file failed to open");
      exit(-1);
   }
   if (!(flammps=fopen("lammps.dump", "w"))) {		// to convert to lammps dump file for lammps2pdb
      printf("dumphst: lammps dump file failed to open.\n");
      exit(-1);
   }

   while (Dump_GetNext(fp)) {			// read in next record/configuration, including init lattice

#ifdef CELL_LIST
      if (init) {
         CL_Init();	 init	=	0;
      }
      else
	 CL_Destroy();
      CL_Build();
#endif
      CalcV();					// need cell list

      if (dynvar==3)
         SampleP2All();				// calculate global P2 and local p2
      SampleSpherical();
      Find_Nuclei();

printf("%d %f %f %d\n", timestep, v[0].tot/NSITES, NSITES/BOX[0].vol, nmax[0][0]);

      // calculate bias for each record/configuration 

      if (dynvar==1)
         var	=	MAXSIZE[0];
      else if (dynvar==3)
	 var	=	P2M[0];

      bias	=	exp(0.5*kP2*(var-P2middle)*(var-P2middle));

      // collect probability distribution of cluster size 
 
      nnmax[MAXSIZE[0]]	++;
      for (i=1; i<=MAXSIZE[0]; i++) {
	 pn[0][i]	+=	bias * sizedist[i];		// probability distribution of cluster size
	 n[0][i]	+=	sizedist[i];

         pn[MAXSIZE[0]][i]	+=	bias * sizedist[i];
         n[MAXSIZE[0]][i]	+=	sizedist[i];
      }
      if (MAXSIZE[0] > max)
	 max	=	MAXSIZE[0];
      if (MAXSIZE[0] < min) 
         min	=	MAXSIZE[0];

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
      } 
      if (!(timestep%IRADIAL))
         radial("sample");		// sample rdf

      correlation();			// calculate correlation function

      //Write_Conf(timestep);		// convert to ASCII dump file

      hst2lammpsdump(flammps);		// convert to lammps dump format to use lammps2pdb in the future

      timestep	+=	ITAPE;		// because the binary hist file was stored every ITAPE cycles

      if (Nblock!=-1 && timestep%Nblock==0) {		// summarize the results for this block and output
         for (i=1; i<=max; i++) 
            Gn[i]	=	-log(pn[0][i]);
         for (i=min; i<=max; i++) {
            for (j=1; j<=i; j++) {
               G[i]	+=	n[i][j] * Gn[j];
            }
         }

	 printout(fsize);

         // clear up the counters
         for (i=0; i<MAXNMOLS; i++) {
	    nnmax[i]	=	0;
	    Gn[i]	=	0.0;
	    G[i]	=	0.0;
	 }
         for (i=0; i<100; i++) { 
            for (j=0; j<MAXNMOLS; j++) {
               pn[i][j]	=	0.0;
	       n[i][j]	=	0;
            }
         }
         max	=	0;
	 min	=	MAXNMOLS;
      }
   }

   /* normalize correlation functions and print out */

   corr_norm();
   corr_print();
   radial("print");

   // summarize results for all records and print out size distribution

   if (Nblock==-1) {
      for (i=1; i<=max; i++) 
         Gn[i]	=	-log(pn[0][i]);			// free energy of nucleus
      for (i=min; i<=max; i++) {
         for (j=1; j<=i; j++) {
            G[i]	+=	n[i][j] * Gn[j];	// free energy of the system
         }
      }
      printout(fsize);
   }

   fclose(flammps);
   fclose(fp);
   fclose(fsize);
   return EXIT_SUCCESS;
}


void printout(FILE *fPtr)
{
   long		i, j;

   fprintf(fPtr, "\nOverall nuclei size distribution P(n):\n\n");
   fprintf(fPtr, "n    N(n)   fixed p(n)   -log(fixed pn[i])\n\n");
   for (i=1; i<=max; i++) 
      fprintf(fPtr, "%4d   %8d   %14.4f   %f\n", i, n[0][i], pn[0][i], Gn[i]);

   fprintf(fPtr, "P(n) for copy purpose...\n");
   for (i=1; i<=max; i++)
      fprintf(fPtr, "%f\n", Gn[i]);

   fprintf(fPtr, "\nDistribution of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, "n_max    N(n_max)\n");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%4d  %8d\n", i, nnmax[i]);

   fprintf(fPtr, "\nSize distribution N(n) for different value of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, " n_max");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%8d  ", i);
   fprintf(fPtr, "\n   n \n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%8d  ", n[j][i]);
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n  Normalized to N(n_max)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", (nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n  -Logrithm of N(n)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", -log(nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
}


void hst2conf()
{
   long		i;
   molstruct	*moli;
   static FILE	*fPtr;
   static long	init =1;

   vector	rA, rB, rAB;
   double	r2, L;
   long		id, print=0;

   if (init) {
      if (!(fPtr=fopen("hst.conf", "w")) )
          Exit("hst", "main", "open conf file failed.");
      init	=	0;
   }
/*
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (sizeofnucl[moli->nuclid[0]]==MAXSIZE[0]) {
         id	=	moli->nuclid[0];
	 rA	=	CenterofMass(moli);
         break;
      }
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (moli->nuclid[0] == id) {
         rB	=	CenterofMass(moli);
	 rAB	=	V_Subtr(&rA, &rB);
         r2	=	V_Dot(&rAB, &rAB);
         L	=	MIN(BOX[0].lx, MIN(BOX[0].ly, BOX[0].lz));
         if (r2 > L * L/2) {
            print	=	1;
	    break;
         }
      }

   if (!print)
      return;
*/
   fprintf(fPtr, "!TIMESTEP %d Nmax = %d 2ndNmax = %d\n", timestep, MAXSIZE[0], secondNmax[0]);
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);

   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx*unit.LENGTH, BOX[i].ly*unit.LENGTH, BOX[i].lz*unit.LENGTH);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
   }
}

void hst2lammpsdump(FILE *fPtr)
{
   long		i;
   molstruct	*moli;

   fprintf(fPtr, "ITEM: TIMESTEP\n");	
   fprintf(fPtr, "%d\n", timestep);		
   fprintf(fPtr, "ITEM: NUMBER OF ATOMS\n");
   fprintf(fPtr, "%d\n", NSITES);
   fprintf(fPtr, "ITEM: BOX BOUNDS\n");
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lx*unit.LENGTH, 0.5*BOX[0].lx*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].ly*unit.LENGTH, 0.5*BOX[0].ly*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lz*unit.LENGTH, 0.5*BOX[0].lz*unit.LENGTH); 
   fprintf(fPtr, "ITEM: ATOMS\n");
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         fprintf(fPtr, "%d %d %f %f %f %f %f %f %d %d %d\n", (moli-mol)*moli->nsites+i, moli->type[i], 
              moli->p[i].x*unit.LENGTH, moli->p[i].y*unit.LENGTH, moli->p[i].z*unit.LENGTH, 0, 0, 0, 0, 0, 0);
}

