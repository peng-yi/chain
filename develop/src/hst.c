/*
    program:    hst.c
    author:     Peng Yi at MIT
    date:       January 12, 2008
    purpose:    get info from binary history files and do analysis
    note:	April 15, 2010 add biggest nucleus shape analysis
*/

#define __MAIN_PROGRAM
#include "header.h"

#include "correlation.h"

long	nrecords,			// # of records
	min=MAXNMOLS, 
	max=0,
	nnmax[MAXNMOLS],		// nnmax[i] is the counts of conf. with nmax=i
	n[1000][MAXNMOLS];		// n[i][j] is the counts of nuclei of size j when nmax=i
double	pn[1000][MAXNMOLS],		// pn[i][j] is the prob. of nuclei of size j when nmax=i
	Gn[MAXNMOLS],			// free energy of one nucleus
	G[MAXNMOLS];			// free energy of the system
long	lencnt[MAXNSYSTEMS][MAXNMOLSITES+1];	// number of chains of different length
long	BIN;				// binsize for nucleus size statistics

// some functions declared here
void printout(FILE *);			// print out probability distribution
void hst2conf();			// convert to configuration file
void hst2lammpsdump(FILE *fPtr);	// convert to dump file in lammps format

int main(int argc, char * argv[])
{
   FILE		*fhist, *fp, *fp1, *fsize, *fconf, *flammps;
   char		file_hst1[256],
		comments[1024];
   long		addflag=0, splitflag=0, normalflag=0, 
		Nsplit, 
		Nblock = -1;		// # of records in each block, could be set by command line
   long		version, system;
   long		i, j, k, m, nuclid; 
   double	dummy;
   double	var, bias; 		// bias variable
   static long	init = 1;
   molstruct	*moli;
   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];	// contain the beads in the biggest nucleus

   if (argc<3) {
      printf("hst (c) 2010 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\thst -normal history_file BINSIZE [N]\n\n");
      printf("\thst -add history_file1 history_file2\n\n");
      printf("\thst -split history_file N\n\n");
      printf("Notes:\n");
      printf("\t* require setup file to get information e.g. Temp for analysis\n\n");
      printf("\t* -normal: normal history file analysis\n");
      printf("\t\t BINSIZE is the binsize for nucleus increment\n\n");
      printf("\t\t if N is given, then do analysis in blocks each with N records, otherwise as one big block\n\n");
      printf("\t* -add: add history_file2 to the end of history_file1\n\n");
      printf("\t* -split: split history_file to binsplit1 and binsplit2 and the former contains N records\n\n");
      exit(1);
   }

printf("start\n");
   if (!strcmp(argv[1], "-add")) {
      addflag	=	1;
      sprintf(file_hst1, argv[2]);
      sprintf(file_hst, argv[3]);
   }
   else if (!strcmp(argv[1], "-split")) {
      splitflag	=	1;
      sprintf(file_hst, argv[2]);
      Nsplit		=	atol(argv[3]);
   }
   else if (!strcmp(argv[1], "-normal")) {
      normalflag	=	1;
      sprintf(file_hst, argv[2]);
      BIN	=	atol(argv[3]);
      if (argc == 5)					// optional
	 Nblock		=	atol(argv[4]);
   }

   InitMols(MAXNMOLS, MAXNMOLS);// allocate max. memory for molecules
   GetSetup(argv);		// most of the variables read from setup file
				// will be updated by the reading from history file
   InitUnits();			// calculate the transformation b/w SI units and system units	
   InitForcefield();		// initialize forcefield
   InitSample(argv);		// initialize sampling

   if (normalflag) {		// open output file for distribution output
      if (!(fsize = fopen("Psize","w"))) 
         Exit("hst", "main","Size prob. distribution file failed to open");
   }
   else if (addflag) {
      if (!(fp1=fopen(file_hst1, "a")))
         Exit("hst", "main", "add target history file fails to open");
   }

/*
   if (!(fdump=fopen("dump", "w"))) 		// to convert binary hist file to dump file through Write_Conf
      Exit("hst", "main", "dump file failed to open");
*/

   if (!(flammps=fopen("lammps.dump", "w")))	// to convert to lammps dump file for lammps2pdb
      Exit("hst", "main", "lammps dump file failed to open");

   if (!(fhist=fopen(file_hst, "r")))
      Exit("hst", "main", "file_hst failed to open.");

   Read_MultiConf(fhist);			// read the first conf., no analysis for this one

   nrecords	=	0;

   while (Read_MultiConf(fhist)) {		// read next conf., start analysis

      printf("read in %d configuration...\n", nrecords);

      /////////////////////////////////////////
      /* PRE-PROCESSING of data for analysis */
      /////////////////////////////////////////
 
      CoorSI2System();				// convert from SI to system units

      for (i=0; i<NSYSTEMS; i++) {
         NMols[i]	=	0;
         NSites[i]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++)
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;			// total # of mols in certain system
            NSites[system]	+=	moli->nsites;	// total # of sites in certain system
            lencnt[system][moli->nsites]	++;	// polydispersity
         }
      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox	=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol	=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].pres 	=	P;
         BOX[i].temp 	=	kT;
         BOX[i].drmax 	=	DRMAX;
         BOX[i].dlmax 	=	DLMAX;
         BOX[i].damax	=	DAMAX;
         BOX[i].rc	=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb	=	Rb;
         BOX[i].rv	=	Rv;
      } 
      for (i=0; i<NMOLS; i++) { 			// note NMOLS might change due to MPI
         for (j=0; j<mol[i].nsites; j++)  {
            mol[i].flags[j]	=	1;		// activate all the sites on this processor
            mol[i].parent[j]	=	j-1;		// initialize parent sites
         }
         mol[i].flip		=	0;		// flip to the original direction
         mol[i].origin		=	CenterofMass(mol+i);
      }
#ifdef CELL_LIST				// cell list needed for energy calc.
      if (init) {
         CL_Init();				// Cell list should be ready before Verlet list
         init	=	0;
      }
      else
         CL_Destroy();
      CL_Build();
#endif	/* CELL_LIST */

      ///////////////////////
      /*  Perform Analysis */
      ///////////////////////

      printf("start performing analysis...\n");

      CalcV();					// need cell list
      if (V_VIRIAL)	Sample_Pressure();
      SampleP2All();				// calculate global P2 and local p2
      Dist_p2();				// collect distribution of local p2
      SampleSpherical();			// sample spherical coordinates
      Dist_Spherical();				// collect distribution

      //////////////////////////////////////////////////
      /* Calculate Bias for each record/configuration */
      //////////////////////////////////////////////////

      printf("calculate bias ...\n");

      Find_Nuclei_p2(1);

//      Find_Nuclei(1);				// sample nuclei for bias calculation
//printf("%d %d %d ", nrecords, nmax[0][0], nmax[0][1]);

      if (dynvar==1)
         var	=	MAXSIZE[0];
      else if (dynvar==3)
	 var	=	P2M[0];
      else if (dynvar==6)
	 var	=	MAXSIZE[0];
    
      bias	=	exp(0.5*kP2*(var-P2middle)*(var-P2middle));

      ////////////////////////////////////////////
      /* Find the nuclid of the biggest nucleus */
      ////////////////////////////////////////////

      sizeofnucl	=	sizeofnuclp2;	// Find_Nuclei_p2() uses sizeofnucleip2
      sizedist		=	sizedistp2;

      system	=	0;		// for now only one system, 4/16/2010
      i	= 1;				// nucleus id starts from 1
      while (sizeofnucl[i] != nmax[system][0]) {	// find the biggest nucleus
         i ++;
      }
      nuclid	=	i;		// store the id of the biggest nucleus
    
      ////////////////////////////////////////////////////////////////////////////
      /* Group the segments belong to the biggest nucleus based on sizeofnucl[] */
      ////////////////////////////////////////////////////////////////////////////

      m	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if (system == moli->box) {
            for (i=0; i<moli->nsites; i++) {
               if (moli->nuclid[i] == nuclid) {
                  nucleus[m].moli	=	moli;
	          nucleus[m].site	=	i;
                  m	++;
	       }
            }
         }
      }

      ////////////////////////////////
      /* Nucleus structure analysis */
      ////////////////////////////////

//      cylindershape(nucleus, nmax[system][0], nuclid); 	// calc. length and radius
							// assuming cylindrical shape

      //////////////////////////////////////////////////////
      /* Collect probability distribution of cluster size */
      //////////////////////////////////////////////////////

      printf("collect histogram ...\n");

      //Find_Nuclei_p2(1);			// Needed if the nucleus definition used
						// in biasing is different from the nucleus
						// definition used to calculate the nucleus 
						// size distribution

      //printf("%d %d\n ", nmax[0][0], nmax[0][1]);

      nnmax[MAXSIZE[0]]	++;
      for (i=1; i<=MAXSIZE[0]; i++) {		// index starts from 1
         if (BIN==1) {
	   pn[0][i] += bias * sizedistp2[i];	// probability distribution of cluster size
	   n[0][i]  += sizedistp2[i];
	 }
         else {
           pn[0][i/BIN+1] += bias * sizedistp2[i];
	   n[0][i/BIN+1]  += sizedistp2[i];
         }
         pn[MAXSIZE[0]][i] += bias * sizedistp2[i];
         n[MAXSIZE[0]][i]  += sizedistp2[i];
      }
      if (MAXSIZE[0] > max)
	 max	=	MAXSIZE[0];
      if (MAXSIZE[0] < min) 
         min	=	MAXSIZE[0];

      ////////////////////////////
      /* Calculate correlations */
      ////////////////////////////

      printf("calculate correlation ...\n");

      correlation();			// calculate correlation function
      if (!((nrecords*ITAPE)%IRADIAL))
         radial("sample");		// sample radial distribution function

      ///////////////////////////////////////
      /* Output to other trajactory format */
      ///////////////////////////////////////

//      hst2lammpsdump(flammps);	// convert to lammps dump format to use lammps2pdb 

//      Write_Conf(TIMESTEP);		// convert to ASCII dump file
//      if (normalflag && MAXSIZE[0]==17 && secondNmax[0]<=6) 
//      if (normalflag && MAXSIZE[0] >=15 && MAXSIZE[0] <=17 )
	// hst2conf(fconf);		// convert to configuration file    

      /////////////////////////////////////////////////////
      /* Summarize the results for this block and output */
      /////////////////////////////////////////////////////

      if (Nblock!=-1 && TIMESTEP%Nblock==0) {
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
      nrecords	++;
   }

   /////////////////////////////////////////////////
   /* Print out sampling results and correlations */
   /////////////////////////////////////////////////

   printf("print out sampling results and correlations...\n");

   S_PrintAll();
   corr_norm();
   corr_print();
   radial("print");

   ///////////////////////////////////////////////////////////////////////
   /* Summarize results for all records and print out size distribution */
   ///////////////////////////////////////////////////////////////////////
 
   if (Nblock==-1) {
      for (i=1; i<=max; i++) 
         Gn[i]	=	-log(pn[0][i]);			// free energy of nuclei
      for (i=min; i<=max; i++) {
         for (j=1; j<=i; j++) {
            G[i]	+=	n[i][j] * Gn[j];	// free energy of system
         }
      }
      printout(fsize);
   }
   if (normalflag) {
      fclose(fsize);
   }
   fclose(fhist);
   return EXIT_SUCCESS;
}


void printout(FILE *fPtr)		// output main information and distribution
{
   long		i, j;

   fprintf(fPtr, "ITAPE = %d (dump file written every ITAPE MC cycles)\n", ITAPE); 
   fprintf(fPtr, "Number of records: %d\n\n", nrecords);
/*
   fprintf(fPtr, "bias nucleus definition = ");
   switch (clusdef) {
      case	1:	fprintf(fPtr, "%s\n", );	break;
      case	1:	fprintf(fPtr, "%s\n", );	break;
      case	1:	fprintf(fPtr, "%s\n", );	break;
      default:		break;
   }
   fprintf(fPtr, "kP2 = %f\t P2middle = %f\n",kP2, P2middle);
   fprintf(fPtr, "distribution nucleus definition = ");
   switch (dynvar) {
      case	1:	fprintf(fPtr, "%s\n", );	break;
      case	1:	fprintf(fPtr, "%s\n", );	break;
      case	1:	fprintf(fPtr, "%s\n", );	break;
      default:		break;
   }
*/  

   /* Output Polydispersity */

   fprintf(fPtr, "\n#Polydispersity\n\n");
   for (i=0; i<NSYSTEMS; i++) {
      fprintf(fPtr, "length      number of chains (SYSTEM #%d)\n\n", i);
         for (j=0; j<=MAXNMOLSITES; j++)
            if (lencnt[i][j]!=0)
               fprintf(fPtr, "%d    %d\n", j, lencnt[i][j]/nrecords);
   }

   /* Output nuclei size distribution */

   fprintf(fPtr, "\n#Overall nuclei size distribution P(n):\n\n");
   fprintf(fPtr, "n    N(n)   fixed p(n)   -log(fixed pn[i])\n\n");
   for (i=1; i<=max/BIN; i++) 
      fprintf(fPtr, "%4d   %8d   %14.4f   %f\n", i*BIN, n[0][i], pn[0][i], Gn[i]-Gn[1]);

   fprintf(fPtr, "\n\n#Throw away poorly sampled nucleus size\n\n");
   fprintf(fPtr, "n    N(n)   fixed p(n)   -log(fixed pn[i])\n\n");
   for (i=1; i<=max/BIN; i++) {
      if (n[0][i] <10) {			// poorly sampled nucleus size
         fprintf(fPtr, "%4d   %8d   nan      nan\n", i*BIN, n[0][i]);
      }
      else {
         fprintf(fPtr, "%4d   %8d   %14.4f   %f\n", i*BIN, n[0][i], pn[0][i], Gn[i]-Gn[1]);
      }
   } 

   fprintf(fPtr, "\n#Nucleus Free energy G(n) for copy purpose...\n");
   for (i=1; i<=max/BIN; i++) {
      if (n[0][i]<10)
         fprintf(fPtr, "nan\n");		// poorly sampled nucleus size
      else
         fprintf(fPtr, "%f\n", Gn[i]-Gn[1]);
   }

   fprintf(fPtr, "\n#Distribution of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, "n_max    N(n_max)\n");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%4d  %8d\n", i, nnmax[i]);

   fprintf(fPtr, "\n#Size distribution N(n) for different value of n_max in the range of [%-3d, %3d]\n\n", min, max);
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
   fprintf(fPtr, "\n#Normalized to N(n_max)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", (nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n#-Logrithm of N(n)\n\n");
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
   fprintf(fPtr, "!TIMESTEP %d Nmax = %d 2ndNmax = %d\n", TIMESTEP, MAXSIZE[0], secondNmax[0]);
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
   long		i, accum;
   molstruct	*moli;

   fprintf(fPtr, "ITEM: TIMESTEP\n");	
   fprintf(fPtr, "%d\n", TIMESTEP);		
   fprintf(fPtr, "ITEM: NUMBER OF ATOMS\n");
   fprintf(fPtr, "%d\n", NSITES);
   fprintf(fPtr, "ITEM: BOX BOUNDS\n");
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lx*unit.LENGTH, 0.5*BOX[0].lx*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].ly*unit.LENGTH, 0.5*BOX[0].ly*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lz*unit.LENGTH, 0.5*BOX[0].lz*unit.LENGTH); 
   fprintf(fPtr, "ITEM: ATOMS\n");
   
   accum	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         fprintf(fPtr, "%d %d %d ", accum+i+1, moli-mol+1, moli->type[i]); 
         fprintf(fPtr, "%f %f %f ", moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
         fprintf(fPtr, "%f %f %f %d %d %d\n", 0.0, 0.0, 0.0, 0, 0, 0);
      }
      accum	+=	moli->nsites;
   }
   return;
}

 
