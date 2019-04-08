/*
    program:    hstcomb.c
    author:     Peng Yi at MIT
    date:       May 1, 2008
    purpose:    combine two binary history files
*/

#define __MAIN_PROGRAM
#include "header.h"

#include "correlation.h"

long	timestep = 0,
	min=MAXNMOLS, max=0,
	nnmax[MAXNMOLS],		// nnmax[i] is the counts of conf. with nmax=i
	n[100][MAXNMOLS];		// n[i][j] is the counts of nuclei of size j when nmax=i
double	pn[100][MAXNMOLS],		// pn[i][j] is the prob. of nuclei of size j when nmax=i
	Gn[MAXNMOLS],			// free energy of one nucleus
	G[MAXNMOLS];			// free energy of the system

// some functions declared here
void printout(FILE *);
void hst2conf(FILE *fPtr);
void hst2dump(FILE *fPtr);

int main(int argc, char * argv[])
{
   FILE		*fp, *fp1, *fsize, *fconf, *fcomb;
   long		addflag=0, splitflag=0, normalflag=0, 
		Nsplit, 
		Nblock = -1;		// # of records in each block, could be set by command line
   long		version;
   char		file_hst1[256], file_hst2[256],
		comments[1024];
   char		data;
   double	dummy;
   long		i, j, k; 
   molstruct	*moli;
   double	var, bias; 

   if (argc<2) {
      printf("hstcombine (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\thstcombine history_file1 history_file2\n\n");
      printf("\t\tcombine two history files and generate a combined one binhstcomb\n\n");
      exit(1);
   }

   sprintf(file_hst1, argv[1]);
   sprintf(file_hst2, argv[2]);
   sprintf(file_hst, "binhstcomb");

/*
   InitMols(MAXNMOLS, MAXNMOLS);	// allocate max. memory for molecules
   GetSetup(argv);			// most of the variables read from setup 
					// will be updated by the reading from history file
   CalcUnits(ConvertUnits);		// so far the binary file is in system units (4/26/08)
*/

   if (!(fp = fread_hist(file_hst1, &version)))		// open the first history file
      Exit("hstcombine", "main", "first history file fails to open.");

   if (!(fcomb=fcreate_hist(file_hst)))			// open the combo file, must be after one history file
      Exit("hstcombine", "main", "history combo file fails to open.");

   data = fgetc(fp);
   while(!feof(fp)) {					// copy 1st history file to combo file
      fputc((int)data, fcomb);
      data = fgetc(fp);
   }
   fclose(fp);

   if (!(fp = fread_hist(file_hst2, &version)))		// open the second history file
      Exit("hstcombine", "main", "second history file fails to open.");
   
   data = fgetc(fp);
   while(!feof(fp)) {					// copy 2nd history file to combo file
      fputc((int)data, fcomb);
      data = fgetc(fp);
   }
   fclose(fp);
   fclose(fcomb);


/*
   if (!(fp = H_GetFirst(file_hst, &version, FALSE))) {
      printf("File format not recognized.\n\n");
   }
   InitSample(argv);		// InitSample will need some commands read in from the first record
   fclose(fp);			// Or we can move the pointer to the beginning of the history file

   if (!(fp = H_GetFirst(file_hst, &version, FALSE))) {	// reopen history file
      printf("File format not recognized.\n\n");
   }
   while (!H_GetNext(fp, version)) {			// read in next record/configuration

      if (dynvar==3)
         SampleP2All();					// calculate global P2 and local p2
      SampleSpherical();
      Find_Nuclei();

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

      correlation();			// calculate correlation function

      // if (normalflag) hst2conf(fconf);		// convert to configuration file
      // hst2dump();			// convert to lammps dump format to use lammps2pdb in the future

      timestep	++;

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

   // normalize correlation functions and print out

   corr_norm();
   corr_print();

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

   fclose(fp);

   if (normalflag) {
      fclose(fsize);
      fclose(fconf);
   }
*/
   return EXIT_SUCCESS;
}

/*
void printout(FILE *fPtr)
{
   long		i, j;

   fprintf(fPtr, "\nOverall nuclei size distribution\n\n");
   fprintf(fPtr, "n    N(n)   fixed p(n)   -log(fixed pn[i])\n\n");
   for (i=1; i<=max; i++) 
      fprintf(fPtr, "%4d   %8d   %14.4f   %f\n", i, n[0][i], pn[0][i], Gn[i]);

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


void hst2conf(FILE *fPtr)
{
   long		i;
   molstruct	*moli;

   fprintf(fPtr, "!TIMESTEP %d\n", timestep);
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


void hst2dump(FILE *fPtr)
{
   long		i;
   molstruct	*moli;

   printf("ITEM: TIMESTEP\n");	
   printf("%d\n", timestep);		
   printf("ITEM: NUMBER OF ATOMS\n");
   printf("%d\n", NSITES);
   printf("ITEM: BOX BOUNDS\n");
   printf("%f %f\n", -BOX[0].lx, BOX[0].lx); 
   printf("%f %f\n", -BOX[0].ly, BOX[0].ly); 
   printf("%f %f\n", -BOX[0].lz, BOX[0].lz); 
   printf("ITEM: ATOMS\n");
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         printf("%d %d %f %f %f %f %f %f %d %d %d\n", (moli-mol)*moli->nsites+i, moli->type[i], 
              moli->p[i].x, moli->p[i].y, moli->p[i].z, 0, 0, 0, 0, 0, 0);
}
*/
