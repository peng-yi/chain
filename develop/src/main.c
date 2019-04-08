/*
    program:    main.c
    author:     Peng Yi at MIT
    date:       October 23, 2006
    purpose:    Single processor NPT ensemble simulation program
    note:	Nov. 5, 2009, prepare to output more variables, and use another
			program to view a selected list, like g_energy in 
			Gromacs does.
*/

#define __MAIN_PROGRAM
#include "header.h"

void Print_hst(FILE *fPtr)
{
   long		i;
   static long	init_hst=1;

   if (init_hst && !FLAG_restart) {			// Print out item name line
      init_hst	=	0;
      fprintf(fPtr, "counter  ");
      for (i=0; i<NBOX; i++) {
         fprintf(fPtr, "vtot[%ld]  ", i);		// 1. total potential energy
         fprintf(fPtr, "vol[%ld]  ", i); 		// 6. volume
         fprintf(fPtr, "prestot[%ld]  ", i);		// 7. pressure
         fprintf(fPtr, "x[%ld]  y[%ld]  z[%ld]  ", i, i, i);
         fprintf(fPtr, "P2[%ld]  P2M[%ld]  P2z[%ld]  ", i, i, i);
         fprintf(fPtr, "transfrac[%ld]  ", i);
         fprintf(fPtr, "Xtal[%ld]  ", i);
         fprintf(fPtr, "realXtal[%ld]  ", i);
         fprintf(fPtr, "Nnucl[%ld] ", i);
         fprintf(fPtr, "Nmax[%ld]  ", i);
         fprintf(fPtr, "2ndNmax[%ld]  ", i);
         fprintf(fPtr, "3rdNmax[%ld]  ", i);
         fprintf(fPtr, "4thNmax[%ld]  ", i);
         fprintf(fPtr, "vbnd[%ld]  ", i);		// 2. bond energy
         fprintf(fPtr, "vang[%ld]  ", i);		// 3. angle energy
         fprintf(fPtr, "vtors[%ld]  ", i);		// 4. torsion energy
         fprintf(fPtr, "vnbnd[%ld]", i);		// 5. nonbonded energy
      }
      fprintf(fPtr, "\n");
   }

   fprintf(fPtr, "%-6ld  ", counter);				// print out value
   for (i=0; i<NBOX; i++) {
      fprintf(fPtr, "%8.3f  ", v[i].tot);	
      fprintf(fPtr, "%8.3f  ", BOX[i].vol);
      fprintf(fPtr, "%8.3f  ", BOX[i].pres);
      fprintf(fPtr, "%6.4f  %6.4f  %6.4f  ", BOX[i].lx, BOX[i].ly, BOX[i].lz);
      fprintf(fPtr, "%6.4f  %6.4f  %6.4f  ", P2[i], P2M[i], P2z[i]);
      fprintf(fPtr, "%6.4f  ", transfrac[i]);
      fprintf(fPtr, "%4ld  ", Xtal[i]);
      fprintf(fPtr, "%4ld  ", realXtal[i]);
      fprintf(fPtr, "%4ld  ", Nnucl[i]);
      fprintf(fPtr, "%4ld  ", nmax[i][0]);
      fprintf(fPtr, "%4ld  ", nmax[i][1]);
      fprintf(fPtr, "%4ld  ", nmax[i][2]);
      fprintf(fPtr, "%4ld  ", nmax[i][4]);
      fprintf(fPtr, "%8.3f  ", v[i].stretch);
      fprintf(fPtr, "%8.3f  ", v[i].bending);
      fprintf(fPtr, "%8.3f  ", v[i].torsion);
      fprintf(fPtr, "%8.3f", v[i].nonbonded);
   }
   fprintf(fPtr, "\n");
   fflush(fPtr);
}

//////////////////////////
//	Main Program	//
//////////////////////////
int main(int argc, char * argv[])
{
   char		par[80];
   molstruct	*moli;
   long		i, ii, k;
   time_t	start, end;
   
   if (argc<2) {
      printf("main (c) 2012 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tmain [-option]\n\n");
      printf("Notes:\n");
      printf("\t* -option = \n");
      printf("\t* -new: a new run\n");
      printf("\t* -restart: restart from an unfinished run\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   FLAG_restart	=	0;		// whether restart from a unfinished run
   for (i=1; i<argc; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-restart")) {
	 FLAG_restart	=	1;
      }
   }
   
   start	=	time(NULL);	// program start time
   
   tim=(int *)malloc(sizeof(int));     	// random number generator
   seed=(long *)malloc(sizeof(long));
   *tim=(int)time(NULL);
   *seed= -1*(*tim);           		// seed must start out as a negative long

   InitAll(argv);			// initialization, initial configuration

   if (!FLAG_restart) 
      TIME_restart	=	0;
   else {
      TIME_restart	=	TIMESTEP;	// TIMESTEP read from initconf
      fprintf(foutput, "\nRESTART an unfinished run, TIME_restart = %ld\n", TIME_restart);
   }

#ifdef MPI				
   // to be added ...
#else
   //StartNewHistoryFile(file_hst, FALSE);	// open binary history file
#endif

   /* do the first sampling */

   for (i=0; i<NMOLS; i++)
      mol[i].origin	=	CenterofMass(mol+i);
         
   SampleP2All();				// calculate global P2 and local p2
   SampleSpherical();
   //Find_Nuclei(dynvar);
   Dist_p2();					// because p2 is used in biasing
   Dist_Spherical();				// because Ptrans is used in biasing
   //S_PrintAll();
   //H_StoreCurrent(file_hst);			// write to binary history file
   Find_Nuclei_p2(1);

   Sample_All();				// initial sampling, energy, pressure, etc
   
   //if (!FLAG_restart) {				// if not a restart run
      //Write_Conf(0);				// dump first configuration
      dump_conf(0);				// dump first configuration
      Print_hst(fhst);				// output variables to histogram file
   //}

   //***************************//
   //	Perform MC cycles	//
   //***************************//

   for (ii=0; ii<2; ii++) {
      if (ii==0) {
	 NCYCLE		=	Nequil;			// Equilibration run
         fprintf(foutput, "\nStart equilibration ......\n");
      }
      else if (ii==1) {
	 NCYCLE		=	Nprod;			// length of production run
         fprintf(foutput, "\nStart production .......\n");
         InitSample();					// clear the distributions to zero
      }
      InitEnsembles();					// reset acceptance, etc

      for (counter=TIME_restart+1; counter<=NCYCLE; counter++) {	// perform NCYCLE MC cycles

         if (counter==10) {				// if initial configuration far from
	    CalcV();					// equil. and has very large energy,
	 }						// the energy update in each MC move
							// might be too noisy, so recalculate 
							// energy here to make sure energy 
							// calculation is correct
/*
         if (!mod(counter+1, NCYCLE/8)) {
            type[2].EPSILON	=	type[0].EPSILON * (2-1/8-((double) (counter+1))/NCYCLE);
            type[3].EPSILON	=	type[1].EPSILON * (2-1/8-((double) (counter+1))/NCYCLE);
            fprintf(foutput, "\ttype[2].EPSILON = %f\ttype[3].EPSILON = %f\n", type[2].EPSILON,
		type[3].EPSILON);
            InitForcefield();
         }
*/

         Cycle();					// perform MC cycle and sampling
         Print_hst(fhst);				// output variables realtime

   //SampleN2NSQ();				// sample distributions
   if (fabs(kP2) > ZERO) {
      if (dynvar==3)
         Dist_p2();				// because p2 is used in bias

      Dist_Spherical();				// because Ptrans is used in bias
      //Find_Nuclei();				// must after p2 calculation
   }

         //*****************************//
         //	Dumping configurations	//
         //*****************************//

         if (mod(counter, ITAPE)==0) {		// dump every ITAPE sequences
            //H_StoreCurrent(file_hst);			// dump to binary file
	    //Write_Conf(counter);
            dump_conf(counter);
	 }

         //*****************************//
         //	Creating restart file	//
         //*****************************//
         /* 
	 if (mod(counter+1, ICONF)==0 && ii==1)	{
	    Write_Conf(counter+1);
	 }
         */
      }	// cycles done

      if (counter!=0) {					// sanity check: energy and virial
         EnergyCheck();      
      }
      if (ii<=1) {
         //S_PrintAll();					// output all distributions
      }
      if (ii==1) {
         Sample_Done();
         Write_Conf(-1);				// output last conf.
         //Visualize(1);
      }
   }

   Print_Summary();
   
   end	=	time(NULL);				// program end time
   fprintf(foutput, "\nTotal running time was %f seconds.\n", difftime(end, start));

   CloseFile();

   return 0;
}
