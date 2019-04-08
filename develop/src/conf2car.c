/*
	program:	conf2car.c
	author:		Peng Yi at MIT
	date:		November 19, 2007
	purpose:	convert configuration files to .car files
	note:		require setup file
			shiftbox() added June 25, 2009
*/

#define __MAIN_PROGRAM
#include "header.h"

#define SIZECAP	1

void nucleus_carfile(long, beadstruct *, long, vector *);

char *Element(long t, char *s)
{
   double	mass = type[t].M;

   s[0]		=	0;

   if (fabs(mass-14) <= 1e-6)		strcpy(s, "C");
   else if (fabs(mass-15) <= 1e-6)	strcpy(s, "C");
   else if (fabs(mass-1.01) <= 1e-6)	strcpy(s, "H");
   else if (fabs(mass-28.086) <= 1e-6)	strcpy(s, "Si");
   else if (fabs(mass-26.982) <= 1e-6)	strcpy(s, "Al");
   else if (fabs(mass-16) <= 1e-6)	strcpy(s, "O");

   strcpy(s, "C");
   return	s;
}


void doublesize()		// double the system size
{ 		
   long		i;
   molstruct	*moli;
			
   printf("%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      printf("%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      printf("%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
         
      for (i=0; i<moli->nsites; i++) 
         printf("%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      printf("%d\t%d\t%d\n", moli-mol+NMOLS, moli->box, moli->nsites);

      for (i=0; i<moli->nsites; i++) 
         printf("%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z + BOX[0].lz);
   }
}


void types()		// temporary code to prepare a system with 4 types
{
   long		i;
   molstruct	*moli;

   printf("%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      printf("%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

   for (moli=0; moli<mol+NMOLS; moli++) {
      printf("%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
         
      for (i=0; i<moli->nsites; i++) {
         if ( mod((moli-mol)/6, 2)) 
            printf("%d\t%f\t%f\t%f\n", moli->type[i]+2, moli->p[i].x, moli->p[i].y, moli->p[i].z);
         else
            printf("%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
      }
   }
}

//============================================================================//

void shift()				// move everything to the central box and 
{					// make the largest nucleus at the center  (5/20/08)
   FILE		*fPtr;
   vector	center, com, 
		rA,			// center of one chain in nucleus
		rO,			// center of nucleus 
		rBA, 			// rB-rA, B and A belong to same nucleus
		rOA;			// rO-rA
   molstruct	*moli;
   long		i, system=0, id, n;

   if (!(fPtr=fopen("shifted", "w"))) {
      printf("shifted file failed to open\n");
      exit(1);
   }

   // step1
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         if (sizeofnucl[moli->nuclid[i]]==MAXSIZE[0]){ 	// remember that sizeofnucl
            id	=	moli->nuclid[i];		// points to sizeofnuclp2 now
            rA	=	moli->p[i];
            break;
         }
      }
   }

   //step 2
   n=0;			// # of sites in this nucleus
   V_Null(&rBA);
   V_Null(&rOA);
   
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         if (moli->nuclid[i] == id) {
            n++;
            com	=	moli->p[i];
            rBA	=	V_Subtr(&com, &rA);
            rOA	=	V_Add(&rOA, &rBA);
         }}}
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);

   //step3
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->p[i]	=	V_Subtr(moli->p+i, &rO);
      }
   }

   // step 4
   for (moli=mol; moli<mol+NMOLS; moli++) {
      MolInBox2(moli);
   }

   // output shifted configuration file

   fprintf(fPtr, "TIMESTEP\t%d\n", TIMESTEP); 
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx * unit.LENGTH, BOX[i].ly * unit.LENGTH, BOX[i].lz * unit.LENGTH);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
   }
   fclose(fPtr);
}

//============================================================================//

void shift1()				// move everything to the central box and 
{					// make the largest nucleus at the center  (5/20/08)
   FILE		*fPtr;
   vector	center, com, 
		rA,			// center of one chain in nucleus
		rO,			// center of nucleus 
		rBA, 			// rB-rA, B and A belong to same nucleus
		rOA;			// rO-rA
   molstruct	*moli;
   long		i, system=0, id, n;

   double	lx[MAXNSYSTEMS], 
		ly[MAXNSYSTEMS], 
		lz[MAXNSYSTEMS];// shrink box to just fit the biggest nucleus in
   long		inthebox, nmol, nsites;
   char		dummy[80];
   
   if (!(fPtr=fopen("shifted", "w")))
      Exit("conf2car", "shift", "shifted failed to open");


   // Step 1: find one chain A belonging to the largest nucleus as a reference

   for (moli=mol; moli<mol+NMOLS; moli++)
      if (sizeofnucl[moli->nuclid[0]]==MAXSIZE[0]) {
	  id	=	moli->nuclid[0];		// nucleus id
	  rA	=	CenterofMass(moli); 		// take one chain as reference point
          break;
      }

   // Step 2: calc. the shift of center of nucleus to this chain

   n	=	0;					// # of chains in this nucleus
   V_Null(&rBA);
   V_Null(&rOA);

   for (moli=mol; moli<mol+NMOLS; moli++) 
      if (moli->nuclid[0] == id) {
	 n	++;
	 com	=	CenterofMass(moli);
	 rBA	=	V_Subtr(&com, &rA);
	 rBA	=	MapInBox2(&rBA, PBC, system);
	 rOA	=	V_Add(&rOA, &rBA);
      }

   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);			// center of nucleus

/*
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (sizeofnucl[moli->nuclid[0]]==MAXSIZE[0]) {
         id	=	moli->nuclid[0];
         rO	=	CenterofNucleus(moli->nuclid[0], moli);
	 break;
      }
*/
   // Step 3: every particle shifts rO so that the center of nucleus becomes center of box

   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         moli->p[i]	=	V_Subtr(moli->p+i, &rO);

   // Step 4: map all in central box

   for (moli=mol; moli<mol+NMOLS; moli++) 
      MolInBox2(moli);

// output the minimum box that contains the biggest nucleus 11/13/08

   for (n=0; n<NSYSTEMS; n++) {
      lx[n]	=	0;
      ly[n]	=	0;
      lz[n]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->nuclid[0]==id) {
	 for (i=0; i<moli->nsites; i++) {
            if (fabs(moli->p[i].x) > lx[moli->box])
	       lx[moli->box]	=	fabs(moli->p[i].x);
            if (fabs(moli->p[i].y) > ly[moli->box])
	       ly[moli->box]	=	fabs(moli->p[i].y);
            if (fabs(moli->p[i].z) > lz[moli->box])
	       lz[moli->box]	=	fabs(moli->p[i].z);
         }
      }
   }

   for (i=0; i<40; i++)
      fprintf(fPtr, " ");
   fprintf(fPtr, "\n");
   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", lx[i] *2 * unit.LENGTH, ly[i]*2 * unit.LENGTH, lz[i]*2 * unit.LENGTH);

   n=0;   nsites=0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      inthebox	=	1;
/*
      com	=	CenterofMass(moli);
      if (fabs(com.x) > lx[0] || fabs(com.y) > ly[0] || fabs(com.z) > lz[0]) 
	 inthebox	=	0;
*/
      for (i=0; i<moli->nsites; i++) {
         if (fabs(moli->p[i].x) > lx[0] || fabs(moli->p[i].y) > ly[0] || fabs(moli->p[i].z) > lz[0]) {
            inthebox	=	0;
            break;
	 }
      }

      if (inthebox) {
         fprintf(fPtr, "%d\t%d\t%d\n", n, moli->box, moli->nsites);
         for (i=0; i<moli->nsites; i++) 
            fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
         n	++;
         nsites	+=	moli->nsites;
      }
   }
   rewind(fPtr);
   sprintf(dummy, "%d\t%d\t%d", NSYSTEMS, n, nsites);
   fprintf(fPtr, "%s", dummy);
   for (i=0; i<40-strlen(dummy); i++)
      fprintf(fPtr, " ");

   fseek(fPtr, 0, SEEK_END);			// go to the end of the file

   /* output shifted conf file */
 
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx * unit.LENGTH, BOX[i].ly * unit.LENGTH, BOX[i].lz * unit.LENGTH);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
   }
   fclose(fPtr);
}

//======================//
//	Main Program	//
//======================//

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   FILE		*fPtr, *fcar, *fpdb, *fxyz;
   molstruct	*moli;
   char		atomname;
   char		s[255], ff[255], filename1[255], filename2[255];
   long		i, j, k, n, system, flag;
   long		nsite, nbond, nangle, ndihedral;
   vector	com;		// center of mass

long	drawmol;
matrix	Mj, Mk;
vector eigj, eigk, dp, nucleusshape;
double	Rg2;

   beadstruct	*nucleus;		// for grouping beads in biggest nucleus
   vector	rbead[MAXNMOLS*MAXNMOLSITES];
   molstruct 	*molid[MAXNMOLS];	// mols that in the largest nuclei
   long		nmol, nuclid;

   tim=(int *)malloc(sizeof(int));     	//random number generator
   seed=(long *)malloc(sizeof(long));
   *tim=(int)time(NULL);
   *seed= -1*(*tim);           		//seed must start out as a negative long

   //------MPI initialization------//
   MPI_numprocs	=	1;			// default value
   MPI_myid	=	0;			// default value
#ifdef myMPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &MPI_numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myid);
#endif

   if (argc<2) {
      printf("conf2car (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tconf2car filename\n\n");
      printf("Notes:\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   strcpy(filename1, argv[1]);
   strcat(filename1, ".car");
   if (! (fPtr = fopen(filename1, "w")))      exit(1);

   strcpy(filename2, argv[1]);
   strcat(filename2, ".pdb");
   if (! (fpdb = fopen(filename2, "w")))      exit(1);
   
   strcpy(filename2, argv[1]);
   strcat(filename2, ".xyz");
   if (! (fxyz = fopen(filename2, "w")))      exit(1);

   freopen("conf2car.out", "w", stdout);	// redirect standard output stream to a file

   InitMols(MAXNMOLS, MAXNMOLS);
   GetSetup(argv);

   InitUnits();
   GetCoordinates(argv[1]);		// include unit conversion
   //Read_Conf(argv[1]);		// read in configuration file
   //CoorSI2System();			// convert coordinates from SI to system units
   InitForcefield();

#ifdef CELL_LIST
   CL_Init();				// Cell list should be ready before Verlet list
   CL_Build();				// needed for Calc_Qlm
#endif	/* CELL_LIST */

   for (moli=mol; moli<mol+NMOLS; moli++) { 		
      for (i=0; i<moli->nsites; i++)  {
         moli->flags[i]		=	1;		// activate all the sites on this processor
         moli->parent[i]	=	i-1;		// initialize parent site
      }
      moli->flip		=	0;		// flip to the original direction
      moli->origin		=	CenterofMass(moli);
   }

   nucleus	=	(beadstruct *) calloc(NSITES, sizeof(beadstruct));
   if (nucleus == NULL) {
      printf("nucleus allocation failed!\n");
      exit(-1);
   }

   //---------------------------//
   //	Perform Analysis	//
   //---------------------------//
   CalcV();

   InitSample();
   SampleP2All();

   SampleSpherical();
   Dist_Spherical();
   //Calc_Qlm(l_of_Ylm);

   for (i=0; i<4; i++) {	// smooth crystal beads until convergence
      xtal_smooth();
   }
   Find_Nuclei_p2(1);

   S_PrintAll();

   for (i=0; i<argc; i++) {
      printf("%s ", argv[i]);
   }
   printf("\n");

   printf("# System Info.\n");
   for (system=0; system <NSYSTEMS; system++) {
      printf("System %d\n", system);
      printf("Vtot/NSITES=%f\n", v[system].tot/NSites[system]);
      printf("Rp = %f\n", Rp);
      printf("Rconn = %f\n", Rconn);
      printf("P2 = %f\n", P2[system]);
      printf("P2m = %f\n", P2M[system]);
      printf("P2z = %f\n", P2z[system]);
      printf("nmax = %d\n", nmax[system][0]);
      printf("2ndnmax = %d\n", nmax[system][1]);
      printf("3rdnmax = %d\n", nmax[system][2]);
      printf("realXtal = %d\n", realXtal[system]);
      printf("Q6 = %f\n", Q6[system]);
      printf("Q4 = %f\n", Q4[system]);
      printf("\n");
   }

   //------Grouping the beads belonging to the Biggest nucleus------//

   system	=	0;		// for now only one system, 4/16/2010
   i	= 1;
   while (sizeofnucl[i] != nmax[system][0]) {
      i ++;
   }
   nuclid	=	i;

   n	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system == moli->box) {
         for (i=0; i<moli->nsites; i++) {
            if (moli->nuclid[i] == nuclid) {
               nucleus[n].moli	=	moli;
	       nucleus[n].site	=	i;
	       n	++;
	    }
         }
      }
   }
   MapInNucleus(system, nucleus, n, rbead);	// move nucleus to the center box
   Find_segments();				// Find all segments
   Seg_Type(nuclid);				// identify the type of segments associated with nuclid
   nucleus_carfile(system, nucleus, n, rbead);	// visualize nucleus

   /* Find out how many chains participate in the biggest nucleus and store them in an array */

   if (n>0) {
      molid[0]	=	nucleus[0].moli;
      nucleus[0].moli->p[nucleus[0].site]	=	rbead[0];
      unfoldchain(nucleus[0].moli, nucleus[0].site);
      nmol	=	1;
      
      for (i=1; i<n; i++) {			// search through nucleus
	 flag	=	0;
	 for (j=0; j<nmol; j++) {
	    if (molid[j]==nucleus[i].moli) {
	       flag	=	1;
	       break;
	    }
	 }
	 if (flag==0) {
	    molid[nmol]		=	nucleus[i].moli;
	    nucleus[i].moli->p[nucleus[i].site]	=	rbead[i];
	    unfoldchain(nucleus[i].moli, nucleus[i].site);
	    nmol	++;
	 }
      }
   }

   //---------------------------------------------------//
   //	Write chains in the biggest nucleus to xyz file	//
   //---------------------------------------------------//
   fprintf(fxyz, "%ld\n", nmol*NSITES/NMOLS);
   fprintf(fxyz, "comments\n");
   for (i=0; i<nmol; i++) {
      k = molid[i]-mol;

      for (j=0; j<molid[i]->nsites; j++) {
	 if (molid[i]->nuclid[j] == nuclid || k==2 || k==5 || k==38 || k==53 || k==58 || k==1 || 
	    k==4 || k==24 || k==28 || k==30 || k==41 || k==46 || k==54 || k==55 || k==56) {
	    fprintf(fxyz, "O  %lf %lf %lf\n", 
	        molid[i]->p[j].x * unit.LENGTH, molid[i]->p[j].y * unit.LENGTH, molid[i]->p[j].z * unit.LENGTH);
	 }
	 else {
	    //fprintf(fxyz, "N  %lf %lf %lf\n", 
	    //    molid[i]->p[j].x * unit.LENGTH, molid[i]->p[j].y * unit.LENGTH, molid[i]->p[j].z * unit.LENGTH);
	 }
      }
   }
   fclose(fxyz);

   //***** create data file for Greg for POVRAY visualization *****//
   /* 
   fprintf(fxyz, "0.0, 0.0, 0.0,\n");		// center of display
   fprintf(fxyz, "%lf, %lf, %lf,\n", BOX[0].lx * unit.LENGTH, BOX[0].ly * unit.LENGTH, BOX[0].lz * unit.LENGTH);
   fprintf(fxyz, "%ld,\n", nmol);
   fprintf(fxyz, "%ld,\n", nmol*NSITES/NMOLS);
   for (i=0; i<nmol; i++) {
      k = molid[i]-mol;

      for (j=0; j<molid[i]->nsites; j++) {
	 fprintf(fxyz, "%lf, %lf, %lf, %ld, ", 
	    molid[i]->p[j].x * unit.LENGTH, molid[i]->p[j].y * unit.LENGTH, molid[i]->p[j].z * unit.LENGTH, k);

	 if (molid[i]->nuclid[j]==nuclid) {
	    fprintf(fxyz, "2,\n");
	 }
	 else if (k==2 || k==5 || k==38 || k==53 || k==58) {
	    fprintf(fxyz, "1,\n");
	 }
	 else {
	    fprintf(fxyz, "0,\n");
	 }
      }
   }
   fclose(fxyz);
   */
   
   // doublesize();
   // types();

   system	=	0;			// only one system for now

   //-------------------------------------------//
   //	Write the whole system to .pdb file	//
   //-------------------------------------------//
   fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
   fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, BOX[system].lz * unit.LENGTH,
	90.0, 90.0, 90.0);

   n		=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {
         //MolInBox2(moli);

         for (i=0; i<moli->nsites; i++) {
            if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])	// nuclid index starts from 1
               atomname	=	'N';				// blue in VMD
	    else if (sizeofnucl[moli->nuclid[i]] >= SIZECAP)
               atomname	=	'O';				// red in VMD
            else
               atomname	=	'C';				// cyan in VMD

	    n	++;

	    fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
            fprintf(fpdb, "%5d ", n);		// atom number
	    fprintf(fpdb, " %c  ", atomname);	// atom name
            fprintf(fpdb, " ");			// alternate location indiator
	    fprintf(fpdb, "   ");		// residue name
	    fprintf(fpdb, " ");			// column 21
            fprintf(fpdb, " ");			// chain identifier, column 22
	    fprintf(fpdb, "    ");		// residue sequence number, 23-26
	    fprintf(fpdb, " ");			// code for insertion of residues, 27
	    fprintf(fpdb, "   ");		// column 28-30
            fprintf(fpdb, "%8.3f%8.3f%8.3f", 	// position
		moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
		moli->p[i].z * unit.LENGTH);
	    fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
	    fprintf(fpdb, "\n"); 
         } 
      }   
   }
   fprintf(fpdb, "END\n");

   //-------------------------------------------//
   //	Write the whole system to .car file	//
   //-------------------------------------------//
   fprintf(fPtr, "!BIOSYM archive 3\n");
   if (1==PBC) {
      fprintf(fPtr, "PBC=ON\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
      fprintf(fPtr, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, 
	BOX[system].lz * unit.LENGTH, 90.0, 90.0, 90.0);
   }
   else if (0==PBC) {
      fprintf(fPtr, "PBC=OFF\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
   }

   n		=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {

         //MolInBox2(moli);
         for (i=0; i<moli->nsites; i++) {
            if (moli->nuclid[i] > 0) {

               if (sizeofnucl[moli->nuclid[i]] == nmax[system][0])
                  sprintf(s, "N%d", n++);
	       else
		  sprintf(s, "C%d", n++);
            }
            else {
	       switch(moli->segtype[i]) {
		  case '2':	sprintf(s, "H%d", n++);	break;	// tail
		  case '3':	sprintf(s, "S%d", n++);	break;	// loop
		  case '4':	sprintf(s, "P%d", n++); break;	// brdg
		  case '5':	sprintf(s, "Z%d", n++); break;	// pbcbrdg
		  case '7':	sprintf(s, "O%d", n++); break;	// dfct
		  default:	break;
	       }
	    }
            //if (sizeofnucl[moli->nuclid[i]] > SIZECAP) {
            fprintf(fPtr, "%-6.6s ", s);
            sprintf(s, "M%d", moli-mol);
            fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
			moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
			moli->p[i].z * unit.LENGTH);
            strcpy(ff, "O");
            fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[i], s));
            //}
         } 
      }   
   }
   fprintf(fPtr, "end\nend\n");

   //-------------------------------------------------------------------//
   //	Write chains containing the biggest nucleus to .car file	//
   //-------------------------------------------------------------------//
   fprintf(fPtr, "**********The Chains Containing Biggest Nucleus Only**********\n");
   fprintf(fPtr, "!BIOSYM archive 3\n");
   if (1==PBC) {
      fprintf(fPtr, "PBC=ON\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
      fprintf(fPtr, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, 
	BOX[system].lz * unit.LENGTH, 90.0, 90.0, 90.0);
   }
   else if (0==PBC) {
      fprintf(fPtr, "PBC=OFF\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
   }
   
   n	=	0;

   for (moli=mol; moli<mol+NMOLS; moli++) {		// for nucleus def based on p2
      if (system==moli->box) {
         //MolInBox2(moli);

         drawmol	=	0;
         for (i=0; i<moli->nsites; i++) {
            if (sizeofnucl[moli->nuclid[i]] == nmax[system][0]) {
            //  if (sizeofnucl[moli->nuclid[i]] == nmax[system][0] 
            //  || sizeofnucl[moli->nuclid[i]] == nmax[system][1]) {
               drawmol	=	1;
               break;
            }
         }
         if (drawmol) {
	    for (i=0; i<moli->nsites; i++) {
	       k	=	(moli-mol)*NSITES/NMOLS+i;
	       
               if (sizeofnucl[moli->nuclid[i]] == nmax[system][0])
                  sprintf(s, "N%d", k);
               // else if (sizeofnucl[moli->nuclid[i]] == nmax[system][1])
               //    sprintf(s, "O%d", n++);
               else
                  sprintf(s, "C%d", k);

               fprintf(fPtr, "%-6.6s ", s);
               sprintf(s, "M%d", moli-mol);

               //moli->p[i]	=	MapInBox2(moli->p+i, PBC, system);

               fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
			moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
			moli->p[i].z * unit.LENGTH);
               strcpy(ff, "O");
               fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, 
			Element(moli->type[i], s));
            }
         }
      }
   }
   fprintf(fPtr, "end\nend\n");

   //-----------------------------------------------------------//
   //	Write the beads in the biggest nucleus to .car file	//
   //-----------------------------------------------------------//
   fprintf(fPtr, "**********The Biggest Nucleus Only**********\n");
   fprintf(fPtr, "!BIOSYM archive 3\n");
   if (1==PBC) {
      fprintf(fPtr, "PBC=ON\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
      fprintf(fPtr, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, 
	BOX[system].lz * unit.LENGTH, 90.0, 90.0, 90.0);
   }
   else if (0==PBC) {
      fprintf(fPtr, "PBC=OFF\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
   }
   n	=	0;

   for (moli=mol; moli<mol+NMOLS; moli++) {		// visualize biggest nucleus
      if (system==moli->box) {
         MolInBox2(moli);

	 for (i=0; i<moli->nsites; i++) {
	    k	=	(moli-mol)*NSITES/NMOLS+i;
	    
            if (sizeofnucl[moli->nuclid[i]] == nmax[system][0]) {
               sprintf(s, "N%d", k);

               fprintf(fPtr, "%-6.6s ", s);
               sprintf(s, "M%d", moli-mol);
               fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
			moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
			moli->p[i].z * unit.LENGTH);
               strcpy(ff, "O");
               fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, 
			Element(moli->type[i], s));
	    }
         }
      }
   }
   /*
   for (moli=mol; moli<mol+NMOLS; moli++) {		// visualize 2nd biggest nucleus
      if (system==moli->box) {
         MolInBox2(moli);

	 for (i=0; i<moli->nsites; i++) {
	    k	=	(moli-mol)*NSITES/NMOLS+i;
	    
            if (sizeofnucl[moli->nuclid[i]] == nmax[system][1]) {
               sprintf(s, "O%d", k);

               fprintf(fPtr, "%-6.6s ", s);
               sprintf(s, "M%d", moli-mol);
               fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
			moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
			moli->p[i].z * unit.LENGTH);
               strcpy(ff, "O");
               fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, 
			Element(moli->type[i], s));
	    }
         }
      }
   }
   */
   fprintf(fPtr, "end\nend\n");


   ////////////////////////////////////////////////////////////
   /* print out center of mass of chains in .car file format */
   ////////////////////////////////////////////////////////////
   
   /*
   n		=	0;
   printf("**********Center of mass of chains**********\n");
   printf("!BIOSYM archive 3\n");
   if (1==PBC) {
      printf("PBC=ON\n\n");
      printf("!DATE %s", asctime(localtime(&t)));
      printf("PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, BOX[system].lz * unit.LENGTH, 90.0, 90.0, 90.0);
   }
   else {
      printf("PBC=OFF\n\n");
      printf("!DATE %s", asctime(localtime(&t)));
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {

         if (sizeofnucl[moli->nuclid[0]] > SIZECAP) {

         com	=	CenterofMass(moli);
         sprintf(s, "M%d", n++);
         printf("%-5.5s ", s);
         sprintf(s, "M%d", moli-mol);
         printf("%14.8g %14.8g %14.8g ", com.x * unit.LENGTH, com.y * unit.LENGTH, com.z * unit.LENGTH);
         strcpy(ff, "O");
         printf("%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[0], s));

         }
      }   
   }
   printf("end\nend\n");
   */

   ////////////////////////////////////////////////////////
   /* Shift the largest nucleus to the center of the box */
   ////////////////////////////////////////////////////////a

   //shift();

   //-------------------//
   //	Close files	//
   //-------------------//
   fclose(fPtr);
   fflush(stdout);
   fclose(stdout);
   free(nucleus);
#ifdef	myMPI
   MPI_Finalize();
#endif
   return	0;
}

//===========================================================================//

void shiftbox(long system, beadstruct *nucleus, long nsites)	// move the biggest nucleus
{								// to the center of the box
   molstruct	*moli;
   long		i, n, site;
   vector	rA, rO, rBA, rOA;

   // step 1: find one segment that belongs to the biggest nucleus
   
   rA	=	nucleus[0].moli->p[nucleus[0].site];
  
   // step 2: calc. the shift of com of nucleus to this segment
   
   V_Null(&rBA);
   V_Null(&rOA);

   for (n=0; n<nsites; n++) { 
      moli	=	nucleus[n].moli;
      site	=	nucleus[n].site;
      rBA	=	V_Subtr(moli->p+site, &rA);
      rBA	=	MapInBox2(&rBA, PBC, system);
      rOA	=	V_Add(&rOA, &rBA);
   }
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);

   // step 3: every segment shift and move to central box
  
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (site=0; site<moli->nsites; site++) {
         moli->p[site]	=	V_Subtr(moli->p+site, &rO);
         moli->p[site]	=	MapInBox2(moli->p+site, PBC, system);
      } 
   }
}

//======================================//
//	Visualize one nucleus (4/2012)	//	
//======================================//
void nucleus_carfile(long system, beadstruct *nucleus, long nsites, vector *rbead)
{
   FILE		*fPtr;
   time_t	t	=	time(NULL);
   char		s[255], ff[255];
   long		i;

   fPtr	=	fopen("nucleus.car", "w");

   fprintf(fPtr, "!BIOSYM archive 3\n");
   if (1==PBC) {
      fprintf(fPtr, "PBC=ON\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
      fprintf(fPtr, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, 
	BOX[system].lz * unit.LENGTH, 90.0, 90.0, 90.0);
   }
   else if (0==PBC) {
      fprintf(fPtr, "PBC=OFF\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
   }

   for (i=0; i<nsites; i++) {
      sprintf(s, "N%d", i+1);		// N: blue color in VMD (O: red)

      fprintf(fPtr, "%-5.5s ", s);
      sprintf(s, "M%d", nucleus[i].moli-mol);
     
      fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
		rbead[i].x * unit.LENGTH, rbead[i].y * unit.LENGTH, rbead[i].z * unit.LENGTH);
      strcpy(ff, "O");
      fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", 
		ff, nucleus[i].moli-mol, Element(nucleus[i].moli->type[nucleus[i].site], s));
   }
   fprintf(fPtr, "end\nend\n");
   fclose(fPtr);
   return;
}
