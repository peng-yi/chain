/*
	program:	lmp_conf.c
	author:		Peng Yi at MIT
	date:		March 2, 2011
	purpose:	Read lammps dump file, do quick analysis related to chain conformation
			or box shape, etc.  It is meant to do quick analysis, so complicated 
			analysis will still be done using lammps2hst.
	note:		
			require setup file
*/

#define __MAIN_PROGRAM
#include "header.h"

#define DEBUG		0
#define NEWVERSION	1
#define Rcut2		9	// 3sigma squared

#include "correlation.h"

long		timestep;
double		rshift2;		// shift of the biggest nucleus
long		nsitesp2;	// # of sites with p2 greater than a threshold
long		nmolsp2;
double		Rg2, Ro2, DRg2, DRo2;		// radius of gyration, n2n distance 
long		nRg2, nRo2;
double		nchains;

void Print_hst(FILE *fPtr)
{
   long		i;
   static long	init = 1;

   if (init) {
      fprintf(fPtr, "timestep ");
      fprintf(fPtr, "volume ");
      fprintf(fPtr, "lx ");
      fprintf(fPtr, "ly ");
      fprintf(fPtr, "lz ");
      fprintf(fPtr, "Rg2 ");		// average radius of gyration squared
      fprintf(fPtr, "Ro2 ");		// average n2n distance squared
      fprintf(fPtr, "nchains");
      fprintf(fPtr, "\n");
      init	=	0;
   }
   fprintf(fPtr, "%-6d ", timestep);
   fprintf(fPtr, "%8.4f ", BOX[0].vol);
   fprintf(fPtr, "%8.4f ", BOX[0].lx);
   fprintf(fPtr, "%8.4f ", BOX[0].ly);
   fprintf(fPtr, "%8.4f ", BOX[0].lz);
   fprintf(fPtr, "%6.4f ", Rg2);
   fprintf(fPtr, "%6.4f ", Ro2);
   fprintf(fPtr, "%6.4f", nchains);
   fprintf(fPtr, "\n");
   fflush(fPtr);
}

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   long		i, j, k, system, m, n, nuclid;
   long		nsites, record;
   long		nx, ny, nz, id, siteid, molid, type;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo;			// box upper boundary
   molstruct	*moli, *molj;
   double	temp1, temp2, temp3;		// temperory variables
   double	R2;

   char		filein[255], s[80], ff[80], par[80], dummy[255];
   FILE		*fin, *fhst, *fout, *fconf, *fpdb, *fdat;
   long		LENGTH, accum, confflag=0, carflag=0, pdbflag=0, polydisperse=0, drawmol;
   char		atomname;
   static long	init=1;

   vector		con;	// center of nucleus
   static vector	rO;	// original position of nucleus

   vector	chainface;
   double	orient;
   long		orientdist[180];	// chain orientation distribution

   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];

   vector	com[MAXNMOLS], temp;		// average center of mass
   long		ncom[MAXNMOLS];

   double	imagex, imagey, imagez;		// variables for creating pbc images
   long		imagen;

   if (argc<2) {
      printf("lmp_conf (c) 2011 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tlmps2hst [-option] [x= 0.1 y= 0.1 z= 0.1 n= 1] lammpsdumpfile\n\n");
      printf("Notes:\n");
      printf("\t* -option = \n");
      printf("\t* -poly: polydisperse system, has molid in dump file\n");
      printf("\t* -conf: configuration file output\n");
      printf("\t* x= y= z=: duplicate the system and shift unit vector\n");
      printf("\t* n=: multiple of shift vector\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-poly"))
         polydisperse	=	1;
      else if (samestr(par, "-conf")) 
         confflag	=	1;
      else if (samestr(par, "x=")) 
         imagex	=	atof(argv[i+1]);
      else if (samestr(par, "y=")) 
         imagey	=	atof(argv[i+1]);
      else if (samestr(par, "z=")) 
         imagez	=	atof(argv[i+1]);
      else if (samestr(par, "n=")) 
         imagen	=	atol(argv[i+1]);
   }
   strcpy(filein, argv[argc-1]);		// input file

   fin		=	fopen(filein, "r");
   fhst		=	fopen("lmp_conf.hst", "w");
   if (confflag)
      fconf	=	fopen("lammps.conf", "w");

   freopen("lammps.out", "w", stdout);	// redirect standard output stream to a file

   if (DEBUG)	printf("Initialization: start ... \n");

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   InitUnits();				// initialize units for calculation
   InitForcefield();			// initialize Lennard-Jones potential mixing rule

   system	=	0;		// for now 2/14/08, only one system

   InitSample();			// initialize sampling

   for (i=0; i<180; i++)		// initialize chain rotation distribution
      orientdist[i]	=	0;

   while (!feof(fin)) {

      if (DEBUG)	printf("Read one configuration: start ...\n");

      /* Read in one configuration from dump file */

      if (!fgets(dummy, sizeof(dummy), fin))	break;			// end of file
      fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &xlo, &xhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &ylo, &yhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &zlo, &zhi);	fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);

      BOX[system].lx	=	xhi-xlo;
      BOX[system].ly	=	yhi-ylo;
      BOX[system].lz	=	zhi-zlo;

      LENGTH		=	NSITES/NMOLS;	// monodisperse system for now (4/26/2008)
      accum		=	0;

      for (i=0; i<nsites; i++) {
         fscanf(fin, "%ld", &id);
         fscanf(fin, "%ld", &molid);
         fscanf(fin, "%ld", &type);
         fscanf(fin, "%lf%lf%lf %lf%lf%lf %ld%ld%ld", &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
         fgets(dummy, sizeof(dummy), fin);

         molid	--;				// note that lammps id starts from 1, not 0
         siteid	=	(id-1) % LENGTH;

         mol[molid].box		=	system;		// for now, only one system
         mol[molid].nsites	=	LENGTH;		// for now, Jan/10/2010
         mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx);
         mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly);
         mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);
         mol[molid].type[siteid]=	type - 1;	// -1 because lammps index starts from 1
      }

      for (system=0; system<NSYSTEMS; system++) {
         NMols[system]	=	0;
         NSites[system]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;			// total # of mols in certain system
            NSites[system]	+=	moli->nsites;	// total # of sites in certain system
         }
      }

      for (moli=mol; moli<mol+NMOLS; moli++) { 		
         for (i=0; i<moli->nsites; i++)  {
            moli->flags[i]	=	1;		// activate all the sites on this processor
            moli->parent[i]	=	i-1;		// initialize parent site
         }
         moli->flip		=	0;		// flip to the original direction
         moli->origin		=	CenterofMass(moli);
      }

      ///////////////////////////////////////////////////////////////////////////
      /* Start: Calculate the average position of center of mass of each chain */
      ///////////////////////////////////////////////////////////////////////////
      for (moli=mol; moli<mol+NMOLS; moli++) {
      	  temp	=	CenterofMass(moli);
          if (fabs(temp.x) < 0.45 * BOX[moli->box].lx && 
		fabs(temp.y) < 0.45 * BOX[moli->box].ly &&
		fabs(temp.z) < 0.45 * BOX[moli->box].lz) {
	     com[moli-mol]	=	V_Add(&temp, com+(moli-mol));
             ncom[moli-mol]	++;
          }
      }

      /////////////////////////////////////////////////////////////
      /* Convert coordinates and box size from SI to system unit */
      /////////////////////////////////////////////////////////////
      CoorSI2System();
      //_________________________________________________________//

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
      } 

      /////////////////////
      /* Build cell list */
      /////////////////////
#ifdef CELL_LIST	
      if (init) {
         CL_Init();		
         init	=	0;
      } 
      else 
         CL_Destroy();

      CL_Build();
#endif	/* CELL_LIST */


      //////////////////////
      /* Perform Analysis */
      //////////////////////
      if (DEBUG)	printf("Analysis: start ...\n");
/*
      CalcV();				// calc energy and virial
      if (V_VIRIAL)	Sample_Pressure();

      if (DEBUG)	printf("...sampling p2\n");

      SampleP2All();			// sample P2 and P2m and local p2
      Dist_p2();			// put local p2 into distribution

      if (DEBUG)	printf("...sampling spherical coordinates\n");

      SampleSpherical();		// sample spherical coordinates
      Dist_Spherical();			// put spherical coord. into distribution
*/
      /////////////////////////////////////////////////////////////////////
      /* Compute chain orientation distribution, test the rotator phase. */
      /////////////////////////////////////////////////////////////////////

      if (DEBUG)	printf("...computer chain orientation\n");

      for (moli=mol; moli<mol+NMOLS; moli++) {
         V_Null(&chainface);

         for (i=0; i<moli->nsites; i++) {
            if (mod(i,2)) {
	       chainface	=	V_Add(&chainface, moli->p+i);
            }
	    else {
	       chainface	=	V_Subtr(&chainface, moli->p+i);
	    }
	 }
	 orient	=	atan2(chainface.y, chainface.x)	+ M_PI;
         //temp=CenterofMass(moli);
         //if (temp.z < -0.25 * BOX[system].lz && temp.z <0) 
	 orientdist[(int)(orient/(2*M_PI)*36)]	++;
      }

      //////////////////////////////////////////////
      /* Chain conformation statistics (3/2/2011) */
      //////////////////////////////////////////////

      Rg2 	= 	0.0;	DRg2	=	0.0;      nRg2	=	0.0;
      Ro2	=	0.0;    DRo2	=	0.0;      nRo2	=	0.0;

      for (moli=mol; moli<mol+NMOLS; moli++) {
         Rg2	+=	R2_gyration(moli);
	 Ro2	+=	R2_n2n(moli);
	 nRg2	++;
	 nRo2	++;
      }
      Rg2	/=	nRg2;
      Ro2	/=	nRo2;

      /////////////////////////////////////////////////////////////////////////////////
      /* calculate on average how many chains appear in the neighborhood of one bead */
      /////////////////////////////////////////////////////////////////////////////////
 
      nchains	=	0.0;
      system	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         for (i=0; i<moli->nsites; i++) {

             for (molj=mol; molj<mol+NMOLS; molj++) {
                if (moli==molj)	continue;

                for (j=0; j<molj->nsites; j++) {

                   R2	=	DistSQ(moli->p[i], molj->p[j], system);
                   if (R2 < Rcut2) {
                      nchains	+=	1;
                      break;
                   }
                }
             }
         }
      }
      nchains	/=	NSITES;

      /////////////////////////////////////////////////////////////////////
      /* calculate the # of sites with local p2 greater than a threshold */
      /////////////////////////////////////////////////////////////////////

      nsitesp2	=	0;
      nmolsp2	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
//       nsitesp2	=	0;
         for (i=0; i<moli->nsites; i++) {
            if (moli->p2[i] > 0.4)
               nsitesp2	++;
         }
//            if (nsitesp2>1 && nsitesp2<6)
//	       nmolsp2	++;
      }

      /////////////////////////////////////////
      /* Correlation and Radial distribution */
      /////////////////////////////////////////

      if (DEBUG)	printf("Calculate correlation: start...\n");

      correlation();			// calculate correlation, in system units
      if (!(timestep%IRADIAL)) {
         radial("sample");		// sample radial distribution function
      }

//      shiftbox(system, nucleus, n);

      /////////////////////////////
      /* Output analysis results */
      /////////////////////////////
      
      if (DEBUG)	printf("Output: start ...\n");

      Print_hst(fhst);		// print out histgram
      CoorSystem2SI();		// convert coordinates and box size back to SI units

      ////////////////////////////////////////////
      /* OUTPUT .car file for VMD visualization */
      ////////////////////////////////////////////
      
      if (carflag) {
         fprintf(fout, "!BIOSYM archive 3\n");
         fprintf(fout, "PBC=ON\n");
         fprintf(fout, "!TIMESTEP %d\n", timestep);
         fprintf(fout, "!DATE %s", asctime(localtime(&t)));
         fprintf(fout, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

	 n	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {

               MolInBox2(moli);
               for (i=0; i<moli->nsites; i++) {

                  if (moli->nuclid[i]>0)	// crystal like particle
                     sprintf(s, "N%d", n++);	// N: blue color in VMD
                  else
	             sprintf(s, "O%d", n++);	// O: red color in VMD

/*
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])		// note nuclid index starts from 1
                     sprintf(s, "M%d", n++);
                  else if (sizeofnucl[moli->nuclid[i]] -MAXSIZE[system] > -3 && MAXSIZE[system]>=10)
                     sprintf(s, "C%d", n++);
	          else if (moli->nuclid[i] >= 1)
                     sprintf(s, "O%d", n++);
                  else
                     sprintf(s, "H%d", n++);
*/
                  fprintf(fout, "%-5.5s ", s);
                  sprintf(s, "M%d", moli-mol);
                  fprintf(fout, "%14.8g %14.8g %14.8g ", moli->p[i].x, moli->p[i].y, moli->p[i].z);
                  strcpy(ff, "O");
                  fprintf(fout, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[i], s));
               } 
            }   
         }
         fprintf(fout, "end\nend\n");
      }
      fflush(fout);
      //___________OUTPUT .car file_____________//


      ///////////////////////////////////////////
      // OUTPUT .pdb file for further analysis */
      ///////////////////////////////////////////

      if (pdbflag) {
         fprintf(fpdb, "HEADER: file created from %s on %s", argv[1], asctime(localtime(&t)));
         fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
        	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

	 m		=	0;	// molecule sequence number
         n		=	0;	// atom sequence number
#define SIZECAP	5
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {
               //MolInBox2(moli);

               drawmol	=	0;
               for (i=0; i<moli->nsites; i++) {
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {
	             drawmol	=	1;		// participate the biggest nucleus
		     break;
                  }
	       }
//temp=CenterofMass(moli);
//if (temp.z < -0.25 * BOX[system].lz) {
               m	++; 
               for (i=0; i<moli->nsites; i++) {
                  if (drawmol) {
                     if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {	// nuclid index starts from 1
                        atomname	=	'N';	// N: blue color in VMD
		        fprintf(fdat, " 10");
                     }
                     else {
                        atomname	=	'O';	// O: red color in VMD
                        fprintf(fdat, " 0");
                     }
                  }
                  else {
		     atomname	=	'C';		// C: cyan color in VMD
		     fprintf(fdat," -1");
	          }

                  n	++;
	          fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
                  fprintf(fpdb, "%5d ", n);	// atom number
                  fprintf(fpdb, " %c  ", atomname);	// atom name
                  fprintf(fpdb, " ");		// alternate location indiator
  	          fprintf(fpdb, " C8");		// residue name
	          fprintf(fpdb, " ");		// column 21
                  fprintf(fpdb, " ");		// chain identifier, column 22
	          fprintf(fpdb, "%4d", m);	// residue sequence number, 23-26
	          fprintf(fpdb, " ");		// code for insertion of residues, 27
                  fprintf(fpdb, "   ");		// column 28-30
                  fprintf(fpdb, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z);
                  fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                  fprintf(fpdb, "%5.5s", "");
                  fprintf(fpdb, "\n"); 

                  if (imagen) {			// for image box
                     n	++;
                     fprintf(fpdb, "ATOM  ");
                     fprintf(fpdb, "%5d ", n);
                     fprintf(fpdb, " %c  ", atomname);
                     fprintf(fpdb, " ");
                     fprintf(fpdb, " C8");
	             fprintf(fpdb, " ");
                     fprintf(fpdb, " ");
	             fprintf(fpdb, "%4d", m);
                     fprintf(fpdb, " ");
                     fprintf(fpdb, "   ");
                     fprintf(fpdb, "%8.3f%8.3f%8.3f", 
			moli->p[i].x + imagex, moli->p[i].y+imagey, moli->p[i].z+imagez);
                     fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                     fprintf(fpdb, "%5.5s", "");
                     fprintf(fpdb, "\n"); 
                  }
               } 
//}
            }   
         }
         fprintf(fpdb, "END\n");
         fflush(fpdb);
      }	//pdbflag
      //____________OUTPUT .pdb file___________//

      ////////////////////////////////////////////////////
      /* OUTPUT configuration file for further analysis */
      ////////////////////////////////////////////////////

      if (confflag) { 
         fprintf(fconf, "!TIMESTEP %d\n", timestep);
         fprintf(fconf, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
         for (i=0; i<NSYSTEMS; i++)
            fprintf(fconf, "%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

         for (moli=mol; moli<mol+NMOLS; moli++) {
            fprintf(fconf, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
            //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
            //MolInBox(moli);
            for (i=0; i<moli->nsites; i++) 
               fprintf(fconf, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
         }
         fflush(fconf);
      }
      //____________OUTPUT configuration file___________//
   } 

   /////////////////////////////////////////////
   /* output average center of mass of chains */
   /////////////////////////////////////////////

   if (pdbflag) {
      fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
      fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
 	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

      n		=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         i	=	moli-mol;
         if (ncom[i]>0) {
//if( fabs(com[i].x/ncom[i])<3 && fabs(com[i].y/ncom[i])<3 && fabs(com[i].z/ncom[i])<7 ) {
            n	++;
            atomname	=	'N';
            fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
            fprintf(fpdb, "%5d ", n);		// atom number
            fprintf(fpdb, " %c  ", atomname);	// atom name
            fprintf(fpdb, " ");		// alternate location indiator
            fprintf(fpdb, "   ");	// residue name
            fprintf(fpdb, " ");		// column 21
            fprintf(fpdb, " ");		// chain identifier, column 22
            fprintf(fpdb, "    ");	// residue sequence number, 23-26
            fprintf(fpdb, " ");		// code for insertion of residues, 27
            fprintf(fpdb, "   ");	// column 28-30
            fprintf(fpdb, "%8.3f%8.3f%8.3f", com[i].x/ncom[i], com[i].y/ncom[i], com[i].z/ncom[i]);
            fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
            fprintf(fpdb, "\n"); 
         }
//}
      }
      fprintf(fpdb, "END\n");
      fflush(fpdb);
   }
   //______OUTPUT com of chains to .pdb file__//

   S_PrintAll();		// print out distributions
   corr_norm();			// normalize correlation function
   corr_print();		// print out correlation function
   //radial("print");		// print out radial distribution function

   for (i=0; i<180; i++) {
      printf("%d %d\n", i, orientdist[i]);
   }

   if (DEBUG)	printf("Closing files ...\n");

   fclose(fin);		// close files
   fclose(fhst);
   if (carflag)		fclose(fout);
   if (confflag)    	fclose(fconf);
   if (pdbflag)	     {	fclose(fpdb); fclose(fdat);}

   fflush(stdout);
   fclose(stdout);

   return	0;
}
