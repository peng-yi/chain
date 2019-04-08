/*
	program:	conf2DL_POLY.c
	author:		Peng Yi at MIT
	date:		July 30, 2009
	purpose:	convert configuration files to DL_POLY input
	note:		require setup file
*/

#define __MAIN_PROGRAM
#include "header.h"

#define VERSION		"Jul_30_2009"

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   long		i, system, n, j;
   char		s[80], filename1[80], filename2[80], filename3[80];
   molstruct	*moli;
   vector	com;		// center of mass
   FILE		*fcfg, *ffld, *fpdb;
   long		nsite, nbond, nangle, ndihedral;
   long		realXtal;
   char		atomname;

   char		typename[16][80];	//max 15 types
   long		ntype[16];

   double	sig, epsilon, kl, ktheta, l0, theta0, k0, k1, k2, k3;

   if (argc<2) {
      printf("confDL_POLY (c) 2009 by Peng Yi at MIT\n\n");
      printf("convert configuration file to DL_POLY input files.\n");
      printf("Usage:\n");
      printf("\tconf2DL_POLY filename\n\n");
      printf("Notes:\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   strcpy(filename1, argv[1]);
   strcat(filename1, ".CONFIG");
   if (! (fcfg = fopen(filename1, "w")))      exit(1);

   strcpy(filename2, argv[1]);
   strcat(filename2, ".FIELD");
   if (! (ffld = fopen(filename2, "w")))      exit(1);

   strcpy(filename3, argv[1]);
   strcat(filename3, ".pdb");
   if (! (fpdb = fopen(filename3, "w")))      exit(1);

   InitMols(MAXNMOLS, MAXNMOLS);
   GetSetup(argv);

   Read_Conf(argv[1]);
   system	=	0;

//   for (i=1; i<NTYPES+1; i++)
//      sprintf(typename+i, "C%d", i);
   sprintf(typename+1, "C1");		// temporary
   sprintf(typename+2, "C2");
   sprintf(typename+3, "C3");
   sprintf(typename+4, "C4");
   for (i=0; i<16; i++)
      ntype[i]	=	0;

   // generate DL_POLY CONFIG file, fixed-formatted

   fprintf(fcfg, "DL_POLY CONFIG file created from %s on %s", argv[1], asctime(localtime(&t)));
   fprintf(fcfg, "%10d%10d%10d%20f\n", 2, 3, NSITES, 0);
   fprintf(fcfg, "%20f%20f%20f\n", BOX[system].lx, 0.0, 0.0); 
   fprintf(fcfg, "%20f%20f%20f\n", 0.0, BOX[system].ly, 0.0); 
   fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, BOX[system].lz); 

   n	=	1;					// index starts from 1
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->type[0]==0) {				// printout all molecules of type 1
         ntype[1]	++;
         for (i=0; i<moli->nsites; i++) {
            fprintf(fcfg, "%-8s%10d\n", typename[moli->type[i]+1], n);
            fprintf(fcfg, "%20f%20f%20f\n", moli->p[i].x, moli->p[i].y, moli->p[i].z);
            fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, 0.0);
            fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, 0.0);
            n	++;
         }
      }
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->type[0]==2) {				// printout all molecules of type 2
         ntype[2]	++;
         for (i=0; i<moli->nsites; i++) {
            fprintf(fcfg, "%-8s%10d\n", typename[moli->type[i]+1], n);
            fprintf(fcfg, "%20f%20f%20f\n", moli->p[i].x, moli->p[i].y, moli->p[i].z);
            fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, 0.0);
            fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, 0.0);
            n	++;
         }
      }
   }
   fclose(fcfg);

   // generate DL_POLY FIELD file that contains force field parameters

   sig		=	4.01;		// in A
   epsilon	=	0.112;		// in kcal/mol
   kl		=	700;		// in kcal/mol/A-2
   l0		=	1.53;		// in A
   ktheta	=	120;
   theta0	=	109.5;
   k0		=	0;
   k1		=	1.6;
   k2		=	-0.867;
   k3		=	3.24;

   fprintf(ffld, "DL_POLY FIELD file created from %s on %s", argv[1], asctime(localtime(&t)));

   fprintf(ffld, "\nunits	kcal\n\n");		// i/o energy unit kcal/mol
   fprintf(ffld, "molecular type %d\n", 2);		// # of molecular types

   /* Register molecule type 1 */

   fprintf(ffld, "octane1\n");				// name of type 1
   fprintf(ffld, "NUMMOLS %d\n", ntype[1]);		// # of type 1 molecules

   fprintf(ffld, "ATOMS %d\n", NSITES/NMOLS);		// # of atoms in each type 1 molecule
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[1], 15.0, 0.0, 1);	 	// head bead
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[2], 14.0, 0.0, NSITES/NMOLS-2); 	// middle beads
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[1], 15.0, 0.0, 1); 		// tail bead

   fprintf(ffld, "BONDS	%d\n", NSITES/NMOLS-1);		// # of bonds in each type 1 molecule
   for (i=1; i<=NSITES/NMOLS-1; i++) {
      fprintf(ffld, "harm	%d %d %f %f\n", i, i+1, kl, l0);
   }

   fprintf(ffld, "ANGLES %d\n", NSITES/NMOLS-2);	// # of angles in each type 1 molecule
   for (i=1; i<=NSITES/NMOLS-2; i++) {
      fprintf(ffld, "harm	%d %d %d %f %f\n", i, i+1, i+2, ktheta, theta0);
   }
   
   fprintf(ffld, "DIHEDRALS %d\n", NSITES/NMOLS-3);	// # of dihes in each type 1 molecule
   for (i=1; i<=NSITES/NMOLS-3; i++) {
      fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", i, i+1, i+2, i+3, k0, k1, k2, k3);
   }

   fprintf(ffld, "FINISH\n");				// finish type 1 molecule

   /* Register molecule type 2 */

   fprintf(ffld, "octane2\n");				// name of type 2
   fprintf(ffld, "NUMMOLS %d\n", ntype[2]);		// # of type 2 molecules

   fprintf(ffld, "ATOMS %d\n", NSITES/NMOLS);		// # of atoms in each molecule
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[3], 15.0, 0.0, 1);	 
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[4], 14.0, 0.0, NSITES/NMOLS-2); 	// 6 middle beads
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[3], 15.0, 0.0, 1); 
   fprintf(ffld, "BONDS	%d\n", NSITES/NMOLS-1);		// # of bonds in each molecule
   for (i=1; i<=NSITES/NMOLS-1; i++) {
      fprintf(ffld, "harm	%d %d %f %f\n", i, i+1, kl, l0);
   }

   fprintf(ffld, "ANGLES %d\n", NSITES/NMOLS-2);	// # of angles in each molecule
   for (i=1; i<=NSITES/NMOLS-2; i++) {
      fprintf(ffld, "harm	%d %d %d %f %f\n", i, i+1, i+2, ktheta, theta0);
   }

   fprintf(ffld, "DIHEDRALS %d\n", NSITES/NMOLS-3);	// # of dihes in each molecule
   for (i=1; i<=NSITES/NMOLS-3; i++) {
      fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", i, i+1, i+2, i+3, k0, k1, k2, k3);
   }

   fprintf(ffld, "FINISH\n");				// finish type 2 molecule

   /* Summary of VDW interaction for all types of molecules */

   fprintf(ffld, "VDW	%d\n", 10);			// nonbonded interaction
   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[1], typename[1], "lj", epsilon, sig);
   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[2], typename[2], "lj", epsilon, sig);
   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[1], typename[2], "lj", epsilon, sig);

   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[3], typename[3], "lj", 2*epsilon, sig);
   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[4], typename[4], "lj", 2*epsilon, sig);
   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[3], typename[4], "lj", 2*epsilon, sig);

   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[1], typename[3], "lj", sqrt(2)*epsilon, sig);
   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[1], typename[4], "lj", sqrt(2)*epsilon, sig);
   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[2], typename[3], "lj", sqrt(2)*epsilon, sig);
   fprintf(ffld, "%-8s%-8s%-8s%8.3f%8.3f\n", typename[2], typename[4], "lj", sqrt(2)*epsilon, sig);

   fprintf(ffld, "CLOSE\n");
   fclose(ffld);

   /* Print out .pdb file for the whole system */

   fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
   fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
	BOX[system].lx, BOX[system].ly, BOX[system].lz,
	90.0, 90.0, 90.0);

   n		=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {
         for (i=0; i<moli->nsites; i++) {
	    n	++;
            atomname	=	'O';

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
		moli->p[i].x, moli->p[i].y, moli->p[i].z);
	    fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
	    fprintf(fpdb, "\n"); 
         } 
      }   
   }
   fprintf(fpdb, "END\n");
   fclose(fpdb);

   return	0;
}
