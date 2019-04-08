/*
	program:	conf2gro.c
	author:		Peng Yi at MIT
	date:		Sept 11, 2009
	purpose:	convert configuration files to GROMACS input
	note:		require setup file
*/

#define __MAIN_PROGRAM
#include "header.h"

#define VERSION		"Sept_11_2009"

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
      printf("conf2GROMACS (c) 2009 by Peng Yi at MIT\n\n");
      printf("convert configuration file to GROMACS input files.\n");
      printf("Usage:\n");
      printf("\tconf2gro filename\n\n");
      printf("Notes:\n");
      printf("\t* require setup file\n");
      printf("\t* length unit is nm, NOT Angstrom\n\n");
      exit(1);
   }

   strcpy(filename1, argv[1]);
   strcat(filename1, ".gro");
   if (! (fcfg = fopen(filename1, "w")))      exit(1);

   strcpy(filename2, argv[1]);
   strcat(filename2, ".top");
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

   // generate GROMACS .gro file

   fprintf(fcfg, "GROMACS coordinate file created on %s", asctime(localtime(&t)));
   fprintf(fcfg, "%5d\n", NSITES);
   n	=	1;				// index starts from 1
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         fprintf(fcfg, "%5d", moli-mol+1);
         fprintf(fcfg, "%-5s", "C8");
         fprintf(fcfg, "%5s", (i==0 ||i==7 ? "CH3" : "CH2"));
         fprintf(fcfg, "%5d", n);
         fprintf(fcfg, "%8.3f%8.3f%8.3f", moli->p[i].x/10, moli->p[i].y/10, moli->p[i].z/10);
         fprintf(fcfg, "%8.3f%8.3f%8.3f", 0.0, 0.0, 0.0);	// velocity
         fprintf(fcfg, "\n");
         n	++;
      }
   }
   fprintf(fcfg, "%10.5f%10.5f%10.5f", BOX[system].lx/10, BOX[system].ly/10, BOX[system].lz/10);
   fprintf(fcfg, "%10.5f%10.5f%10.5f", 0.0, 0.0, 0.0);
   fprintf(fcfg, "%10.5f%10.5f%10.5f", 0.0, 0.0, 0.0);
   fprintf(fcfg, "\n");
   fclose(fcfg);

   // generate GROMACS .top file

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
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[1], 15.0, 0.0, 1);	 
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[2], 14.0, 0.0, 6); 	// 6 middle beads
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[1], 15.0, 0.0, 1); 

   fprintf(ffld, "BONDS	%d\n", NSITES/NMOLS-1);		// # of bonds in each type 1 molecule
   fprintf(ffld, "harm	%d %d %f %f\n", 1, 2, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 2, 3, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 3, 4, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 4, 5, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 5, 6, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 6, 7, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 7, 8, kl, l0);

   fprintf(ffld, "ANGLES %d\n", NSITES/NMOLS-2);	// # of angles in each type 1 molecule
   fprintf(ffld, "harm	%d %d %d %f %f\n", 1, 2, 3, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 2, 3, 4, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 3, 4, 5, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 4, 5, 6, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 5, 6, 7, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 6, 7, 8, ktheta, theta0);
   
   fprintf(ffld, "DIHEDRALS %d\n", NSITES/NMOLS-3);	// # of dihes in each type 1 molecule
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 1, 2, 3, 4, k0, k1, k2, k3);
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 2, 3, 4, 5, k0, k1, k2, k3);
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 3, 4, 5, 6, k0, k1, k2, k3);
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 4, 5, 6, 7, k0, k1, k2, k3);
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 5, 6, 7, 8, k0, k1, k2, k3);

   fprintf(ffld, "FINISH\n");				// finish type 1 molecule

   /* Register molecule type 2 */

   fprintf(ffld, "octane2\n");				// name of type 2
   fprintf(ffld, "NUMMOLS %d\n", ntype[2]);		// # of type 2 molecules

   fprintf(ffld, "ATOMS %d\n", NSITES/NMOLS);		// # of atoms in each molecule
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[3], 15.0, 0.0, 1);	 
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[4], 14.0, 0.0, 6); 	// 6 middle beads
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[3], 15.0, 0.0, 1); 
   fprintf(ffld, "BONDS	%d\n", NSITES/NMOLS-1);		// # of bonds in each molecule
   fprintf(ffld, "harm	%d %d %f %f\n", 1, 2, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 2, 3, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 3, 4, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 4, 5, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 5, 6, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 6, 7, kl, l0);
   fprintf(ffld, "harm	%d %d %f %f\n", 7, 8, kl, l0);
   fprintf(ffld, "ANGLES %d\n", NSITES/NMOLS-2);	// # of angles in each molecule
   fprintf(ffld, "harm	%d %d %d %f %f\n", 1, 2, 3, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 2, 3, 4, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 3, 4, 5, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 4, 5, 6, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 5, 6, 7, ktheta, theta0);
   fprintf(ffld, "harm	%d %d %d %f %f\n", 6, 7, 8, ktheta, theta0);
   fprintf(ffld, "DIHEDRALS %d\n", NSITES/NMOLS-3);	// # of dihes in each molecule
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 1, 2, 3, 4, k0, k1, k2, k3);
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 2, 3, 4, 5, k0, k1, k2, k3);
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 3, 4, 5, 6, k0, k1, k2, k3);
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 4, 5, 6, 7, k0, k1, k2, k3);
   fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", 5, 6, 7, 8, k0, k1, k2, k3);
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
