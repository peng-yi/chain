/*
	program:	conf2lammps.c
	author:		Peng Yi at MIT
	date:		Feb. 14, 2008
	purpose:	convert configuration files to lammps data files
	note:		need corresponding setup file
*/

#define __MAIN_PROGRAM
#include "header.h"


void findatomtypes(long * natomtypes)		// find the number of atom types
{						// unfinished ...
   molstruct	*moli;
   long		i;
}


int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   long		i, system = 0, n = 0, j;
   char		s[255], ff[255], filename[255];
   molstruct	*moli;
   FILE		*fPtr;
   long		nsite, nbond, nangle, ndihedral;
   long		natomtype, nbondtype, nangletype, ndihedraltype;
   vector	r;

   if (argc<3) {
      printf("conf2lammps (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tconf2lammps atom_style conf.file\n\n");
      printf("\tatom_style = lj or chain\n\n");
      printf("Notes:\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   InitMols(MAXNMOLS, MAXNMOLS);
   GetSetup(argv);

   Read_Conf(argv[2]);				// read in configuration file

   //sprintf(filename, argv[1]);
   //strcat(filename, ".lammpsdata");

   if ( (fPtr=fopen("lammps.data", "w"))==NULL )
      Exit("conf2lammps", "main", "Open output lammps data file failed");
   else {
      fprintf(fPtr, "LAMMPS Description\n\n");

      if (!strcmp(argv[1], "chain")) {
         fprintf(fPtr, "%d atoms\n", NSITES);
         fprintf(fPtr, "%d bonds\n", NSITES - NMOLS);
         fprintf(fPtr, "%d angles\n", NSITES - 2*NMOLS);
         fprintf(fPtr, "%d dihedrals\n", NSITES - 3*NMOLS);
         fprintf(fPtr, "\n");
         fprintf(fPtr, "%d atom types\n", 4);
         fprintf(fPtr, "%d bond types\n", 1);
         fprintf(fPtr, "%d angle types\n", 1);
         fprintf(fPtr, "%d dihedral types\n", 1);
         fprintf(fPtr, "%f %f xlo xhi\n", -0.5*BOX[0].lx, 0.5*BOX[0].lx);
         fprintf(fPtr, "%f %f ylo yhi\n", -0.5*BOX[0].ly, 0.5*BOX[0].ly);
         fprintf(fPtr, "%f %f zlo zhi\n", -0.5*BOX[0].lz, 0.5*BOX[0].lz);
      }
      else if (!strcmp(argv[1], "lj")) {
         fprintf(fPtr, "%d atoms\n", NSITES);
         fprintf(fPtr, "%d atom types\n", 2);
         fprintf(fPtr, "%f %f xlo xhi\n", -0.5*BOX[0].lx, 0.5*BOX[0].lx);
         fprintf(fPtr, "%f %f ylo yhi\n", -0.5*BOX[0].ly, 0.5*BOX[0].ly);
         fprintf(fPtr, "%f %f zlo zhi\n", -0.5*BOX[0].lz, 0.5*BOX[0].lz);
      }

      if (!strcmp(argv[1], "chain")) {
         fprintf(fPtr, "\nAtoms\n\n");
         nsite	=	0;
         system	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            for (i=0; i<moli->nsites; i++) {
               nsite	++;

               //r	=	MapInBox2(moli->p+i, PBC, system);
               r	=	moli->p[i];

               fprintf(fPtr, "%ld %ld %ld %lf %lf %lf\n", nsite, moli-mol+1, moli->type[i]+1, r.x, r.y, r.z);
//               fprintf(fPtr, "%d %d %d %f %f %f\n", nsite, moli-mol+1, moli->type[i]+1, moli->p[i].z, moli->p[i].x, moli->p[i].y);	// rotate coordinate system
            }
         }

         fprintf(fPtr, "\nBonds\n\n");
         nbond	=	0;
         nsite	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            for (i=0; i<moli->nsites-1; i++) {
               nbond	++;
               fprintf(fPtr, "%d %d %d %d\n", nbond, 1, nsite+i+1, nsite+i+2);
            }
            nsite	+=	moli->nsites;
         }

         fprintf(fPtr, "\nAngles\n\n");
         nangle	=	0;
         nsite	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            for (i=0; i<moli->nsites-2; i++) {
               nangle	++;
               fprintf(fPtr, "%d %d %d %d %d\n", nangle, 1, nsite+i+1, nsite+i+2, nsite+i+3);
            }
            nsite	+=	moli->nsites;
         }

         fprintf(fPtr, "\nDihedrals\n\n");
         ndihedral	=	0;
         nsite	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            for (i=0; i<moli->nsites-3; i++) {
               ndihedral	++;
               fprintf(fPtr, "%d %d %d %d %d %d\n", ndihedral, 1, nsite+i+1, nsite+i+2, nsite+i+3, nsite+i+4);
            }
            nsite	+=	moli->nsites;
         }
      }
      else if (!strcmp(argv[1], "lj")) {
         fprintf(fPtr, "\nAtoms\n\n");
         nsite	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            for (i=0; i<moli->nsites; i++) {
               nsite	++;
               fprintf(fPtr, "%d %d %f %f %f\n", nsite, moli->type[i]+1, moli->p[i].x, moli->p[i].y, moli->p[i].z);
            }
         }
      }
      fclose(fPtr);
   }

   return	0;
}
