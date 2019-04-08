/*
	program:	lammps2car.c
	author:		Peng Yi at MIT
	date:		Feb. 11, 2008
	purpose:	convert lammps dump file to .car files
*/

#define __MAIN_PROGRAM
#include "header.h"

char *Element(long t, char *s)
{
   double	mass = type[t].M;

   s[0]		=	0;
/*
   if (14==mass)		strcpy(s, "C");
   else if (15==mass)		strcpy(s, "C");
   else if (1.01==mass)		strcpy(s, "H");
   else if (28.086==mass)	strcpy(s, "Si");
   else if (26.982==mass)	strcpy(s, "Al");
   else if (16==mass)		strcpy(s, "O");
*/
   strcpy(s, "C");
   return	s;
}

#define LENGTH	8		// monodisperse chain length

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   long		i, system = 0, n, j;
   char		s[255], ff[255];
   molstruct	*moli;
   vector	com;		// center of mass
   FILE		*fPtr;

   char		filein[255];
   FILE		*fin, *fout, *fhst;
   double	x, y, z, vx, vy, vz;
   double	xhi, xlo, yhi, ylo, zhi, zlo;
   long		timestep, nsites, record;
   long		nx, ny, nz, id, siteid, molid, type;
   char		a[255], b[255], c[255], d[255], dummy[255];

   InitMols(MAXNMOLS, MAXNMOLS);
   GetSetup(argv);
   InitUnits();

   strcpy(filein, argv[1]);

   if ( (fin=fopen(filein, "r"))==NULL )
      Exit("lammps2car", "main", "open input file failed.");
   else if ( (fout=fopen("lammps.car", "w"))==NULL )
      Exit("lammps2car", "main", "open output file failed.");
   else if ( (fhst=fopen("hst", "w"))==NULL )
      Exit("lammps2car", "main", "open hst file failed.");
   else {

      while (!feof(fin)) {

	 /* read in */

         if (!fgets(a, sizeof(a), fin))
	    break;
         fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);
         fgets(b, sizeof(b), fin);
         fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);
         fgets(c, sizeof(c), fin);
         fscanf(fin, "%lf%lf", &xlo, &xhi); 	fgets(dummy, sizeof(dummy), fin);
         fscanf(fin, "%lf%lf", &ylo, &yhi); 	fgets(dummy, sizeof(dummy), fin);
         fscanf(fin, "%lf%lf", &zlo, &zhi);	fgets(dummy, sizeof(dummy), fin);
         fgets(d, sizeof(d), fin);

	 system		=	0;		// for now (2/12/08), only one system
      	 BOX[system].lx	=	xhi-xlo;
         BOX[system].ly	=	yhi-ylo;
         BOX[system].lz	=	zhi-zlo;

         for (i=0; i<nsites; i++) {
            fscanf(fin, "%ld%ld %lf%lf%lf %lf%lf%lf %ld%ld%ld", 
			&id, &type, &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
            fgets(dummy, sizeof(dummy), fin);

            id		--;				// because lammps starts from 1 not 0
            type	--;
            molid			=	(long) (id/LENGTH);
	    siteid			=	mod(id, LENGTH);
            mol[molid].box		=	system;		// for now, only one system
            mol[molid].nsites		=	LENGTH;
            mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx);
            mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly);
            mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);
            mol[molid].type[siteid]	=	type;
         }

         /* analysis */
/* 
         InitSample();
         SampleP2All();
         Dist_p2();
         Find_Nuclei();
         S_PrintAll();
	 fprintf(fhst, "%f %d %d\n", P2[0], MAXSIZE[0], Xtal[0]);
*/
	 /* output */

         fprintf(fout, "!BIOSYM archive 3\n");
         fprintf(fout, "PBC=ON\n\n");
         fprintf(fout, "!DATE %s", asctime(localtime(&t)));
         fprintf(fout, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

	 n	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {

               //MolInBox2(moli);
               for (i=0; i<moli->nsites; i++) {
/*
                  if (moli->p2[i]>critp2)
                     sprintf(s, "M%d", n++);
                  else
	             sprintf(s, "C%d", n++);
*/
                  sprintf(s, "M%d", n++);
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
      fclose(fout);
      fclose(fin);
      fclose(fhst);
   }   

   /* Print out coordinates for .car file */
#ifdef TEST
   printf("!BIOSYM archive 3\n");
   if (1==PBC) {
      printf("PBC=ON\n\n");
      printf("!DATE %s", asctime(localtime(&t)));
      printf("PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);
   }
   else if (0==PBC) {
      printf("PBC=OFF\n\n");
      printf("!DATE %s", asctime(localtime(&t)));
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {

         MolInBox2(moli);
         for (i=0; i<moli->nsites; i++) {

            if (moli->p2[i]>critp2)
               sprintf(s, "M%d", n++);
	    else
	       sprintf(s, "C%d", n++);

//	    sprintf(s, "M%d", n++);
            printf("%-5.5s ", s);
            sprintf(s, "M%d", moli-mol);
            printf("%14.8g %14.8g %14.8g ", moli->p[i].x, moli->p[i].y, moli->p[i].z);
            strcpy(ff, "O");
            printf("%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[i], s));
         } 
      }   
   }
   printf("end\nend\n");

   /* print out center of mass of chains in .car file format */

   n	=	0;
   printf("!BIOSYM archive 3\n");
   if (1==PBC) {
      printf("PBC=ON\n\n");
      printf("!DATE %s", asctime(localtime(&t)));
      printf("PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);
   }
   else {
      printf("PBC=OFF\n\n");
      printf("!DATE %s", asctime(localtime(&t)));
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {
         com	=	CenterofMass(moli);
         sprintf(s, "M%d", n++);
         printf("%-5.5s ", s);
         sprintf(s, "M%d", moli-mol);
         printf("%14.8g %14.8g %14.8g ", com.x, com.y, com.z);
         strcpy(ff, "O");
         printf("%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[0], s));
      }   
   }
   printf("end\nend\n");
#endif
   return	0;
}
