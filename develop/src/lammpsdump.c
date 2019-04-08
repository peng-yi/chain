/*
	program:	lammpsdump.c
	author:		Peng Yi at MIT
	date:		July 7, 2008
	purpose:	reduce the size of a lammps dump file by only keeping 
			one of every dn configurations
*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

long mod (long numb, long divisor)
{
  while (numb<0) {
    numb+=divisor;
  }
  return numb%divisor;
}

int main(int argc, char *argv[])
{
   FILE		*fin, *fout;
   long		i, nsites, dn, counter, timestep;
   char		filein[255], fileout[255];
   char		a[255], b[255], c[255], d[255], dummy[255];
   char		dimx[255], dimy[255], dimz[255];

   if (argc<3) {
      printf("lammpsdump (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tlammpsdump dumpfile dn\n\n");
      printf("Notes:\n");
      printf("\t* dn is an integer, this command will keep one out of \n");
      printf("\t* every dn configurations in the original dump file,\n\n");
      exit(1);
   }

   strcpy(filein, argv[1]);
   strcpy(fileout, filein);
   strcat(fileout, ".");
   strcat(fileout, argv[2]);
   dn	=	atol(argv[2]);

   if ( (fin=fopen(filein, "r"))==NULL )
      exit(1);
   else if ( (fout=fopen(fileout, "w"))==NULL )
      exit(1);
   else {
      counter	=	0;

      while (!feof(fin)) {

	 /* Read in one snapshot from dump file */

         if (!fgets(a, sizeof(a), fin))		break;
         fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);
         fgets(b, sizeof(b), fin);
         fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);
         fgets(c, sizeof(c), fin);
         fgets(dimx, sizeof(dimx), fin);
         fgets(dimy, sizeof(dimy), fin);
         fgets(dimz, sizeof(dimz), fin);
         fgets(d, sizeof(d), fin);

         if (mod(counter, dn)==0) {			// output every dn configurations
            fprintf(fout, "ITEM: TIMESTEP\n");
	    fprintf(fout, "%d\n", timestep);
	    fprintf(fout, "ITEM: NUMBER OF ATOMS\n");
            fprintf(fout, "%d\n", nsites);
	    fprintf(fout, "ITEM: BOX BOUNDS\n");
	    fputs(dimx, fout);
	    fputs(dimy, fout);
	    fputs(dimz, fout);
	    fputs(d, fout);
         } 

         for (i=0; i<nsites; i++) {
            fgets(dummy, sizeof(dummy), fin);

	    if (mod(counter, dn)==0) {
	       fputs(dummy, fout);
            }
         }
	 counter	++;
      } 

      fclose(fin);
      fclose(fout);
   }   

   return	0;
}
