/*
	program:	dumpreduce.c
	author:		Peng Yi at MIT
	date:		Feb. 15, 2009
	purpose:	reduce the size of a dump file by only keeping 
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
   long		i, nsystems, nmols, nsites, dn, counter, timestep;
   char		filein[255], fileout[255];
   char		dummy[255], eol[255];

   if (argc<3) {
      printf("dumpreduce (c) 2009 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tdumpreduce dumpfile dn\n\n");
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

         fscanf(fin, "%s%ld", dummy, &timestep);			// read in timestep
         fgets(eol, sizeof(eol), fin);					// get end of line

         if (strcmp(dummy, "TIMESTEP"))	break;				// if hit end of file, exit

         fscanf(fin, "%ld%ld%ld", &nsystems, &nmols, &nsites);		// # of systems, mols, and sites
         fgets(eol, sizeof(eol), fin);

         if (mod(counter, dn)==0) {					// output every dn configurations
	    fprintf(fout, "TIMESTEP	%d\n", timestep);
	    fprintf(fout, "%d\t%d\t%d\n", nsystems, nmols, nsites);
         }
         for (i=0; i<nsystems+nmols+nsites; i++) {
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
