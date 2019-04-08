/*
	program:	dumpreorder.c
	author:		Peng Yi at MIT
	date:		Feb. 28, 2009
	purpose:	when two dump files combined, reorder the timestep
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
      printf("dumpreorder (c) 2009 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tdumpreorder dumpfile ITAPE\n\n");
      printf("Notes:\n");
      printf("\t* dumpfile is a combined dump file\n");
      printf("\t* ITAPE is the dumping frequency\n\n");
      exit(1);
   }

   strcpy(filein, argv[1]);
   strcpy(fileout, filein);
   strcat(fileout, ".reorder");
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

	 fprintf(fout, "TIMESTEP	%d\n", counter);
	 fprintf(fout, "%d\t%d\t%d\n", nsystems, nmols, nsites);

         for (i=0; i<nsystems+nmols+nsites; i++) {
            fgets(dummy, sizeof(dummy), fin);
	    fputs(dummy, fout);
         }
	 counter	+=	dn;
      } 

      fclose(fin);
      fclose(fout);
   }   

   return	0;
}
