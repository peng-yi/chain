#include "header.h"

int main(int argc, char * argv[])
{
   double	LBOX;
   long		i, nparts;
   molstruct	* particle;
   FILE 	* fconf;

   argc		=	2;
   nparts	=	0;
/*
   if ( (particle=calloc(20000, sizeof(molstruct)))==NULL)
      printf("Memory Allocation Fail.\n");
   
   if ( (fconf=fopen(argv[1],"r"))==NULL) 
      printf("Open configuration file failed!\n");
   else {
      fread(&LBOX, sizeof(double), 1, fconf);
      i	=	0;
      fread(&(particle[i].p), sizeof(vector), 1, fconf);
      while (!feof(fconf)) {
         i	++;
         fread(&(particle[i].p), sizeof(vector), 1, fconf);
      }
      nparts	=	i;		//read in nparts particles
      fclose(fconf);
   }
*/
   printf("%d\n", nparts);

   realloc(particle, nparts*sizeof(molstruct));		//reallocate the memory for particles
/*
   New_Vlist();
   New_Clist();
   Find_Nuclei();
   Print_Nuclei();
*/
   return 0;
}

