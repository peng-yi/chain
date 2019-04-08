/*
    program:    histogram.c
    author:     Peng Yi at MIT
    date:       Aug 8, 2007
    purpose:    histogram file analysis program
*/

#define __MAIN_PROGRAM
#include "header.h"

int main(int argc, char * argv[])
{
   long		i, j;
   FILE		*fhist, *fp, *foutput;
   double	*newpQ, *newpN;
   double	*GN, *GQ;
   char		comments[1024];

   Read_Setup(argv);

   if (dynvar==1) {
      newpN	=	calloc(NMAXbins, sizeof(double));
      GN	=	calloc(NMAXbins, sizeof(double));
      for (j=0; j<NMAXbins; j++)
         newpN[j]	=	0;
   }
   else if (dynvar==2) {
      newpQ	=	calloc(Qlbins, sizeof(double));
      GQ	=	calloc(Qlbins, sizeof(double));
      for (j=0; j<Qlbins; j++)
         newpQ[j]	=	0;
   }

   if (! (fhist=fopen("histogram.dat","r"))) {
      printf("histogram.dat cannot open.\n");
   }
   else {
      if (! (fp=fopen("p.dat", "w"))){
	 printf("output file cannot open.\n");
      }
      else {
   
         printf("start scanning\n");

         for (i=0; i<NCYCLE*(TRIALRUN+PROD); i++) {
            fscanf(fhist, "%ld%lf%lf%ld%ld%ld", &counter, &VSYSTEM.tot, &Ql, &MAXSIZE, &Xtal, &Nnucl);
	    fscanf(fhist, "%lf%lf%lf%lf", &LBOX, &Rv, &DRMAX, &DLMAX);
	    fgets(comments, sizeof(comments), fhist);

            if (i>=NCYCLE*TRIALRUN) { 
               if (dynvar==1)
                  newpN[(int) ((MAXSIZE - NMAXmiddle + NMAXbins/2.0 * NMAXbinsize)/NMAXbinsize)] +=	1.0;
               else if (dynvar==2) 
                  newpQ[(int) ((Ql - Qlmiddle + Qlbins/2.0 * Qlbinsize)/Qlbinsize)] += exp(-etaQl(Ql));	
            }
         }
         
         if (dynvar==1) {
            for (i=0; i<NMAXbins; i++) {
	       GN[i]	=	-log(newpN[i]);
               GN[i]	+=	etaNMAX( i*NMAXbinsize + NMAXmiddle -(int)(NMAXbins/2.0)*NMAXbinsize );
            }
         }
         else if (dynvar==2) {
            for (i=0; i<Qlbins; i++) {
               GQ[i]	=	-log(newpQ[i]);
            }
         }
      }

      if (dynvar==1) {
         for (i=0; i<NMAXbins; i++) {
            fprintf(fp, "%d\t%f\t%f\n", i, newpN[i], GN[i]);
         }
      }
      else if (dynvar==2) {
         for (i=0; i<Qlbins; i++) {
            fprintf(fp, "%d\t%f\t%f\n", i, newpQ[i], GQ[i]);
         }
      }
      fclose(fp);
   }
   fclose(fhist);    
   return	0;
}


