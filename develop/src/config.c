/*
    program:    config.c (Analysis program)
    author:     Peng Yi at MIT
    date:       Jul 23, 2007
    purpose:    Transfer mechine-readable configuration file(s) to human-readable one(s) 
		and calculate cluster information.
*/

#define __MAIN_PROGRAM
#include "header.h"

int main(int argc, char * argv[])	//read in three parameters: start file, end file, and particle number
{

   FILE 	*fout;
   long		i, start, end, file;
   double	scale, L;
   vector	p;
   char		color;
   int		a, b, c;
   long		secondsize;

   if (argc!=5) {
      printf("Usage: config startfile endfile partilenumber PBC\n");
      return;
   }

   start	=	atol(argv[1]);
   end		=	atol(argv[2]);
   NPARTS	=	atol(argv[3]);
   PBC		=	atol(argv[4]);

   part	=	(molstruct *) calloc(NPARTS, sizeof(molstruct));
   for (i=0; i<NPARTS; i++) {						//allocate space for connected neighbor lists
      part[i].clist	=	calloc(MAXCONNEIGH, sizeof(long));	//15 = 4pi/3 * Rb^3 * density
   }
   sizeofnucl	=	calloc(NPARTS+1, sizeof(long));			//allocate memory for cluster info.
   sizedist	=	calloc(NPARTS+1, sizeof(long));

   for (file=start; file<=end; file++) {
      sprintf(tempname,"%d",file);
      strncpy(vsulname,"00000000",8-strlen(tempname));
      vsulname[8-strlen(tempname)]='\0';
      strcat(vsulname,tempname);
      printf("%s\n",vsulname);

      if (!Read_Conf(vsulname))
	 return;

      scale	=	1.0;
      Rc	=	Rc0 * scale;
      Rb	=	Rb0 * scale;
      Rv	=	Rv0 * scale;

      printf("Rc=%f\tRb=%f\tRv=%f\n", Rc, Rb, Rv);

#ifdef CELL_LIST
      New_CL();
#endif
      New_Qlm(l_of_Ylm);
      New_Clist();
      Find_Nuclei();

      fout	=	fopen(strcat(vsulname,"img.pdb"), "w");
      fprintf(fout, "#MAXSIZE=%d\t# of biggest nuclei=%d\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);

      L	=	LBOX;      
      for (i=0; i<NPARTS; i++) {
	   if (part[i].nuclid!=-1) {				//crystal-like particle
              if (sizeofnucl[part[i].nuclid] == MAXSIZE) {	//the biggest nucleus (nuclei), blue color	
		 color	=	'N';
	      }
/*	   else {						//green color
		 color	=	'F';
	      }
	   }
	   else {						//liquid-like particle, red color
	      color	=	'O';
	   }
*/
	   if (PBC==1){		//cubic pbc images
	      for (a=-1; a<=1; a++) {
	         for (b=-1; b<=1; b++) {
	            for (c=-1; c<=1; c++) {
		       p.x	=	part[i].p.x + a * L;
		       p.y	=	part[i].p.y + b * L;
		       p.z	=	part[i].p.z + c * L;
                       fprintf(fout,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
		    }
	         }
	      }	
	   }
	   else if (PBC==2) {
              p	=	part[i].p;
              fprintf(fout,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);

	      for (a=-1; a<=1; a+=2) {			//octahedron images
		 for (b=-1; b<=1; b+=2) {
		    for (c=-1; c<=1; c+=2) {
		       p.x	=	part[i].p.x + 0.5 * a * L;
		       p.y	=	part[i].p.y + 0.5 * b * L;
		       p.z	=	part[i].p.z + 0.5 * c * L;
                       fprintf(fout,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
		    }
                 }
              }
	      for (a=-1; a<=1; a+=2) {			//cubic images
		 p.x	=	part[i].p.x + a * L;
		 p.y	=	part[i].p.y;
		 p.z	=	part[i].p.z;
                 fprintf(fout,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
	      }	
	      for (b=-1; b<=1; b+=2) {
		 p.x	=	part[i].p.x;
		 p.y	=	part[i].p.y + b * L;
		 p.z	=	part[i].p.z;
                 fprintf(fout,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
	      }	
	      for (c=-1; c<=1; c+=2) {
		 p.x	=	part[i].p.x;
		 p.y	=	part[i].p.y;
		 p.z	=	part[i].p.z + c * L;
                 fprintf(fout,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
	      }	
	   }
         }
      }
      fclose(fout);

      fout	=	fopen(strcat(vsulname,".pdb"), "w");
      fprintf(fout, "#MAXSIZE=%d\t# of biggest nuclei=%d\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);

      L	=	LBOX;      
      for (i=0; i<NPARTS; i++) {
	   if (part[i].nuclid!=-1) {				//crystal-like particle
              if (sizeofnucl[part[i].nuclid] == MAXSIZE) {	//the biggest nucleus (nuclei), blue color	
		 color	=	'N';
	      }
              p	=	part[i].p;
              fprintf(fout,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
           }
      } 
      fclose(fout);

      secondsize	=	0;
      if (sizedist[MAXSIZE]>1)
	 secondsize	=	MAXSIZE;
      else {
         for (i=MAXSIZE-1; i>=1; i--) {
	    if (sizedist[i] != 0) {
	       secondsize	=	i;
	       break;
	    }
         }
      }
      printf("dynvar\t%f\t%d\t%d\t%d\t%d\n", Ql, MAXSIZE, secondsize, Xtal, Nnucl);
      for (i=1; i<=MAXSIZE; i++) {
	 printf("%d\t",sizedist[i]);
      }
      printf("\n");
   }
   return;
}
