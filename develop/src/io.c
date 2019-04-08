/*
    program:	io.c
    author:     Peng Yi at MIT
    date:       October 24, 2006
    purpose:    Integral module for input/output handling
*/
#define __IO_MODULE
#include "io.h"

char		title1[256], title2[256], comments[256];

void Print_Header(FILE * filehandle)
{
  curtime	=	time(NULL);
  loctime	=	localtime(&curtime);
  fprintf(filehandle, "%s\n", asctime(loctime));
  return;
}

int samestr(char *string1, char *string2)	// compare two strings, insensitive to case
{						// return 1 if same
   int	i;
   if (strlen(string1) != strlen(string2))
      return	0;
   for (i=0; i<strlen(string1); i++)
      if (toupper(string1[i]) != toupper(string2[i]))
         return	0;
   return 1;
}

void Exit(char *module, char *procedure, char *error)
{
  fprintf(foutput, "%s::%s: %s\n", module, procedure, error);
#ifdef MPI
  fprintf(foutput, "[%ld] %s::%s: %s\n", CommRank, module, procedure, error);
  if (CommSize>1) MPI_Abort(MPI_COMM_WORLD, -1);
  else
#endif
    fprintf(foutput, "Program aborted\n\n");
  fclose(foutput);
//  close(stderr);
  exit(-1);
}

void GetLVar(FILE *fPtr, long maxn, long *L)		// read-in long type variables
{							// only work for maxn==1 so far
   long		i;
   char		c[255], dummy[255];
   for (i=0; i<maxn; i++)
      fscanf(fPtr, "%s%ld", c, L+i);
   fgets(dummy, sizeof(dummy), fPtr);
}

void GetDVar(FILE *fPtr, long maxn, double *d)		// read-in double type variables
{							// so far only work for maxn=1
   long		i;
   char		c[255], a[255], dummy[255];
 
   if (maxn == 1)
      fscanf(fPtr, "%s%lf", c, d);
   else {
      fscanf(fPtr, "%s%s", c, a);
      for (i=0; i<maxn; i++)
         fscanf(fPtr, "%lf", d+i);
   }
   fgets(dummy, sizeof(dummy), fPtr);
}

void GetSetup(char * argv[])
{
   FILE		*fPtr;
   long		i;
   char		c[256], name[256], a[256];
   char 	DYNVAR[256];

//   if ( (fPtr=fopen(argv[1], "r"))==NULL )
   if ( (fPtr=fopen("setup","r"))==NULL )
      Exit("io", "GetSetup", "Open setup file failed");
   else {
      fgets(title1, sizeof(title1), fPtr);
      fgets(title2, sizeof(title2), fPtr);
      fgets(comments, sizeof(comments), fPtr);

      fscanf(fPtr, "%s", name);
      while(strcmp(name, "END_SETUP")) {

	 /* SYSTEM PROPERTY */

	 if (!strcmp(name, "moltype")) {		
	    fscanf(fPtr, "%s%s", c, moltype);	   	fgets(comments,sizeof(comments),fPtr);
	 }
         else if (!strcmp(name, "ConvertUnits"))	GetLVar(fPtr, 1, &ConvertUnits);
	 else if (!strcmp(name, "NPARTS"))		GetLVar(fPtr, 1, &NPARTS);
         else if (!strcmp(name, "NMOLS"))		GetLVar(fPtr, 1, &NMOLS);
         else if (!strcmp(name, "NSITES"))		GetLVar(fPtr, 1, &NSITES);
         else if (!strcmp(name, "NBOX")) {		GetLVar(fPtr, 1, &NBOX);
	    NSYSTEMS=NBOX; 				// temporary
	 }
         else if (samestr(name, "Dpoly"))		GetDVar(fPtr, 1, &Dpoly);
         else if (!strcmp(name, "Rho"))			GetDVar(fPtr, 1, &Rho);
         else if (!strcmp(name, "kT")) 			GetDVar(fPtr, 1, &kT);
         else if (!strcmp(name, "P")) 			GetDVar(fPtr, 1, &P);

         else if (!strcmp(name, "V_LJ"))		GetLVar(fPtr, 1, &V_LJ); 
         else if (!strcmp(name, "V_HS")) 		GetLVar(fPtr, 1, &V_HS);
         else if (!strcmp(name, "V_RPL")) 		GetLVar(fPtr, 1, &V_RPL);
         else if (!strcmp(name, "V_LJSHIFT")) 		GetLVar(fPtr, 1, &V_LJSHIFT);
         else if (!strcmp(name, "V_LJLRC")) 		GetLVar(fPtr, 1, &V_LJLRC);
         else if (!strcmp(name, "V_STRETCH")) 		GetLVar(fPtr, 1, &V_STRETCH);
         else if (!strcmp(name, "V_BENDING")) 		GetLVar(fPtr, 1, &V_BENDING);
         else if (!strcmp(name, "V_TORSION")) 		GetLVar(fPtr, 1, &V_TORSION);
         else if (!strcmp(name, "V_VIRIAL")) 		GetLVar(fPtr, 1, &V_VIRIAL);

         else if (!strcmp(name, "Rc"))			GetDVar(fPtr, 1, &Rc); 
         else if (!strcmp(name, "Rclow")) 		GetDVar(fPtr, 1, &Rclow);
         else if (!strcmp(name, "Rv"))			GetDVar(fPtr, 1, &Rv);
         else if (!strcmp(name, "Rb")) 			GetDVar(fPtr, 1, &Rb);
         else if (!strcmp(name, "Rp")) 			GetDVar(fPtr, 1, &Rp);
         else if (!strcmp(name, "Rconn")) 		GetDVar(fPtr, 1, &Rconn);
         else if (!strcmp(name, "critconnect")) 	GetLVar(fPtr, 1, &critconnect);
         else if (!strcmp(name, "SCALECUTOFF")) 	GetLVar(fPtr, 1, &V_SCALECUTOFF);
	
	 /* SIMULATION SPECIFICS */
	 
         else if (!strcmp(name, "INITCONF")) { 
	    fscanf(fPtr, "%s%s", c, INITCONF);	    	fgets(comments,sizeof(comments),fPtr);
	 }
         else if (!strcmp(name, "Stage"))		GetLVar(fPtr, 1, &Stage);
         else if (!strcmp(name, "PBC")) 		GetLVar(fPtr, 1, &PBC);
         else if (!strcmp(name, "DRMAX")) 		GetDVar(fPtr, 1, &DRMAX);
         else if (!strcmp(name, "DLMAX")) {		GetDVar(fPtr, 1, &DLMAX);
            GDLMAX	=	DLMAX;			// Gibbs volume change step
	 }
         else if (!strcmp(name, "DAMAX")) 		GetDVar(fPtr, 1, &DAMAX);
         else if (!strcmp(name, "E_NVT")) 		GetLVar(fPtr, 1, &E_NVT);
         else if (!strcmp(name, "E_NPT")) 		GetLVar(fPtr, 1, &E_NPT);
         else if (!strcmp(name, "E_GIBBS")) 		GetLVar(fPtr, 1, &E_GIBBS);
         else if (!strcmp(name, "NDISPLACE")) 		GetLVar(fPtr, 1, &NDISPLACE);
         else if (!strcmp(name, "NREPTATION")) 		GetLVar(fPtr, 1, &NREPTATION);
         else if (!strcmp(name, "NENDROT")) 		GetLVar(fPtr, 1, &NENDROT);
         else if (!strcmp(name, "NVOLCHANGE")) 		GetLVar(fPtr, 1, &NVOLCHANGE);
         else if (!strcmp(name, "NGIBBSVOL")) 		GetLVar(fPtr, 1, &NGIBBSVOL);
         else if (!strcmp(name, "NSWAP")) 		GetLVar(fPtr, 1, &NSWAP);
         else if (!strcmp(name, "NCBMC")) 		GetLVar(fPtr, 1, &NCBMC);
	 else if (samestr(name, "NENDBR"))		GetLVar(fPtr, 1, &NENDBR);
	 else if (samestr(name, "NREBR"))		GetLVar(fPtr, 1, &NREBR);
	 else if (samestr(name, "NFLIP"))		GetLVar(fPtr, 1, &NFLIP);
	 else if (samestr(name, "NDB"))			GetLVar(fPtr, 1, &NDB);
	 else if (samestr(name, "NIDR"))		GetLVar(fPtr, 1, &NIDR);

         else if (!strcmp(name, "NTRIALCONF")) 		GetLVar(fPtr, 1, &NTRIALCONF);
         else if (!strcmp(name, "NTRIALFIRSTBEAD")) 	GetLVar(fPtr, 1, &NTRIALFIRSTBEAD);
         else if (!strcmp(name, "SUCC_DISP")) 		GetDVar(fPtr, 1, &SUCC_DISP);
         else if (!strcmp(name, "SUCC_VOL")) 		GetDVar(fPtr, 1, &SUCC_VOL);
         else if (!strcmp(name, "Nequil")) 		GetLVar(fPtr, 1, &Nequil);
         else if (!strcmp(name, "Nprod")) 		GetLVar(fPtr, 1, &Nprod);
         else if (!strcmp(name, "SEQUENCE")) 		GetLVar(fPtr, 1, &SEQUENCE);
         else if (!strcmp(name, "TRIALRUN")) 		GetLVar(fPtr, 1, &TRIALRUN);
         else if (!strcmp(name, "PROD")) 		GetLVar(fPtr, 1, &PROD);
         else if (!strcmp(name, "ITAPE")) 		GetLVar(fPtr, 1, &ITAPE);
         else if (!strcmp(name, "ICONF")) 		GetLVar(fPtr, 1, &ICONF);
         else if (!strcmp(name, "NBLOCKS")) 		GetLVar(fPtr, 1, &NBLOCKS);
         else if (!strcmp(name, "NGSAMPLE")) 		GetLVar(fPtr, 1, &NGSAMPLE);
         else if (!strcmp(name, "IRADIAL")) 		GetLVar(fPtr, 1, &IRADIAL);
         else if (!strcmp(name, "NGRBINS")) 		GetLVar(fPtr, 1, &NGRBINS);
         else if (!strcmp(name, "Alpha")) 		GetDVar(fPtr, 1, &Alpha);
         else if (!strcmp(name, "CRIT")) 		GetDVar(fPtr, 1, &CRIT);
 
	 /* SAMPLING PARAMETERS */
         
         else if (!strcmp(name, "D_DENSITY")) 		GetLVar(fPtr, 1, &D_DENSITY);
         else if (!strcmp(name, "D_ENERGY")) 		GetLVar(fPtr, 1, &D_ENERGY);
         else if (!strcmp(name, "D_PRESSURE"))		GetLVar(fPtr, 1, &D_PRESSURE);
         else if (!strcmp(name, "D_DRIFT"))		GetLVar(fPtr, 1, &D_DRIFT);
         else if (!strcmp(name, "D_TORSION"))		GetLVar(fPtr, 1, &D_TORSION);
         else if (!strcmp(name, "D_BONDA"))		GetLVar(fPtr, 1, &D_BONDA);
         else if (!strcmp(name, "D_BONDL"))		GetLVar(fPtr, 1, &D_BONDL);
         else if (!strcmp(name, "D_RADIAL"))		GetLVar(fPtr, 1, &D_RADIAL);
         else if (!strcmp(name, "D_LOCALP2"))		GetLVar(fPtr, 1, &D_LOCALP2);
         else if (!strcmp(name, "D_XTALSIZE"))		GetLVar(fPtr, 1, &D_XTALSIZE);

	 else if (!strcmp(name, "DYNVAR")) {
	    fscanf(fPtr, "%s%s", c, DYNVAR);	   	fgets(comments,sizeof(comments),fPtr);
	 }
         else if (!strcmp(name, "kP2"))			GetDVar(fPtr, 1, &kP2);
         else if (!strcmp(name, "P2middle"))		GetDVar(fPtr, 1, &P2middle);
         else if (!strcmp(name, "NMAXmiddle"))		GetLVar(fPtr, 1, &NMAXmiddle);
         else if (!strcmp(name, "NMAXbinsize"))		GetLVar(fPtr, 1, &NMAXbinsize);
         else if (!strcmp(name, "NMAXbins")) 		GetLVar(fPtr, 1, &NMAXbins);
         else if (!strcmp(name, "kN")) 			GetDVar(fPtr, 1, &kN);
         else if (!strcmp(name, "Qlmiddle")) 		GetDVar(fPtr, 1, &Qlmiddle);
         else if (!strcmp(name, "Qlbinsize")) 		GetDVar(fPtr, 1, &Qlbinsize);
         else if (!strcmp(name, "Qlbins")) 		GetLVar(fPtr, 1, &Qlbins);
         else if (!strcmp(name, "kQ")) 			GetDVar(fPtr, 1, &kQ);
         else if (!strcmp(name, "critqlproduct"))	GetDVar(fPtr, 1, &critqlproduct);
         else if (!strcmp(name, "critconnect")) 	GetLVar(fPtr, 1, &critconnect);
         else if (!strcmp(name, "critp2")) 		GetDVar(fPtr, 1, &critp2);
         else
            fgets(comments, sizeof(comments),fPtr);
	 
         fscanf(fPtr, "%s", name);
      }
      fgets(comments,sizeof(comments),fPtr);

      /* Initialize memory, some of which, e.g., eta[] is necessayr for next step read-in */
    /* 
      sizeofnucl	=	calloc(NPARTS+1, sizeof(long));
      sizedist		=	calloc(NPARTS+1, sizeof(long));

      sizeofnucl2	=	calloc(NPARTS+1, sizeof(long));
      sizedist2		=	calloc(NPARTS+1, sizeof(long));

      if (sizeofnucl==NULL || sizedist==NULL)
         printf("size distribution calloc failed.\n");
    */
      if (!strcmp(DYNVAR, "NMAX")) {		//NMAX as dynamic variable
	 dynvar	=	1;
/*
         p	=	calloc(NMAXbins, sizeof(double));
         eta	=	calloc(NMAXbins, sizeof(double));
         pt	=	calloc(NMAXbins, sizeof(double));
         G 	=	calloc(NMAXbins, sizeof(double));
         if (p==NULL || eta==NULL || pt==NULL || G==NULL)  printf("p and eta calloc failed.\n");
*/
      }
      else if (!strcmp(DYNVAR, "Q6")) {		//Ql as dynamic variable
	 dynvar	=	2;	
/*
         pQ	=	calloc(Qlbins, sizeof(double));
         etaQ	=	calloc(Qlbins, sizeof(double));
         ptQ	=	calloc(Qlbins, sizeof(double));
         GQ 	=	calloc(Qlbins, sizeof(double));
         if (pQ==NULL || etaQ==NULL || ptQ==NULL || GQ==NULL)  printf("pQ and etaQ calloc failed.\n");
*/
      }
      else if (!strcmp(DYNVAR, "P2M")) {
         dynvar	=	3;
      }
      else if (!strcmp(DYNVAR, "P2")) {		// P2 as umbrella sampling order parameter
         dynvar	=	4;
      }
      else if (!strcmp(DYNVAR, "RHO")) {	// density as umbrella sampling order parameter
         dynvar	=	5;
      }
      else if (samestr(DYNVAR, "NMAXp2")) {	// nmax defined based on local p2
         dynvar	=	6;
      }

      /* Read in the pre-set eta */
/*
      rewind(fPtr);     //return to the beginning of setup file

      fscanf(fPtr, "%s", name);
      while(strcmp(name, "END_SETUP")) {
         if (!strcmp(name, "etaN")) {           // two if's cannot be combined, especially cannot write
            if (dynvar==1) {                    // if (dynvar==1 && !strcmp(name,"etaN")
               fscanf(fPtr, "%s%s", c, a);
               for (i=0; i<NMAXbins; i++)
                  fscanf(fPtr, "%lf", &eta[i]);
            }
            fgets(comments, sizeof(comments), fPtr);
         }
         else if(!strcmp(name, "etaQ")) {
            if (dynvar==2) {
               fscanf(fPtr, "%s%s", c, a);
               for (i=0; i<Qlbins; i++)
                  fscanf(fPtr, "%lf", &etaQ[i]);
            }
            fgets(comments, sizeof(comments), fPtr);
         }
         else
            fgets(comments, sizeof(comments), fPtr);

         fscanf(fPtr, "%s", name);
      }
      fgets(comments, sizeof(comments), fPtr);
*/
      /* READ IN FORCEFIELD MODEL */	
	
      rewind(fPtr);			// return to the beginning of the setup file
      fscanf(fPtr, "%s", name);
      while(strcmp(name, "END_SETUP")) {
	 if (!strcmp(name, "NTYPES")) 			GetLVar(fPtr, 1, &NTYPES);
         else if (!strcmp(name, "DLJ")) 		GetLVar(fPtr, 1, &DLJ);
         else if (!strcmp(name, "TYPEMASS")) {  
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++) 
               fscanf(fPtr, "%lf", &(type[i].M));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "SIGMA")) {  
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].SIGMA));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "EPSILON")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].EPSILON));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "KSTRETCH")) { 
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].KSTRETCH));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "LSTRETCH")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].LSTRETCH));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "FIXBONDA")) 			GetLVar(fPtr, 1, &FIXBONDA);
         else if (!strcmp(name, "FIXBONDL")) 			GetLVar(fPtr, 1, &FIXBONDL);
         else if (!strcmp(name, "KBENDING")) { 
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].KBENDING));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "THETA")) { 
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].THETA));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "HS")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].HS));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (samestr(name, "TORTYPE")) {
	    fscanf(fPtr, "%s%s", c, TORTYPE);
	    fgets(comments,sizeof(comments),fPtr);
	 }
         else if (!strcmp(name, "TORSION0")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[0]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION1")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[1]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION2")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[2]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION3")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[3]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION4")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[4]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION5")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[5]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else {
            fgets(comments, sizeof(comments), fPtr);
         }
         fscanf(fPtr, "%s", name);
      }
      fgets(comments, sizeof(comments), fPtr);
      fclose(fPtr);     //read in finished
   }
   return;
}


void PrintSetup()		//print out setup variables
{
   long		i, j;
   char		name[256];

   fprintf(foutput, "%s\n", title1);
   fprintf(foutput, "%s\n", title2);

   fprintf(foutput, "\n/* Basic system property. */\n\n");

      fprintf(foutput, "Molecule Type\t=\t%s\n", moltype);
      fprintf(foutput, "ConvertUnits\t=\t%s\n", (ConvertUnits ? "Yes":"No"));
      fprintf(foutput, "Rho 	\t=\t%f\n", Rho);
      fprintf(foutput, "Dpoly   \t=\t%f\n", Dpoly);
      fprintf(foutput, "kT	\t=\t%f\n", kT);
      fprintf(foutput, "P 	\t=\t%f\n", P);
      fprintf(foutput, "NBOX	\t=\t%ld\n", NBOX);
      fprintf(foutput, "NMOLS    \t=\t%ld\n", NMOLS);
      fprintf(foutput, "NSITES    \t=\t%ld\n", NSITES);
      fprintf(foutput, "NPARTS    \t=\t%ld\n", NPARTS);

   fprintf(foutput, "\n/* Basic simulation setup */\n\n");
      
      fprintf(foutput, "INITCONF  \t=\t%s\n", INITCONF);
      fprintf(foutput, "Stage 	\t=\t%ld\n", Stage);
      fprintf(foutput, "PBC 	\t=\t%ld\n", PBC);
#ifdef CELL_LIST
      fprintf(foutput, "CELL_LIST \t=\t%s\n", "Yes");
#else
      fprintf(foutput, "CELL_LIST \t=\t%s\n", "No");
#endif
      fprintf(foutput, "Nequil 	\t=\t%ld\n", Nequil);
      fprintf(foutput, "Nprod 	\t=\t%ld\n", Nprod);
      fprintf(foutput, "SEQUENCE  \t=\t%ld\n", SEQUENCE);
      fprintf(foutput, "TRIALRUN  \t=\t%ld\n", TRIALRUN);
      fprintf(foutput, "PROD 	\t=\t%ld\n", PROD);
      fprintf(foutput, "ITAPE 	\t=\t%ld\n", ITAPE);
      fprintf(foutput, "ICONF 	\t=\t%ld\n", ICONF);
      fprintf(foutput, "NBLOCKS   \t=\t%ld\n", NBLOCKS);
      fprintf(foutput, "NGSAMPLE  \t=\t%ld\n", NGSAMPLE);
      fprintf(foutput, "IRADIAL   \t=\t%ld\n", IRADIAL);
      fprintf(foutput, "NGRBINS   \t=\t%ld\n", NGRBINS);
      fprintf(foutput, "CRIT 	\t=\t%f\n", CRIT);
      fprintf(foutput, "Alpha 	\t=\t%f\n", Alpha);

   fprintf(foutput, "\n/* Forcefield model */\n\n");

      fprintf(foutput, "NTYPES	\t=\t%ld\n", NTYPES);
      fprintf(foutput, "DLJ	\t=\t%ld\n", DLJ);
      fprintf(foutput, "TYPEMASS\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(foutput, "%f\t", type[i].M);
      fprintf(foutput, "}\n");
      fprintf(foutput, "SIGMA	\t=\t{ ");
      for (i=0; i<NTYPES; i++)	        fprintf(foutput, "%f\t", type[i].SIGMA);
      fprintf(foutput, "}\n");
      fprintf(foutput, "EPSILON	\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(foutput, "%f\t", type[i].EPSILON);
      fprintf(foutput, "}\n");
      fprintf(foutput, "KSTRETCH\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(foutput, "%f\t", type[i].KSTRETCH);
      fprintf(foutput, "}\n");
      fprintf(foutput, "LSTRETCH\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(foutput, "%f\t", type[i].LSTRETCH);
      fprintf(foutput, "}\n");
      fprintf(foutput, "KBENDING\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(foutput, "%f\t", type[i].KBENDING);
      fprintf(foutput, "}\n");
      fprintf(foutput, "THETA	\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(foutput, "%f\t", type[i].THETA);
      fprintf(foutput, "}\n");
      fprintf(foutput, "HS	\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(foutput, "%f\t", type[i].HS);
      fprintf(foutput, "}\n");

      fprintf(foutput, "TORTYPE	\t=\t%s\n",TORTYPE);

      for (j=0; j<6; j++) {
         sprintf(name, "TORSION%c", '0'+j);
         fprintf(foutput, "%s\t=\t{ ", name);
         for (i=0; i<NTYPES; i++)	fprintf(foutput, "%f\t", type[i].TORSION[j]);
         fprintf(foutput, "}\n");
      }
      fprintf(foutput, "V_LJ	\t=\t%s\n", (V_LJ ? "Yes" : "No"));
      fprintf(foutput, "V_LJSHIFT \t=\t%s\n", (V_LJSHIFT ? "Yes" : "No"));
      fprintf(foutput, "V_LJLRC	\t=\t%s\n", (V_LJLRC ? "Yes" : "No"));
      fprintf(foutput, "V_HS	\t=\t%s\n", (V_HS ? "Yes" : "No"));
      fprintf(foutput, "V_RPL	\t=\t%s\n", (V_RPL ? "Yes" : "No"));
      fprintf(foutput, "V_STRETCH \t=\t%s\n", (V_STRETCH ? "Yes" : "No"));
      fprintf(foutput, "V_BENDING \t=\t%s\n", (V_BENDING ? "Yes" : "No"));
      fprintf(foutput, "V_TORSION \t=\t%s\n", (V_TORSION ? "Yes" : "No"));
      fprintf(foutput, "FIXBONDL  \t=\t%s\n", (FIXBONDL ? "Yes" : "No"));
      fprintf(foutput, "FIXBONDA  \t=\t%s\n", (FIXBONDA ? "Yes" : "No"));      

      fprintf(foutput, "V_VIRIAL  \t=\t%s\n", (V_VIRIAL ? "Yes" : "No"));
      fprintf(foutput, "Rc	\t=\t%f\n", Rc);
      fprintf(foutput, "Rclow	\t=\t%f\n", Rclow);
      fprintf(foutput, "Rv	\t=\t%f\n", Rv);
      fprintf(foutput, "Rb	\t=\t%f\n", Rb);
      fprintf(foutput, "Rp	\t=\t%f\n", Rp);
      fprintf(foutput, "Rconn	\t=\t%f\n", Rconn);
      fprintf(foutput, "critconn\t=\t%ld\n", critconnect);
      fprintf(foutput, "SCALECUTOFF \t=\t%s\n", (V_SCALECUTOFF ? "Yes" : "No"));

   fprintf(foutput, "\n/* Ensemble setup */\n\n");
      
      fprintf(foutput, "E_NVT	\t=\t%s\n", (E_NVT ? "Yes" : "No"));
      fprintf(foutput, "E_NPT	\t=\t%s\n", (E_NPT ? "Yes" : "No"));
      fprintf(foutput, "E_GIBBS \t=\t%s\n", (E_GIBBS ? "Yes" : "No"));
      fprintf(foutput, "NDISPLACE \t=\t%ld\n", NDISPLACE);
      fprintf(foutput, "NSWAP	\t=\t%ld\n", NSWAP);
      fprintf(foutput, "NREPTATION\t=\t%ld\n", NREPTATION);
      fprintf(foutput, "NENDROT   \t=\t%ld\n", NENDROT);
      fprintf(foutput, "NCBMC	\t=\t%ld\n", NCBMC);
      fprintf(foutput, "NENDBR	\t=\t%ld\n", NENDBR);
      fprintf(foutput, "NREBR	\t=\t%ld\n", NREBR);
      fprintf(foutput, "NFLIP	\t=\t%ld\n", NFLIP);
      fprintf(foutput, "NDB  	\t=\t%ld\n", NDB);
      fprintf(foutput, "NIDR 	\t=\t%ld\n", NIDR);
      fprintf(foutput, "NVOLCHANGE\t=\t%ld\n", NVOLCHANGE);
      fprintf(foutput, "NGIBBSVOL \t=\t%ld\n", NGIBBSVOL);
      fprintf(foutput, "NTRIALCONF\t=\t%ld\n", NTRIALCONF);
      fprintf(foutput, "NTRIALFIRSTBEAD\t=\t%ld\n", NTRIALFIRSTBEAD);
      fprintf(foutput, "SUCC_DISP \t=\t%f\n", SUCC_DISP);
      fprintf(foutput, "SUCC_VOL \t=\t%f\n", SUCC_VOL);
      fprintf(foutput, "DRMAX 	\t=\t%f\n", DRMAX);
      fprintf(foutput, "DLMAX 	\t=\t%f\n", DLMAX);
      fprintf(foutput, "DAMAX 	\t=\t%f\n", DAMAX);
      
   fprintf(foutput, "\n/* Sampling parameters */\n\n");
      
      fprintf(foutput, "D_DENSITY\t=\t%s\n", (D_DENSITY ? "Yes" : "No"));
      fprintf(foutput, "D_ENERGY \t=\t%s\n", (D_ENERGY ? "Yes" : "No"));
      fprintf(foutput, "D_PRESSURE\t=\t%s\n", (D_PRESSURE ? "Yes" : "No"));
      fprintf(foutput, "D_DRIFT	 \t=\t%s\n", (D_DRIFT ? "Yes" : "No"));
      fprintf(foutput, "D_TORSION\t=\t%s\n", (D_TORSION ? "Yes" : "No"));
      fprintf(foutput, "D_BONDA	 \t=\t%s\n", (D_BONDA ? "Yes" : "No"));
      fprintf(foutput, "D_BONDL	 \t=\t%s\n", (D_BONDL ? "Yes" : "No"));
      fprintf(foutput, "D_RADIAL \t=\t%s\n", (D_RADIAL ? "Yes" : "No"));
      fprintf(foutput, "D_LOCALP2\t=\t%s\n", (D_LOCALP2 ? "Yes" : "No"));
      fprintf(foutput, "D_XTALSIZE\t=\t%s\n", (D_XTALSIZE ? "Yes" : "No"));
     
      switch(dynvar) {
         case	1:	fprintf(foutput, "DYNVAR  \t=\t%s\n", "NMAX");	break;
	 case	2:	fprintf(foutput, "DYNVAR  \t=\t%s\n", "Q6");	break;
	 case	3:	fprintf(foutput, "DYNVAR  \t=\t%s\n", "P2M");	break;
	 case	4:	fprintf(foutput, "DYNVAR  \t=\t%s\n", "P2");	break;
	 case	5:	fprintf(foutput, "DYNVAR  \t=\t%s\n", "RHO");	break;
	 case	6:	fprintf(foutput, "DYNVAR  \t=\t%s\n", "NMAXp2");break;
         default:	fprintf(foutput, "DYNVAR  \t=\t%s\n", "ERROR");	break;
      }
      fprintf(foutput, "kP2	\t=\t%f\n", kP2);
      fprintf(foutput, "P2middle  \t=\t%f\n", P2middle);
      fprintf(foutput, "NMAXmiddle\t=\t%ld\n", NMAXmiddle);
      fprintf(foutput, "NMAXbinsize\t=\t%ld\n", NMAXbinsize);
      fprintf(foutput, "NMAXbins  \t=\t%ld\n", NMAXbins);
      fprintf(foutput, "kN	\t=\t%f\n", kN);
      fprintf(foutput, "Qlmiddle  \t=\t%f\n", Qlmiddle);
      fprintf(foutput, "Qlbinsize \t=\t%f\n", Qlbinsize);
      fprintf(foutput, "Qlbins	\t=\t%ld\n", Qlbins);
      fprintf(foutput, "kQ 	\t=\t%f\n", kQ);
      fprintf(foutput, "l_of_Ylm \t=\t%d\n", l_of_Ylm);
      fprintf(foutput, "critqlproduct \t=\t%f\n", critqlproduct);
      fprintf(foutput, "critconnect   \t=\t%ld\n", critconnect);
      fprintf(foutput, "critp2        \t=\t%f\n", critp2);
/*
      if (dynvar==1) {
         fprintf(foutput, "etaN	\t=\t{ ");
         for (i=0; i<NMAXbins; i++)
            fprintf(foutput, "%f\t", eta[i]);
         fprintf(foutput, "}\n");
      }
      if (dynvar==2) { 
         fprintf(foutput, "etaQ	\t=\t{ ");
         for (i=0; i<Qlbins; i++)
            fprintf(foutput, "%f\t", etaQ[i]);
         fprintf(foutput, "}\n");
      }
*/
      fprintf(foutput, "\n");
      for (i=0; i<80; i++)
	 fprintf(foutput, "*");
      fprintf(foutput, "\n");
      fflush(foutput);

      return;
}

//////////////////////////////////
//	Open output files	//
//////////////////////////////////
void InitFile()	
{
   long	i;

   if (!FLAG_restart) {
      fhst	=	fopen("histogram.dat", "w");
      fdump	=	fopen("dump", "w");
   }
   else {
      fhst	=	fopen("histogram.dat", "a");
      fdump	=	fopen("dump", "a");
   }
         
   if (fhst==NULL)
      Exit("io", "InitFile", "Histogram file cannot open!");
   if (fdump==NULL)
      Exit("io", "InitFile", "Dump file cannot open!");

   sprintf(file_hst, "binhist");
   frame = 0;

   if ((foutput=fopen("output","a"))==NULL)
      Exit("io", "InitFile", "Output file cannot open!");

   //***** Output command line and time for record *****//

   fprintf(foutput, "////////////////////////////////////////\n");
   fprintf(foutput, "////////////////////////////////////////\n");
   fprintf(foutput, "////////////////////////////////////////\n");

   fprintf(foutput, "\n### BEGINNING OF OUTPUT FILE\n");
   curtime	=	time(NULL);
   loctime	=	localtime(&curtime);
   fprintf(foutput, "%s\n", asctime(loctime));
   fprintf(foutput, "version: %s\n\n", VERSION);

   return;
}

void Print_Summary()			// print final summary
{
   long		i;
 
   for (i=0; i<NBOX; i++) {
      fprintf(foutput, "\n\tBOX[%ld]:\n", i);
      fprintf(foutput, "\t# of mols:\t%ld\n", NMols[i]);
      fprintf(foutput, "\t# of sites:\t%ld\n", NSites[i]);
      fprintf(foutput, "\tbox length (in system unit):\t%f\n", BOX[i].lbox);
      fprintf(foutput, "\tbox length x:\t%f\n", BOX[i].lx);
      fprintf(foutput, "\tbox length y:\t%f\n", BOX[i].ly);
      fprintf(foutput, "\tbox length z:\t%f\n", BOX[i].lz);
      fprintf(foutput, "\tbox volume:\t%f\n", BOX[i].vol);
      fprintf(foutput, "\trc:\t%f\n", BOX[i].rc);
      fprintf(foutput, "\trb:\t%f\n", BOX[i].rb);
      fprintf(foutput, "\trv:\t%f\n", BOX[i].rv);
      fprintf(foutput, "\tdrmax:\t%f\n", BOX[i].drmax);
      fprintf(foutput, "\tdlmax:\t%f\n", BOX[i].dlmax);
   }
   fprintf(foutput, "\n");

   for (i=0; i<NBOX; i++) {
      fprintf(foutput, "\n\tBOX[%ld]:\n", i);
      fprintf(foutput, "\tmove:\t%ld\t acc:\t%f\n", 
			av[i].move, ((double)av[i].acc_move)/av[i].move);
      fprintf(foutput, "\tvol: \t%ld\t acc:\t%f\n", 
			av[i].vol, ((double)av[i].acc_vol)/av[i].vol);
      fprintf(foutput, "\tcbmc: \t%ld\t acc:\t%f\n", 
			av[i].cbmc, ((double)av[i].acc_cbmc)/av[i].cbmc);
      fprintf(foutput, "\trep: \t%ld\t acc:\t%f\n", 
			av[i].rep, ((double)av[i].acc_rep)/av[i].rep);
      fprintf(foutput, "\tendrot:\t%ld\t acc:\t%f\n", 
			av[i].erot, ((double)av[i].acc_erot)/av[i].erot);
      fprintf(foutput, "\tendbr: \t%ld\t acc:\t%f\n", 
			av[i].eb, ((double)av[i].acc_eb)/av[i].eb);
      fprintf(foutput, "\trebr: \t%ld\t acc:\t%f\n", 
			av[i].re, ((double)av[i].acc_re)/av[i].re);
      fprintf(foutput, "\tflip: \t%ld\t acc:\t%f\n", 
			av[i].flip, ((double)av[i].acc_flip)/av[i].flip);
      fprintf(foutput, "\tdb: \t%ld\t acc:\t%f\n", 
			av[i].db, ((double)av[i].acc_db)/av[i].db);
      fprintf(foutput, "\tidr: \t%ld\t acc:\t%f\n", 
			av[i].idr, ((double)av[i].acc_idr)/av[i].idr);
      fprintf(foutput, "\tswap: \t%ld\t acc:\t%f\n", 
			av[i].swap, ((double)av[i].acc_swap)/av[i].swap);
      fprintf(foutput, "\trot: \t%ld\t acc:\t%f\n", 
			av[i].rot, ((double)av[i].acc_rot)/av[i].rot);
      fprintf(foutput, "\tmovemol: \t%ld\t acc:\t%f\n", 
			av[i].movemol, ((double)av[i].acc_movemol)/av[i].movemol);
      fprintf(foutput, "\tseq: \t%ld\t acc:\t%f\n", 
			av[i].seq, ((double)av[i].acc_seq)/av[i].seq);
   }
   if (NCBMC) {
      fprintf(foutput, "\tcbmc statistics:\n");
      for (i=0; i<NSITES/NMOLS; i++)
         fprintf(foutput, "\tcut point=\t%ld\t accepted:\t%ld\n", i, cbmcsucc[i]);
   }

   fprintf(foutput, "### END OF OUTPUT FILE\n");
   curtime	=	time(NULL);
   loctime	=	localtime(&curtime);
   fprintf(foutput, "%s\n", asctime(loctime));

   return;
}

//////////////////////////////////
//	Close output files	//
//////////////////////////////////
void CloseFile()
{
   fclose(foutput);
   fclose(fhst);
   fclose(fdump);
   fflush(stdin);
   return;
}

/****************************************************************************************/
/*	Print_Nuclei()									*/
/*											*/
/*	Print out nuclei information.							*/
/****************************************************************************************/

void Print_Nuclei()		//print out nuclei information
{
  long	i, size, id;

  fprintf(foutput,"Printout nuclei size information......\n");
  /*
  for (id=1; id<NPARTS; id++) {			//nuclei index starts from 1
    if (sizeofnucl[id] != 0) {
      fprintf(foutput,"Size of nucleus #%ld\t=%ld\n", id, sizeofnucl[id]);
    }
  }
  */
  /*
  for (i=0; i<NPARTS; i++) {		//print out nuclei info. of each particle
    fprintf(foutput, "particle #%ld\tbelongs to nuclei #%ld\n", i+1, part[i].nuclid);
  }
  */
  for (i=1; i<NPARTS+1; i++) {		//printout size distribution
    if (sizedist[i]!=0) {
      fprintf(foutput, "Print_Nuclei: Number of nuclei of size\t%ld\t=%ld\n", i, sizedist[i]);
    }
  }
  fprintf(foutput,"Printout nuclei size information finished.\n\n");
  
  return;
}

/*********************************************************************************************************/
#ifdef VERLET_LIST

void Print_Verlet()		//Print out Verlet list
{
  long	i, jj;
  fprintf(foutput,"Printout Verlet list for each particle......\n");
  for (i=0; i<NPARTS; i++) {
    fprintf(foutput,"Particle #%ld has %ld Verlet neighbors:\n",i+1, part[i].nverlet);  
    fprintf(foutput,"They are:\t");
    for (jj=0; jj<part[i].nverlet; jj++) {
      fprintf(foutput,"%ld\t", part[i].vlist[jj]+1);
    }
    fprintf(foutput,"\n");
  }
  fprintf(foutput,"Printout Verlet list for each particle finished.\n\n");
  return;
}

void Print_Verlet2()		//Print out affected particles in Verlet list
{
  long	i;
  fprintf(foutput,"Printout Verlet list affected by a displacement/volume change:\n");
  for (i=0; i<nverletplus; i++) {
   fprintf(foutput,"%ld\t",vlistplus[i]+1); 
  }
  fprintf(foutput,"\nPrintout Verlet list affected by a displacement/volume change finished.\n\n");
  return;
}

#endif
/*********************************************************************************************************/
/*
void Print_Clist()		//Print out the connected neighbors
{
  long	i, jj;
  fprintf(foutput,"Printout Connected neighbors if any...\n");
  for (i=0; i<NPARTS; i++) {
    if (part[i].nconnect >= 1) {
      fprintf(foutput,"Particle #%ld is connected to %ld neighbors. They are:\t", i+1, part[i].nconnect);
      for (jj=0; jj<part[i].nconnect; jj++) {
        fprintf(foutput,"%ld\t", part[i].clist[jj]+1);
      }
      fprintf(foutput,"\n");
    }
  }
  fprintf(foutput,"Printout Connected neighbors finished.\n\n");
  return;
}
*/
/*********************************************************************************************************/
/*
void Print_q()			//Print out q information
{
  long	i;
  int	m, l=l_of_Ylm;
  for (i=0; i<NPARTS; i++) {
    for (m=-l; m<=l; m++) {
      fprintf(foutput,"q(l=%ld,m=%ld,%ld)=\t%+f+\t%+f*i\n", l, m, i+1, creal(part[i].qlm[m+l]), cimag(part[i].qlm[m+l]));
    }
    fprintf(foutput,"q(l=%ld,%ld)=\t%f\n\n", l, i+1, ql(l, i));
  }
  fprintf(foutput,"Printout q finished.\n\n");
  return;
}
*/
/*********************************************************************************************************/
/*
void Print_Q()			//Print out Q information
{
  int	m, l=l_of_Ylm;
  for (m=-l; m<=l; m++) {
    fprintf(foutput, "aveQ(l=%ld, m=%ld)=\t%+f+\t%+f*i\n", l, m, creal(aveQlm(l, m)), cimag(aveQlm(l, m)));
  } 
  fprintf(foutput, "Q(l=%ld)=\t%f\n", l, CalcQl(l));
  fprintf(foutput, "Printout Q finished.\n\n");
  return;
}
*/
/*********************************************************************************************************/
/*
void Print_qproduct()		//Print out qproduct, which is used to judge connectivity
{
  long		i, jj, k;
  complex	term;

  fprintf(foutput, "Printout ql product......\n");
  for (i=0; i<NPARTS; i++) {
    for (jj=0; jj<part[i].nverlet; jj++) {
      k	=	part[i].vlist[jj];
      if (k>i) {
	term	=	qlproductSQ(l_of_Ylm, i, k);
        fprintf(foutput, "%ld -> %ld\t%+f+\t%+f*i\n", i+1, k+1, creal(term), cimag(term));
      }
    }
  }
  fprintf(foutput, "Printout ql product finished.\n\n");
  return;
}
*/

//======================================================================================//
// 	Visualize(int type)								//
//											//
// 	Print out crystal nuclei configuration in pdb format for VMD visualization.	//
//	type = 0 (all particles) or 1 (biggest crystals only)				//
//	Color: O (red), N (blue), F (Green), H (White), C (blue-green)			//
//======================================================================================//

#ifdef TEST
int Read_Conf(char * filename)
{
   long		i;
   FILE *	fconf;

   if ((fconf=fopen(filename, "r"))==NULL) {
      printf("%s not found or failed to open!\n", filename);
      return	0;
   }
   else {
//      fseek(fconf, NPARTS*sizeof(vector)+sizeof(double), SEEK_SET);
      fread(&(LBOX), sizeof(double), 1, fconf);
      for (i=0; i<NPARTS; i++) {
         fread(&(part[i].p), sizeof(vector), 1, fconf);
      }
      fclose(fconf);
      return	1;
   }
}


int Write_Conf()
{
   long		i;
   FILE *	fconf;
   char		confname[256], tempname[256];
   
   frame ++;				//index the visualization output file name
   sprintf(tempname,"%ld",frame);
   strncpy(confname,"00000000",8-strlen(tempname));
   confname[8-strlen(tempname)]='\0';
   strcat(confname,tempname);
   
   if ((fconf=fopen(confname,"w"))==NULL) {		
      printf("Visualization file cannot be opened!\n");
      return	0;
   }
   else {
      fwrite(&(LBOX), sizeof(double), 1, fconf);
      //printf("%-12.8f\n", LBOX);
      for (i=0; i<NPARTS; i++) {
         fwrite(&(part[i].p), sizeof(vector), 1, fconf);
         //printf("%-12.8f%-12.8f%-12.8f\n", part[i].p.x, part[i].p.y, part[i].p.z);
      }
      fclose(fconf);
      return	1;
   }
}


int Visualize(int vtype)
{
   long		i;	
   FILE		*fhere;
   vector	pvisual, p0, p;
   double	L = LBOX;
   char		color;
   long		a, b, c;
   char		vsulname[256], tempname[256];
	
   frame ++;				//index the visualization output file name
   sprintf(tempname,"%ld",frame);
   strncpy(vsulname,"00000000",8-strlen(tempname));
   vsulname[8-strlen(tempname)]='\0';
   strcat(vsulname,tempname);

   switch(vtype) {
      case 1:   strcat(vsulname, ".pdb");       break;
      case 0:   strcat(vsulname, "full.pdb");   break;
      case 2:   strcat(vsulname, "img.pdb");    break;
      case 3:   strcat(vsulname, "cell.pdb");   break;
      default:  break;
   }

   if ((fhere=fopen(vsulname,"w"))==NULL) {
      printf("Visualization file cannot be opened!\n");
      return 0;
   }
   
   if (vtype == 1) {
      fprintf(fhere, "#MAXSIZE=%ld\t# of biggest nuclei=%ld\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);

      for (i=0; i<NPARTS; i++) {
         if (sizeofnucl[part[i].nuclid] == MAXSIZE) {           //print out only the biggest nucleus (nuclei)   

            pvisual     =       part[i].p;

            fprintf(fhere,"ATOM%7d  O   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, i+1, pvisual.x, pvisual.y, pvisual.z, 1.0, 0.0);
         }      //print out only the biggest nucleus/nuclei
      }
   }
   else if (vtype == 0) {                                       //all particles in different colors
/*
      fprintf(fhere, "#MAXSIZE=%ld\t# of biggest nuclei=%ld\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);
      for (i=0; i<NPARTS; i++) {
         if (part[i].nuclid!=-1) {                              //crystal-like particle
            if (sizeofnucl[part[i].nuclid] == MAXSIZE) {        //the biggest nucleus (nuclei), blue color      
                color    =       'N';
            }
            else {                                              //green color
                color    =       'F';
            }
         }
         else {                                                 //liquid-like particle, red color
            color       =       'O';
         }
         p      =       part[i].p;
         fprintf(fhere,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
      }
*/
   }
   else if (vtype == 2) {                                       //all particles with nearest images
      fprintf(fhere, "#MAXSIZE=%ld\t# of biggest nuclei=%ld\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);
      for (i=0; i<NPARTS; i++) {
         if (part[i].nuclid!=-1) {                              //crystal-like particle
            if (sizeofnucl[part[i].nuclid] == MAXSIZE) {        //the biggest nucleus (nuclei), blue color      
                color    =       'N';
            }
            else {                                              //green color
                color    =       'F';
            }
         }
         else {                                         //liquid-like particle, red color
            color       =       'O';
         }

         if (PBC==1){           //cubic pbc images
            for (a=-1; a<=1; a++) {
               for (b=-1; b<=1; b++) {
                  for (c=-1; c<=1; c++) {
                     p.x        =       part[i].p.x + a * L;
                     p.y        =       part[i].p.y + b * L;
                     p.z        =       part[i].p.z + c * L;
                     fprintf(fhere,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
                  }
               }
            }
         }
         else if (PBC==2) {     //truncated octahedral pbc images
            p   =       part[i].p;
            fprintf(fhere,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);

            for (a=-1; a<=1; a+=2) {                    //octahedron images
               for (b=-1; b<=1; b+=2) {
                  for (c=-1; c<=1; c+=2) {
                     p.x        =       part[i].p.x + 0.5 * a * L;
                     p.y        =       part[i].p.y + 0.5 * b * L;
                     p.z        =       part[i].p.z + 0.5 * c * L;
                     fprintf(fhere,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
                  }
               }
            }
            for (a=-1; a<=1; a+=2) {                    //cubic images
               p.x      =       part[i].p.x + a * L;
               p.y      =       part[i].p.y;
               p.z      =       part[i].p.z;
               fprintf(fhere,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
            }
            for (b=-1; b<=1; b+=2) {
	       p.x      =       part[i].p.x;
	       p.y      =       part[i].p.y + b * L;
	       p.z      =       part[i].p.z;
	       fprintf(fhere,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
            }
            for (c=-1; c<=1; c+=2) {
	       p.x      =       part[i].p.x;
	       p.y      =       part[i].p.y;
	       p.z      =       part[i].p.z + c * L;
	       fprintf(fhere,"ATOM%7d  %c   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, color, i+1, p.x, p.y, p.z, 1.0, 0.0);
            }
         }
      }
   }
   else if (vtype==3) {         //Cell list illustration
      fprintf(fhere, "#MAXSIZE=%ld\t# of biggest nuclei=%ld\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);
      for (i=0; i<NPARTS; i++) {
         if (mod(part[i].icell,4)==0)
            fprintf(fhere,"ATOM%7d  N   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", 1, 1, part[i].p.x, part[i].p.y, part[i].p.z, 1.0, 0.0);
         else if (mod(part[i].icell,4)==1)
            fprintf(fhere,"ATOM%7d  F   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", 1, 1, part[i].p.x, part[i].p.y, part[i].p.z, 1.0, 0.0);
         else if (mod(part[i].icell,4)==2)
            fprintf(fhere,"ATOM%7d  O   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", 1, 1, part[i].p.x, part[i].p.y, part[i].p.z, 1.0, 0.0);
         else if (mod(part[i].icell,4)==3)
            fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", 1, 1, part[i].p.x, part[i].p.y, part[i].p.z, 1.0, 0.0);

      }
   }
   
   //print out reference points
   fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", -1, i+1, 0.0, 0.0, 0.0, 1.0, 0.0);
   fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", -1, i+1, L/2, L/2, L/2, 1.0, 0.0);
   fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", -1, i+1, L/2, L/2, -L/2, 1.0, 0.0);
   fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", -1, i+1, L/2, -L/2, L/2, 1.0, 0.0);
   fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", -1, i+1, L/2, -L/2, -L/2, 1.0, 0.0);
   fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", -1, i+1, -L/2, L/2, L/2, 1.0, 0.0);
   fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", -1, i+1, -L/2, L/2, -L/2, 1.0, 0.0);
   fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", -1, i+1, -L/2, -L/2, L/2, 1.0, 0.0);
   fprintf(fhere,"ATOM%7d  H   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", -1, i+1, -L/2, -L/2, -L/2, 1.0, 0.0);
   
   fclose(fhere);
   return 1;
}   
#endif /* TEST */


int Read_Conf(char * filename)		// read in ONE configuration 
{
   FILE *	fconf;
   char		dummy[80];
   long		i, j, id;

   if ((fconf=fopen(filename, "r"))==NULL)
      Exit("io", "Read_Conf", "Input conf. file not found or failed to open!");
   else {
      fscanf(fconf, "%s%ld", dummy, &TIMESTEP);			// read timestep
      fscanf(fconf, "%ld%ld%ld", &NSYSTEMS, &NMOLS, &NSITES);	// # of systems, mols, and sites
      NBOX	=	NSYSTEMS;				// temporary

      for (i=0; i<NSYSTEMS; i++)
         fscanf(fconf, "%lf%lf%lf", &(BOX[i].lx), &(BOX[i].ly), &(BOX[i].lz));	// box dimensions

      for (i=0; i<NMOLS; i++) {
         fscanf(fconf, "%ld%ld%ld", &id, &(mol[i].box), &(mol[i].nsites));	// molecule info.
         //fscanf(fconf, "%ld%ld%ld%ld%ld", &id, &(mol[i].box), &(mol[i].nsites), &(mol[i].fix), &(mol[i].flip));
         if (id!=i)	
            Exit("io", "Read_Conf", "molecule id mismatch");
         for (j=0; j<mol[i].nsites; j++)					// site info.
            fscanf(fconf, "%ld%lf%lf%lf", &(mol[i].type[j]), 
		&(mol[i].p[j].x), &(mol[i].p[j].y), &(mol[i].p[j].z));
      }
   }
   fclose(fconf);
   return	1;
}

int Read_MultiConf(FILE *fPtr)		// read in multiple configurations
{
   char		dummy[80];
   long		i, j, id;

   fscanf(fPtr, "%s%ld", dummy, &TIMESTEP);			// read in timestep
   if (!feof(fPtr)) {
      fscanf(fPtr, "%ld%ld%ld", &NSYSTEMS, &NMOLS, &NSITES);	// # of systems, mols, and sites
      NBOX	=	NSYSTEMS;				// temporary
      for (i=0; i<NSYSTEMS; i++)
         fscanf(fPtr, "%lf%lf%lf", &(BOX[i].lx), &(BOX[i].ly), &(BOX[i].lz));	// box dimension
      for (i=0; i<NMOLS; i++) {
         fscanf(fPtr, "%ld%ld%ld", &id, &(mol[i].box), &(mol[i].nsites));	// molecule info.
         //if (id!=i)	Exit("io", "Read_MultiConf", "molecule id mismatch");
         for (j=0; j<mol[i].nsites; j++)					// site info.
            fscanf(fPtr, "%ld%lf%lf%lf", &(mol[i].type[j]), 
		&(mol[i].p[j].x), &(mol[i].p[j].y), &(mol[i].p[j].z));
      }
   }
   return	(!feof(fPtr));
}

/*
int Read_Conf(char * filename)		// read in LBOX and coord. of NMOLS molecules
{					// from a binary file
   long		i, j;
   FILE *	fconf;
   molstruct	*moli;

   if ((fconf=fopen(filename, "r"))==NULL) {
      Exit("io", "Read_Conf", "Input conf. file not found or failed to open!");
   }
   else {
//      fseek(fconf, NPARTS*sizeof(vector)+sizeof(double), SEEK_SET);
      moli	=	mol;

      fread(&(LBOX), sizeof(double), 1, fconf);

      for (i=0; i<NMOLS; i++) {
         for (j=0; j<moli->nsites; j++) {
            fread(moli->p+j, sizeof(vector), 1, fconf);
         }
         moli	++;
      }
      fclose(fconf);
      return	1;
   }
}
*/

int Write_Conf(long timestep)			// write one configuration
{
   FILE		*fPtr;
   molstruct	*moli;
   long		i, j;
   double	dL = unit.LENGTH;

   if (timestep==-1) {				// open a file for a single conf.
      if ((fPtr=fopen("finalconf", "w"))==NULL)		
         Exit("io", "Write_Conf", "finalconf file cannot open");
   }
   else {
      fPtr	=	fdump;			// goes to dump file for continuous dumping
   }

   fprintf(fPtr, "TIMESTEP	%ld\n", timestep);
   fprintf(fPtr, "%ld\t%ld\t%ld\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx * dL, BOX[i].ly * dL, BOX[i].lz * dL);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%ld\t%ld\t%ld\n", moli-mol, moli->box, moli->nsites);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%ld\t%f\t%f\t%f\n", moli->type[i], 
		moli->p[i].x * dL, moli->p[i].y * dL, moli->p[i].z * dL);
   }

   if (timestep==-1) {
      fclose(fPtr);
   }
   return	1;
}


void dump_conf(long timestep)	// dump conf. in lammps format, added on 8/2/2012
{
   molstruct	*moli;
   long		system = 0, atomid, i;
   long		ix=0, iy=0, iz=0;
   double	vx=0.0, vy=0.0, vz=0.0;
   static FILE	*fPtr;
   static long	init_dump=1;

   if (init_dump) {
      init_dump	=	0;
      if (!FLAG_restart) {
         fPtr	=	fopen("dumpfile","w");
      }
      else {
         fPtr	=	fopen("dumpfile","a");
      }
   }

   fprintf(fPtr, "ITEM: TIMESTEP\n");
   fprintf(fPtr, "%ld\n", timestep);
   fprintf(fPtr, "ITEM: NUMBER OF ATOMS\n");
   fprintf(fPtr, "%ld\n", NSITES);
   fprintf(fPtr, "ITEM: BOX BOUNDS\n");
   fprintf(fPtr, "%lf\t%lf\n", -0.5*BOX[system].lx*unit.LENGTH, 0.5*BOX[system].lx*unit.LENGTH);
   fprintf(fPtr, "%lf\t%lf\n", -0.5*BOX[system].ly*unit.LENGTH, 0.5*BOX[system].ly*unit.LENGTH);
   fprintf(fPtr, "%lf\t%lf\n", -0.5*BOX[system].lz*unit.LENGTH, 0.5*BOX[system].lz*unit.LENGTH);
   fprintf(fPtr, "ITEM: ATOMS id mol type x y z vx vy vz ix iy iz\n");

   atomid	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         fprintf(fPtr, "%ld %ld %ld", atomid+1, moli-mol+1, moli->type[i]+1); 
         fprintf(fPtr, " %lf %lf %lf", 
		moli->p[i].x*unit.LENGTH, moli->p[i].y*unit.LENGTH, moli->p[i].z*unit.LENGTH);
         fprintf(fPtr, " %lf %lf %lf", vx, vy, vz);	
         fprintf(fPtr, " %ld %ld %ld\n", ix, iy, iz); 
         atomid	++;
      }
   }
   return;
}

/*
int Write_Conf(long timestep)		// Write configuration files in SI units in human-readable format
{
   long		i, j;
   FILE *	fconf;
   char		confname[256], tempname[256];
   molstruct	*moli;
   vector	l;

 
//   frame ++;				// index the visualization output file name
//   sprintf(tempname,"%ld",frame);

   if (timestep!=-1) {
      sprintf(tempname, "%ld", timestep);

      strncpy(confname,"conf.00000000",13-strlen(tempname));
      confname[13-strlen(tempname)]='\0';
      strcat(confname,tempname);
   }
   else {
      strcpy(confname, "finalconf");
   }

   if ((fconf=fopen(confname,"w"))==NULL) {
      Exit("io", "Write_Conf", "Visualization file cannot be opened!");
   }
   else {
      fprintf(fconf, "%ld\t%ld\t%ld\n", NSYSTEMS, NMOLS, NSITES);

      for (i=0; i<NSYSTEMS; i++)
         fprintf(fconf, "%f\t%f\t%f\n", BOX[i].lx*unit.LENGTH, BOX[i].ly*unit.LENGTH, BOX[i].lz*unit.LENGTH);

      for (moli=mol; moli<mol+NMOLS; moli++) {
         fprintf(fconf, "%ld\t%ld\t%ld\n", moli-mol, moli->box, moli->nsites);
         //fprintf(fconf, "%ld\t%ld\t%ld\t%ld\t%ld\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
         //MolInBox(moli);
         for (i=0; i<moli->nsites; i++) 
            fprintf(fconf, "%ld\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
      }
      fclose(fconf);
      return	1;
   }
}
*/
/*
int Write_Conf()	// write to binary files
{
   long		i, j;
   FILE *	fconf;
   char		confname[256], tempname[256];
   molstruct	*moli;
   vector	l;
 
   frame ++;				//index the visualization output file name
   sprintf(tempname,"%ld",frame);
   strncpy(confname,"00000000",8-strlen(tempname));
   confname[8-strlen(tempname)]='\0';
   strcat(confname,tempname);
   
   if ((fconf=fopen(confname,"w"))==NULL) {		
      Exit("io", "Write_Conf", "Visualization file cannot be opened!");
   }
   else {
      moli	=	mol;

      fwrite(&(LBOX), sizeof(double), 1, fconf);

      for (i=0; i<NMOLS; i++) {
         for (j=0; j<moli->nsites; j++) {
            l	=	moli->p[j];
            l	=	MapInBox2(&l, PBC, moli->box);
            fwrite(&l, sizeof(vector), 1, fconf);
         }
         moli	++;
      }
      fclose(fconf);
      return	1;
   }
}
*/


int Visualize(int vtype)
{
   long		i, j;	
   FILE		*fhere;
   vector	pvisual, p0, p;
   vector	l;
   double	L = LBOX;
   char		color;
   long		a, b, c;
   char		vsulname[256], tempname[256];
   molstruct	*moli;
	
   frame ++;				//index the visualization output file name
   sprintf(tempname,"%ld",frame);
   strncpy(vsulname,"00000000",8-strlen(tempname));
   vsulname[8-strlen(tempname)]='\0';
   strcat(vsulname,tempname);

   switch(vtype) {
      case 1:   strcat(vsulname, ".pdb");       break;
      case 0:   strcat(vsulname, "full.pdb");   break;
      case 2:   strcat(vsulname, "img.pdb");    break;
      case 3:   strcat(vsulname, "cell.pdb");   break;
      default:  break;
   }

   if ((fhere=fopen(vsulname,"w"))==NULL) {
      Exit("io", "Visualize", "Visualization file cannot be opened!");
   }

   moli		=	mol;
   for (i=0; i<NMOLS; i++) {			// for all chain molecules
      for (j=0; j<moli->nsites; j++) {		// on each chain molecule
         l	=	moli->p[j];
         l	=	MapInBox2(&l, PBC, moli->box);
         fprintf(fhere, "ATOM%7ld  C%12ld    %8.3f%8.3f%8.3f\n", i*(moli->nsites)+j+1, i+1, l.x, l.y, l.z);
      }
      moli	++;
   }
   fclose(fhere);
   return	1;
}


void Printout()			// output everything
{
  long 		i, j;
  vector 	dp;
  double 	r2, costheta, sintheta, phi, alpha;
 
  //fprintf(foutput, "Box length=\t%f\n", LBOX);

  /*	
  fprintf(foutput, "Mol#\tx\t\ty\t\tz\t\tnbond\tql\t\tnconnect\n");
  for (i=0; i<NPARTS; i++) 		//print out particles coordinates 
    fprintf(foutput, "%ld\t%f\t%f\t%f\t%ld\t%f\t%ld\n", i+1, part[i].p.x, part[i].p.y, part[i].p.z, part[i].nbond, part[i].ql,part[i].nconnect);
  for (i=0; i<NPARTS-1; i++) {
    for (j=i+1; j<NPARTS; j++) {
      dp.x	=	part[j].p.x	-	part[i].p.x;
      dp.y	=	part[j].p.y	-	part[i].p.y;
      dp.z	=	part[j].p.z	-	part[i].p.z;
      MapInBox(&dp);
      r2	=	dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;
      fprintf(foutput, "r(%ld,%ld)=%f\n",i+1,j+1,sqrt(r2));
      costheta	=	dp.z/sqrt(r2);
      sintheta	=	sqrt(1-costheta*costheta);
        phi	=	atan2( dp.y, dp.x );
      alpha	=	pow((sqrt(r2)-BONDCUT)/(sigma-BONDCUT),2);
      fprintf(foutput, "%ld->%ld:\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i+1, j+1, dp.x, dp.y, dp.z, costheta, sintheta, phi, alpha);
    }
  }	
  fprintf(foutput, "Energy:\t%f\t%f\t%f\n",VSYSTEM.lj6, VSYSTEM.lj12, VSYSTEM.tot);

  Print_Verlet();
  Print_Verlet2();
  Print_q();
  Print_qproduct();
*/

  //Print_Verlet();
  //Print_Clist(); 
  //Print_Q();
  //Print_Nuclei();
  //Print_Histogram();
  //Visualize(0);  
  return;
}
