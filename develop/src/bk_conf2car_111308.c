src/bk_embed.c                                                                                      0000600 0143352 0000144 00000007135 11107577042 012762  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	embed.c
	author:		Peng Yi at MIT
	date:		May 18, 2008
	purpose:	embed a smaller system to the center of a larger system, do NO analysis
	note:		require setup file
*/

#define __MAIN_PROGRAM
#include "header.h"


int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   char		s[255], ff[255];
   molstruct	*moli;
   FILE		*fPtr;

   vector	com,
		Linner[MAXNSYSTEMS];		// inner BOX dimension
   long		i, j, n, system,
		nsystems,			// for one system only, 5/18/2008
		nmols1, nsites1,		// NMOLS and NSITES for inner box
		nmolstot, nsitestot,		// total NMOLS and NSITES
		id;
   double	lambda;
   long		overlap;

   if (argc<4) {
      printf("embed (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tembed file1 file2 lambda\n\n");
      printf("Notes:\n");
      printf("\t* file1 is the configuration file for inner box, file2 for outer box\n\n");
      printf("\t* lambda <= 1.0 \n\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }
   lambda	=	atof(argv[3]);

   InitMols(MAXNMOLS, 1);	// allocate memory for molecules, MAXNMOLS must > the sum of NMOLS
   GetSetup(argv);		// mainly PBC I think

   // Read in the conf. of the inner system

   Read_Conf(argv[1]);
   for (moli=mol; moli<mol+NMOLS; moli++)	// arrange the inner system
      MolInBox2(moli);

   nsystems	=	NSYSTEMS;		// save inner system information
   nmols1	=	NMOLS;
   nsites1	=	NSITES;
   nmolstot	=	NMOLS;
   nsitestot	=	NSITES;

   for (i=0; i<NSYSTEMS; i++) {
      Linner[i].x	=	BOX[i].lx;
      Linner[i].y	=	BOX[i].ly;
      Linner[i].z	=	BOX[i].lz;
   }

   // Read in the conf. of the outer system

   if (!(fPtr=fopen(argv[2], "r")))		// open the conf. file of the outer system
      Exit("embed", "main", "big system conf. file failed to open.");

   fscanf(fPtr, "%ld%ld%ld", &NSYSTEMS, &NMOLS, &NSITES);		// update NSYSTEMS, NMOLS, NSITES
   if (nsystems != NSYSTEMS) {}			// might do something in the future
   for (i=0; i<NSYSTEMS; i++)
      fscanf(fPtr, "%lf%lf%lf", &(BOX[i].lx), &(BOX[i].ly), &(BOX[i].lz));	// update BOX dimension

   moli	=	mol + nmols1;			// clear, although not necessary
   for (i=0; i<NMOLS; i++) {			// start reading in the molecules in the outer box
      fscanf(fPtr, "%ld%ld%ld", &id, &(moli->box), &(moli->nsites));
      for (j=0; j<moli->nsites; j++)
         fscanf(fPtr, "%ld%lf%lf%lf", &(moli->type[j]), &(moli->p[j].x), &(moli->p[j].y), &(moli->p[j].z));

      MolInBox2(moli);				// map back in outer box

      overlap	=	0;			// is this molecule overlap with inner box?
      for (j=0; j<moli->nsites; j++) {
         if (fabs(moli->p[j].x)*lambda < 0.5*Linner[moli->box].x  || 
		fabs(moli->p[j].y)*lambda < 0.5*Linner[moli->box].y || 
		fabs(moli->p[j].z)*lambda < 0.5*Linner[moli->box].z ) {
            overlap	=	1;
            break;
	 }
      }
      if (overlap==0) {
         nmolstot	++;
	 nsitestot	+=	moli->nsites;
	 moli	++;
      }
/*
      com	=	CenterofMass(moli);	// center of mass of one chain
      com.x	*=	lambda;		// give some gap
      com.y	*=	lambda;
      com.z	*=	lambda;

      if (!(com.x < 0.5*Linner[moli->box].x && com.x >= -0.5*Linner[moli->box].x &&
	   com.y < 0.5*Linner[moli->box].y && com.y >= -0.5*Linner[moli->box].y &&
	   com.z < 0.5*Linner[moli->box].z && com.z >= -0.5*Linner[moli->box].z) ) {	// com NOT in the inner box
         nmolstot	++;
	 nsitestot	+=	moli->nsites;	
	 moli		++;
      }
*/
   }
   NMOLS	=	nmolstot;
   NSITES	=	nsitestot;

   // PRINT OUT CONFIGURATION

   system	=	0;			// only one system for now
   unit.LENGTH	=	1.0;			// not doing any analysis
   Write_Conf(-1);				// output combined configuration

   return	0;
}
                                                                                                                                                                                                                                                                                                                                                                                                                                   src/block.c                                                                                         0000600 0143352 0000144 00000024733 11231704771 012326  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    block.c
    author:     Peng Yi at MIT
    date:       Sept 4, 2007
    purpose:    block average
    note:	Flyvbjerg and Petersen, JCP 91, 461	
*/

#define __MAIN_PROGRAM
#include "header.h"

#define	MAXNVARS	50			// maximum number of variables that need to do block average


typedef struct{
   double	var[MAXNVARS];			// NVARS variables to do block average
} avestruct;


void SetZero(avestruct *a)
{
   long		i;
   for (i=0; i<MAXNVARS; i++) {
      a->var[i]	=	0;
   }
}


avestruct Add(avestruct a, avestruct b)
{
   long		i;
   avestruct	c;
   for (i=0; i<MAXNVARS; i++)
      c.var[i]	=	a.var[i] + b.var[i];

   return	c;
}


avestruct Subtr(avestruct a, avestruct b)
{
   long		i;
   avestruct	c;
   for (i=0; i<MAXNVARS; i++)
      c.var[i]	=	a.var[i] - b.var[i];

   return	c;
}


avestruct Mult(avestruct a, double x)
{
   long		i;
   avestruct	b;
   for (i=0; i<MAXNVARS; i++)
      b.var[i]	=	a.var[i] * x;

   return	b;      
}


avestruct Square(avestruct a)
{
   long		i;
   avestruct	b;
   for (i=0; i<MAXNVARS; i++)
      b.var[i]	=	a.var[i] * a.var[i];
   return	b;
}


avestruct SQRT(avestruct a)
{
   long		i;
   avestruct	b;
   for (i=0; i<MAXNVARS; i++) {
      if (a.var[i] < -ZERO)
	 printf("block average, sqrt to negtive number\n");
      else
         b.var[i]	=	sqrt(a.var[i]);
   }
   return	b;
}


void Print_System(avestruct a)
{
   long		i;
   for (i=0; i<MAXNVARS; i++)
      printf("%f\t", a.var[i]);
   printf("\n");
}


int main(int argc, char * argv[])
{
   FILE			*fhist;
   char			temp[80],
			comments[1024],
			varname[MAXNVARS][256];	// variable names
   double		dummy;
   long			NVAR=0, 		// number of variables in the data file
			flag=1;
   char			current, previous='1';
   long			i, j, k, nrecords, 
			layer,			// identify the way to block the samples
						// each block contains 2^layer samples
			nblock[20],		// block numbers in each way of blocking, maxlayer=19
 			NBLOCK = 16;		// particularly calculate the case with 16 blocks

   avestruct		item;			// read-in			
   avestruct		elem[20],		// average within each block
			ave[20],		// average of elem
			aveSq[20],		// average of sq of elem
			block[16],		// particular case: divide all sampling into 16 blocks
   			var[20],		// variance of the sampling data
			se[20],			// standard deviation of the mean
			see[20];		// standard deviation of the error

   double		biasvar, bias,
   			pNmax[MAXNMOLS],
			dP2m =	0.005,
			pP2m[200],		// P2m has maximum value 1, grid width dP2m
			dQ6  =  0.001,
			pQ6[1000];		// Q6 maximum value 1, grid width dQ6
   long			nNmax[MAXNMOLS],	// # of samples
			nP2m[200],
			nQ6[1000];
   long			minNmax = MAXNMOLS-1, 	// min INDEX that has non-zero count
			maxNmax = 0,		// max INDEX that has count
			minP2m	= 200-1,
			maxP2m	= 0,
			minQ6	= 1000-1,
			maxQ6	= 0;

   if (argc<3) {
      printf("block (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tblock filename Nrecords\n\n");
      printf("Notes:\n");
      printf("\t* Nrecords is the number of records you want to do block statistics\n\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }
   nrecords	=	atol(argv[2]);

   GetSetup(argv);				// Read in setup information
   D_PRESSURE = 1;	// temporary
   D_XTALSIZE = 1;	// temporary
   InitSample();	// temporary

   printf("kP2 = %f  P2middle = %f\n", kP2, P2middle);

   if ( (fhist=fopen(argv[1],"r"))==NULL ) {		// read in data file
      Exit("block", "main", "block average input file cannot open.");
   }
   else { 
      while(flag) {						// determine how many variables
         if ( (current=fgetc(fhist)) == '\n') {
            NVAR	++;
            flag	=	0;
         }
	 else if (isspace(current) && !isspace(previous)) 
            NVAR	++;

	 previous	=	current;
      }

      if (NVAR > MAXNVARS) {
	 Exit("block", "main", "block average too many variables.\n");
      }
      rewind(fhist);

      // Initialization

      for (layer=0; layer<20; layer++) {
         SetZero(&elem[layer]);
         SetZero(&ave[layer]);
         SetZero(&aveSq[layer]);
         nblock[layer]	=	0;
      }
      for (i=0; i<MAXNMOLS; i++) {
         pNmax[i]	=	0.0;
         nNmax[i]	=	0;
      }
      for (i=0; i<200; i++) {
         pP2m[i]	=	0.0;
         nP2m[i]	=	0;
      }
      for (i=0; i<1000; i++) {
         pQ6[i]		=	0.0;
	 nQ6[i]		=	0;
      }

      // Read in data

      for (i=0; i<NVAR; i++) {					// scan in titles on the first line
         fscanf(fhist, "%s", varname[i]);
      } 
      fgets(comments, sizeof(comments), fhist);

      fgets(comments, sizeof(comments), fhist);			// disgard the first record 

      for (j=0; j<nrecords; j++) {      
	 if (feof(fhist))					// if end of data file
	    break;

         for (i=0; i<NVAR; i++) {				// scan in one line of data 
	    fscanf(fhist, "%s", temp);
	    item.var[i]	=	atof(temp);
         }
         fgets(comments, sizeof(comments), fhist);

         if (j>=Nequil) {					// block average for production run
	    k	=	j - Nequil;

            // calculate bias of each record

            if (1==dynvar)					// Nmax as dynamic variable
               biasvar	=	item.var[11];	
	    else if (3==dynvar)					// P2m
	       biasvar	=	item.var[7];
	    else if (2==dynvar)					// Q6
     	       biasvar	=	item.var[16];
 
            bias	=	exp(0.5*kP2*(biasvar-P2middle)*(biasvar-P2middle));

            // collect bias-fixed record count of dynamic variable

            if (dynvar==1) {
               nNmax[(int) item.var[11]]	++;
	       pNmax[(int) item.var[11]]	+=	bias;

               if ((int) item.var[11] > maxNmax)
	          maxNmax	=	(int) item.var[11];
               if ((int) item.var[11] < minNmax)
                  minNmax	=	(int) item.var[11];
            }
            else if (dynvar==2) {
	       nQ6[(int) (item.var[16]/dQ6)]	++;
               pQ6[(int) (item.var[16]/dQ6)]	+=	bias;

               if ((int) (item.var[16]/dQ6) > maxQ6)
                  maxQ6		=	(int) (item.var[16]/dQ6);
               if ((int) (item.var[16]/dQ6) < minQ6)
                  minQ6		=	(int) (item.var[16]/dQ6);
            }
            else if (dynvar==3) {
               nP2m[(int) (item.var[7]/dP2m)]	++;
	       pP2m[(int) (item.var[7]/dP2m)]	+=	bias;

	       if ((int) (item.var[7]/dP2m) > maxP2m)
                  maxP2m	=	(int) (item.var[7]/dP2m);
	       if ((int) (item.var[7]/dP2m) < minP2m)
                  minP2m	=	(int) (item.var[7]/dP2m);
	    }

            /*printf("P2m=%f  pP2m[%d]=%f   Nmax=%f  pNmax[%d]=%f  bias=%f \n", 
		item.var[7], (int) (item.var[7]/0.01), pP2m[(int) (item.var[7]/0.01)], 
		item.var[11], (int) item.var[11], pNmax[(int) item.var[11]], bias);*/

	    // collect biased probability distribution 

	    if (dynvar==1)
	       PutInDistribution(D_Xtalsize, item.var[11], 1.0, 1.0);	//temporary, collect Nmax counts
	    //PutInDistribution(D_Energy, item.var[1], 1.0, 1.0);		//temporary
            else if (dynvar==3)
               PutInDistribution(D_Pressure, item.var[7], 1.0, 1.0);	//temporary, collect P2m counts

	    /* block average */

            for (layer=0; layer<20; layer++) {			// each block has 2^layer samples

               for (i=0; i<NVAR; i++) {         
                  elem[layer].var[i]	+=	item.var[i]/intpow(2, layer);	// ave within each block in each layer
               }

               if ( !mod(k+1, intpow(2, layer)) ) {

                  if (layer==intlog(2, Nprod/NBLOCK)) {
                     block[(k+1)*NBLOCK/(Nprod)-1]	=	elem[layer];
                  }
	
                  ave[layer]	=	Add( ave[layer], elem[layer] );		// ave across the blocks in each layer
                  aveSq[layer]	=	Add( aveSq[layer], Square(elem[layer]) ); // ave2 across the blocks in each layer
                  nblock[layer]	++;
                  SetZero(&elem[layer]);
               }
            }
         }	// if (j>= Nequil)
      }

      S_PrintAll();	//temporary

      if (dynvar==1) {
      printf("\nBias-fixed probability distribution of Nmax.  Index = [%d, %d] \n", minNmax, maxNmax);
      for (i=minNmax; i<=maxNmax; i++)
         printf("%f\n", -log(pNmax[i]));

      printf("\nCount of samples of Nmax. Index = [%d, %d] \n", minNmax, maxNmax);
      for (i=minNmax; i<=maxNmax; i++)
         printf("%d\n", nNmax[i]);
      }

      else if (dynvar==3) {
      printf("\nBias-fixed probability distribution of P2m. Index = [%d, %d], binsize = %f.\n", minP2m, maxP2m, dP2m);
      for (i=minP2m; i<=maxP2m; i++)
         printf("%f\n", -log(pP2m[i]));

      printf("\nCount of samples of P2m. Index = [%d, %d]\n", minP2m, maxP2m);
      for (i=minP2m; i<=maxP2m; i++)
         printf("%d\n", nP2m[i]);
      }

      else if (dynvar==2) {
      printf("\nBias-fixed probability distribution of Q6. Index = [%d, %d], binsize = %f.\n", minQ6, maxQ6, dQ6);
      for (i=minQ6; i<=maxQ6; i++)
         printf("%f\n", -log(pQ6[i]));

      printf("\nCount of samples of Q6.  Index = [%d, %d]\n", minQ6, maxQ6);
      for (i=minQ6; i<=maxQ6; i++)
         printf("%d\n", nQ6[i]);
      }

      for (layer=0; layer<20; layer++) {
         if (nblock[layer] > 1) {
            ave[layer]		=	Mult( ave[layer], 1.0/nblock[layer] );
            aveSq[layer]	=	Mult( aveSq[layer], 1.0/nblock[layer] );
            var[layer]		=	Mult( Subtr(aveSq[layer], Square(ave[layer]) ), 1.0/(nblock[layer]-1) );
            se[layer]		=	SQRT( var[layer] );
            see[layer]		=	Mult( se[layer], 1.0/sqrt(2*(nblock[layer]-1)) );
	 }
      }

      printf("\nStart block average ......\n\n");
      printf("%d blocks ... \n", NBLOCK);
      printf("\t#block");
      for (j=0; j<NVAR; j++) {
         printf("\t %s", varname[j]);
      }

      for (i=0; i<NBLOCK; i++) {			// print out the results with 16 blocks
	 printf("\n\t%d", i+1);

         for (j=0; j<NVAR; j++) {
	    printf("\t%f", block[i].var[j]);
         }
      }
      printf("\n");

      printf("\n\tAve:");
      for (j=0; j<NVAR; j++) {
         printf("\t%f", ave[0].var[j]);	//ave[0]=ave[1]=ave[2]=...
      }
      printf("\n\n");
      printf("\n\tSe:");
      for (j=0; j<NVAR; j++) {
         printf("\t%f", se[intlog(2,(Nprod))-4].var[j]);	
      }
      printf("\n");

      printf("\nBlock average calculation ... \n");
      printf("\ts.p.b.");
      printf("\tstderror");
      printf("\tstdee");
      printf("\n");

      for (layer=0; layer<20; layer++) {
         if (nblock[layer] > 1) {
            printf("\t%d", intpow(2,layer));
	    for (j=0; j<NVAR; j++){
               printf("\t%f\t%f", se[layer].var[j], see[layer].var[j]);
            }
            printf("\n");
         }
      }

      fclose(fhist);    
   }

   return	0;
}


                                     src/builder.c                                                                                       0000600 0143352 0000144 00000060523 11531746063 012662  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    builder.c
    author:     Peng Yi at MIT
    date:       November 15, 2007
    purpose:    Initialize polymer configuration for one single box.
    notes:
		July 25, 2009 enable polydispersity in Amorphous()
*/

#define __MAIN_PROGRAM
#include "header.h"

#define Rcut2	0.8*16.0801		// cutoff^2, unit Angstrom

/* Variables that are only useful in builder.c */

long		NX, NY, NZ;
double		cella, cellb, cellc, angle;
char		title1[256], title2[256], comments[256];

void BuilderReadSetup()
{
   FILE		*fPtr;
   long		i;
   char		c[256], name[256], a[256];

   if ( (fPtr=fopen("in","r"))==NULL )
      printf("Open setup file failed.\n");
   else {
      fgets(title1, sizeof(title1), fPtr);
      fgets(title2, sizeof(title2), fPtr);
      fgets(comments, sizeof(comments), fPtr);

      fscanf(fPtr, "%s", name);
      while(strcmp(name, "END_SETUP")) {

         if (!strcmp(name, "NMOLS"))			GetLVar(fPtr, 1, &NMOLS);
         else if (!strcmp(name, "NSITES"))		GetLVar(fPtr, 1, &NSITES);
         else if (!strcmp(name, "NBOX")) {		GetLVar(fPtr, 1, &NBOX);
            NSYSTEMS	=	NBOX;
	 }
         else if (!strcmp(name, "Dpoly"))		GetDVar(fPtr, 1, &Dpoly);
         else if (!strcmp(name, "Rho"))			GetDVar(fPtr, 1, &Rho);
         else if (!strcmp(name, "kT"))			GetDVar(fPtr, 1, &kT);
         else if (!strcmp(name, "PBC"))			GetLVar(fPtr, 1, &PBC);
         else if (!strcmp(name, "NX"))			GetLVar(fPtr, 1, &NX);
         else if (!strcmp(name, "NY"))			GetLVar(fPtr, 1, &NY);
         else if (!strcmp(name, "NZ"))			GetLVar(fPtr, 1, &NZ);
         else if (!strcmp(name, "a"))			GetDVar(fPtr, 1, &cella);
         else if (!strcmp(name, "b"))			GetDVar(fPtr, 1, &cellb);
         else if (!strcmp(name, "c"))			GetDVar(fPtr, 1, &cellc);
         else if (!strcmp(name, "angle"))		GetDVar(fPtr, 1, &angle);

         /* Forcefield parameters */

	 else if (!strcmp(name, "NTYPES")) 		GetLVar(fPtr, 1, &NTYPES);
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
         else if (!strcmp(name, "FIXBONDL")) {
	    fscanf(fPtr, "%s%ld", c, &(FIXBONDL));	fgets(comments,sizeof(comments),fPtr);
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
         else if (!strcmp(name, "FIXBONDA")) {
	    fscanf(fPtr, "%s%ld", c, &(FIXBONDA));	fgets(comments,sizeof(comments),fPtr);
	 }
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
         else
            fgets(comments, sizeof(comments),fPtr);
	 
         fscanf(fPtr, "%s", name);
      }
      fgets(comments,sizeof(comments),fPtr);
      fclose(fPtr);
   }
   return;
}


void BuilderPrintSetup()
{
   FILE		*fPtr;
   long		i, j;
   char		name[256];
 
   if ( (fPtr=fopen("builder.log","a"))==NULL )
      printf("Builder log file failed.\n");
   else {
      fprintf(fPtr, "%s", title1);
      fprintf(fPtr, "%s", title2);
      fprintf(fPtr, "Program Version: \t%s\n", VERSION);
      for (i=0; i<80; i++)	
	 fprintf(fPtr, "*");
      fprintf(fPtr, "\n");

      fprintf(fPtr, "\n/* Box dimension */\n\n");

      fprintf(fPtr, "NBOX   \t=\t%d\n", NBOX);
      fprintf(fPtr, "NMOLS  \t=\t%d\n", NMOLS);
      fprintf(fPtr, "NSITES \t=\t%d\n", NSITES);
      fprintf(fPtr, "Dpoly  \t=\t%f\n", Dpoly);
      fprintf(fPtr, "Rho    \t=\t%f\n", Rho);
      fprintf(fPtr, "PBC    \t=\t%d\n", PBC);
      fprintf(fPtr, "NX     \t=\t%d\n", NX);
      fprintf(fPtr, "NY     \t=\t%d\n", NY);
      fprintf(fPtr, "NZ     \t=\t%d\n", NZ);
      fprintf(fPtr, "cella  \t=\t%f\n", cella);
      fprintf(fPtr, "cellb  \t=\t%f\n", cellb);
      fprintf(fPtr, "cellc  \t=\t%f\n", cellc);
      fprintf(fPtr, "angle  \t=\t%f\n", angle);
     
      fprintf(fPtr, "\n/* Forcefield parameters */\n\n");
 
      fprintf(fPtr, "NTYPES	\t=\t%d\n", NTYPES);
      fprintf(fPtr, "DLJ	\t=\t%d\n", DLJ);
      fprintf(fPtr, "TYPEMASS\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].M);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "SIGMA	\t=\t{ ");
      for (i=0; i<NTYPES; i++)	        fprintf(fPtr, "%f\t", type[i].SIGMA);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "EPSILON	\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].EPSILON);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "KSTRETCH\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].KSTRETCH);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "LSTRETCH\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].LSTRETCH);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "KBENDING\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].KBENDING);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "THETA	\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].THETA);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "HS	\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].HS);
      fprintf(fPtr, "}\n");
      for (j=0; j<6; j++) {
         sprintf(name, "TORSION%c", '0'+j);
         fprintf(fPtr, "%s\t=\t{ ", name);
         for (i=0; i<NTYPES; i++)	fprintf(fPtr, "%f\t", type[i].TORSION[j]);
         fprintf(fPtr, "}\n");
      }
      fclose(fPtr);
   }
   return; 
}

double tors_angle()
{
   double		angle, cosb, select;
   static double	value[360], accum[360], sum;
   static long		init =1;
   long			i;
   typestruct		t = type[0];

   if (init) {
      sum	=	0;
      for (i=0; i<360; i++) {
	 value[i]	=	0;
      }
      for (i=0; i<360; i++) {
         angle	=	M_PI/180.0*i;    	// [0, 2pi]  
         cosb	=	cos(angle);
         value[i]	=	exp(-OPLS2(cosb, t.TORSION[1], t.TORSION[2], t.TORSION[3])/(kT* N_KB * N_NAV * 1e-3));
	 sum	+=	value[i];
      }
/*
      cosb	=	1.0;			// trans
      value[0]	=	exp(-OPLS2(cosb, t.TORSION[1], t.TORSION[2], t.TORSION[3])/kT);
      cosb	=	-0.5;			// g+ and g-
      value[1]	=	exp(-OPLS2(cosb, t.TORSION[1], t.TORSION[2], t.TORSION[3])/kT);
      value[2]	=	value[1];
      sum	=	value[0] + value[1] + value[2];
*/
      init	=	0; 
   }
   select	=	ran1(seed)*sum;

   i	=	0;
   while (select>value[i]) {
      select	-=	value[i];
      i	++;
   }
   return	M_PI + M_PI/180.0*i;

/*
   if (select < value[0]/sum) 
      return	0;
   else if (select < (value[0]+value[1])/sum)
      return	M_PI * 2.0/3; 
   else
      return	- M_PI * 2.0/3;
*/
}


void Amorphous(long nsystems, long nmols, long nsites)	// generate random chain conf.
{
   molstruct	*moli, *molj;
   long		i, j, n, avensites, ndelta, ibox, system, reject, nreject;
   double	r2, vol, numberdensity;
   double	LX, LY, LZ;
   vector	dp;
   sphere	s;
   typestruct	t = type[0];

  static long		init=1;
  static long		Dnsites, maxnsites, minnsites;		// for polydispersity

  if (init) {
     Dnsites	=	(int) (Dpoly*nsites/nmols);
     maxnsites	=	nsites/nmols + Dnsites;
     minnsites	=	nsites/nmols - Dnsites;
     init	=	0;
  }

   double 	l0 = type[0].LSTRETCH, 		// bond length
		theta0 = type[0].THETA,		// bond angle
                mass = type[0].M;

   numberdensity	=	Rho*0.6022/mass;	// Rho in unit of g/cm^3
   vol	=	(nsites/nsystems)/numberdensity;	// vol in unit of Anstrom^3
   LX	=	pow(vol, 1.0/3);			// cubic
   LY	=	LX;
   LZ	=	LX;

   printf("Random number test (print out 5 random numbers):\n");
   for (i=0; i<5; i++) {
      printf("%f\n", ran1(seed));
   }

   for (i=0; i<nsystems; i++) {
      BOX[i].lx		=	LX;
      BOX[i].ly		=	LY;
      BOX[i].lz		=	LZ;
      BOX[i].temp	=	kT;
   }

   for (moli=mol; moli<mol+nmols; moli++) {
      n			=	(int) (ran1(seed) * (2*Dnsites+1));	// [0, 2Dnsites]
      moli->nsites	=	minnsites + n;		// polydispersity

      moli->type[0]			=	0;	// set the type of end sites
      moli->type[moli->nsites-1]	=	0;
      for (i=1; i<moli->nsites-1; i++) {		// set the type of middle sites
         moli->type[i]			=	1;
      }
      for (i=0; i<moli->nsites; i++)			// set parent site
         moli->parent[i]	=	i-1;
   }

   for (moli=mol; moli<mol+nmols; moli++) {
      system	=	(moli-mol)/(nmols/nsystems);	// set box id
      moli->box	=	system;

      // Place the first atom on a chain
printf("moli=%d \n", moli-mol);
      do {
         reject	=	0;
         moli->p[0].x	=	BOX[system].lx * (ran1(seed)-0.5);
         moli->p[0].y	=	BOX[system].ly * (ran1(seed)-0.5);
         moli->p[0].z	=	BOX[system].lz * (ran1(seed)-0.5);
/*
         moli->p[0].x	=	0;
         moli->p[0].y	=	0;
         moli->p[0].z	=	0;
*/

         for (molj=mol; molj<moli; molj++) {		// overlap test
            for (j=0; j<molj->nsites; j++) {
               r2	=	DistSQ(moli->p[0], molj->p[j], system);
               if (r2 < Rcut2) {			// 1.25sigma ~ 5A
                  reject	=	1;
                  break;
               }
            }
            if (reject)
               break;
         }
      } while (reject);

      for (i=1; i<moli->nsites; i++) {			// subsequent atoms
         s.d		=	l0;			// atom spherical coord.
         s.alpha	=	theta0;
//printf("moli=%d i=%d\n", moli-mol, i);
         nreject	=	0;		// # of rejection for this site
         do {
            s.beta	=	tors_angle();
            moli->p[i]	=	SiteCartesian(moli, i, s);	// cart. coord.

            reject	=	0;

            for (molj=mol; molj<=moli; molj++) {		// overlapping test
               for (j=0; j< (molj==moli ? i-3 : molj->nsites); j++) {
                  r2	=	DistSQ(moli->p[i], molj->p[j], system);
//printf("molj=%d j=%d\n", molj-mol, j);
//printf("r2=%f Rcut2=%f\n", r2, Rcut2);

                  if (r2 < Rcut2) {			// 1.25 sigma ~ 5 Angstrom
                     reject	=	1;
		     nreject	++;
	             break;
                  }
               }
	       if (reject)
		  break;
            }
         } while (nreject >= 1 && nreject <= 5);

         if (nreject > 5) {		// if continuously rejected
            i	-=	2;		// back one atom
            if (i==-1) {
               do {
                  reject	=	0;
                  moli->p[0].x	=	BOX[system].lx * (ran1(seed)-0.5);
                  moli->p[0].y	=	BOX[system].ly * (ran1(seed)-0.5);
                  moli->p[0].z	=	BOX[system].lz * (ran1(seed)-0.5);

                  for (molj=mol; molj<moli; molj++) {		// overlap test
                     for (j=0; j<molj->nsites; j++) {
                        r2	=	DistSQ(moli->p[0], molj->p[j], system);
                        if (r2 < Rcut2) {			// 1.25sigma ~ 5A
                           reject	=	1;
                           break;
                        }
                     }
                     if (reject)               break;
                  }
               } while (reject);
	       i	++;
            }
         }
      }
   }

   return;
}


void Orthorhombic(long nbox, long nmols, long nsites, long nx, long ny, long nz, double a, double b, double c, double angle)
{
   // nx, ny and nz are the number of unit cells on x, y and z direction, respectively
   // a, b and c are the dimension of unit cells on x, y and z direction, respectively
   // angle is the angle between chain surface and x direction

   long		i, j, k;
   double	bondlx, bondly, bondlz, theta;
   double	cosy, siny, angletilt;
   vector	com;

   theta	=	type[0].THETA / 180.0 * M_PI;	// bond angle
   angle	=	angle / 180.0 * M_PI;		// packing angle on x-y plane
//   angletilt	=	15.0 /180.0 * M_PI;		
   angletilt	=	0;

   bondlx	=	type[0].LSTRETCH * sin(0.5 * theta) * fabs(cos(angle)); 
   bondly	=	type[0].LSTRETCH * sin(0.5 * theta) * fabs(sin(angle));
   bondlz	=	type[0].LSTRETCH * cos(0.5 * theta);
//   c		=	bondlz * nsites/nmols;

   /* determine box size */

   for (i=0; i<nbox; i++) {
      BOX[i].lx		=	nx * a;
      BOX[i].ly		=	ny * b;
      BOX[i].lz		=	nz * c;
   }
   for (i=0; i<nmols; i++) {
      mol[i].nsites	=	nsites/nmols;		// monodisperse
   }

   /* Chains on the first z layer */

   // first 2 sites of the 2 chains in the first unit cell in the first z layer

   mol[0].p[0].x	=	0 - 0.5 * bondlx - 0.5*nx*a;
   mol[0].p[0].y	=	0 - 0.5 * bondly - 0.5*ny*b;
   mol[0].p[1].x	=	0 + 0.5 * bondlx - 0.5*nx*a;
   mol[0].p[1].y	=	0 + 0.5 * bondly - 0.5*ny*b;
   mol[1].p[0].x	=	0.5*a - 0.5 * bondlx - 0.5*nx*a;
   mol[1].p[0].y	=	0.5*b + 0.5 * bondly - 0.5*ny*b;
   mol[1].p[1].x	=	0.5*a + 0.5 * bondlx - 0.5*nx*a;
   mol[1].p[1].y	=	0.5*b - 0.5 * bondly - 0.5*ny*b;

   // first 2 sites of all other chains in the first z layer

   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         for (k=0; k<2; k++) {
            mol[(i*ny+j)*2].p[k].x	=	mol[0].p[k].x	+ i * a;
            mol[(i*ny+j)*2].p[k].y	=	mol[0].p[k].y	+ j * b;
            mol[(i*ny+j)*2+1].p[k].x	=	mol[1].p[k].x	+ i * a;
            mol[(i*ny+j)*2+1].p[k].y	=	mol[1].p[k].y	+ j * b;
         }
      }
   }

   // all other sites of all chains in the first z layer
 
   for (i=0; i<nx*ny*2; i++) {
      for (j=2; j<nsites/nmols; j++) {
         mol[i].p[j].x		=	mol[i].p[mod(j,2)].x;
         mol[i].p[j].y		=	mol[i].p[mod(j,2)].y;
      }
      for (j=0; j<nsites/nmols; j++)
         mol[i].p[j].z		=	j * bondlz - 0.5*NZ*c;
   }

   /* Multiple layers along z-direction */

   // chains on other z layers

   for (k=1; k<nz; k++) {
      for (i=0; i<2*nx*ny; i++) {			// one z-layer has 2*NX*NY chains
         for (j=0; j<nsites/nmols; j++) {
            if (!mod(nsites/nmols,2)) {			// even number sites 
               mol[k*2*nx*ny+i].p[j].x	=	mol[i].p[j].x;
               mol[k*2*nx*ny+i].p[j].y	=	mol[i].p[j].y;
               mol[k*2*nx*ny+i].p[j].z	=	mol[i].p[j].z + k*c;
            }
            else {
               mol[k*2*nx*ny+i].p[j].x	=	mol[i].p[(j-1)<0 ? j+1 : j-1].x;
               mol[k*2*nx*ny+i].p[j].y	=	mol[i].p[(j-1)<0 ? j+1 : j-1].y;
               mol[k*2*nx*ny+i].p[j].z	=	mol[i].p[j].z + k*c;
            }
         }
      }
   }

   /* Rotating */


   cosy		=	cos(angletilt);
   siny		=	sin(angletilt);

   for (i=0; i<2*nx*ny*nz; i++) {
      com	=	CenterofMass(mol+i);

      for (j=0; j<mol[i].nsites; j++) {
         mol[i].p[j].x	-=	com.x; 
         mol[i].p[j].y	-=	com.y; 
         mol[i].p[j].z	-=	com.z; 

         // rotate about y axis
         mol[i].p[j].x	=	cosy * mol[i].p[j].x + siny * mol[i].p[j].z;
         mol[i].p[j].z	=      -siny * mol[i].p[j].x + cosy * mol[i].p[j].z;

         mol[i].p[j].x	+=	com.x; 
         mol[i].p[j].y	+=	com.y; 
         mol[i].p[j].z	+=	com.z; 
      }
   } 

   /* Multiple boxes */

   for (i=1; i<nbox; i++)
      for (j=0; j<nmols/nbox; j++)
         mol[i*nmols/nbox+j]	=	mol[j];
   for (i=0; i<nmols; i++)
      mol[i].box	=	i/(nmols/nbox);

   /* molecule site type */

   for (i=0; i<nmols; i++) {
      for (j=0; j<mol[i].nsites; j++)
         mol[i].type[j]		=	0;

      if (NTYPES>=2) {
         for (j=1; j<mol[i].nsites-1; j++) 
            mol[i].type[j]	=	1;
      }
   }
/*
   for (i=0; i<nmols; i++) {			// for melting simulation
      if ( mod(i/(nx*ny), 2)) {
         for (j=0; j<mol[i].nsites; j++)
            mol[i].type[j]	+=	2;
         mol[i].fix	=	1;		// not used yet
      }
      else
         mol[i].fix	=	0;
   }
*/
   for (i=0; i<nmols; i++)
      if (mol[i].p[0].y > 0) {		// interface along x direction
//	if (mol[i].p[0].x > 0) {		// interface along y direction
//      if ( (mol[i].p[0].y+mol[i].p[1].y) >  (mol[i].p[0].x+mol[i].p[1].x) * b/a) {	// interface along y=b/a*x direction
	 for (j=0; j<mol[i].nsites; j++)
	     mol[i].type[j]	+=	2;
      }
/*
   for (i=0; i<nmols; i++) {			// for melting, create a cylinder
      if ( mol[i].p[0].x * mol[i].p[0].x + mol[i].p[0].y * mol[i].p[0].y < 169 ) {
         for (j=0; j<mol[i].nsites; j++) {
               mol[i].type[j]	+=	2;
         }
      }
   }
*/
   return;
}


void sc_lattice(long N, double rho, long PBC)		// total particle number N = NC * NC * NC
{							// box dimension L
   long		i, NC;
   double	L, cell;

   for (i=0; i<N; i++)
      mol[i].nsites	=	1;			// only one atom in one molecule

   if (PBC==1)
      NC	=	(int) rint(pow(N, 1.0/3));

   L	=	pow(N/rho, 1.0/3);			// cubic box length

   BOX[0].lx	=	L;				// only for one box now 4/18/08
   BOX[0].ly	=	L;
   BOX[0].lz	=	L;
   BOX[0].lbox	=	L;

   if (PBC==1 && N!=NC*NC*NC)
      Exit("builder", "sc_lattice", "number of particle incomp. with sc lattice.\n");

   if (PBC==1) {
      cell	=	L/NC;

      for (i=0; i<N; i++) {
         mol[i].p[0].x	=	(mod(i, NC) - NC/2) * cell;			//or (int)mod(i,M)/1
         mol[i].p[0].y	=	((int)(mod(i, NC*NC)/NC) - NC/2) * cell;	//or (int)mod(i,M*M)/M
	 mol[i].p[0].z	=	((int)(i/(NC*NC)) - NC/2) * cell;		//or (int)mod(i,M*M*M)/(M*M)
      }
   }
}


void fcc_lattice(long N, double rho, long PBC)			//reference: F.23 of Allen and Tildesley
{							// now only for one box 4/8/2008
   long		i, j, k, m, ref;
   double	L, cell, cell2;
   long		NC;

   for (i=0; i<N; i++) 
      mol[i].nsites	=	1;				// only one atom on one chain

   if (PBC==1)
      NC	=	(int) rint(pow(N/4, 1.0/3));		// rint: nearest integer value
   else if (PBC==2)
      NC	=	(int) rint(pow(N/16, 1.0/3));

   L		=	pow(N/rho, 1.0/3);			// cubic box length

   BOX[0].lx	=	L;
   BOX[0].ly	=	L;
   BOX[0].lz	=	L;
   BOX[0].lbox	=	L;

   if ( (PBC==1 && N!=4*NC*NC*NC) || (PBC==2 && N!=16*NC*NC*NC))
      Exit("builder", "fcc_lattice", "number of particles incomp. with fcc lattice.\n");

   if (PBC==1) {
      cell	=	L / NC;				// the reason to use plus sign here is to adapt with our
      cell2	=	0.5 * cell;			// p.b.c., that is if X<-LBOX/2, then X+=LBOX
							// if our pbc is if X<=-LBOX/2, then X+=LBOX, then minus sign
							// in another words, our box is [-L, +L)
      mol[0].p[0].x	=	0;	mol[0].p[0].y	=	0;	mol[0].p[0].z	=	0;
      mol[1].p[0].x	=	cell2;	mol[1].p[0].y	=	cell2;	mol[1].p[0].z	=	0;
      mol[2].p[0].x	=	0;	mol[2].p[0].y	=	cell2;	mol[2].p[0].z	=	cell2;
      mol[3].p[0].x	=	cell2;	mol[3].p[0].y	=	0;	mol[3].p[0].z	=	cell2;

      m	=	0;

      for (i=0; i<NC; i++) {
         for (j=0; j<NC; j++) {
	    for (k=0; k<NC; k++) {
	       for (ref=0; ref<4; ref++) {
		  mol[ref+m].p[0].x	=	mol[ref].p[0].x	+	cell*k;
		  mol[ref+m].p[0].y	=	mol[ref].p[0].y	+	cell*j;
		  mol[ref+m].p[0].z	=	mol[ref].p[0].z	+	cell*i;
	       }
	       m	+=	4;
	    }
	 }
      }
      for (i=0; i<N; i++) {
	  mol[i].p[0].x	-=	(double)L/2;		//same reason to use += instead of -= as mentioned above
	  mol[i].p[0].y	-=	(double)L/2;
	  mol[i].p[0].z	-=	(double)L/2;
      }
   }
   if (PBC==2) {					//truncated octahedron periodic boundary condition
      printf("Set up fcc lattice in truncated octahedron box.\n");

      cell	=	L / (2*NC);
      cell2	=	0.5 * cell;

      mol[0].p[0].x	=	0;	mol[0].p[0].y	=	0;	mol[0].p[0].z	=	0;
      mol[1].p[0].x	=	cell2;	mol[1].p[0].y	=	cell2;	mol[1].p[0].z	=	0;
      mol[2].p[0].x	=	0;	mol[2].p[0].y	=	cell2;	mol[2].p[0].z	=	cell2;
      mol[3].p[0].x	=	cell2;	mol[3].p[0].y	=	0;	mol[3].p[0].z	=	cell2;

      m	=	0;

      for (i=0; i<NC; i++) {			// (2*NC, 2*NC, NC)
         for (j=0; j<2*NC; j++) {
	    for (k=0; k<2*NC; k++) {
               for (ref=0; ref<4; ref++) {
		  mol[ref+m].p[0].x	=	mol[ref].p[0].x	+	cell*k;
		  mol[ref+m].p[0].y	=	mol[ref].p[0].y	+	cell*j;
		  mol[ref+m].p[0].z	=	mol[ref].p[0].z	+	cell*i;
	       }
	       m	+=	4;
	    }
	 }
      }
     
      for (i=0; i<N; i++) {
	 mol[i].p[0].x	-=	0.5 * L;
         mol[i].p[0].y	-=	0.5 * L;

	 if ( (fabs(mol[i].p[0].x) + fabs(mol[i].p[0].y) + fabs(mol[i].p[0].z)) > 0.75*L) {
            mol[i].p[0].x	+=	((mol[i].p[0].x >=0) ? -0.5 : 0.5) * L;
            mol[i].p[0].y	+=	((mol[i].p[0].y >=0) ? -0.5 : 0.5) * L;
            mol[i].p[0].z	+=	((mol[i].p[0].z >=0) ? -0.5 : 0.5) * L;
	 }
	 else if ( (fabs(mol[i].p[0].x) + fabs(mol[i].p[0].y) + fabs(mol[i].p[0].z) == 0.75 * L) && mol[i].p[0].z>=0) {
	    mol[i].p[0].x	+=	((mol[i].p[0].x >= 0) ? -0.5 : 0.5) * L;
	    mol[i].p[0].y	+=	((mol[i].p[0].y >= 0) ? -0.5 : 0.5) * L;
	    mol[i].p[0].z	-=	0.5 * L;
	 }
      }
   }
}


void BuilderCheckSetup(char * argv[])
{
   long	error=0;

   if (!strcmp(argv[1], "chain")) {
      if (NMOLS/NSYSTEMS != NX*NY*NZ*2) {
         printf("NMOLS mismatch!\n");
         exit(-1);
      }
   }
   if (NMOLS > MAXNMOLS) 		{printf("NMOLS>MAXNMOLS\n");	error=1;}
   if (NSITES/NMOLS > MAXNMOLSITES) 	{printf("NMOLSITES error\n");	error=1;}

   if (error) {
      printf("setup parameter error\n");
      exit(-1);
   }
}


int main(int argc, char * argv[])
{
   
   if (argc<2) {
      printf("builder (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tbuilder structure\n\n");
      printf("\tstructure:	ljfcc\n");
      printf("\t	 	ljsc\n");
      printf("\t	 	chain\n");
      printf("\t	 	chainrandom\n");
      printf("Notes:\n");
      printf("\t* require in file\n\n");
      exit(1);
   }

   tim=(int *)malloc(sizeof(int));     	//random number generator
   seed=(long *)malloc(sizeof(long));
   *tim=(int)time(NULL);
   *seed= -1*(*tim);           		//seed must start out as a negative long

//   InitMols(MAXNMOLS, MAXNMOLS);	// allocate max. memory for molecules
   InitMols(MAXNMOLS, 0);

   BuilderReadSetup();
   BuilderCheckSetup(argv);
   BuilderPrintSetup();

   if (!strcmp(argv[1], "chain")) 
      Orthorhombic(NBOX, NMOLS, NSITES, NX, NY, NZ, cella, cellb, cellc, angle);	
   else if (!strcmp(argv[1], "chainrandom"))
      Amorphous(NBOX, NMOLS, NSITES);
   else if (!strcmp(argv[1], "ljfcc"))
      fcc_lattice(NMOLS, Rho, PBC);
   else if (!strcmp(argv[1], "ljsc"))
      sc_lattice(NMOLS, Rho, PBC);



   CalcUnits(0);			// unit = 1, unit needed for Write_Conf() 
   //printf("%f  %f  %f\n", mol[0].p[0].x, mol[0].p[0].y, mol[0].p[0].z);
   Write_Conf(-1);
   return 0;
}

                                                                                                                                                                             src/cluster.c                                                                                       0000600 0143352 0000144 00000001570 10714730275 012712  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  #include "header.h"

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

                                                                                                                                        src/conf2car.c                                                                                      0000600 0143352 0000144 00000061106 11602261057 012721  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	conf2car.c
	author:		Peng Yi at MIT
	date:		November 19, 2007
	purpose:	convert configuration files to .car files
	note:		require setup file
			shiftbox() added June 25, 2009
*/

#define __MAIN_PROGRAM
#include "header.h"

#define SIZECAP	1

char *Element(long t, char *s)
{
   double	mass = type[t].M;

   s[0]		=	0;

   if (fabs(mass-14) <= 1e-6)		strcpy(s, "C");
   else if (fabs(mass-15) <= 1e-6)	strcpy(s, "C");
   else if (fabs(mass-1.01) <= 1e-6)	strcpy(s, "H");
   else if (fabs(mass-28.086) <= 1e-6)	strcpy(s, "Si");
   else if (fabs(mass-26.982) <= 1e-6)	strcpy(s, "Al");
   else if (fabs(mass-16) <= 1e-6)	strcpy(s, "O");

   strcpy(s, "C");
   return	s;
}


void doublesize()		// double the system size
{ 		
   long		i;
   molstruct	*moli;
			
   printf("%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      printf("%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      printf("%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
         
      for (i=0; i<moli->nsites; i++) 
         printf("%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      printf("%d\t%d\t%d\n", moli-mol+NMOLS, moli->box, moli->nsites);

      for (i=0; i<moli->nsites; i++) 
         printf("%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z + BOX[0].lz);
   }
}


void types()		// temporary code to prepare a system with 4 types
{
   long		i;
   molstruct	*moli;

   printf("%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      printf("%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

   for (moli=0; moli<mol+NMOLS; moli++) {
      printf("%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
         
      for (i=0; i<moli->nsites; i++) {
         if ( mod((moli-mol)/6, 2)) 
            printf("%d\t%f\t%f\t%f\n", moli->type[i]+2, moli->p[i].x, moli->p[i].y, moli->p[i].z);
         else
            printf("%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
      }
   }
}

/*******************************************************************/

void shift()				// move everything to the central box and 
{					// make the largest nucleus at the center  (5/20/08)
   FILE		*fPtr;
   vector	center, com, 
		rA,			// center of one chain in nucleus
		rO,			// center of nucleus 
		rBA, 			// rB-rA, B and A belong to same nucleus
		rOA;			// rO-rA
   molstruct	*moli;
   long		i, system=0, id, n;

   if (!(fPtr=fopen("shifted", "w"))) {
      printf("shifted file failed to open\n");
      exit(1);
   }

   // step1
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         if (sizeofnucl[moli->nuclid[i]]==MAXSIZE[0]){ 	// remember that sizeofnucl
            id	=	moli->nuclid[i];		// points to sizeofnuclp2 now
            rA	=	moli->p[i];
            break;
         }
      }
   }

   //step 2
   n=0;			// # of sites in this nucleus
   V_Null(&rBA);
   V_Null(&rOA);
   
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         if (moli->nuclid[i] == id) {
            n++;
            com	=	moli->p[i];
            rBA	=	V_Subtr(&com, &rA);
            rOA	=	V_Add(&rOA, &rBA);
         }}}
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);

   //step3
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->p[i]	=	V_Subtr(moli->p+i, &rO);
      }
   }

   // step 4
   for (moli=mol; moli<mol+NMOLS; moli++) {
      MolInBox2(moli);
   }

   // output shifted configuration file

   fprintf(fPtr, "TIMESTEP\t%d\n", TIMESTEP); 
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx * unit.LENGTH, BOX[i].ly * unit.LENGTH, BOX[i].lz * unit.LENGTH);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
   }
   fclose(fPtr);
}

/**********************************************************************/

void shift1()				// move everything to the central box and 
{					// make the largest nucleus at the center  (5/20/08)
   FILE		*fPtr;
   vector	center, com, 
		rA,			// center of one chain in nucleus
		rO,			// center of nucleus 
		rBA, 			// rB-rA, B and A belong to same nucleus
		rOA;			// rO-rA
   molstruct	*moli;
   long		i, system=0, id, n;

   double	lx[MAXNSYSTEMS], 
		ly[MAXNSYSTEMS], 
		lz[MAXNSYSTEMS];// shrink box to just fit the biggest nucleus in
   long		inthebox, nmol, nsites;
   char		dummy[80];
   
   if (!(fPtr=fopen("shifted", "w")))
      Exit("conf2car", "shift", "shifted failed to open");


   // Step 1: find one chain A belonging to the largest nucleus as a reference

   for (moli=mol; moli<mol+NMOLS; moli++)
      if (sizeofnucl[moli->nuclid[0]]==MAXSIZE[0]) {
	  id	=	moli->nuclid[0];		// nucleus id
	  rA	=	CenterofMass(moli); 		// take one chain as reference point
          break;
      }

   // Step 2: calc. the shift of center of nucleus to this chain

   n	=	0;					// # of chains in this nucleus
   V_Null(&rBA);
   V_Null(&rOA);

   for (moli=mol; moli<mol+NMOLS; moli++) 
      if (moli->nuclid[0] == id) {
	 n	++;
	 com	=	CenterofMass(moli);
	 rBA	=	V_Subtr(&com, &rA);
	 rBA	=	MapInBox2(&rBA, PBC, system);
	 rOA	=	V_Add(&rOA, &rBA);
      }

   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);			// center of nucleus

/*
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (sizeofnucl[moli->nuclid[0]]==MAXSIZE[0]) {
         id	=	moli->nuclid[0];
         rO	=	CenterofNucleus(moli->nuclid[0], moli);
	 break;
      }
*/
   // Step 3: every particle shifts rO so that the center of nucleus becomes center of box

   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         moli->p[i]	=	V_Subtr(moli->p+i, &rO);

   // Step 4: map all in central box

   for (moli=mol; moli<mol+NMOLS; moli++) 
      MolInBox2(moli);

// output the minimum box that contains the biggest nucleus 11/13/08

   for (n=0; n<NSYSTEMS; n++) {
      lx[n]	=	0;
      ly[n]	=	0;
      lz[n]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->nuclid[0]==id) {
	 for (i=0; i<moli->nsites; i++) {
            if (fabs(moli->p[i].x) > lx[moli->box])
	       lx[moli->box]	=	fabs(moli->p[i].x);
            if (fabs(moli->p[i].y) > ly[moli->box])
	       ly[moli->box]	=	fabs(moli->p[i].y);
            if (fabs(moli->p[i].z) > lz[moli->box])
	       lz[moli->box]	=	fabs(moli->p[i].z);
         }
      }
   }

   for (i=0; i<40; i++)
      fprintf(fPtr, " ");
   fprintf(fPtr, "\n");
   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", lx[i] *2 * unit.LENGTH, ly[i]*2 * unit.LENGTH, lz[i]*2 * unit.LENGTH);

   n=0;   nsites=0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      inthebox	=	1;
/*
      com	=	CenterofMass(moli);
      if (fabs(com.x) > lx[0] || fabs(com.y) > ly[0] || fabs(com.z) > lz[0]) 
	 inthebox	=	0;
*/
      for (i=0; i<moli->nsites; i++) {
         if (fabs(moli->p[i].x) > lx[0] || fabs(moli->p[i].y) > ly[0] || fabs(moli->p[i].z) > lz[0]) {
            inthebox	=	0;
            break;
	 }
      }

      if (inthebox) {
         fprintf(fPtr, "%d\t%d\t%d\n", n, moli->box, moli->nsites);
         for (i=0; i<moli->nsites; i++) 
            fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
         n	++;
         nsites	+=	moli->nsites;
      }
   }
   rewind(fPtr);
   sprintf(dummy, "%d\t%d\t%d", NSYSTEMS, n, nsites);
   fprintf(fPtr, "%s", dummy);
   for (i=0; i<40-strlen(dummy); i++)
      fprintf(fPtr, " ");

   fseek(fPtr, 0, SEEK_END);			// go to the end of the file

   /* output shifted conf file */
 
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx * unit.LENGTH, BOX[i].ly * unit.LENGTH, BOX[i].lz * unit.LENGTH);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
   }
   fclose(fPtr);
}

/***************************************************************************/

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   long		i, system, n, j;
   char		s[255], ff[255], filename1[255], filename2[255];
   molstruct	*moli;
   vector	com;		// center of mass
   FILE		*fPtr, *fcar, *fpdb;
   long		nsite, nbond, nangle, ndihedral;
   //long		realXtal;
   char		atomname;

long	drawmol;
matrix	Mj, Mk;
vector eigj, eigk, dp, nucleusshape;
double	Rg2;
beadstruct	*nucleus;		// for grouping beads in biggest nucleus

   tim=(int *)malloc(sizeof(int));     	//random number generator
   seed=(long *)malloc(sizeof(long));
   *tim=(int)time(NULL);
   *seed= -1*(*tim);           		//seed must start out as a negative long

   if (argc<2) {
      printf("conf2car (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tconf2car filename\n\n");
      printf("Notes:\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   strcpy(filename1, argv[1]);
   strcat(filename1, ".car");
   if (! (fPtr = fopen(filename1, "w")))      exit(1);

   strcpy(filename2, argv[1]);
   strcat(filename2, ".pdb");
   if (! (fpdb = fopen(filename2, "w")))      exit(1);

   freopen("conf2car.out", "w", stdout);	// redirect standard output stream to a file

   InitMols(MAXNMOLS, MAXNMOLS);
   GetSetup(argv);

/*
   // Direct conversion to a pdb file for visualization
   system=0;
   Read_Conf(argv[1]);			// read in WITHOUT unit conversion
   fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
   fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
	BOX[system].lx, BOX[system].ly, BOX[system].lz,
	90.0, 90.0, 90.0);

   n		=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {

         for (i=0; i<moli->nsites; i++) {
               atomname	=	'M';
	    n	++;

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
		moli->p[i].x , moli->p[i].y , moli->p[i].z);
	    fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
	    fprintf(fpdb, "\n"); 
         } 
      }   
   }
   fprintf(fpdb, "END\n");
*/   

   InitUnits();
   GetCoordinates(argv[1]);		// include unit conversion
   //Read_Conf(argv[1]);		// read in configuration file
   //CoorSI2System();			// convert coordinates from SI to system units
   InitForcefield();

#ifdef CELL_LIST
   CL_Init();				// Cell list should be ready before Verlet list
   CL_Build();				// needed for Calc_Qlm
#endif	/* CELL_LIST */

   for (moli=mol; moli<mol+NMOLS; moli++) { 		
      for (i=0; i<moli->nsites; i++)  {
         moli->flags[i]		=	1;		// activate all the sites on this processor
         moli->parent[i]	=	i-1;		// initialize parent site
      }
      moli->flip		=	0;		// flip to the original direction
      moli->origin		=	CenterofMass(moli);
   }

   nucleus	=	(beadstruct *)calloc(NSITES, sizeof(beadstruct));
   if (nucleus == NULL) {
      printf("nucleus allocation failed!\n");
      exit(-1);
   }

   //////////////////////
   /* Do Some Analysis */
   //////////////////////

   InitSample();
   SampleP2All();

   SampleSpherical();
   Dist_Spherical();

// Find_Nuclei(1);
   Find_Nuclei_p2(1);

   sizeofnucl	=	sizeofnuclp2;
   sizedist	=	sizedistp2;
//shift();

   //Calc_Qlm(l_of_Ylm);
   S_PrintAll();
   CalcV();

   //realXtal	=	Xtal[0];
   //for (i=1; i<=5; i++)
      //realXtal	-=	sizedist[i]*i;

   for (i=0; i<argc; i++) {
      printf("%s ", argv[i]);
   }
   printf("\n");

   printf("**********conf2car Some statistics**********\n");
   for (system=0; system <NSYSTEMS; system++) {
      printf("Vtot/NSITES=%f\n", v[system].tot/NSites[system]);
      printf("Rp = %f\n", Rp);
      printf("Rconn = %f\n", Rconn);
      printf("P2 = %f\n", P2[system]);
      printf("P2m = %f\n", P2M[system]);
      printf("P2z = %f\n", P2z[system]);
      printf("nmax = %d\n", nmax[system][0]);
      printf("2ndnmax = %d\n", nmax[system][1]);
      printf("3rdnmax = %d\n", nmax[system][2]);
      printf("realXtal = %d\n", realXtal[system]);
      printf("Q6 = %f\n", Q6[system]);
      printf("Q4 = %f\n", Q4[system]);
   }
   system	=	0;		// reset system=0

   /* Grouping the segments belonging to the Biggest nucleus */

   n	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system == moli->box) {
         for (i=0; i<moli->nsites; i++) {
            if (sizeofnucl[moli->nuclid[i]] == nmax[system][0]) {
               nucleus[n].moli	=	moli;
	       nucleus[n].site	=	i;
	       n	++;
	    }
         }
      }
   }
// shiftbox(system, nucleus, n);	// shift the biggest nucleus to the box center

   // doublesize();
   // types();

   system	=	0;			// only one system for now

shift();

   //////////////////////////////////////////////
   /* Print out .pdb file for the whole system */
   //////////////////////////////////////////////

   fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
   fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, BOX[system].lz * unit.LENGTH,
	90.0, 90.0, 90.0);

   n		=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {
         //MolInBox2(moli);

         for (i=0; i<moli->nsites; i++) {
            if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])	// nuclid index starts from 1
               atomname	=	'N';				// blue in VMD
	    else if (sizeofnucl[moli->nuclid[i]] >= SIZECAP)
               atomname	=	'O';				// red in VMD
            else
               atomname	=	'C';				// cyan in VMD

	    n	++;

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
		moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
		moli->p[i].z * unit.LENGTH);
	    fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
	    fprintf(fpdb, "\n"); 
         } 
      }   
   }
   fprintf(fpdb, "END\n");

   //////////////////////////////////////////////
   /* Print out .car file for the WHOLE SYSTEM */
   //////////////////////////////////////////////

   fprintf(fPtr, "!BIOSYM archive 3\n");
   if (1==PBC) {
      fprintf(fPtr, "PBC=ON\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
      fprintf(fPtr, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, 
	BOX[system].lz * unit.LENGTH, 90.0, 90.0, 90.0);
   }
   else if (0==PBC) {
      fprintf(fPtr, "PBC=OFF\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
   }

   n		=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {

         //MolInBox2(moli);
         for (i=0; i<moli->nsites; i++) {
/*
            if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])	// nuclid index starts from 1
               sprintf(s, "N%d", n++);
	    else if (sizeofnucl[moli->nuclid[i]] >= SIZECAP)
	       sprintf(s, "O%d", n++);
            else
	       sprintf(s, "C%d", n++);
*/
            if (moli->nuclid[i] >= 0) {
               if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])// nuclid index starts from 1
                  sprintf(s, "N%d", n++);
     	       else if (sizeofnucl[moli->nuclid[i]] >= SIZECAP)
	          sprintf(s, "O%d", n++);
            }
            else
	       sprintf(s, "C%d", n++);


            //if (sizeofnucl[moli->nuclid[i]] > SIZECAP) {
            fprintf(fPtr, "%-6.6s ", s);
            sprintf(s, "M%d", moli-mol);
            fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
			moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
			moli->p[i].z * unit.LENGTH);
            strcpy(ff, "O");
            fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[i], s));
            //}
         } 
      }   
   }
   fprintf(fPtr, "end\nend\n");

   //////////////////////////////////////////////////////////////////////////////
   /* Print out .car file for the Chains that contain the Biggest Nucleus only */
   //////////////////////////////////////////////////////////////////////////////

   fprintf(fPtr, "**********The Chains Containing Biggest Nucleus Only**********\n");
   fprintf(fPtr, "!BIOSYM archive 3\n");
   if (1==PBC) {
      fprintf(fPtr, "PBC=ON\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
      fprintf(fPtr, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, 
	BOX[system].lz * unit.LENGTH, 90.0, 90.0, 90.0);
   }
   else if (0==PBC) {
      fprintf(fPtr, "PBC=OFF\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
   }
   n	=	0;
/*
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {

         MolInBox2(moli);
         for (i=0; i<moli->nsites; i++) {
            if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {	// note nuclid index starts from 1
               sprintf(s, "O%d", n++);

               fprintf(fPtr, "%-6.6s ", s);
               sprintf(s, "M%d", moli-mol);
               fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
			moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
			moli->p[i].z * unit.LENGTH);
               strcpy(ff, "O");
               fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[i], s));
            }
         } 
      }   
   }
*/

   for (moli=mol; moli<mol+NMOLS; moli++) {		// for nucleus def based on p2
      if (system==moli->box) {
         MolInBox2(moli);

         drawmol	=	0;
         for (i=0; i<moli->nsites; i++) {
            if (sizeofnucl[moli->nuclid[i]] == nmax[system][0]) {
//            if (sizeofnucl[moli->nuclid[i]] == nmax[system][0] 
//		|| sizeofnucl[moli->nuclid[i]] == nmax[system][1]) {
               drawmol	=	1;
               break;
            }
         }
         if (drawmol) {
	    for (i=0; i<moli->nsites; i++) {
               if (sizeofnucl[moli->nuclid[i]] == nmax[system][0])
                  sprintf(s, "N%d", n++);
//               else if (sizeofnucl[moli->nuclid[i]] == nmax[system][1])
//                  sprintf(s, "O%d", n++);
               else
                  sprintf(s, "C%d", n++);

               fprintf(fPtr, "%-6.6s ", s);
               sprintf(s, "M%d", moli-mol);

               //moli->p[i]	=	MapInBox2(moli->p+i, PBC, system);

               fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
			moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
			moli->p[i].z * unit.LENGTH);
               strcpy(ff, "O");
               fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, 
			Element(moli->type[i], s));
            }
         }
      }
   }
   fprintf(fPtr, "end\nend\n");

   //////////////////////////////////////////////////////
   /* Print out .car file for the Biggest Nucleus only */
   //////////////////////////////////////////////////////

   fprintf(fPtr, "**********The Biggest Nucleus Only**********\n");
   fprintf(fPtr, "!BIOSYM archive 3\n");
   if (1==PBC) {
      fprintf(fPtr, "PBC=ON\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
      fprintf(fPtr, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, 
	BOX[system].lz * unit.LENGTH, 90.0, 90.0, 90.0);
   }
   else if (0==PBC) {
      fprintf(fPtr, "PBC=OFF\n\n");
      fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
   }
   n	=	0;

   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {
         MolInBox2(moli);

	 for (i=0; i<moli->nsites; i++) {
            if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {
               sprintf(s, "N%d", n++);

               fprintf(fPtr, "%-6.6s ", s);
               sprintf(s, "M%d", moli-mol);
               fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
			moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
			moli->p[i].z * unit.LENGTH);
               strcpy(ff, "O");
               fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, 
			Element(moli->type[i], s));
	    }
         }
      }
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {
         MolInBox2(moli);

	 for (i=0; i<moli->nsites; i++) {
            if (sizeofnucl[moli->nuclid[i]] == nmax[system][1]) {
               sprintf(s, "O%d", n++);

               fprintf(fPtr, "%-6.6s ", s);
               sprintf(s, "M%d", moli-mol);
               fprintf(fPtr, "%14.8g %14.8g %14.8g ", 
			moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, 
			moli->p[i].z * unit.LENGTH);
               strcpy(ff, "O");
               fprintf(fPtr, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, 
			Element(moli->type[i], s));
	    }
         }
      }
   }
   fprintf(fPtr, "end\nend\n");

   ////////////////////////////////////////////////////////////
   /* print out center of mass of chains in .car file format */
   ////////////////////////////////////////////////////////////
   
/*
   n		=	0;
   printf("**********Center of mass of chains**********\n");
   printf("!BIOSYM archive 3\n");
   if (1==PBC) {
      printf("PBC=ON\n\n");
      printf("!DATE %s", asctime(localtime(&t)));
      printf("PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
	BOX[system].lx * unit.LENGTH, BOX[system].ly * unit.LENGTH, BOX[system].lz * unit.LENGTH, 90.0, 90.0, 90.0);
   }
   else {
      printf("PBC=OFF\n\n");
      printf("!DATE %s", asctime(localtime(&t)));
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {

         if (sizeofnucl[moli->nuclid[0]] > SIZECAP) {

         com	=	CenterofMass(moli);
         sprintf(s, "M%d", n++);
         printf("%-5.5s ", s);
         sprintf(s, "M%d", moli-mol);
         printf("%14.8g %14.8g %14.8g ", com.x * unit.LENGTH, com.y * unit.LENGTH, com.z * unit.LENGTH);
         strcpy(ff, "O");
         printf("%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[0], s));

         }
      }   
   }
   printf("end\nend\n");
*/

   ////////////////////////////////////////////////////////
   /* Shift the largest nucleus to the center of the box */
   ////////////////////////////////////////////////////////a

   //shift();

   /////////////////
   /* close files */
   /////////////////

   fclose(fPtr);
   fflush(stdout);
   fclose(stdout);
   free(nucleus);
   return	0;
}


void shiftbox(long system, beadstruct *nucleus, long nsites)	// move the biggest nucleus
{								// to the center of the box
   molstruct	*moli;
   long		i, n, site;
   vector	rA, rO, rBA, rOA;

   // step 1: find one segment that belongs to the biggest nucleus
   
   rA	=	nucleus[0].moli->p[nucleus[0].site];
  
   // step 2: calc. the shift of com of nucleus to this segment
   
   V_Null(&rBA);
   V_Null(&rOA);

   for (n=0; n<nsites; n++) { 
      moli	=	nucleus[n].moli;
      site	=	nucleus[n].site;
      rBA	=	V_Subtr(moli->p+site, &rA);
      rBA	=	MapInBox2(&rBA, PBC, system);
      rOA	=	V_Add(&rOA, &rBA);
   }
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);

   // step 3: every segment shift and move to central box
  
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (site=0; site<moli->nsites; site++) {
         moli->p[site]	=	V_Subtr(moli->p+site, &rO);
         moli->p[site]	=	MapInBox2(moli->p+site, PBC, system);
      } 
   }
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                          src/conf2DL_POLY.c                                                                                  0000600 0143352 0000144 00000020370 11320471376 013320  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	conf2DL_POLY.c
	author:		Peng Yi at MIT
	date:		July 30, 2009
	purpose:	convert configuration files to DL_POLY input
	note:		require setup file
*/

#define __MAIN_PROGRAM
#include "header.h"

#define VERSION		"Jul_30_2009"

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
      printf("confDL_POLY (c) 2009 by Peng Yi at MIT\n\n");
      printf("convert configuration file to DL_POLY input files.\n");
      printf("Usage:\n");
      printf("\tconf2DL_POLY filename\n\n");
      printf("Notes:\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   strcpy(filename1, argv[1]);
   strcat(filename1, ".CONFIG");
   if (! (fcfg = fopen(filename1, "w")))      exit(1);

   strcpy(filename2, argv[1]);
   strcat(filename2, ".FIELD");
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

   // generate DL_POLY CONFIG file, fixed-formatted

   fprintf(fcfg, "DL_POLY CONFIG file created from %s on %s", argv[1], asctime(localtime(&t)));
   fprintf(fcfg, "%10d%10d%10d%20f\n", 2, 3, NSITES, 0);
   fprintf(fcfg, "%20f%20f%20f\n", BOX[system].lx, 0.0, 0.0); 
   fprintf(fcfg, "%20f%20f%20f\n", 0.0, BOX[system].ly, 0.0); 
   fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, BOX[system].lz); 

   n	=	1;					// index starts from 1
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->type[0]==0) {				// printout all molecules of type 1
         ntype[1]	++;
         for (i=0; i<moli->nsites; i++) {
            fprintf(fcfg, "%-8s%10d\n", typename[moli->type[i]+1], n);
            fprintf(fcfg, "%20f%20f%20f\n", moli->p[i].x, moli->p[i].y, moli->p[i].z);
            fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, 0.0);
            fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, 0.0);
            n	++;
         }
      }
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->type[0]==2) {				// printout all molecules of type 2
         ntype[2]	++;
         for (i=0; i<moli->nsites; i++) {
            fprintf(fcfg, "%-8s%10d\n", typename[moli->type[i]+1], n);
            fprintf(fcfg, "%20f%20f%20f\n", moli->p[i].x, moli->p[i].y, moli->p[i].z);
            fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, 0.0);
            fprintf(fcfg, "%20f%20f%20f\n", 0.0, 0.0, 0.0);
            n	++;
         }
      }
   }
   fclose(fcfg);

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
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[1], 15.0, 0.0, 1);	 	// head bead
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[2], 14.0, 0.0, NSITES/NMOLS-2); 	// middle beads
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[1], 15.0, 0.0, 1); 		// tail bead

   fprintf(ffld, "BONDS	%d\n", NSITES/NMOLS-1);		// # of bonds in each type 1 molecule
   for (i=1; i<=NSITES/NMOLS-1; i++) {
      fprintf(ffld, "harm	%d %d %f %f\n", i, i+1, kl, l0);
   }

   fprintf(ffld, "ANGLES %d\n", NSITES/NMOLS-2);	// # of angles in each type 1 molecule
   for (i=1; i<=NSITES/NMOLS-2; i++) {
      fprintf(ffld, "harm	%d %d %d %f %f\n", i, i+1, i+2, ktheta, theta0);
   }
   
   fprintf(ffld, "DIHEDRALS %d\n", NSITES/NMOLS-3);	// # of dihes in each type 1 molecule
   for (i=1; i<=NSITES/NMOLS-3; i++) {
      fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", i, i+1, i+2, i+3, k0, k1, k2, k3);
   }

   fprintf(ffld, "FINISH\n");				// finish type 1 molecule

   /* Register molecule type 2 */

   fprintf(ffld, "octane2\n");				// name of type 2
   fprintf(ffld, "NUMMOLS %d\n", ntype[2]);		// # of type 2 molecules

   fprintf(ffld, "ATOMS %d\n", NSITES/NMOLS);		// # of atoms in each molecule
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[3], 15.0, 0.0, 1);	 
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[4], 14.0, 0.0, NSITES/NMOLS-2); 	// 6 middle beads
   fprintf(ffld, "%-8s%12f%12f%8d\n", typename[3], 15.0, 0.0, 1); 
   fprintf(ffld, "BONDS	%d\n", NSITES/NMOLS-1);		// # of bonds in each molecule
   for (i=1; i<=NSITES/NMOLS-1; i++) {
      fprintf(ffld, "harm	%d %d %f %f\n", i, i+1, kl, l0);
   }

   fprintf(ffld, "ANGLES %d\n", NSITES/NMOLS-2);	// # of angles in each molecule
   for (i=1; i<=NSITES/NMOLS-2; i++) {
      fprintf(ffld, "harm	%d %d %d %f %f\n", i, i+1, i+2, ktheta, theta0);
   }

   fprintf(ffld, "DIHEDRALS %d\n", NSITES/NMOLS-3);	// # of dihes in each molecule
   for (i=1; i<=NSITES/NMOLS-3; i++) {
      fprintf(ffld, "opls	%d %d %d %d %f %f %f %f\n", i, i+1, i+2, i+3, k0, k1, k2, k3);
   }

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
                                                                                                                                                                                                                                                                        src/conf2gro.c                                                                                      0000600 0143352 0000144 00000022600 11255725214 012743  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
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
                                                                                                                                src/conf2lammps.c                                                                                   0000600 0143352 0000144 00000010452 11531771745 013456  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
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
         for (moli=mol; moli<mol+NMOLS; moli++) {
            for (i=0; i<moli->nsites; i++) {
               nsite	++;
               fprintf(fPtr, "%d %d %d %f %f %f\n", nsite, moli-mol+1, moli->type[i]+1, moli->p[i].x, moli->p[i].y, moli->p[i].z);
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
                                                                                                                                                                                                                      src/config.c                                                                                        0000600 0143352 0000144 00000011711 10714730275 012474  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
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
                                                       src/distributions.c                                                                                 0000600 0143352 0000144 00000031247 10745445712 014142  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	distributions.c
    author:	Pieter J. in 't Veld for UT at Austin
    date:	January 27, June 11, 1998, April 1, 1999.
    purpose:	Definition of a multi-dimensional self-defining distribution
		formalism

    Notes:
      19980127	Creation date
      19980611	- Adaption to pointers
		- Alteration of energy functions
		- Alteration of distribution storage functions
		- Transfer to module style
      19990401	Change in finding minimum: Newton-like root finding algorithm
		through local linearization of the scaled force (see notes).
      20040306	Addition of multi-dimensionality

      20071014  Perfect multi-dimensionality, basically change InitDist() --Peng
*/
#define __DISTRIBUTIONS_MODULE
#include "distributions.h"

void D_Allocate(diststruct *d, long n)
{
  long			i, i0 = 0, dnbins = 0;
  long			*newbin;
  double		*newdata, *newcweight;
  diststruct		*newdist;

  if (!d->nbins) 
    d->startbin		= n;
  n			-= d->startbin;
  if (!(dnbins = n<0 ? i0 = -n : n>=d->nbins ? n+1-d->nbins : 0)) 
    return;
  if (d->level)
  {
    if (!(newdist = (diststruct *) calloc(d->nbins+dnbins, sizeof(diststruct))))
      Exit("distributions", "D_AllocateDist", "calloc error");
    for (i=0; i<d->nbins; ++i)
      newdist[i+i0]	= d->dist[i];
    i0			= n<0 ? 0 : d->nbins;
    for (i=i0; i<i0+dnbins; ++i)
    {
      newdist[i].binsize= d->binsize+1;
      newdist[i].level	= d->level-1;
    }
    if (d->dist) free(d->dist);	
    d->dist		= newdist;
  }
  else
  {
    if ((!d->average)&&
	  !(d->average = (double *) calloc(D_NAVERAGE, sizeof(double))))
      Exit("distributions", "D_AllocateData", "calloc error");
    if (!((newbin = (long *) calloc(d->nbins+dnbins, sizeof(long)))&&
	  (newdata = (double *) calloc(d->nbins+dnbins, sizeof(double)))&&
	  (newcweight = (double *) calloc(d->nbins+dnbins, sizeof(double)))))
      Exit("distributions", "D_AllocateData", "calloc error");
    for (i=0; i<d->nbins; ++i)
    {
      newbin[i+i0]	= d->bin[i];
      newdata[i+i0]	= d->data[i];
      newcweight[i+i0]	= d->cweight[i];
    }
    if (d->bin) free(d->bin);
    if (d->data) free(d->data);
    if (d->cweight) free(d->cweight);
    d->bin		= newbin;
    d->data		= newdata;
    d->cweight		= newcweight;
  }
  d->nbins		+= dnbins;
  if (n<0)
    d->startbin		+= n;
}


void D_Reset(diststruct *d)
{
  long			i;

  if (d->level)
  {
    if (d->dist)
    {
      for (i=0; i<d->nbins; ++i)
	D_Reset(d->dist+i);
      free(d->dist);
    }
  }
  else
  {
    if (d->bin) free(d->bin);
    if (d->data) free(d->data);
    if (d->cweight) free(d->cweight);
    if (d->average) free(d->average);
  }
  d->n			= 0;
  d->nbins		= 0;
  d->ncount		= 0;
  d->startbin		= 0;
  d->dist		= NULL;
  d->bin		= NULL;
  d->data		= NULL;
  d->cweight		= NULL;
  d->average		= NULL;
}


void D_CopyRecursive(
  diststruct *dest, diststruct *src, double *binsize, int minus)
{
  long			i;
  
  dest->nbins		= src->nbins;
  dest->startbin	= src->startbin;
  dest->ncount		= minus ? -src->ncount : src->ncount;
  dest->n		= minus ? -src->n : src->n;
  if (!binsize)
  {
    if (!(dest->binsize||
	 (dest->binsize = (double *) calloc(src->level+1, sizeof(double)))))
      Exit("distributions", "D_CopyRecursive", "binsize calloc error");
    for (i=0; i<=src->level; ++i)
      dest->binsize[i] = src->binsize[i];
  }
  if (dest->level = src->level)			// distribution fork
  {
    if (src->nbins)
    {
      if (!(dest->dist = (diststruct *) calloc(src->nbins, sizeof(diststruct))))
        Exit("distributions", "D_CopyRecursive", "dist fork calloc error");
      for (i=0; i<src->nbins; ++i)
        D_CopyRecursive(dest->dist+i, src->dist+i, dest->binsize+1, minus);
    }
  }
  else if (src->nbins)				// data fork
  {
    if (!((dest->bin = (long *) calloc(src->nbins, sizeof(long)))&&
	  (dest->data = (double *) calloc(src->nbins, sizeof(double)))&&
	  (dest->cweight = (double *) calloc(src->nbins, sizeof(double)))&&
	  (dest->average = (double *) calloc(D_NAVERAGE, sizeof(double)))))
      Exit("distributions", "D_CopyRecursive", "data fork calloc error");
    for (i=0; i<D_NAVERAGE; ++i)
      dest->average[i]	= minus ? -src->average[i] : src->average[i];
    for (i=0; i<src->nbins; ++i)
    {
      dest->bin[i]	= minus ? -src->bin[i] : src->bin[i];
      dest->data[i]	= minus ? -src->data[i] : src->data[i];
      dest->cweight[i]	= minus ? -src->cweight[i] : src->cweight[i];
    }
  }
}


void D_Copy(diststruct *dest, diststruct *src)
{
  D_Reset(dest);
  D_CopyRecursive(dest, src, NULL, FALSE);
}


void D_Add(diststruct *dest, diststruct *src)
{
  long			i, n;

  if (!dest->nbins)
  {
    D_Reset(dest);
    D_CopyRecursive(dest, src, NULL, FALSE);
    return;
  }
  for (i=0; i<src->nbins; ++i)
  {
    D_Allocate(dest, i+src->startbin);
    n			= i+src->startbin-dest->startbin;
    if (src->level)
      D_Add(dest->dist+n, src->dist+i);
    else
    {
      dest->bin[n]	+= src->bin[i];
      dest->data[n]	+= src->data[i];
      dest->cweight[n]	+= src->cweight[i];
    }
  }
  if (dest->average&&src->average)
    for (i=0; i<D_NAVERAGE; ++i)
      dest->average[i]+= src->average[i];
  dest->n		+= src->n;
  dest->ncount		+= src->ncount;
}


void D_Subtr(diststruct *dest, diststruct *src)
{
  long			i, n;

  if (!dest->nbins)
  {
    D_Reset(dest);
    D_CopyRecursive(dest, src, NULL, TRUE);
    return;
  }
  for (i=0; i<src->nbins; ++i)
  {
    D_Allocate(dest, i+src->startbin);
    n			= i+src->startbin-dest->startbin;
    if (src->level)
      D_Add(dest->dist+n, src->dist+i);
    else
    {
      dest->bin[n]	-= src->bin[i];
      dest->data[n]	-= src->data[i];
      dest->cweight[n]	-= src->cweight[i];
    }
  }
  if (dest->average&&src->average)
    for (i=0; i<D_NAVERAGE; ++i)
      dest->average[i]-= src->average[i];
  dest->n		-= src->n;
  dest->ncount		-= src->ncount;
}


void D_Submit(diststruct *d, double *x, double *y, double *weight)
{
  register long		n = floor(*x/d->binsize[0])-d->startbin;
  register double	*avg = (double *) d->level;
  double		z;

  if ((n<0)||(n>=d->nbins))
  {
    D_Allocate(d, n += d->startbin);
    n			-= d->startbin;
  }
  while (d->level)				// advance to last level
    if (((n = floor(*(++x)/(d = d->dist+n)->binsize[0])-d->startbin)<0)||
	 (n>=d->nbins))
    {
      D_Allocate(d, n += d->startbin);
      n			-= d->startbin;
    }
  d->data[n]		+= *y * *weight;
  d->cweight[n]		+= *weight;
  ++(d->bin[n]);
  ++d->n;
  if ((!avg)&&(avg = d->average))		// only average for 1D
  {
    z			= *x;
    for (n=0; n<D_NAVERAGE; ++n)
    {
      *(avg++)		+= z;
      z			*= *x;
    }
  }
}


typedef
  struct {
    double		avg[D_NAVERAGE];
    long		ntotal, nbins;
  } distheaderstruct;

//#define DISTTEST

#ifdef DISTTEST
void D_Header(
  diststruct *d, int dist_type, distheaderstruct *header, long ncount)
{
  long			i, j;
  double		x, z;

  if (!d) return;
  if (d==D_Density3D)
    fprintf(stderr, "");
  if (dist_type&&(dist_type<3))
  {
    if (d->level)
    {
      for (i=0; i<d->nbins; ++i)
	D_Header(d->dist+i, dist_type, header, ncount);
      return;
    }
    header->nbins	+= d->nbins;
    header->ntotal	+= d->n;
    for (i=0; i<d->nbins; ++i)
    {
      switch (dist_type) {
        case 1: x = z = d->cweight[i] ? d->data[i]/d->cweight[i] : 0.0; break;
        case 2: x = z = ncount ? d->data[i]/ncount : 0.0; break;
      }         
      for (j=0; j<D_NAVERAGE; ++j)
      {
	header->avg[j]	+= z;
	z		*= x;
      }
    }	 
    return;
  }
  if (!d->average) return;
  for (i=0; i<D_NAVERAGE; ++i)
    header->avg[i]	= d->average[i]*d->nbins/(d->n ? d->n : 1.0);
  header->nbins		= d->nbins;
  header->ntotal	= d->n;
}


void D_PrintHeader(diststruct *d, int dist_type, double *L)
{
  char			s[40];
  long			i, n = d->level+1;
  double		x[n], stddev = 0;
  distheaderstruct	header;

  if (!d) return;
  if (d==D_Density3D)
    fprintf(stderr, "");
  for (i=0; i<D_NAVERAGE; ++i)
    header.avg[i]	= 0.0;
  header.nbins		= 0;
  header.ntotal		= 0;
  D_Header(d, dist_type, &header, d->ncount);
  printf("/* %s */\n\n", d->header);
  for (i=0; i<n; ++i)
    x[i]		= d->binsize[i]* (L ? *L : 1.0);
  PrintDVar("BINSIZE", n, x);
  PrintLVar("NCOUNT", 1, &d->ncount);
  PrintLVar("NBINS", 1, &header.nbins);
  PrintLVar("NTOTAL", 1, &header.ntotal);
  if (header.nbins)
  {
    for (i=0; i<D_NAVERAGE; ++i)
      header.avg[i]	/= header.nbins;
    stddev		= sqrt(header.avg[1]-header.avg[0]*header.avg[0]);
  }
  PrintDVar("STDDEV", 1, &stddev);
  PrintDVar("AVERAGE", D_NAVERAGE, header.avg);
  printf("\n");
  for (i=0; i<=d->level; ++i)
  {
    sprintf(s, "x%d", i);
    if (i) printf(", %12.12s", s);
    else printf("%12.12s", s);
  }
  printf(", %12.12s\n", "y");
  PrintLine(14*(d->level+2)-2);
}
#endif /* DISTTEST */


void D_PrintData(
  diststruct *d, int dist_type, double *x, double *L, long ncount, long nlevels)
{
  long			i, j;
  double		y, binsize = *(d->binsize)*(L ? *L : 1.0);

  for (i=0; i<d->nbins; ++i)
  {
/*    printf("nbins=%d\n", d->nbins);
    printf("startbin=%d\n", d->startbin);
    printf("level=%d\n", d->level);
    printf("i+d->startbin=%d\n", i+d->startbin);
    printf("binsize=%f\n", binsize);
    printf("binsize*(i+d->startbin)=%f\n", binsize*((double)(i+d->startbin)));
    printf("%d\t%f\n", nlevels-1-d->level, binsize*(i+(d->startbin)));
*/
    x[nlevels-1-d->level] = binsize*(i+(d->startbin));

    if (d->level)
      D_PrintData(d->dist+i, dist_type, x, L, ncount, nlevels);
    else
    {
      for (j=0; j<nlevels; ++j)
	printf("%s%12.6g", j ? ", " : "", x[j]);
      if (dist_type)
      {
	switch (dist_type) {
	  case 1: y = d->cweight[i] ? d->data[i]/d->cweight[i] : 0.0; break;
          case 2: y = ncount ? d->data[i]/ncount : 0.0; break;
	  case 3: y = d->data[i]; break;
	}
        printf(", %12.6g\n", y);
      }
      else printf(", %12d\n", d->bin[i]);
    }
  }
}


void D_Print(diststruct *d, int dist_type, double *L)
{
  double		*x;

  if (!d) return;
//  D_PrintHeader(d, dist_type, L);
  if (!(x = (double *) calloc(d->level+1, sizeof(double))))
    Exit("distributions", "D_Print", "calloc error");
  D_PrintData(d, dist_type, x, L, d->ncount, d->level+1);
  printf("\n");
  free(x);
}

#ifdef DISTTEST
void D_PrintMath(diststruct *d, int dist_type, double *L)
{
  char			s[40], mem[256];
  long			i, n = d->n ? d->n : 1;
  double		*avg = d->average;
  
  sprintf(mem, "{\nnbins -> %d, binsize -> %s, startbin -> %d, ntotal -> %d, ",
    d->nbins, DToMath(s, *(d->binsize)* (L ? *L : 1.0)), d->startbin, d->n);
  PrintToMath(mem, FALSE);
  sprintf(mem, "ncount -> %d", d->ncount); PrintToMath(mem, FALSE);
  if (d->average)
  {
    //D_Average(d, dist_type);
    for (i=0; i<D_NAVERAGE; ++i)
    {
      if (i) sprintf(mem, ", avg%d -> %s", i+1, DToMath(s, avg[i]/n));
      else sprintf(mem, ", avg -> %s", DToMath(s, avg[i]/n));
      PrintToMath(mem, FALSE);
    }
  }
  if (d->level)
  {
    PrintToMath(", dist -> {", FALSE);
    for (i=0; i<d->nbins; ++i)
      D_PrintMath(d->dist+i, dist_type, L);
    PrintToMath("}", FALSE);
  }
  else
  {
    PrintToMath(", data -> \n{\n", FALSE);
    for (i=0; i<d->nbins; ++i)
    {
      if (i) PrintToMath(", ", FALSE);
      sprintf(mem,"{%s", DToMath(s, d->data[i])); PrintToMath(mem, FALSE);
      sprintf(mem,", %s}", DToMath(s, d->cweight[i])); PrintToMath(mem, FALSE);
    }
    PrintToMath("\n}", FALSE);
  }
  PrintToMath("\n}\n", FALSE);
}
#endif /* DISTTEST */

// Only create a distribution once up initialization;  double creation leads
// to memory leaks;  distributions should only be created through D_Init.

void D_Create(diststruct *d, char *header, long nlevels, double *binsize)
{
  long			i;

  D_Reset(d);
  if (nlevels<1) 
    nlevels		= 1;
  d->level		= nlevels-1;
  if (!(d->header = (char *) calloc(strlen(header)+1, sizeof(char))))
    Exit("distributions", "D_Create", "header calloc error");
  strcpy(d->header, header);
  if (!(d->binsize = (double *) calloc(nlevels, sizeof(double))))
    Exit("distributions", "D_Create", "header/binsize calloc error");
  for (i=0; i<nlevels; ++i)
    d->binsize[i]	= binsize[i];
}


void D_Init(diststruct **dist, char *header, long nlevels, double *binsize)
{
  char			s[80];
  long			system;

  if (*dist) 					// if existed already, then reset
  { 
    for (system=0; system<NSYSTEMS; ++system) 
      D_Reset(*dist+system); 
    return; 
  }
  if (!(*dist=(diststruct *) calloc(NSYSTEMS, sizeof(diststruct))))
    Exit("distributions", "D_Init", "calloc error");
  for (system=0; system<NSYSTEMS; ++system)
  {
    sprintf(s, "%s distribution of system %d", header, system);
    D_Create(*dist+system, s, nlevels, binsize);
  }
}


void PutInDistribution(diststruct *d, double x, double y, double weight)
{
  D_Submit(d, &x, &y, &weight);
}


void PrintDistribution(diststruct *d)
{
  D_Print(d, 0, NULL);
}


void InitDist(long flag, diststruct **dist, char *header, long nlevels, double *binsize)
{
  if (flag)
    D_Init(dist, header, nlevels, binsize);
}
                                                                                                                                                                                                                                                                                                                                                         src/dumphst.c                                                                                       0000600 0143352 0000144 00000027277 11047110557 012724  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    dumphst.c
    author:     Peng Yi at MIT
    date:       January 12, 2008
    purpose:    get info from dump files and do analysis
    note:	previously I have a program hst.c to deal with the binary history file
		produced by Pieter's code.  However later I realize that I didn't need
		that much information to be stored in the history file.  Even it is 
		in binary format, it still take a huge amount of disk space.
		I also noticed that lammps dumps in human-readable format.  Of course 
		as an MD code, lammps can dump coordinates as velocity and force and 
		etc.  So I decide to only dump coordinates in my MC code and in human-
		readable format.  The analysis program is thus change to dumphst.c
*/

#define __MAIN_PROGRAM
#include "header.h"

#include "correlation.h"

long	timestep = 0,
	min=MAXNMOLS, 
	max=0,
	nnmax[MAXNMOLS],		// nnmax[i] is the counts of conf. with nmax=i
	n[1000][MAXNMOLS];		// n[i][j] is the counts of nuclei of size j when nmax=i
double	pn[1000][MAXNMOLS],		// pn[i][j] is the prob. of nuclei of size j when nmax=i
	Gn[MAXNMOLS],			// free energy of one nucleus
	G[MAXNMOLS];			// free energy of the system

// some functions declared here
void printout(FILE *);
void hst2conf();
void hst2lammpsdump(FILE *fPtr);

long Dump_GetNext(FILE *fPtr)		// get next configuration from dump file, and primary processing
{					// very much like GetCoordinate() function in io.c
   char			name[80], dummy[80];
   long			i, j, id, system;
   molstruct		*moli;

   // read in data

   fscanf(fPtr, "%s", name);
   fscanf(fPtr, "%ld", &timestep);			// read in timestep
   fscanf(fPtr, "%ld%ld%ld", &NSYSTEMS, &NMOLS, &NSITES);
   NBOX	=	NSYSTEMS;		// temporary
   for (i=0; i<NSYSTEMS; i++)
      fscanf(fPtr, "%lf%lf%lf", &(BOX[i].lx), &(BOX[i].ly), &(BOX[i].lz));
   for (i=0; i<NMOLS; i++) {
      fscanf(fPtr, "%ld%ld%ld", &id, &(mol[i].box), &(mol[i].nsites));
      //fscanf(fconf, "%ld%ld%ld%ld%ld", &id, &(mol[i].box), &(mol[i].nsites), &(mol[i].fix), &(mol[i].flip));
      for (j=0; j<mol[i].nsites; j++)
         fscanf(fPtr, "%ld%lf%lf%lf", &(mol[i].type[j]), &(mol[i].p[j].x), &(mol[i].p[j].y), &(mol[i].p[j].z));
   }
   fgets(dummy, sizeof(dummy), fPtr);

   // some primary data processing

   CoorSI2System();				// the dump file is in real units (8/7/08)
   for (i=0; i<NSYSTEMS; i++) {
      NMols[i]	=	0;
      NSites[i]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if ( (system=moli->box) >= 0) {
         NMols[system]	++;				// total # of mols in certain system
         NSites[system]	+=	moli->nsites;		// total # of sites in certain system
      }
   }
   for (i=0; i<NSYSTEMS; i++) {
      BOX[i].lbox	=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
      BOX[i].vol	=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
      BOX[i].pres 	=	P;
      BOX[i].temp 	=	kT;
      BOX[i].drmax 	=	DRMAX;
      BOX[i].dlmax 	=	DLMAX;
      BOX[i].damax	=	DAMAX;
      BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
      BOX[i].rb		=	Rb;
      BOX[i].rv		=	Rv;
   } 
   for (i=0; i<NMOLS; i++) { 				// note NMOLS might change due to MPI
      for (j=0; j<mol[i].nsites; j++)  {
         mol[i].flags[j]	=	1;		// activate all the sites on this processor
         mol[i].parent[j]	=	j-1;		// initialize parent sites
      }

      mol[i].flip		=	0;		// flip to the original direction
      mol[i].origin		=	CenterofMass(mol+i);
   }

   // return 

   return feof(fPtr)==0;
}


int main(int argc, char * argv[])
{
   FILE		*fhist, *fp, *fp1, *fsize, *fconf, *frdf, *flammps;
   long		Nsplit, 
		Nblock = -1;		// # of records in each block, could be set by command line
   long		version, system;
   char		file_hst1[256],
		comments[1024];
   double	dummy;
   long		i, j, k; 
   molstruct	*moli;
   double	var, bias; 
   static long	init = 1;

   if (argc<2) {
      printf("dumphst (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tdumphst dumpfile [N]\n\n");
      printf("Notes:\n");
      printf("\t* require setup file to get e.g. Temp for analysis\n\n");
      printf("\t\t if N is given, then do analysis in blocks each with N records, otherwise as one big block\n\n");
      exit(-1);
   }

   if (!(fp=fopen(argv[1], "r"))) {			// open dump file
      printf("dumphst: dump file cannot open!\n");
      exit(-1);
   }
   if (argc == 3)
      Nblock		=	atol(argv[2]);		// otherwise = -1

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate max. memory for molecules
   GetSetup(argv);			// most of the variables read from setup 
					// will be updated by the reading from history file
   InitUnits();				// calculate the units needed for system<->SI conversion
					// and convert setup variable to system units 
   InitForcefield();
 
   if (!(fsize = fopen("Psize","w"))) {			// open output file for size distribution information
      printf("dumphst: Size prob. distribution file failed to open");
      exit(-1);
   }
   if (!(flammps=fopen("lammps.dump", "w"))) {		// to convert to lammps dump file for lammps2pdb
      printf("dumphst: lammps dump file failed to open.\n");
      exit(-1);
   }

   while (Dump_GetNext(fp)) {			// read in next record/configuration, including init lattice

#ifdef CELL_LIST
      if (init) {
         CL_Init();	 init	=	0;
      }
      else
	 CL_Destroy();
      CL_Build();
#endif
      CalcV();					// need cell list

      if (dynvar==3)
         SampleP2All();				// calculate global P2 and local p2
      SampleSpherical();
      Find_Nuclei();

printf("%d %f %f %d\n", timestep, v[0].tot/NSITES, NSITES/BOX[0].vol, nmax[0][0]);

      // calculate bias for each record/configuration 

      if (dynvar==1)
         var	=	MAXSIZE[0];
      else if (dynvar==3)
	 var	=	P2M[0];

      bias	=	exp(0.5*kP2*(var-P2middle)*(var-P2middle));

      // collect probability distribution of cluster size 
 
      nnmax[MAXSIZE[0]]	++;
      for (i=1; i<=MAXSIZE[0]; i++) {
	 pn[0][i]	+=	bias * sizedist[i];		// probability distribution of cluster size
	 n[0][i]	+=	sizedist[i];

         pn[MAXSIZE[0]][i]	+=	bias * sizedist[i];
         n[MAXSIZE[0]][i]	+=	sizedist[i];
      }
      if (MAXSIZE[0] > max)
	 max	=	MAXSIZE[0];
      if (MAXSIZE[0] < min) 
         min	=	MAXSIZE[0];

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
      } 
      if (!(timestep%IRADIAL))
         radial("sample");		// sample rdf

      correlation();			// calculate correlation function

      //Write_Conf(timestep);		// convert to ASCII dump file

      hst2lammpsdump(flammps);		// convert to lammps dump format to use lammps2pdb in the future

      timestep	+=	ITAPE;		// because the binary hist file was stored every ITAPE cycles

      if (Nblock!=-1 && timestep%Nblock==0) {		// summarize the results for this block and output
         for (i=1; i<=max; i++) 
            Gn[i]	=	-log(pn[0][i]);
         for (i=min; i<=max; i++) {
            for (j=1; j<=i; j++) {
               G[i]	+=	n[i][j] * Gn[j];
            }
         }

	 printout(fsize);

         // clear up the counters
         for (i=0; i<MAXNMOLS; i++) {
	    nnmax[i]	=	0;
	    Gn[i]	=	0.0;
	    G[i]	=	0.0;
	 }
         for (i=0; i<100; i++) { 
            for (j=0; j<MAXNMOLS; j++) {
               pn[i][j]	=	0.0;
	       n[i][j]	=	0;
            }
         }
         max	=	0;
	 min	=	MAXNMOLS;
      }
   }

   /* normalize correlation functions and print out */

   corr_norm();
   corr_print();
   radial("print");

   // summarize results for all records and print out size distribution

   if (Nblock==-1) {
      for (i=1; i<=max; i++) 
         Gn[i]	=	-log(pn[0][i]);			// free energy of nucleus
      for (i=min; i<=max; i++) {
         for (j=1; j<=i; j++) {
            G[i]	+=	n[i][j] * Gn[j];	// free energy of the system
         }
      }
      printout(fsize);
   }

   fclose(flammps);
   fclose(fp);
   fclose(fsize);
   return EXIT_SUCCESS;
}


void printout(FILE *fPtr)
{
   long		i, j;

   fprintf(fPtr, "\nOverall nuclei size distribution P(n):\n\n");
   fprintf(fPtr, "n    N(n)   fixed p(n)   -log(fixed pn[i])\n\n");
   for (i=1; i<=max; i++) 
      fprintf(fPtr, "%4d   %8d   %14.4f   %f\n", i, n[0][i], pn[0][i], Gn[i]);

   fprintf(fPtr, "P(n) for copy purpose...\n");
   for (i=1; i<=max; i++)
      fprintf(fPtr, "%f\n", Gn[i]);

   fprintf(fPtr, "\nDistribution of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, "n_max    N(n_max)\n");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%4d  %8d\n", i, nnmax[i]);

   fprintf(fPtr, "\nSize distribution N(n) for different value of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, " n_max");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%8d  ", i);
   fprintf(fPtr, "\n   n \n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%8d  ", n[j][i]);
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n  Normalized to N(n_max)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", (nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n  -Logrithm of N(n)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", -log(nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
}


void hst2conf()
{
   long		i;
   molstruct	*moli;
   static FILE	*fPtr;
   static long	init =1;

   vector	rA, rB, rAB;
   double	r2, L;
   long		id, print=0;

   if (init) {
      if (!(fPtr=fopen("hst.conf", "w")) )
          Exit("hst", "main", "open conf file failed.");
      init	=	0;
   }
/*
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (sizeofnucl[moli->nuclid[0]]==MAXSIZE[0]) {
         id	=	moli->nuclid[0];
	 rA	=	CenterofMass(moli);
         break;
      }
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (moli->nuclid[0] == id) {
         rB	=	CenterofMass(moli);
	 rAB	=	V_Subtr(&rA, &rB);
         r2	=	V_Dot(&rAB, &rAB);
         L	=	MIN(BOX[0].lx, MIN(BOX[0].ly, BOX[0].lz));
         if (r2 > L * L/2) {
            print	=	1;
	    break;
         }
      }

   if (!print)
      return;
*/
   fprintf(fPtr, "!TIMESTEP %d Nmax = %d 2ndNmax = %d\n", timestep, MAXSIZE[0], secondNmax[0]);
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);

   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx*unit.LENGTH, BOX[i].ly*unit.LENGTH, BOX[i].lz*unit.LENGTH);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
   }
}

void hst2lammpsdump(FILE *fPtr)
{
   long		i;
   molstruct	*moli;

   fprintf(fPtr, "ITEM: TIMESTEP\n");	
   fprintf(fPtr, "%d\n", timestep);		
   fprintf(fPtr, "ITEM: NUMBER OF ATOMS\n");
   fprintf(fPtr, "%d\n", NSITES);
   fprintf(fPtr, "ITEM: BOX BOUNDS\n");
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lx*unit.LENGTH, 0.5*BOX[0].lx*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].ly*unit.LENGTH, 0.5*BOX[0].ly*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lz*unit.LENGTH, 0.5*BOX[0].lz*unit.LENGTH); 
   fprintf(fPtr, "ITEM: ATOMS\n");
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         fprintf(fPtr, "%d %d %f %f %f %f %f %f %d %d %d\n", (moli-mol)*moli->nsites+i, moli->type[i], 
              moli->p[i].x*unit.LENGTH, moli->p[i].y*unit.LENGTH, moli->p[i].z*unit.LENGTH, 0, 0, 0, 0, 0, 0);
}

                                                                                                                                                                                                                                                                                                                                 src/dumpreduce.c                                                                                    0000600 0143352 0000144 00000003762 11146140254 013363  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
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
              src/dumpreorder.c                                                                                   0000600 0143352 0000144 00000003422 11152246351 013551  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
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
                                                                                                                                                                                                                                              src/embed.c                                                                                         0000600 0143352 0000144 00000010415 11110121660 012262  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	embed.c
	author:		Peng Yi at MIT
	date:		May 18, 2008
	purpose:	embed a smaller system to the center of a larger system, do NO analysis
	note:		require setup file
*/

#define __MAIN_PROGRAM
#include "header.h"


int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   char		s[255], ff[255], dummy[80];
   molstruct	*moli, *mol_in;
   FILE		*fPtr;

   vector	com,
		Linner[MAXNSYSTEMS];		// inner BOX dimension
   long		i, j, site_in, n, system=0,
		nsystems,			// for one system only, 5/18/2008
		nmols1, nsites1,		// NMOLS and NSITES for inner box
		nmolstot, nsitestot,		// total NMOLS and NSITES
		id;
   double	lambda;
   long		overlap;

   if (argc<4) {
      printf("embed (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tembed file1 file2 lambda\n\n");
      printf("Notes:\n");
      printf("\t* file1 is the configuration file for inner box, file2 for outer box\n\n");
      printf("\t* lambda <= 1.0 is the overlap tuning parameter\n\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }
   lambda	=	atof(argv[3]);

   InitMols(MAXNMOLS, 1);	// allocate memory for molecules, MAXNMOLS must > the sum of NMOLS
   GetSetup(argv);		// mainly PBC I think

   // Read in the conf. of the inner system

   Read_Conf(argv[1]);
   for (moli=mol; moli<mol+NMOLS; moli++)	// arrange the inner system
      MolInBox2(moli);

   nsystems	=	NSYSTEMS;		// save inner system information
   nmols1	=	NMOLS;
   nsites1	=	NSITES;
   nmolstot	=	NMOLS;
   nsitestot	=	NSITES;

   for (i=0; i<NSYSTEMS; i++) {
      Linner[i].x	=	BOX[i].lx;
      Linner[i].y	=	BOX[i].ly;
      Linner[i].z	=	BOX[i].lz;
   }

   // Read in the conf. of the outer system

   if (!(fPtr=fopen(argv[2], "r")))		// open the conf. file of the outer system
      Exit("embed", "main", "big system conf. file failed to open.");

   fscanf(fPtr, "%s%ld", dummy, &TIMESTEP);			// read in timestep
   fscanf(fPtr, "%ld%ld%ld", &NSYSTEMS, &NMOLS, &NSITES);	// update NSYSTEMS, NMOLS, NSITES
   if (nsystems != NSYSTEMS) {}					// might do something in the future
   for (i=0; i<NSYSTEMS; i++)
      fscanf(fPtr, "%lf%lf%lf", &(BOX[i].lx), &(BOX[i].ly), &(BOX[i].lz));	// update BOX dimension

   moli	=	mol + nmols1;			// clear statement, although not necessary
   for (i=0; i<NMOLS; i++) {			// start reading in the molecules in the outer box
      fscanf(fPtr, "%ld%ld%ld", &id, &(moli->box), &(moli->nsites));
      for (j=0; j<moli->nsites; j++)
         fscanf(fPtr, "%ld%lf%lf%lf", &(moli->type[j]), &(moli->p[j].x), &(moli->p[j].y), &(moli->p[j].z));

      MolInBox2(moli);				// map back in outer box

      overlap	=	0;
      for (j=0; j<moli->nsites; j++) {		// check each bead in the outer box
         for (mol_in=mol; mol_in<mol+nmols1; mol_in++) {
            for (site_in=0; site_in<mol_in->nsites; site_in++) {
	       if (DistSQ(moli->p[j], mol_in->p[site_in], system) < type[0].SIGMA * type[0].SIGMA * lambda) {
		  overlap	=	1;
		  break;
	       }
	    }
	    if (overlap)
		break;
	 }
	 if (overlap)
	    break;
      }
      if (!overlap) {
         nmolstot	++;
	 nsitestot	+=	moli->nsites;
	 moli	++;
      }

/*
      overlap	=	0;			// is this molecule overlap with inner box?
      for (j=0; j<moli->nsites; j++) {		// check every bead
         if (fabs(moli->p[j].x)*lambda < 0.5*Linner[moli->box].x  && 
		fabs(moli->p[j].y)*lambda < 0.5*Linner[moli->box].y && 
		fabs(moli->p[j].z)*lambda < 0.5*Linner[moli->box].z ) {
            overlap	=	1;
            break;
	 }
      }
      if (!overlap) {
         nmolstot	++;
	 nsitestot	+=	moli->nsites;
	 moli	++;
      }
*/
/*
      com	=	CenterofMass(moli);	// center of mass of one chain
      com.x	*=	lambda;		// give some gap
      com.y	*=	lambda;
      com.z	*=	lambda;

      if (!(com.x < 0.5*Linner[moli->box].x && com.x >= -0.5*Linner[moli->box].x &&
	   com.y < 0.5*Linner[moli->box].y && com.y >= -0.5*Linner[moli->box].y &&
	   com.z < 0.5*Linner[moli->box].z && com.z >= -0.5*Linner[moli->box].z) ) {	// com NOT in the inner box
         nmolstot	++;
	 nsitestot	+=	moli->nsites;	
	 moli		++;
      }
*/
   }
   NMOLS	=	nmolstot;
   NSITES	=	nsitestot;

   // PRINT OUT CONFIGURATION

   system	=	0;			// only one system for now
   unit.LENGTH	=	1.0;			// not doing any analysis
   Write_Conf(-1);				// output combined configuration

   return	0;
}
                                                                                                                                                                                                                                                   src/ensembles.c                                                                                     0000600 0143352 0000144 00000207165 11317503016 013205  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /* 
	program:	ensembles.c
	purpose:	Collection of subroutines for Monte Carlo moves
			in a NPT ensemble, using Verlet list
	author: 	Peng Yi at MIT
	date:		October 19, 2006
	notes:		
			April 28, 2009.	Add flip move
			July 25, 2009.	Add NnPTmu ensemble in endbridge() 
					for polydisperse system.
                        December, 2009, Add double bridging (DB), 
					and intramolecular double rebridging (IDR)
*/

#define __ENSEMBLES_MODULE
#include "ensembles.h"

/*
	General procedure to make a Monte Carlo move

	First record all the old configuration properties that will be 
	changed during the move, e.g., particle coordinates, energy, 
	system volume, etc.  Then perform the random move and calculate 
	the new properties and keep them as the system variables.  
	If accept the move, then we simply increase the accept counter 
	by one, otherwise, we recover all the old configuration properties.

	The reason to update the system variables is because the calculate 
	will depend on them, so better let them go with any possible move 
	and we just store the old information.

	Notes:
		Apr, 2007:	Verlet neigh list creation.
		Jun 4, 2007:	Cell neigh list creation, program can choose from no neighbor list,
				Verlet list only, Cell list only, or Verlet+Cell.
		Sept. 2007:	Enable multiple boxes.  Add Gibbs ensemble implementation.
		Oct. 2007:	Enable chain molecule simulation.
		Nov. 8, 2007:	rotation move
*/

boxstruct	oldBOX[MAXNBOX];
long		oldNSites[MAXNBOX], oldNMols[MAXNBOX];
vstruct		oldv[MAXNBOX];
wstruct		oldvir[MAXNBOX];

double		oldP2[MAXNSYSTEMS];     	// global orientation order
double		oldP2M[MAXNSYSTEMS];   		// global orientation order
double		oldP2z[MAXNSYSTEMS];    	// P2 with respect to z-axis
double		oldQ6[MAXNSYSTEMS];		// Q6
double		oldtransfrac[MAXNSYSTEMS]; 	// trans state fraction
long		oldMAXSIZE[MAXNSYSTEMS];
long		oldNnucl[MAXNSYSTEMS];
long		oldXtal[MAXNSYSTEMS];
long		oldrealXtal[MAXNSYSTEMS];
long		oldsecondNmax[MAXNSYSTEMS];

void StoreSystem()				//store system information
{
   long		i;

   for (i=0; i<NBOX; i++) {
      oldv[i]		=	v[i];
      oldvir[i]		=	vir[i];
      oldBOX[i]		=	BOX[i];
      oldNSites[i]	=	NSites[i];
      oldNMols[i]	=	NMols[i];

      oldP2[i]		=	P2[i];
      oldP2M[i]		=	P2M[i];
      oldP2z[i]		=	P2z[i];
      oldQ6[i]		=	Q6[i];
      oldtransfrac[i]	=	transfrac[i];      
      oldMAXSIZE[i]	=	MAXSIZE[i];
      oldNnucl[i]	=	Nnucl[i];
      oldXtal[i]	=	Xtal[i];
      oldrealXtal[i]	=	realXtal[i];
      oldsecondNmax[i]	=	secondNmax[i];
   }
}


void RestoreSystem()				//restore system information
{
   long		i;

   for (i=0; i<NBOX; i++) {
      v[i]	=	oldv[i];
      vir[i]	=	oldvir[i];
      BOX[i]	=	oldBOX[i];
      NSites[i]	=	oldNSites[i];
      NMols[i]	=	oldNMols[i];

      P2[i]		=	oldP2[i];
      P2M[i]		=	oldP2M[i];
      P2z[i]		=	oldP2z[i];
      Q6[i]		=	oldQ6[i];
      transfrac[i]	=	oldtransfrac[i];      
      MAXSIZE[i]	=	oldMAXSIZE[i];
      Nnucl[i]		=	oldNnucl[i];
      Xtal[i]		=	oldXtal[i];
      realXtal[i]	=	oldrealXtal[i];
      secondNmax[i]	=	oldsecondNmax[i];
   }
}


void StoreOneMol(long n)	//store system and single molecule information
{
   oldmol[0]	=	mol[n];
   StoreSystem();
}


void RestoreOneMol(long n)
{
   mol[n]	=	oldmol[0];
   RestoreSystem();
}


void StoreMols()
{
   long		i;
   for (i=0; i<NMOLS; i++) {
      oldmol[i]	=	mol[i];
   }
   StoreSystem();
}


void RestoreMols()
{
   long		i;
   for (i=0; i<NMOLS; i++) {
      mol[i]	=	oldmol[i];
   }
   RestoreSystem();
}


void ResetAcceptance()
{
   long		ib, i;

   for (ib=0; ib<NBOX; ib++) {
      av[ib].move		=	0;	// single site displacement
      av[ib].acc_move		=	0;
      av[ib].vol		=	0;	// volume change
      av[ib].acc_vol		=	0;
      av[ib].gibbsvol		=	0;	// gibbs volume change
      av[ib].acc_gibbsvol	=	0;
      av[ib].swap		=	0;	// gibbs swap
      av[ib].acc_swap		=	0;
      av[ib].cbmc		=	0;	// cbmc
      av[ib].acc_cbmc		=	0;
      av[ib].rep		=	0;	// end reptation
      av[ib].acc_rep		=	0;
      av[ib].erot		=	0;	// end rotation
      av[ib].acc_erot		=	0;
      av[ib].flip		=	0;	// flip
      av[ib].acc_flip		=	0;
      av[ib].rot		=	0;
      av[ib].acc_rot		=	0;
      av[ib].eb			=	0;	// end bridging
      av[ib].acc_eb		=	0;
      av[ib].re			=	0;	// rebridging
      av[ib].acc_re		=	0;
      av[ib].db			=	0;	// double bridging
      av[ib].acc_db		=	0;
      av[ib].idr		=	0;	// intramolecular double rebridging
      av[ib].acc_idr		=	0;
      av[ib].movemol		=	0;
      av[ib].acc_movemol	=	0;
      av[ib].seq		=	0;
      av[ib].acc_seq		=	0;

      av_past[ib]	=	av[ib];
   }
   for (i=0; i<NSITES/NMOLS; i++)
      cbmcsucc[i]	=	0;
}


void Adjust_Stepsize()
{
   long		ib;
   double	frac, dold;
   double	drbond;

   for (ib=0; ib<NBOX; ib++) {

      // adjust displacement step size

      if (av[ib].move == 0 || av_past[ib].move > av[ib].move) {
         av_past[ib].move	=	av[ib].move;
         av_past[ib].acc_move	=	av[ib].acc_move;
      }
      else {
         frac	=	((double) (av[ib].acc_move - av_past[ib].acc_move))/(av[ib].move-av_past[ib].move);
         dold	=	BOX[ib].drmax;
         BOX[ib].drmax	*=	fabs(frac/SUCC_DISP);

         if (BOX[ib].drmax / dold > 1.5)	BOX[ib].drmax = dold * 1.5;
         if (BOX[ib].drmax / dold < 0.5)	BOX[ib].drmax = dold * 0.5;

#ifdef VERLET_LIST
	 drbond		=	2.0 * (BOX[ib].rv - BOX[ib].rc);
#else
	 drbond		=	0.25 * BOX[ib].lbox;
#endif	/* VERLET_LIST */

         if (BOX[ib].drmax > drbond)	BOX[ib].drmax = drbond;
         
         av_past[ib].move	=	av[ib].move;
         av_past[ib].acc_move	=	av[ib].acc_move;
      }

      // adjust cbmc first atom displacement step size

      if (av[ib].cbmc == 0 || av_past[ib].cbmc > av[ib].cbmc) {
         av_past[ib].cbmc	=	av[ib].cbmc;
         av_past[ib].acc_cbmc	=	av[ib].acc_cbmc;
      }
      else {
         frac	=	((double) (av[ib].acc_cbmc - av_past[ib].acc_cbmc))/(av[ib].cbmc-av_past[ib].cbmc);
         dold	=	BOX[ib].drmax;
         BOX[ib].drmax	*=	fabs(frac/SUCC_DISP);

         if (BOX[ib].drmax / dold > 1.5)	BOX[ib].drmax = dold * 1.5;
         if (BOX[ib].drmax / dold < 0.5)	BOX[ib].drmax = dold * 0.5;

#ifdef VERLET_LIST
	 drbond		=	2.0 * (BOX[ib].rv - BOX[ib].rc);
#else
	 drbond		=	0.25 * BOX[ib].lbox;
#endif	/* VERLET_LIST */

         if (BOX[ib].drmax > drbond)	BOX[ib].drmax = drbond;
         
         av_past[ib].cbmc	=	av[ib].cbmc;
         av_past[ib].acc_cbmc	=	av[ib].acc_cbmc;
      }

      // adjust volume change step size
      
      if (av[ib].vol == 0 || av_past[ib].vol > av[ib].vol ) {
         av_past[ib].vol	=	av[ib].vol;
         av_past[ib].acc_vol	=	av[ib].acc_vol;
      }
      else {
         frac	=	((double) (av[ib].acc_vol-av_past[ib].acc_vol))/(av[ib].vol-av_past[ib].vol);
         dold	=	BOX[ib].dlmax;
         BOX[ib].dlmax	*=	fabs(frac/SUCC_VOL);
        
         if (BOX[ib].dlmax / dold > 1.5)	BOX[ib].dlmax = dold * 1.5;
         if (BOX[ib].dlmax / dold < 0.5)	BOX[ib].dlmax = dold * 0.5;
 
         av_past[ib].vol	=	av[ib].vol;
         av_past[ib].acc_vol	=	av[ib].acc_vol;
      }
   }
}


void MolInBox(molstruct *molm)		// map a molecule back into central box
{					// by mapping its first site into central box
   long		i;
   vector	dp;

   dp		=	molm->p[0];
   molm->p[0]	=	MapInBox2(molm->p, PBC, molm->box);
   dp.x		-=	molm->p->x;
   dp.y		-=	molm->p->y;
   dp.z		-=	molm->p->z;
   if ( (fabs(dp.x)>1e-8) || (fabs(dp.y)>1e-8) || (fabs(dp.z)>1e-8) ) {
      for (i=1; i<molm->nsites; i++) {
         molm->p[i].x	-=	dp.x;
         molm->p[i].y	-=	dp.y;
         molm->p[i].z	-=	dp.z;
      }
      molm->origin.x	-=	dp.x;	// make sure that drift doesn't change
      molm->origin.y	-=	dp.y;
      molm->origin.z	-=	dp.z;
   }
}

void MolInBox2(molstruct *molm)		// map a molecule back into central box
{					// based on its center of mass
   long		i;
   vector	cm, dp;
   
   cm		=	CenterofMass(molm);
   dp		=	MapInBox2(&cm, PBC, molm->box);
   dp.x		-=	cm.x;
   dp.y		-=	cm.y;
   dp.z		-=	cm.z;
   if ( (fabs(dp.x)>1e-8) || (fabs(dp.y)>1e-8) || (fabs(dp.z)>1e-8) ) {
      for (i=0; i<molm->nsites; i++) {
         molm->p[i].x	+=	dp.x;
         molm->p[i].y	+=	dp.y;
         molm->p[i].z	+=	dp.z;
      }
      molm->origin.x	+=	dp.x;	// make sure that drift doesn't change
      molm->origin.y	+=	dp.y;
      molm->origin.z	+=	dp.z;
   }
}
   	

long SiteSelect(molstruct **moli)	// randomly select an atom from all atoms in all boxes
{					// return the molecule and site id of this atom
   long		site;

   do {
      *moli	=	mol;
      site	=	ran1(seed) * NSITES;	// NSITES is total # of sites in all mols in all boxes

      while ( site >= (*moli)->nsites ) {
         site	-=	(*moli)->nsites;
         (*moli)	++;
      }
   } while ( (*moli)->flags[site] == 0 );	// if site is not active
   return	site;
}

long molcheck(molstruct *moli)
{
   long		i;
   for (i=0; i<moli->nsites; i++) {
     if (moli->flags[i] != 1)
        printf("molcheck flag error, molid %d siteid %d\n", moli-mol, i);
     if (moli->parent[i] != i-1)
        printf("molcheck parent error, molid %d siteid %d\n", moli-mol, i);
   } 

   return	i;
}

long mcmove()
{
   long		ibox, site, i;
   double	dV, arg;
   molstruct	*molm, mol_old;
   vstruct	v_old;
   wstruct	vir_old;
   vector	dp;

   site		=	SiteSelect(&molm);	// select an active site in random
   ibox		=	molm->box;

   if (molm->fix)				// if molecule is fixed
      return	ibox;

   mol_old	=	*molm;
   v_old	=	v[ibox];
   vir_old	=	vir[ibox];

   dV		=	VDeleteSites(molm, site, site);
   dp.x		=	BOX[ibox].drmax * (ran1(seed)-0.5);
   dp.y		=	BOX[ibox].drmax * (ran1(seed)-0.5);
   dp.z		=	BOX[ibox].drmax * (ran1(seed)-0.5);

   molm->p[site]	=	V_Add(molm->p+site, &dp);

   // do not do MapInBox to change particle position, just let the particle move.  
   // The PBC effect will be taken care of in energy calc. and other places
//printf("mcmove del dV = %f\n", dV);
   dV		+=	VAddSites(molm, site, site);
   arg		=	dV/(BOX[ibox].temp);

/*printf("mcmove add dV = %f\n", dV);
if (fabs(dV-(v[ibox].tot-v_old.tot)) > 0.000001)  {printf("mcmove_err\n"); exit(1);}
*/
   if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {
      av[ibox].acc_move	++; 
/*printf("%d %d\n", molm-mol, site);
printf("%f %f %f %f %f %f\n", mol_old.p[site].x, mol_old.p[site].y, mol_old.p[site].z, molm->p[site].x, molm->p[site].y, molm->p[site].z);
printf("displacement accept  old energy =%f  new energy =%f\n", v_old.tot, v[ibox].tot);
*/
/*
      for (i=0; i<4; i++)		// update spherical coord.
         if (site+i < molm->nsites)
            molm->s[site+i]	=	SiteSpherical(molm, site+i);
*/
   }
   else {
      v[ibox]	=	v_old;
      vir[ibox]	=	vir_old;
#ifdef CELL_LIST
      CL_Delete(molm, site);			// Unregister trial site
#endif
      *molm	=	mol_old;
#ifdef CELL_LIST
      CL_Add(molm, site);			// Register old site
#endif
   }
   av[ibox].move	++;
   return	ibox;
}


long displacemol()		// displace one molecule without changing its internal conf.
{
   long		i, system;
   molstruct	*moli;
}


long rotation(long n)				// rotate the end n (n<=3) sites, 
{						// bond length and angle unchanged
   long		nsites, i, system;
   double	dV, arg;
   sphere	s[3];
   molstruct	*molm, mol_old;
   vstruct	v_old;
   wstruct	vir_old;

   molm		=	mol + (long) (NMOLS * ran1(seed));
   if (ran1(seed) < 0.5)
      MolFlip(molm);
   system	=	molm->box;

   n	=	n<1 ? 1 : (n>3 ? 3 : n);
   for (nsites=0; nsites<n; nsites++)
      if ( (nsites>=molm->nsites) || (0==molm->flags[molm->nsites-nsites-1]) )
         break;

   if (!nsites)
      return system;

   mol_old	=	*molm;
   v_old	=	v[system];
   vir_old	=	vir[system];

   for (i=0; i<nsites; i++)
      s[i]	=	SiteSpherical(molm, molm->nsites-nsites+i);
   VDeleteSites(molm, molm->nsites-nsites, molm->nsites-1);
   for (i=0; i<nsites; i++) {
      s[i].beta	=	AdjustAngle( s[i].beta+(ran1(seed)-0.5)*BOX[system].damax );
      molm->p[molm->nsites-nsites+i]
		=	SiteCartesian(molm, molm->nsites-nsites+i, s[i]);
   }
   VAddSites(molm, molm->nsites-nsites, molm->nsites-1);
   arg 		=	(v[system].tot - v_old.tot)/BOX[system].temp;

   if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {
      av[system].acc_rot	++;
      //MolInBox(molm);				// no mapping into central box in moves
   }
   else {
      v[system]		=	v_old;
      vir[system]	=	vir_old;
#ifdef CELL_LIST
      for (i=molm->nsites-nsites; i<molm->nsites; i++)
         CL_Delete(molm, i);
#endif
      *molm		=	mol_old;
#ifdef CELL_LIST
      for (i=molm->nsites-nsites; i<molm->nsites; i++)
         CL_Add(molm, i);
#endif
   }
   av[system].rot	++;
   return	system;  
}


long movemol()				// move a molecule with bond length fixed
{
   long		i, system;
   molstruct	*molm, mol_old;
   vstruct	v_old;
   wstruct	vir_old;
   double	arg;

   molm		=	mol + (long) (NMOLS * ran1(seed));	// pick up a mol in random
   if (ran1(seed) > 0.5)
      MolFlip(molm);
  
   MolSpherical(molm);		// can be remove if spherical coord. is updated after every move
  
   system	=	molm->box;
   mol_old	=	*molm;
   v_old	=	v[system];
   vir_old	=	vir[system];

   VDeleteSites(molm, 0, molm->nsites-1);
   
   /* Perturb the first site using Cartesian coordinates */

   molm->p->x	+=	(ran1(seed) - 0.5) * BOX[system].drmax;
   molm->p->y	+=	(ran1(seed) - 0.5) * BOX[system].drmax;
   molm->p->z	+=	(ran1(seed) - 0.5) * BOX[system].drmax;

   if (molm->nsites > 1) {
      for (i=1; i<molm->nsites; i++) {
         molm->s[i].alpha	=	AdjustAngle(molm->s[i].alpha + (ran1(seed)-0.5) * BOX[system].damax);
         molm->s[i].beta	=	AdjustAngle(molm->s[i].beta  + (ran1(seed)-0.5) * BOX[system].damax);
      }
   }
   MolCartesian(molm);			// calculate Cartesian coordinates of the molecule
   VAddSites(molm, 0, molm->nsites-1);

   arg	=	(v[system].tot - v_old.tot)/BOX[system].temp;

   if (( arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1 )) {
      av[system].acc_movemol	++;
   }
   else {
      v[system]		=	v_old;
      vir[system]	=	vir_old;
#ifdef CELL_LIST
      for (i=0; i<molm->nsites; i++)
         CL_Delete(molm, i);
#endif
      *molm	=	mol_old;
#ifdef CELL_LIST
      for (i=0; i<molm->nsites; i++)
         CL_Add(molm, i);
#endif
   }
   av[system].movemol	++;
   return	system;
}


double bondl_g(double L0, double kT)	// p(l) prop. to l^2 exp(-0.5*kv/kT*(l-L0)^2)
{
   double	std, r, l, kv;
   long		ready;

   if (FIXBONDL)			// bond length fixed
      return	L0;

   kv	=	type[0].KSTRETCH;
   std	=	sqrt(kT/kv);
   r	=	(L0 + 3.0 * std) * (L0 + 3.0 * std);
   ready	=	0;

   while (!ready) {
      l	=	gauss(std, L0);
      if (ran1(seed) <= l*l/r)
         ready	=	1;
   }
   return	l;
}


void bonda_tors(double theta0, double kT, double *theta, double *phi)
{
   double	std, k, bonda, tors, V;
   long		ready = 0;

   k	=	type[0].KBENDING;
   std	=	sqrt(kT/k);

   while (!ready) {
      bonda	=	gauss(std, theta0);
      if (ran1(seed) <= sin(bonda)) {			// p(theta) prop. to sin(theta)*exp(-0.5*k/kT*(t-t0)^2)
         while (!ready) {
            tors	=	(ran1(seed)-0.5) * 2 * M_PI;
            V	=	OPLS(tors+M_PI, type[0].TORSION[1], type[0].TORSION[2], type[0].TORSION[3]);

            if (ran1(seed) < exp(-V/kT))			// p(phi) prop. to exp(-Vtors/kT)
               ready	=	1;
         }
      }
   }
   *theta	=	bonda;
   *phi		=	tors;
   return; 
}


double bonda_g(double theta0, double kT)		// p(theta) prop. to sin(theta)*exp(-0.5*k/kT*(t-t0)^2)
{
   double	std, k, theta;
   long		ready = 0;

   k	=	type[0].KBENDING;
   std	=	sqrt(kT/k);

   while (!ready) {
      theta	=	gauss(std, theta0);
      if (ran1(seed) <= sin(theta))
	 ready	=	1;
   }
   return	theta; 
}


double tors(double kT)			// generate torsional angle phi, p(phi) prop. to exp(-Vtors/kT)
{
   double	tors, V;
   long		ready = 0;

   while(!ready) {
      tors	=	(ran1(seed) - 0.5) * 2 * M_PI;		// (-pi, pi) 
      V		=	OPLS(tors+M_PI, type[0].TORSION[1], type[0].TORSION[2], type[0].TORSION[3]);

      if (ran1(seed) < exp(-V/kT))
         ready	=	1;
   }
   return	tors;   
}


vector tors_bonda(molstruct *molm, long site)
{
   long		ready = 0;
   vector	dp;
   double	V;

   if ( (!V_BENDING && !V_TORSION) || (site<=1) )	// bending and torsion interaction off
      return	ranor();

   while (!ready) {
      dp		=	ranor();
      molm->p[site]	=	V_Add(molm->p+site-1, &dp);	// for energy calc.
      V			=	VTorsionSite(molm, site) + VBendingSite(molm, site);

      if (ran1(seed) < exp(-V/kT) ) 
	 ready		=	1;
   }

   return	dp;
}


long mccbmc()					// conf. biased Monte Carlo move
{
   long		site, ib, i;
   molstruct	*molm, mol_old;
   vstruct	v_old;				// dV is the energy of the regrown part
   wstruct	vir_old;			// dVIR is the virial of the regrown part
   double	Wold, Wnew;			// W is Rosenbluth factor for old and new conf.

   /* Pick up one molecule and determine the cut point */

   site	=	SiteSelect(&molm);		// pick up molm, longer mol w/ higher prob.
						// cut point will be redetermined below
   if (molm->fix)				// if molecule is fixed, end this move
      return	molm->box;

   site	=	(long) ( molm->nsites * (1-0.8*ran1(seed)) );	// cut length < 0.8 Lchain

   if (ran1(seed)>0.5) 				// pick up either end w/ equal prob.
      MolFlip(molm);

/*
   if (site < 0.5 * molm->nsites) {		// do not regrow more than half
      MolFlip(molm);
      site	=	molm->nsites-1 - site;
   }
*/

   /* Store information */
   ib		=	molm->box;
   mol_old	=	*molm;
   v_old	=	v[ib];
   vir_old	=	vir[ib];

   /* Calc. Rosenbluth factor, update energy and virial */

   Wold		=	grow("old", molm, site);		// system energy and virial updated
								// and calculate Rosenbluth factor
   Wnew		=	grow("new", molm, site);

   /* Determine acceptance */

   if (ran1(seed) < exp(Wnew-Wold)) {
      av[ib].acc_cbmc	++;			// coord. cell list have been updated in grow()
      cbmcsucc[site]	++;
/*
      for (i=site; i<molm->nsites; i++)		// update spherical coord.
         molm->s[i]	=	SiteSpherical(molm, i);
*/
   }
   else {
      v[ib]	=	v_old;			// restore energy and virial
      vir[ib]	=	vir_old;
#ifdef CELL_LIST
      for (i=site; i<molm->nsites; i++)		// get ready to recover original cell list
         CL_Delete(molm, i);
#endif
      *molm	=	mol_old;

#ifdef CELL_LIST
      for (i=site; i<molm->nsites; i++)		// recover original cell list, must after
         CL_Add(molm, i);			// restoring the atom positions
#endif
   }
   av[ib].cbmc	++;
   return	ib;
}


long mcvol()					// volume change move
{
   long		i, j, system, axis;
   molstruct	*moli;
   vstruct	v_old;
   wstruct	vir_old;
   vector	scale;
//   double	scale;
   double	volscale, volold, dV, arg;
 
   i		=	(int) (ran1(seed) * NMOLS);	// pick up one mol. in random
   system	=	(mol + i)->box;			// get the system id

   v_old	=	v[system];
   vir_old	=	vir[system];
   volold	=	BOX[system].vol;

   // each dimention changes independently
   scale.x	=	exp( BOX[system].dlmax * (ran1(seed)-0.5) );
   scale.y	=	exp( BOX[system].dlmax * (ran1(seed)-0.5) );
   scale.z	=	exp( BOX[system].dlmax * (ran1(seed)-0.5) );

   ChangeAxis(system, scale);	// scale coordinates and box sizes

   //scale	=	exp( BOX[system].dlmax * (ran1(seed)-0.5) );
   //ChangeVolume(system, scale);
   
   //CalcV();
   volscale	=	scale.x * scale.y * scale.z;
   CalcV_mcvol(volscale);
 
   dV	=	v[system].tot - v_old.tot;
   arg	=	(dV + P * (BOX[system].vol - volold)) / BOX[system].temp - (NMols[system]+1)*log(volscale);
   //arg	=	(dV + P * (BOX[system].vol - volold)) / BOX[system].temp - (NMols[system]+1)*3*log(scale);

   if (( arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1 )) {
      av[system].acc_vol	++;
   }
   else {
      v[system]		=	v_old;
      vir[system]	=	vir_old;
      scale.x		=	1.0/scale.x;
      scale.y		=	1.0/scale.y;
      scale.z		=	1.0/scale.z;
      ChangeAxis(system, scale);
      //ChangeVolume(system, 1.0/scale);
   }
   av[system].vol	++;
   return	system;
}


/* flip: Mavrantzas Macromolecules v31, 6310 (1998) */
/* added: April 28, 2009 */

long flip()
{
   long		ib, i, site, acceptmove;
   molstruct	*molm, mol_old;
   vstruct	v_old;
   wstruct	vir_old;
   vector	pL, pR, p0, u1, u2, u3, r1, r2, rp, rR;
   double	costheta, sintheta, phi, cosphi, sinphi, d1, d2;
   double	dV, arg;

   do {
      site	=	SiteSelect(&molm);		// select one site in random
   } while (site <= 0 || site >= molm->nsites-1);	// except the end sites

   ib		=	molm->box;
   mol_old	=	*molm;
   v_old	=	v[ib];
   vir_old	=	vir[ib];

   pR	=	molm->p[site+1];		// pR, pL and p0 are reference pts
   pL	=	molm->p[site-1];
   if (site-2 >= 0)	
      p0	=	molm->p[site-2];
   else
      p0	=	molm->p[site+2];

   if( Frame(&pL, &pR, &p0, &u1, &u2, &u3) )	return	0;	// build a frame

   r1	=	V_Subtr(&pR, &pL);
   r2	=	V_Subtr(molm->p+site, &pL);
   d1	=	sqrt(V_Dot(&r1, &r1));  
   d2	=	sqrt(V_Dot(&r2, &r2));  

   costheta	=	V_Dot(&r1, &r2)/(d1 * d2);
   sintheta	=	sqrt(1-costheta*costheta);

   rp.x		=	pL.x + d2 * costheta * u1.x;
   rp.y		=	pL.y + d2 * costheta * u1.y;
   rp.z		=	pL.z + d2 * costheta * u1.z;

   rR		=	V_Subtr(molm->p+site, &rp);
   phi		=	atan2(V_Dot(&rR, &u3), V_Dot(&rR, &u2));

   VDeleteSites(molm, site, site);

   phi		+=	(ran1(seed)-0.5) * 2 * pi/18;	// +- 10 degrees
   cosphi	=	cos(phi);
   sinphi	=	sin(phi);

   molm->p[site].x =	rp.x + d2 * sintheta * (cosphi * u2.x + sinphi * u3.x);
   molm->p[site].y =	rp.y + d2 * sintheta * (cosphi * u2.y + sinphi * u3.y);
   molm->p[site].z =	rp.z + d2 * sintheta * (cosphi * u2.z + sinphi * u3.z);

   VAddSites(molm, site, site);

   dV	=	v[ib].tot - v_old.tot;
   arg	=	dV/BOX[ib].temp;

   if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {	// accept
      av[ib].acc_flip	++;  
   }
   else {
      v[ib]	=	v_old;
      vir[ib]	=	vir_old;
#ifdef CELL_LIST
      CL_Delete(molm, site);
#endif
      *molm	=	mol_old;
#ifdef CELL_LIST
      CL_Add(molm, site);
      CL_Relink(molm);
#endif
   }
   av[ib].flip	++;
   return	ib;
}


/* End-rotation: change the last torsional angle randomly */
/* for chains of length at least 4 */

long end_rotation()
{
   long		ib, i, site, acceptmove = 0;
   molstruct	*molm, mol_old, mol_rev;
   vstruct	v_old;
   wstruct	vir_old;
   sphere	s;
   double	dV, arg;
   
   molm	=	mol + (int) (NMOLS * ran1(seed));	// select a mol in random
   if (molm->fix)
      return	molm->box;

   ib		=	molm->box;
   mol_old	=	*molm;
   v_old	=	v[ib];
   vir_old	=	vir[ib];
 
   if (ran1(seed) > 0.5) {				// pick up either end
      site	=	molm->nsites - 1;		// tail end
      s		=	SiteSpherical(molm, site);
   }
   else {
      site	=	0;				// head end
      for (i=0; i<4; i++) {
         mol_rev.p[i]		=	molm->p[site + 3 - i];
         mol_rev.parent[i]	=	i - 1;
      }
      s		=	SiteSpherical(&mol_rev, 3);
   }

   for (i=0; i<4; i++) { 				// try at most 4 times
      VDeleteSites(molm, site, site); 			// remove this site

      s.beta	=	(ran1(seed)-0.5) * 2 * M_PI;	// random torsional angle
      if (site==molm->nsites - 1) {			// tail end
         molm->p[site]	=	SiteCartesian(molm, site, s);
         molm->s[site].beta	=	s.beta;		// update s coordinates
      }
      else {						// head end
         mol_rev.p[3]	=	SiteCartesian(&mol_rev, 3, s);
         molm->p[0]	=	mol_rev.p[3];
         molm->s[3].beta	=	s.beta;		// update s coordinates 
      }

      VAddSites(molm, site, site);

      dV	=	v[ib].tot - v_old.tot;
      arg	=	dV/BOX[ib].temp;

      if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {	// accept
         acceptmove	=	1;
	 break;
      }
   }
/*
printf("counter = %d chainid = %d site = %d oldenergy = %f newenergy = %f \n", counter, molm-mol, site, v_old.tot, v[ib].tot);
*/
   if (acceptmove) {
//printf("accept endrotation\n");
      av[ib].acc_erot	++;  
   }
   else {
      v[ib]	=	v_old;
      vir[ib]	=	vir_old;
#ifdef CELL_LIST
      CL_Delete(molm, site);
#endif
      *molm	=	mol_old;
#ifdef CELL_LIST
      CL_Add(molm, site);
      CL_Relink(molm);
#endif
   }
   av[ib].erot	++;
   return	ib;
}

/* Reptation: move one site from one end to the other end of the same molecule */
/* To add: take care of difference b/w CH3 and CH2 */

long reptation()
{
   long		i_m, ib, flag, i, acceptmove=0;
   molstruct	*molm, mol_old;
   vstruct	v_old;
   wstruct	vir_old;
   sphere	s;
   double	dV, arg;
   vector	dp;

   av[ib].rep	++;
   molm	=	mol + (int) (NMOLS * ran1(seed));	// pick up a mol in random

   if (molm->fix)		return	molm->box;	// if mol is immobile 
   if (ran1(seed) > 0.5)	MolFlip(molm);		// pick up either end in random
/*
printf("reptation trial\n");
V_Print(molm->p[0]);
V_Print(molm->p[molm->nsites-2]);
V_Print(molm->p[molm->nsites-1]);
*/
   ib		=	molm->box;
   mol_old	=	*molm;
   v_old	=	v[ib];
   vir_old	=	vir[ib];

   i_m		=	molm->nsites-1;			// the atom on the tail
   if (molm->type[i_m-1] == molm->type[0])
      flag	=	0;
   else	
      flag	=	1;

   if (1==NTYPES && 1==flag)
      Exit("ensemble.c", "reptation", "NTYPES = 1 && flag = 1");

   for (i=0; i<4; i++) {			// try at most 4 different torsional angles

      VDeleteSites(molm, i_m-flag, i_m);	// substract the energy by the tail atom(s)
      if (flag)					// Because the change of type will also change
         VDeleteSites(molm, 0, 0); 		// the energy

      molm->nsites	--;
      MolFlip(molm);				// Flip mol., including CL_Relink()
      molm->nsites	++;
   
      // last MolFlip() might copy something undefined to molm site i_m because the 
      // molecule length changes, so we need to take special care here
      SiteCopy(molm, i_m, &mol_old, i_m, 0);	// copy p, s, flags, parent and cell
      molm->flags[i_m]	=	0;		// because the SiteCopy changed flags[i_m] to 1
      if (flag) {					
         molm->type[0]	=	mol_old.type[0];	// Swap types
         molm->type[i_m-1]	=	mol_old.type[i_m-1];
      }

      s		=	SiteSpherical(&mol_old, i_m);	// use spherical coord. to update coord.
      //s.d		+=	(ran1(seed) - 0.5) * BOX[ib].drmax;
      //s.alpha	+=	(ran1(seed) - 0.5) * 0.1 * pi;
      s.beta	=	(ran1(seed) - 0.5) * 2 * pi;
      molm->s[i_m]	=	s;
      molm->p[i_m]	=	SiteCartesian(molm, i_m, s);
/*
V_Print(molm->p[0]);
V_Print(molm->p[molm->nsites-2]);
V_Print(molm->p[molm->nsites-1]);
*/
      VAddSites(molm, i_m-flag, i_m);
      if (flag)
         VAddSites(molm, 0, 0);

      dV	=	v[ib].tot - v_old.tot;
      arg	=	dV/BOX[ib].temp;	
   /*
      printf("s.d=%f\ts.alpha=%f\ts.beta=%f\tdV=%f\n", s.d, s.alpha, s.beta, dV);
      if (dV>100) {
         for (i=0; i<molm->nsites; i++)
            printf("\t%f\t%f\t%f\n", (molm->p+i)->x, (molm->p+i)->y, (molm->p+i)->z);
         printf("LBOX = %f\n", BOX[ib].lbox);
      }
   */
      if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {
	 acceptmove	=	1;
	 break;
      }
   }	// up to 4 trials
   if (acceptmove) {
      av[ib].acc_rep	++;  
//printf("reptation accept, oldenergy= %f  newenergy = %f\n", v_old.tot, v[ib].tot);

      //MolInBox(molm);					// NO mapping in central box in moves
   }
   else {
      v[ib]	=	v_old;
      vir[ib]	=	vir_old;
#ifdef CELL_LIST
      CL_Delete(molm, i_m);
      if (flag) {
         CL_Delete(molm, i_m-1);
         CL_Delete(molm, 0);
      }
#endif
      *molm	=	mol_old;
#ifdef CELL_LIST
      CL_Add(molm, i_m);
      if (flag) {
         CL_Add(molm, i_m-1);
         CL_Add(molm, 0);
      }
      CL_Relink(molm);
#endif
   }
//molcheck(molm);
   return	ib;
}


long endbridge()	// end bridging move, copied from Pieter's code on 4/27/09
{
  vstruct		v_old;
  wstruct		vir_old;
  molstruct		*molm, *moln, mol1, mol2;
  neighborlist		list_forward, list_reverse;
  vector		p[7], p_reverse[7], dr;
  sphere		s[7];
  double		dvbeta, p_bias, J_forward, J_reverse;
  long			n, n_forward, n_reverse, site, reverse, 
  			system, flag, i, old_site, length;
  long			short_loop;

  static long		init=1;
  static long		Dnsites, maxnsites, minnsites;		// for polydispersity

  if (init) {
     Dnsites	=	(int) (Dpoly*NSITES/NMOLS);
     maxnsites	=	NSITES/NMOLS + Dnsites;
     minnsites	=	NSITES/NMOLS - Dnsites;
     init	=	0;
  }

  // Select tail
  
  molm	=	mol + (int) (NMOLS * ran1(seed));	// select one mol
  if (ran1(seed) > 0.5)					// select one end
     MolFlip(molm);

  system		= molm->box;
  ++av[system].eb;
//  ++av[system].CC0;
  
  // Determine target candidates
  if (NeighborList(molm, &list_forward))	// if any neighbor exists
  {
    n 		= (long) (list_forward.n*ran1(seed));	// choose one neighbor
    moln	= list_forward.mol[n];
    dr		= list_forward.dr[n];
    site	= list_forward.site[n];
    old_site	= site;
    reverse	= list_forward.reverse[n];
//    short_loop		= (moln->nsites<=17);
//    if (short_loop) ++av[system].CC1;
    length	=	(reverse ? moln->nsites-1-site : site)-3;
    if (length < minnsites || moln->nsites+molm->nsites-length >maxnsites)// polydispersity control
    {						// For control only
//      ++av[system].CC02;			// Filtered out in NeighborList
//      if (short_loop) ++av[system].CC12;
      return -1;
    }

    // Store old molecules
    oldmol[0]	= *molm;
    oldmol[1]	= *moln;
    v_old	= v[system];
    vir_old	= vir[system];

    // Create new molecules

    mol1		= *molm;		// Define acceptor
    if (reverse)
    {
      MolFlip(moln);				// Reverse donor
      site		= moln->nsites-1-site;
    }
    for (i=site-3; i<moln->nsites; ++i) 	// Sever moln
      SiteCopy(&mol2, i-(site-3), moln, i, -(site-3));
    mol2.box   		= moln->box;
//    mol2.fix		= moln->fix&(-1^1);
    mol2.nsites		= moln->nsites-(site-3);
    if (mol1.nsites+mol2.nsites>=MAXNMOLSITES)
    { 						// Exit on error
//      ++av[system].CCXX;
      return -1;
    }
    for (i=0; i<mol2.nsites; ++i) 		// Translate mol2 to mol1 frame
      mol2.p[i]		= V_Add(mol2.p+i, &dr);
    mol1		= MolAdd(&mol1, &mol2); // Combine mol1 and mol2
//    if (mol1.nsites<17)
//    {
//      fprintf(stderr, "MOVE = %d\n");
//      Exit("ensembles", "NextEndBridge", "nsites<17!");
//    }
    mol2		= *moln; 		// Copy remainder
    mol2.nsites		= site-3;
//    mol2.fix		&= -1^2;
    mol1.type[molm->nsites-1]			// Move terminator type
    			= mol2.type[site-4];
    mol2.type[site-4]	= molm->type[molm->nsites-1];
    
    RebridgeSetup(moln, site+1, 0, p_reverse, s);
    for (i=2; i<7; ++i)
      p[i]		= V_Add(p_reverse+i, &dr);
    p[0]		= molm->p[molm->nsites-2];
    p[1]		= molm->p[molm->nsites-1];
    
    if (n_forward = Rebridge(p, s))
    {						// Rebridge successful
      J_forward		= Jacobian(p, s);
      flag		= (molm->type[molm->nsites-1]!=
      				molm->type[molm->nsites-2]) ? 1 : 0;
/*
printf("another endbridging move, site = %d\n", site);
printf("Old %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
printf("molm\n");
for (i=molm->nsites-4; i<molm->nsites; i++)
   printf("%d %d %d %f %f %f\n", i, molm->parent[i], molm->flags[i], molm->p[i].x, molm->p[i].y, molm->p[i].z);      
printf("moln\n");
for (i=site-7; i<=site+2; i++)
   printf("%d %d %d %f %f %f\n", i, moln->parent[i], moln->flags[i], moln->p[i].x, moln->p[i].y, moln->p[i].z); 
*/
      VDeleteSites(moln, site-3-flag, site-1);	// Delete changed sites
      if (flag)
        VDeleteSites(molm, molm->nsites-1, molm->nsites-1);
/*
printf("Del done %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
printf("molm\n");
for (i=molm->nsites-4; i<molm->nsites; i++)
   printf("%d %d %d %f %f %f\n", i, molm->parent[i], molm->flags[i], molm->p[i].x, molm->p[i].y, molm->p[i].z);      
printf("moln\n");
for (i=site-7; i<=site+2; i++)
   printf("%d %d %d %f %f %f\n", i, moln->parent[i], moln->flags[i], moln->p[i].x, moln->p[i].y, moln->p[i].z); 
*/
      site		= molm->nsites;
      *molm		= mol1;
      *moln		= mol2;
//      MolConnect(molm);				// Repair connectivity and
//      MolConnect(moln);				// mol.nactive
      for (i=site-flag; i<=site+2; ++i)
      {
	//molm->flags[i]	^= -1;			// Switch changed sites off
        molm->flags[i]	=	0;		// Switch changed sites off because VDelete
        molm->p[i]	= p[i+2-site];		// Transcribe rebridge solution
      }
      if (flag) {
        //moln->flags[moln->nsites-1]		^= -1;
        moln->flags[moln->nsites-1]	=	0;
      }

/*
printf("New molecules:\n");
printf("molm\n");
for (i=site-4; i<=site+5; i++)
   printf("%d %d %d %f %f %f\n", i, molm->parent[i], molm->flags[i], molm->p[i].x, molm->p[i].y, molm->p[i].z);      
printf("moln\n");
for (i=moln->nsites-4; i<moln->nsites; i++)
   printf("%d %d %d %f %f %f\n", i, moln->parent[i], moln->flags[i], moln->p[i].x, moln->p[i].y, moln->p[i].z); 
*/

#ifdef CELL_LIST
      CL_Relink(molm);				// Relink lists
      CL_Relink(moln);
#endif
      VAddSites(molm, site-flag, site+2);	// Add changed sites
/*      if (flag) {
         VAddSites(molm, site-flag, site-flag);	// Add changed sites
printf("Add molm site 1 %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
      }
      VAddSites(molm, site, site);	// Add changed sites
printf("Add molm site 2 %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
      VAddSites(molm, site+1, site+1);	// Add changed sites
printf("Add molm site 3 %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
      VAddSites(molm, site+2, site+2);	// Add changed sites
printf("Add molm site 4 %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
*/

      if (flag) {
        VAddSites(moln, moln->nsites-1, moln->nsites-1);
//printf("Add moln end bead done %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
      }
/*
printf("Add all done %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
printf("molm\n");
for (i=0; i<molm->nsites; i++)
   printf("%d %d %d %f %f %f\n", i, molm->parent[i], molm->flags[i], molm->p[i].x, molm->p[i].y, molm->p[i].z);      
printf("moln\n");
for (i=0; i<moln->nsites; i++)
   printf("%d %d %d %f %f %f\n", i, moln->parent[i], moln->flags[i], moln->p[i].x, moln->p[i].y, moln->p[i].z); 
*/
      NeighborList(moln, &list_reverse);
      J_reverse		= Jacobian(p_reverse, s);
      n_reverse		= Rebridge(p_reverse, s)+1;
      list_reverse.n	= list_reverse.n ? list_reverse.n : 1;
						// Round off error warranty
      p_bias		= ((double) n_forward*list_forward.n*J_forward)/
			    ((double) n_reverse*list_reverse.n*J_reverse);
      dvbeta		= (v[system].tot-v_old.tot)/BOX[system].temp;

      // Check acceptance

      if (p_bias*exp(-dvbeta)>=ran1(seed))
      { 					// Accept
	MolInBox(molm);
        MolInBox(moln);
	++av[system].acc_eb;
//printf("accept %d\n",counter);
//exit(1);
      }
      else 
      { 					// Reject and restore
#ifdef CELL_LIST
	for (i=site-flag; i<=site+2; ++i)
	  CL_Delete(molm, i);			// Unregister new sites
	if (flag)
	  CL_Delete(moln, moln->nsites-1);
#endif
	v[system]	= v_old;
	vir[system]	= vir_old;
	*molm		= oldmol[0];
	*moln		= oldmol[1];
#ifdef CELL_LIST
        if (reverse)				// Register old sites
	  for (i=old_site+1; i<=old_site+3+flag; ++i)
	    CL_Add(moln, i);
	else
	  for (i=old_site-3-flag; i<old_site; ++i)
	    CL_Add(moln, i);
	if (flag)
	  CL_Add(molm, molm->nsites-1);
	CL_Relink(molm);
	CL_Relink(moln);
#endif
//	++av[system].CC03; 			// Energetic rejection
//	if (short_loop) ++av[system].CC13;
      }
    }
    else
    { 						// Restore
//      ++av[system].CC04; 			// Rebridging failed
//      if (short_loop) ++av[system].CC14;
      if (reverse)
      {
	*moln		= oldmol[1];
#ifdef CELL_LIST
	CL_Relink(moln);
#endif
      }
    }
  }
//  else 						// No bridging candidates
//    ++av[system].CC05; 
  return system;
}

///////////////////////////////////////////////
/* Double bridging move, added on 12/06/2009 */
///////////////////////////////////////////////

long doublebridge()
{
  vstruct		v_old;
  wstruct		vir_old;
  molstruct		*moli, *molj, mol1, mol2;
  neighborlist		list_forward, list_reverse;
  vector		p[7], pj[7], pi_rev[7],  pj_rev[7], dr;
  sphere		si[7], sj[7];
  double		dvbeta, p_bias, J_forward, J_reverse;
  long			n, ni_forward, ni_reverse, nj_forward, nj_reverse, 
			sitei, sitej, sitei2, sitej2,
			reverse, system, flag, i, old_site, length;
  long			N;

  do {
    sitei	=	SiteSelect(&moli);
  } while (sitei < 2 || sitei > moli->nsites-3) ;	// >= two bonds away from the end

  system	= moli->box;
  ++av[system].db;
  N	=	moli->nsites; 	// length, monodisperse

  if (DB_NeighborList(moli, sitei, &list_forward)) {	// if any neighbor exists

    // pick up a neighbor in random
    n		= (long) (list_forward.n * ran1(seed));
    molj	= list_forward.mol[n];
    sitej	= list_forward.site[n];
    dr		= list_forward.dr[n];			// vector pointing from molj to moli
    reverse	= list_forward.reverse[n];		// (0-3) four possibilities
    
    switch (reverse) {
      case 0:	sitei2	= sitei+4;	sitej2	= sitej+4;	break;
      case 1:	sitei2	= sitei-4;	sitej2	= sitej+4;	break;
      case 2:	sitei2	= sitei+4; 	sitej2	= sitej-4;	break;
      case 3:	sitei2	= sitei-4;	sitej2	= sitej-4;	break;
      default:	break;
    }

    // store old information
    oldmol[0]	=	*moli;
    oldmol[1]	=	*molj;
    v_old	=	v[system];
    vir_old	=	vir[system];

    // bridging setup, pi, pj for forward, pi_rev, pj_rev for reverse

    RebridgeSetup(moli, MAX(sitei, sitei2)+1, 0, pi_rev, si);
    RebridgeSetup(molj, MAX(sitej, sitej2)+1, 0, pj_rev, sj);

    for (i=2; i<7; i++) {
       p[i]	=	pi_rev[i];	// do not use pi[] because pi=M_PI
       pj[i]	=	pj_rev[i];
    }
    if (reverse==0 || reverse==3) {
      p[6]	=	V_Add(pj_rev+0, &dr);
      p[5]	=	V_Add(pj_rev+1, &dr);
      pj[0]	=	V_Subtr(pi_rev+6, &dr);
      pj[1]	=	V_Subtr(pi_rev+5, &dr);
    }
    else if (reverse==1 || reverse==2) {
      p[6]	=	V_Add(pj_rev+6, &dr);
      p[5]	=	V_Add(pj_rev+5, &dr);
      pj[6]	=	V_Subtr(pi_rev+6, &dr);
      pj[5]	=	V_Subtr(pi_rev+5, &dr);
    }

    // do bridging 
    if ((ni_forward=Rebridge(p, si)) && (nj_forward=Rebridge(pj, si))) {// bridging successful 
      J_forward	= new_Jacobian(p, si) * new_Jacobian(pj, sj);
      
      // delete changed sites
      VDeleteSites(moli, MIN(sitei, sitei2)+1, MAX(sitei, sitei2)-1);
      VDeleteSites(molj, MIN(sitej, sitej2)+1, MAX(sitej, sitej2)-1);

      // update molecules
      if (reverse==1 || reverse==2) {
         for (i=0; i<N-MAX(sitei, sitei2); i++)
            SiteCopy(&mol1, i, moli, MAX(sitei, sitei2)+i, -MAX(sitei, sitei2));
         mol1.box	=	moli->box;
         mol1.nsites	=	N-MAX(sitei, sitei2);

         for (i=0; i<N-MAX(sitej, sitej2); i++)
            SiteCopy(&mol2, i, molj, MAX(sitej, sitej2)+i, -MAX(sitej, sitej2));
         mol2.box	=	molj->box;
         mol2.nsites	=	N-MAX(sitej, sitej2);

         for (i=0; i<mol2.nsites; i++) {
            mol2.p[i]	=	V_Add(mol2.p+i, &dr);
            SiteCopy(moli, MAX(sitei, sitei2)+i, &mol2, i, MAX(sitei, sitei2));
         }

         for (i=0; i<mol1.nsites; i++) {
            mol1.p[i]	=	V_Subtr(mol1.p+i, &dr);
            SiteCopy(molj, MAX(sitej, sitej2)+i, &mol1, i, MAX(sitej, sitej2));
         }
      }
      else if (reverse==0 || reverse==3) {
         for (i=0; i<N-MAX(sitei, sitei2); i++) 
            SiteCopy(&mol1, i, moli, MAX(sitei, sitei2)+i, -MAX(sitei, sitei2));
         mol1.box	=	moli->box;
         mol1.nsites	=	N-MAX(sitei, sitei2);
         MolFlip(&mol1);

         for (i=0; i<MIN(sitej, sitej2)+1; i++)
            SiteCopy(&mol2, i, molj, i, 0);
         mol2.box	=	molj->box;
         mol2.nsites	=	MIN(sitej, sitej2)+1;
         MolFlip(&mol2);

         for (i=0; i<mol1.nsites; i++) {
            mol1.p[i]	=	V_Subtr(mol1.p+i, &dr);
            SiteCopy(molj, i, &mol1, i, 0);
         }
         for (i=0; i<mol2.nsites; i++) {
            mol2.p[i]	=	V_Add(mol2.p+i, &dr);
            SiteCopy(moli, MAX(sitei, sitei2)+i, &mol2, i, MAX(sitei, sitei2));
         }
      }
      
      moli->p[MIN(sitei, sitei2)+1]	=	p[2];
      moli->p[MIN(sitei, sitei2)+2]	=	p[3];
      moli->p[MIN(sitei, sitei2)+3]	=	p[4];
      if (reverse==0 || reverse==3) {
         molj->p[N-MAX(sitei, sitei2)]		=	pj[2];
         molj->p[N-MAX(sitei, sitei2)+1]	=	pj[3];
         molj->p[N-MAX(sitei, sitei2)+2]	=	pj[4];
      }
      else if (reverse==1 || reverse==2) {
         molj->p[MIN(sitej, sitej2)+1]	=	pj[2];
         molj->p[MIN(sitej, sitej2)+2]	=	pj[3];
         molj->p[MIN(sitej, sitej2)+3]	=	pj[4];
      }

#ifdef CELL_LIST
      CL_Relink(moli);		// relink cell list because of the chain and site
      CL_Relink(molj);          // identity change
#endif

      // Add changed sites

      VAddSites(moli, MIN(sitei, sitei2)+1, MAX(sitei, sitei2)-1);
      VAddSites(molj, MIN(sitej, sitej2)+1, MAX(sitej, sitej2)-1);

      DB_NeighborList(molj, sitej2, &list_reverse);		// new conf. neighbor list
      J_reverse	 = new_Jacobian(pi_rev, si) * new_Jacobian(pj_rev, sj);
      ni_reverse = Rebridge(pi_rev, si) + 1;
      nj_reverse = Rebridge(pj_rev, sj) + 1;

     // if (list_reverse.n==0) { printf("reverse neighbor==0, exit\n"), exit(1); }
      list_reverse.n = list_reverse.n ? list_reverse.n : 1;	// in case =0, ?
       
      p_bias	=  (double) (ni_forward * nj_forward * list_forward.n * J_forward);
      p_bias	/= (double) (nj_reverse * nj_reverse * list_reverse.n * J_reverse);
      dvbeta	=  (v[system].tot - v_old.tot)/BOX[system].temp;

      // Check acceptance
      if (p_bias*exp(-dvbeta) >= ran1(seed)) {	// accept
        ++av[system].acc_db;
      }
      else {					// reject and restore

#ifdef CELL_LIST
        // unregister new sites
        if (reverse==0 && reverse==3) {
           for (i=1; i<=3; i++) {
              CL_Delete(moli, MIN(sitei, sitei2)+i);
              CL_Delete(molj, N-MIN(sitei, sitei2)-1+i);
           }
        }
        else if (reverse==1 && reverse==2) {
           for (i=1; i<=3; i++) {
              CL_Delete(moli, MIN(sitei, sitei2)+i);
              CL_Delete(molj, MIN(sitej, sitej2)+i);
           }
        }
#endif        
        v[system]	= v_old;
        vir[system]	= vir_old;
        *moli		= oldmol[0];
        *molj		= oldmol[1];

#ifdef CELL_LIST
        // register old sites
        for (i=1; i<=3; i++) {
           CL_Add(moli, MIN(sitei, sitei2)+i);
           CL_Add(molj, MIN(sitej, sitej2)+i);
        }
        CL_Relink(moli);
        CL_Relink(molj);
#endif
      }
    }	// bridge successful
  }	// any neighbor exist
  return 	system;
}


long idr()	// intramolecular double rebridging move
{
  vstruct		v_old;
  wstruct		vir_old;
  molstruct		*moli, *molj, mol1, mol2;
  neighborlist		list_forward, list_reverse;
  vector		p[7], pj[7], pi_rev[7],  pj_rev[7], dr;
  sphere		si[7], sj[7];
  double		dvbeta, p_bias, J_forward, J_reverse;
  long			n, ni_forward, ni_reverse, nj_forward, nj_reverse, 
			sitei, sitej, sitei2, sitej2,
			reverse, system, flag, i, old_site, length,
                        head, tail;
  long			N;

  do {
    sitei	=	SiteSelect(&moli);
  } while (sitei < 2 || sitei > moli->nsites-3) ;	// >= two bonds away from the end

  system	= moli->box;
  ++av[system].idr;
  N	=	moli->nsites; 	// length, monodisperse

  if (IDR_NeighborList(moli, sitei, &list_forward)) {	// if any neighbor exists

    // pick up a neighbor in random
    n		= (long) (list_forward.n * ran1(seed));
    molj	= list_forward.mol[n];
    if (molj!=moli)	{printf("idr_error not the same molecule.\n"); exit(1);}
    sitej	= list_forward.site[n];
    dr		= list_forward.dr[n];
    reverse	= list_forward.reverse[n];		// (0-1) two possibilities
    sitei2	=	sitei + 4 - reverse * 8;	// if reverse=0, site2 = site +4
    sitej2	=	sitej + 4 - reverse * 8;	// if reverse=1, site2 = site -4

    // store old information
    oldmol[0]	=	*moli;
    v_old	=	v[system];
    vir_old	=	vir[system];

    // bridging setup, pi, pj for forward, pi_rev, pj_rev for reverse
    head	=	MIN(MIN(sitei, sitei2), MIN(sitej, sitej2));
    tail	=	MAX(MAX(sitei, sitei2), MAX(sitej, sitej2));

    RebridgeSetup(moli, head+5, 0, pi_rev, si);
    RebridgeSetup(moli, tail+1, 0, pj_rev, sj);
    
    for (i=2; i<7; i++) {
        p[i]	=	pi_rev[i];
        pj[i]	=	pj_rev[i];
    }
    p[6]	=	pj_rev[0];
    p[5]	=	pj_rev[1];
    pj[0]	=	pi_rev[6];
    pj[1]	=	pi_rev[5];

    // do bridging 
    if ((ni_forward=Rebridge(p, si)) && (nj_forward=Rebridge(pj, si))) {// bridging successful 
      J_forward	= new_Jacobian(p, si) * new_Jacobian(pj, sj);
      
      // delete changed sites
      VDeleteSites(moli, head+1, head+3);
      VDeleteSites(molj, tail-3, tail-1);

      // update molecules
      for (i=head+1; i<=tail-1; i++) {
        SiteCopy(&mol1, i-head-1, moli, i, -(head+1));     
      }
      mol1.box		=	moli->box;
      mol1.nsites	=	tail-head-1;
      MolFlip(&mol1);

      for (i=head+1; i<=tail-1; i++) {
        SiteCopy(moli, i, &mol1, i-head-1, head+1);
      }
      moli->p[head+1]	=	p[2];
      moli->p[head+2]	=	p[3];
      moli->p[head+3]	=	p[4];
      moli->p[tail-3]	=	pj[2];
      moli->p[tail-2]	=	pj[3];
      moli->p[tail-1]	=	pj[4];

#ifdef CELL_LIST
      CL_Relink(moli);		// relink cell list because of site identity change
#endif
      // Add changed sites

      VAddSites(moli, head+1, head+3);
      VAddSites(moli, tail-3, tail-1);

      IDR_NeighborList(moli, sitej2, &list_reverse);		// new conf. neighbor list
      J_reverse	 = new_Jacobian(pi_rev, si) * new_Jacobian(pj_rev, sj);
      ni_reverse = Rebridge(pi_rev, si) + 1;
      nj_reverse = Rebridge(pj_rev, sj) + 1;
      list_reverse.n = list_reverse.n ? list_reverse.n : 1;

      p_bias	=  (double) (ni_forward * nj_forward * list_forward.n * J_forward);
      p_bias	/= (double) (nj_reverse * nj_reverse * list_reverse.n * J_reverse);
      dvbeta	=  (v[system].tot - v_old.tot)/BOX[system].temp;

      // Check acceptance
      if (p_bias*exp(-dvbeta) >= ran1(seed)) {	// accept
        ++av[system].acc_db;
      }
      else {					// reject and restore

#ifdef CELL_LIST
        // unregister new sites
        for (i=1; i<=3; i++) {
           CL_Delete(moli, head+i);
           CL_Delete(moli, tail-i);
        }
#endif        
        v[system]	= v_old;
        vir[system]	= vir_old;
        *moli		= oldmol[0];

#ifdef CELL_LIST
        // register old sites
        for (i=1; i<=3; i++) {
           CL_Add(moli, head+i);
           CL_Add(moli, tail-i);
        }
        CL_Relink(moli);
#endif
      }
    }	// bridge successful
  }	// any neighbor exist
  return 	system;
}


////////////////////////////////////////////////////////////////////
/* ConRot (rebridging) move, copied from Pieter's code on 4/20/09 */
////////////////////////////////////////////////////////////////////

long rebridge()	// 
{
  vstruct		v_old;
  wstruct		vir_old;
  molstruct		*molm;
  vector		p[7], p_old[7];
  sphere		s[7];
  double		dvbeta, J_forward, J_reverse, p_bias;
  long			i, ib, i_start, nr, system, flag;
  
  do {
     i_start	=	SiteSelect(&molm);	// select mol and starting site
						// [i_start, i_start+6]
  } while (i_start<0 || i_start+6>molm->nsites-1);
  ib	=	molm->box;

  if (nr = !RebridgeSetup(molm, i_start+6, 0, p, s))
  {

/*for debug
    for (i=0; i<7; i++){
       p_old[i]	=	p[i];
    }
*/
    J_reverse		= Jacobian(p, s);
    nr			= Rebridge(p, s);
  }
  if (nr)
  { 						// Successful rebridging
    J_forward		= Jacobian(p, s);
    oldmol[0]		= *molm;
    v_old		= v[ib];
    vir_old		= vir[ib];
    VDeleteSites(molm, i_start+2, i_start+4);
    
    for (i=i_start+2; i<=i_start+4; ++i)
      molm->p[i]	= p[i-i_start];		// Transcribe new positions

    VAddSites(molm, i_start+2, i_start+4);
   
    p_bias		= J_forward/J_reverse;	// Jacobian to keep detailed balance
    dvbeta		= (v[ib].tot - v_old.tot) / BOX[ib].temp;
  
    // Check acceptance

/*printf("rebridge trial\n");
printf("J_forward = %f J_reverse = %f dvbeta = %f\n", J_forward, J_reverse, dvbeta);     
printf("counter = %d chainid = %d i_start = %d oldenergy = %f newenergy = %f \n", counter, molm-mol, i_start, v_old.tot, v[ib].tot);
for (i=0; i<7; i++) {
   V_Print(p_old[i]);
   V_Print(p[i]);
}
*/
    if (p_bias*exp(-dvbeta)>=ran1(seed)) {	// Accept
       ++av[ib].acc_re;
//printf("rebridge accept! old energy = %f  new energy = %f\n", v_old.tot, v[ib].tot); 
    }
    else {					// Reject
#ifdef CELL_LIST
      for (i=2; i<=4; ++i)			// Unregister new sites
        CL_Delete(molm, i_start+i);
#endif 
      v[ib]		= v_old;
      vir[ib]		= vir_old;
      *molm		= oldmol[0];
#ifdef CELL_LIST
      for (i=2; i<=4; ++i)			// Register old sites
        CL_Add(molm, i_start+i);
#endif
    }
  }
//molcheck(molm);
  ++av[ib].re;
  return ib;
}


/* Gibbsmove() includes volume exchange move and molecule swap move */

void Gibbsmove()		//10/28/2007
{
   long		i, system, in, out, m, n;
   long		choice = NGIBBSVOL + NSWAP;
   systemstruct	BOX_old[2];
   vstruct	v_old[2];
   wstruct	vir_old[2];
   molstruct	*molm, mol_old;

   double	vo1, vo2, vn1, vn2, lnvn, voltot, scale[2];
   double	dVcorr, Wold, Wnew, arg;

   /* Pick up two systems in random */

   out		=	(int) (NSYSTEMS * ran1(seed));
   if (2==NSYSTEMS)
      in	=	1 - out;
   else if (NSYSTEMS > 2)
      while ( (in = (int) (NSYSTEMS * ran1(seed))) == out );

   /* Store system information */
   
   for (i=0; i<2; i++) {
      system		=	(0==i ? in : out);		// in -> 0, out -> 1
      v_old[i]		=	v[system];
      vir_old[i]	=	vir[system];
   }
   
   /* Volume exchange  OR  Swap molecule */

   if ((int) (choice * ran1(seed)) < NGIBBSVOL) {		// Volume exchange move

      vo1	=	BOX[in].vol;		// old volume of in 
      vo2	=	BOX[out].vol;		// old volume of out
      voltot	=	vo1 + vo2;		// voltot remain unchanged
      lnvn	=	log(vo1/vo2) + (ran1(seed)-0.5) * 3 * BOX[in].dlmax;	// random walk in ln(vol1/vol2)
      
      vn1	=	voltot * exp(lnvn) / (1.0+exp(lnvn));
      vn2	=	voltot - vn1;
  
      scale[0]	=	pow(vn1/vo1, 1.0/3.0);
      scale[1]	=	pow(vn2/vo2, 1.0/3.0);
      ChangeVolume(in, scale[0]);		// box dimension and mol. coord.
      ChangeVolume(out, scale[1]);		// also include cell list update

      CalcV();					// recalculate ALL energy

      arg	=	((v[in].tot - v_old[0].tot) + (v[out].tot - v_old[1].tot)) / BOX[in].temp;
      arg	-=	(NMols[in]+1) * 3 * log(scale[0]) + (NMols[out]+1) * 3 * log(scale[1]);

      if ((arg>0 ? (arg<75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {
         av[in].acc_vol		++;
         av[out].acc_vol	++;
      }
      else {
         v[in]		=	v_old[0];	// restore energy
         v[out]		=	v_old[1];
         vir[in]	=	vir_old[0];
         vir[out]	=	vir_old[1];
         ChangeVolume(in, 1.0/scale[0]);
         ChangeVolume(out, 1.0/scale[1]);
      }
      av[in].vol	++;
      av[out].vol	++;
      return;
   }
   else {					// Molecule swap move
      /* Check whether the out box is empty */

      if (NMols[out]<=0) {
         av[in].swap	++;
	 return;
      }

//   printf("counter=%d\tNMols[out]=%d\tNSites[out]=%d\n", counter, NMols[out], NSites[out]);
//   fflush(stdout);

      /* Pick up a molecule from the out box */

      n		=	(long) (NMols[out] * ran1(seed) +1);
      m		=	-1;
      while (m<NMOLS && n) {
         m	++;
         if (mol[m].box == out)
            n	--;
      }
      if (m<0 || m>=NMOLS) 
	 Exit("ensemble", "Gibbsmove", "molecule to swap invalid");

      molm	=	mol + m;

      if (ran1(seed)>0.5)			// because W factor depends on 
	 MolFlip(molm); 			// the retrace direction
      mol_old	=	*molm;			// store molecule info.

      /* Update system energy and virial of out, calc. Rosenbluth factor */

      Wold	=	grow("old", molm, 0);
      molm->box	=	in;			// move molecule from out to in
      Wnew	=	grow("new", molm, 0);

      /* Update NMols and NSites */
      
      NSites[out]	-=	molm->nsites;
      NMols[out]	--;
      NSites[in]	+=	molm->nsites;
      NMols[in]		++;

      /* Update long range correction */

      if (V_LJLRC) {
         CalcVCorr();
         v[in].nonbonded	+=	v[in].corr - v_old[0].corr;
         v[in].tot		+=	v[in].corr - v_old[0].corr;
         v[out].nonbonded	+=	v[out].corr - v_old[1].corr;
         v[out].tot		+=	v[out].corr - v_old[1].corr;
      }
      dVcorr	=	v[in].corr + v[out].corr - v_old[0].corr - v_old[1].corr;

      arg	=	exp(Wnew-Wold-dVcorr/BOX[in].temp);
      arg	*=	BOX[in].vol * (NMols[out]+1) / (BOX[out].vol * NMols[in]);

      if ( ran1(seed) < arg ) {
         av[in].acc_swap	++;
      }
      else {
         v[in]		=	v_old[0];
         v[out]		=	v_old[1];
         vir[in]	=	vir_old[0];
         vir[out]	=	vir_old[1];
         NSites[out]	+=	molm->nsites;
         NMols[out]	++;
         NSites[in]	-=	molm->nsites;
         NMols[in]	--;
#ifdef CELL_LIST
	 for (i=0; i<molm->nsites; i++)
            CL_Delete(molm, i);
#endif
         *molm		=	mol_old;	// restore molm's identity as well
#ifdef CELL_LIST
         for (i=0; i<molm->nsites; i++)
 	    CL_Add(molm, i);
#endif
      }
      av[in].swap	++;
      return;
   }
}	/* Gibbsmove() */


void ParentCheck()	// check parent relation
{
   long		i;
   molstruct	*moli;

   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         if (moli->parent[i] != i-1)
            fprintf(foutput, "mol[%d] site[%d] parent site = %d\n", moli-mol, i, moli->parent[i]);
}



#define NVFIX		5 	// first NVFIX moves with NO volume change
				// in case init conf. far away from equilibrium volume
#define E_MAXCHOICE	9	// total types of moves except volume change
void NextMove()
{
   long		nptchoice = (long) (ran1(seed) * (NSITES + NVOLCHANGE));
   long		choice, total=0, system=-1;
   long		n[E_MAXCHOICE] = {NDISPLACE, NREPTATION, NENDROT, NCBMC, NENDBR, NREBR, 
					NFLIP, NDB, NIDR};
   // remember to update the value of macro E_MAXCHOICE 

   if (E_NPT && !NVOLCHANGE) {
      printf("NPT but no volume change move!\n");
      exit(1);
   }
   if (nptchoice > NSITES && counter>NVFIX) {
      mcvol();
   }
   else {

   // 	determine the move type in NVT ensemble
   //	
   //	0: local site displacement
   //	1: reptation
   //	2: end-mer rotation
   //	3: configuration biased Monte Carlo
   //	4: end-bridging move
   //	5: rebridging move
   //	6: flip
   //   7: double bridging
   //   8: intramolecular double rebridging

      for (choice=0; choice<E_MAXCHOICE; choice++)
         total	+=	n[choice];
      total	=	(long) (ran1(seed) * total);
      choice	=	0;
      while ( (total>=n[choice]) && (choice<E_MAXCHOICE) )
         total	-=	n[choice++];

      switch (choice) {
         case 0:	system	=	mcmove();	break;
         case 1:	system	=	reptation();	break;
	 case 2:	system	=	end_rotation();	break;
         case 3:	system	=	mccbmc();	break;
         case 4:	system	=	endbridge();	break;
	 case 5:	system	=	rebridge();	break;
	 case 6:	system	=	flip();		break;
         case 7:	system  =	doublebridge();	break;
         case 8:	system	=	idr();		break;
         default:	break;
      }
   }
}


void Cycle()		// one MC cycle contains one move per particle on average
{
   long		i, j, k;
   double	arg=0.0;
   double	choice;
   long		nmoves;
   double	ReCoor[MAXNSYSTEMS], oldReCoor[MAXNSYSTEMS];	// reaction coordinates

   StoreMols();		//store system and molecule info. before a sequence of MC moves

   for (k=0; k<SEQUENCE; k++) {		//a sequence of trial MC moves, containing SEQUENCE MC cycles
      for (i=0; i<NSITES; i++) {
         NextMove();
      }
   }

   /* Should I put Normal samplings after each sequence here or in main.c ?*/

   /* BIAS SAMPLING */

   if (fabs(kP2) > ZERO) {

      /* Store old reaction coordinate value */

      for (i=0; i<NSYSTEMS; i++) {
	 switch (dynvar) {
	    case 1:	oldReCoor[i]	=	MAXSIZE[i];	// MAXSIZE as reaction coordinate
			break;
	    case 2:	oldReCoor[i]	=	Q6[i];		// Q6 as reaction coordinate
			break;						
	    case 3:	oldReCoor[i]	=	P2M[i];		// P2 modified
			break;
	    case 4:	oldReCoor[i]	=	P2[i];		// P2
			break;
            case 5:	oldReCoor[i]	=	oldNSites[i]/oldBOX[i].vol;	// Rho
			break;
	    case 6:	oldReCoor[i]	=	MAXSIZE[i];	// nmaxp2 as reaction coordinate
			break;
	    default:	break;
         }
      }

      /* Sample the reaction coordinate */

      switch (dynvar) {
	 case 1:	Find_Nuclei(dynvar);
			break;
	 case 2:	Calc_Qlm(6);				// calc. Q6
			break;
	 case 3:	SampleP2All();				// calc. global P2, P2M and local p2
			break;
	 case 4:	SampleP2();
			break;
	 case 5:	break;
	 case 6:	SampleP2All();
			Find_Nuclei_p2(1);
			break;
	 default:	break;
      }

      SampleSpherical();

      for (i=0; i<NSYSTEMS; i++) {
	 switch (dynvar) {
	    case 1:	ReCoor[i]	=	MAXSIZE[i];
			break;
	    case 2:	ReCoor[i]	=	Q6[i];
			break;
	    case 3:	ReCoor[i]	=	P2M[i];
			break;
	    case 4:	ReCoor[i]	=	P2[i];
			break;
	    case 5:	ReCoor[i]	=	NSites[i]/BOX[i].vol;
			break;
	    case 6:	ReCoor[i]	=	MAXSIZE[i];
			break;
	    default:	break;
	 }
      }

      /* calculate bias potential */

      arg	=	(oldReCoor[0]-P2middle) * (oldReCoor[0]-P2middle) - (ReCoor[0]-P2middle) * (ReCoor[0]-P2middle);
      arg	*=	(-0.5 * kP2);

/*
      for (i=0; i<NSYSTEMS; i++)		// added 5/1/08, allow only one nucleus in the system at any time
         if (dynvar==1)
	    if (sizedist[MAXSIZE[i]]>1 || secondNmax[i]>0)
	       arg	=	10000;
*/
   }
   else {
      SampleP2();					// still do some sampling
      Find_Nuclei(dynvar);
      arg	=	0.0;				// always accept sequence if NO bias
   }

   /* Determine acceptance of the sequence */

   if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {
      av[0].acc_seq	++;
   }
   else {							// reject this sequence
      for (i=0; i<NSYSTEMS; i++) {
	 switch (dynvar) {
	    case 1:	MAXSIZE[i]	=	oldReCoor[i];
			break;
	    case 2:	Q6[i]		=	oldReCoor[i];
			break;
	    case 3:	P2M[i]		=	oldReCoor[i];
			break;
	    case 4:	P2[i]		=	oldReCoor[i];
			break;
	    case 5:	break;
	    case 6:	MAXSIZE[i]	=	oldReCoor[i];
			break;
	    default:	break;
	 }
      }
      RestoreMols();
#ifdef CELL_LIST
      CL_Destroy();
      CL_Build();
#endif
   }
   av[0].seq	++;

   //ParentCheck();

   if (!mod(counter, 1024))
      Adjust_Stepsize();				// Adjust MC move step sizes

   return;
}


void InitEnsembles()
{
   ResetAcceptance();
}


#ifdef TEST
void gibbsvol()				// Gibbs ensemble change volume move, F&S algorithm 18
{
   systemstruct	BOX_old[2];			// two boxes are involved in Gibbs volume change
   vstruct	v_old[2], V;
   wstruct	vir_old[2], VIR;

   long		i, in, out, ibox;
   double	vo1, vo2, lnvn, voltot;
   double	scale[2], scale3, scale6;
   double	dV, arg;

   if (NBOX==2) {
      in	=	0;
      out	=	1;
   }
   else if (NBOX>2) {
      in	=	(int) (NBOX * ran1(seed));
      while ( (out = (int) (NBOX*ran1(seed))) == in );
   }

   for (i=0; i<2; i++) {				// Store information
      ibox	=	(i==0 ? in : out);		// in -> 0; out -> 1

      BOX_old[i]	=	BOX[ibox];
      v_old[i]		=	v[ibox];
      vir_old[i]	=	vir[ibox];      
   }

   vo1		=	BOX[in].vol;					// old volume
   vo2		=	BOX[out].vol;
   voltot	=	vo1 + vo2;					// total volume unchanged
   lnvn		=	log(vo1/vo2) + (ran1(seed)-0.5) * 3 * GDLMAX;	// random walk in ln vol1/vol2	
 
   BOX[in].vol	=	voltot * exp(lnvn)/ (1.0+exp(lnvn));		// assign new volume
   BOX[out].vol	=	voltot - BOX[in].vol;				// total volume remains the same

   BOX[in].lbox		=	pow(BOX[in].vol, 1.0/3);		// calculate LBOX
   BOX[out].lbox	=	pow(BOX[out].vol, 1.0/3);	
   scale[0]		=	BOX[in].lbox / BOX_old[0].lbox;
   scale[1]		=	BOX[out].lbox / BOX_old[1].lbox;

   for (i=0; i<NPARTS; i++) {						// scale particle coordinates
      ibox		=	part[i].box;

      if (ibox == in || ibox == out) {
         part[i].p	=	V_Mult(scale[(ibox==in ? 0 : 1)], &(part[i].p));
#ifdef VERLET_LIST
         part[i].pv	=	V_Mult(scale[(ibox==in ? 0 : 1)], &(part[i].pv));	
#endif
      }
   }

   if (SCALECUTOFF) {						// scale cutoff radii with box dimension
      for (i=0; i<2; i++) {

//	 ibox	=	in * (1-i) + out * i;			// i=0 -> ibox=in,  i=1 -> ibox=out
	 ibox	=	(i == 0 ? in : out);

	 BOX[ibox].rc		*=	scale[i];
	 BOX[ibox].rv		*=	scale[i];
	 BOX[ibox].rb		*=	scale[i];
	 BOX[ibox].drmax	*=	scale[i];
     
	 scale3	=	scale[i] * scale[i] * scale[i];
         scale6	=	scale3 * scale3;

         V	=	v[ibox];	 			// update potential energy
         VIR	=	vir[ibox];				// update virial
         if (V_LJ) {
            V.tot	-=	(V.lj6 + V.lj12);
            V.lj6	*=	1.0/ scale6;
            V.lj12	*=	1.0/ (scale6 * scale6);
            V.tot	+=	(V.lj6 + V.lj12);

            if (V_VIRIAL) {
	       VIR.tot	-=	(VIR.lj6 + VIR.lj12);
	       VIR.lj6	*=	1.0 / scale6;
               VIR.lj12	*=	1.0 / (scale6 * scale6);
	       VIR.tot	+=	(VIR.lj6 + VIR.lj12);
            }   
         }
         if (V_RPL) {	
	    V.tot	-=	V.rpl;
	    V.rpl	*=	1.0 / (scale6 * scale6);
            V.tot	+=	V.rpl;
         }
         if (V_LJLRC) {
            V.tot	-=	N[ibox] * Vtailco(BOX_old[i].rc, N[ibox]/BOX_old[i].vol);
            V.tot	+=	N[ibox] * Vtailco(BOX[ibox].rc, N[ibox]/BOX[ibox].vol);
         }
         v[ibox]	=	V;
         vir[ibox]	=	VIR;
      }
   }
   else {
#ifdef CELL_LIST
      New_CL();
#endif
      Vtotal();
   }

   dV	=	0;
   dV	+=	v[in].tot - v_old[0].tot;
   dV	+=	v[out].tot - v_old[1].tot;

   arg	=	dV/kT;
   arg	-=	(N[in]+1) * 3 * log(scale[0]);
   arg	-=	(N[out]+1) * 3 * log(scale[1]);

   if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {
      acc_gibbsvol	++;
   }
   else {
      rjc_gibbsvol	++;

      for (i=0; i<2; i++) {
         ibox	=	(i==0 ? in : out);

         BOX[ibox]	=	BOX_old[i];
         v[ibox]	=	v_old[i];
         vir[ibox]	=	vir_old[i];
         scale[i]	=	1.0/scale[i];
      }

      for (i=0; i<NPARTS; i++) {			// scale part. coord. back
         ibox		=	part[i].box;

         if (ibox == in || ibox == out) {
            part[i].p	=	V_Mult(scale[(ibox==in ? 0 : 1)], &(part[i].p));
#ifdef VERLET_LIST
            part[i].pv	=	V_Mult(scale[(ibox==in ? 0 : 1)], &(part[i].pv));	
#endif
         }
      }

#ifdef CELL_LIST
      if (SCALECUTOFF==0)
         New_CL();
#endif
   }  
   return;    
}
#endif /* TEST */


#ifdef TEST
void mcswap()				// Gibbs ensemble Swap move, F&S algorithm 19
{
   vstruct	v_old[2];
   wstruct	vir_old[2];

   long		i, n, in, out, ibox;
   molstruct	ghost;
   molstruct	mol_old;
   vector	p;
   double	dV, arg, Vghost;
   long		NCELLS;

   out	=	(int) (NBOX * ran1(seed));		// pick up two boxes
   if (NBOX==2) {					// even if N[out]==0, we don't reject it right away
      in	=	1 - out;			// 	because we can still calculate 
   }							//	chemical potential
   else if (NBOX > 2) {
      while ( (in = (int) (NBOX * ran1(seed))) == out );
   }

   for (i=0; i<2; i++) {				// Store information
      ibox		=	(i==0 ? in : out);
      v_old[i]		=	v[ibox];
      vir_old[i]	=	vir[ibox];
   }

   if (PBC==1) {					// new particle at a random position
      p.x	=	(ran1(seed)-0.5) * BOX[in].lbox;
      p.y	=	(ran1(seed)-0.5) * BOX[in].lbox;
      p.z	=	(ran1(seed)-0.5) * BOX[in].lbox;
   }
   ghost.p	=	p;				// create a test particle
   ghost.box	=	in;

#ifdef CELL_LIST
   ghost.icell	=	CL_Findcell(ghost.p, in, PBC);	// find a cell for this particle
#endif

   dV		=	0.0;
   dV		+=	VAddTestMol(&ghost, in);	// system energy update due to ghost
   N[in]	++;

   Vghost	=	dV;				// calculate single energy of ghost particle
   if (V_LJLRC) {					// Note!  N[in] has been updated
      Vghost	-=	(N[in]) * Vtailco(BOX[in].rc, (N[in])/BOX[in].vol);
      Vghost	+=	(N[in]-1) * Vtailco(BOX[in].rc, (N[in]-1)/BOX[in].vol);
      Vghost	+=	2 * Vtailco(BOX[in].rc, (N[in])/BOX[in].vol);
   }
   
   chp[in]	=	BOX[in].vol * exp(-Vghost/kT) / (N[in]);	// note N[in] has been updated
   cchp[in]	+=	chp[in];
   if (N[in]==NPARTS+1) {
      cchp[in]	+=	chp[in];
   }
   ichp[in]	++;					// sampling accumulator

   if (N[out]==0) {					// if box out is empty
      rjc_swap	++;

      N[in]	--;
      v[in]	=	v_old[0];			// box out hasn't been changed at this point
      vir[in]	=	vir_old[0];
      return;
   }

   ibox	=	-1;					// find a particle to be removed
   while (ibox != out) {
      n		=	(int) (NPARTS*ran1(seed));
      ibox	=	part[n].box;
   }
   dV	+=	VDeleteMol(n);
   N[out]	--;


   arg	=	dV/kT + log( BOX[out].vol * N[in] / (BOX[in].vol * (N[out]+1)) );	// formula is a little 

			// different from that given in the book because N[ibox] has been updated here

   if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {
      acc_swap	++;
#ifdef CELL_LIST
      CL_Update(n, part[n].icell, ghost.icell);
#endif
      part[n]	=	ghost;				// officially swap the particle
   }
   else {
      rjc_swap	++;
      N[in]	--;
      N[out]	++;
      for (i=0; i<2; i++) {
         ibox		=	(i==0 ? in : out);

         v[ibox]	=	v_old[i];
         vir[ibox]	=	vir_old[i];
      }
   }
   return;
}
#endif /* TEST */



/****************************************************************************************/
/*	Update_Eta (char *, int etaswitch)						*/
/*											*/
/*	Reweighting scheme for specified dynamics variable.				*/
/*	if etaswitch = 0, then initialize eta and p for no reweighting run		*/
/*	if etaswitch = 1, then initialize eta and p for reweighting run, and set up	*/
/*		sampling window								*/
/*	if etaswitch = 2, then update eta according to reweighting technique		*/
/****************************************************************************************/
/*
void Update_Eta(int etaswitch)
{
   int		i;
   double	dpQ[Qlbins], dp[NMAXbins];
   int		maxdpid, maxpid;
   double	maxdp, maxp;
   double	damping;
   double	total;
 
   if (etaswitch==0) {
      if (dynvar==2) {			//Ql as dynamic variable
         for (i=0; i<Qlbins; i++) {
            etaQ[i]	=	0;
         }
      }
      else if (dynvar==1) {
         for (i=0; i<NMAXbins; i++) {	//NMAX as dynamics variable
            eta[i]	=	0;
         }
      } 
   }
   else if (etaswitch==1) {		//update eta
      maxdp	=	0;
      maxp	=	0;         
 
      if (dynvar==2) {
	 total		=	0;
	 for (i=0; i<Qlbins; i++) {
	    total	+=	pQ[i];
	 }
         for (i=0; i<Qlbins; i++) {
	    pQ[i]	/=	total;

            ptQ[i]	=	1.0/Qlbins;
            dpQ[i]	=	fabs(pQ[i]/ptQ[i] - 1);
            if (dpQ[i] > maxdp) {
               maxdpid	=	i;
               maxdp	=	dpQ[i];
            }
            if (pQ[i] > maxp) {
               maxpid	=	i;
               maxp	=	pQ[i];
            }
         }
      }
      else if (dynvar==1) {
         for (i=0; i<NMAXbins; i++) {
            pt[i]	=	1.0/NMAXbins;
            dp[i]	=	fabs(p[i]/pt[i] - 1);
            if (dp[i] > maxdp) {
               maxdpid	=	i;
               maxdp	=	dp[i];
            }
            if (p[i] > maxp) {
               maxpid	=	i;
               maxp	=	p[i];
            }
         }
      }

      if (maxdp < CRIT) 	 trial=TRIALRUN-1;	//if already good enough, finish trial runs
 
      damping	=	Alpha;


      if (dynvar==2) {
         etaQ[maxpid]	-=	damping * log( pQ[maxpid]/ptQ[maxpid] );	//update the most sampled bin
         for (i = maxpid+1; i<Qlbins; i++) {					//update other bins		
            if (pQ[i] > ZERO) {
               etaQ[i]	-=	damping * log( pQ[i]/ptQ[i] );
            }
            else if( pQ[i-1]>ZERO) {				
               etaQ[i]	=	etaQ[i-1] + (etaQ[i-1] - etaQ[i-2]);
            }				
            else {
               etaQ[i]	=	etaQ[i-1];
            }
         }
         for (i = maxpid-1; i>=0; i--) {
            if (pQ[i] > ZERO) {
               etaQ[i]	-=	damping * log( pQ[i]/ptQ[i] );
            }
            else if ( pQ[i+1] > ZERO) {				
               etaQ[i]	=	etaQ[i+1] + (etaQ[i+1] - etaQ[i+2]);
            }	
            else {
               etaQ[i]	=	etaQ[i+1];
            }
         }
      }
      else if (dynvar==1) {
         eta[maxpid]	-=	damping * log( p[maxpid]/pt[maxpid] );	//update the most sampled bin
         for (i = maxpid+1; i<NMAXbins; i++) {					//update other bins		
            if (p[i] > ZERO) {
               eta[i]	-=	damping * log( p[i]/pt[i] );
            }
            else if( p[i-1]>ZERO) {				
               eta[i]	=	eta[i-1] + (eta[i-1] - eta[i-2]);
            }				
	    else {
	       eta[i]	=	eta[i-1];
            }
         }
         for (i = maxpid-1; i>=0; i--) {
            if (p[i] > ZERO) {
               eta[i]	-=	damping * log( p[i]/pt[i] );
            }
            else if ( p[i+1] > ZERO) {				
               eta[i]	=	eta[i+1] + (eta[i+1] - eta[i+2]);
	    }	
            else {
	       eta[i]	=	eta[i+1];
	    }
         }
      }
   }

   if (dynvar==2) {
      for (i=0; i<Qlbins; i++) {		//initialize p, both at the very beginning and after eta update
         pQ[i]	=	0;
      } 
   }
   else if (dynvar==1) {
      for (i=0; i<NMAXbins; i++) {
         p[i]	=	0;
      }
   }
   return;  
}
*/
                                                                                                                                                                                                                                                                                                                                                                                                           src/extractlammps.c                                                                                 0000600 0143352 0000144 00000041225 11565016324 014113  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	extractlammps.c
	author:		Peng Yi at MIT
	date:		May 18, 2011
	purpose:	read lammps dump file, extract the configuration of one chain
			from different frames and output one single vmd-readable file
	note:		
			require setup file

			2011/5/18 the box dimension output is not exactly correct, 
			the box size of the last snapshot is output in .car file.
*/

#define __MAIN_PROGRAM
#include "header.h"

#define CAP	5		// largest nuclei considered still the melt
#define p2threshold	0.4
#define DEBUG	0

#include "correlation.h"

long		timestep;
double		rshift2;	// shift of the biggest nucleus
long		nsitesp2;	// # of sites with p2 greater than a threshold
long		nmolsp2;
long		nmaxp2_1[10];	// nmax using p2 nucleus definition
long		nmaxp2_2[10], nmaxp2_3[10], nmaxp2_4[10];

char *Element(long t, char *s)
{
   double	mass = type[t].M;

   s[0]		=	0;
/*
   if (14==mass || 15==mass)	strcpy(s, "C");
   else if (1.01==mass)		strcpy(s, "H");
   else if (28.086==mass)	strcpy(s, "Si");
   else if (26.982==mass)	strcpy(s, "Al");
   else if (16==mass)		strcpy(s, "O");
*/
   strcpy(s, "C");
   return	s;
}

void shiftbox(long system, beadstruct *nucleus, long nsites) //copied from conf2car.c (07/13/09)
{
   molstruct	*moli;
   long		i, n, site;
   vector	rA, rO, rBA, rOA;

   // step 1: find one segment that belongs to the biggest nucleus
   
   rA	=	nucleus[0].moli->p[nucleus[0].site];
  
   // step 2: calc. the shift of com of nucleus to this segment
   
   V_Null(&rBA);
   V_Null(&rOA);

   for (n=0; n<nsites; n++) { 
      moli	=	nucleus[n].moli;
      site	=	nucleus[n].site;
      rBA	=	V_Subtr(moli->p+site, &rA);
      rBA	=	MapInBox2(&rBA, PBC, system);
      rOA	=	V_Add(&rOA, &rBA);
   }
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);

   // step 3: every segment shift and move to central box
  
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (site=0; site<moli->nsites; site++) {
         moli->p[site]	=	V_Subtr(moli->p+site, &rO);
         moli->p[site]	=	MapInBox2(moli->p+site, PBC, system);
      }
   }
} 

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   molstruct	*moli;
   molstruct	*moln;
   long		i, j, k, system, m, n, nuclid, nsites;
   long		previous_timestep=-1;
   long		nx, ny, nz, id, siteid, molid, type;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo;			// box upper boundary
   double	temp1, temp2, temp3;		// dummy variables

   char		infile[8][255], filein[80], filename[80];
   char		s[80], ff[80], par[80], dummy[255];
   FILE		*fin, *fhst, *fout, *fconf, *fpdb, *fdat;
   long		LENGTH, accum, confflag=0, carflag=0, pdbflag=0, polydisperse=0, drawmol;
   long		nframe, dnframe=1;		// analyze only every dnframe
   long		nfiles=1, ifile;		// number of input files 
   long		nshots=1, ishots=0;		// number of snapshots
   long		chainid=0;
   char		atomname;
   static long	init=1;	

   vector		con;	// center of nucleus
   static vector	rO;	// original position of nucleus

   // chain rotational angle distribution
   vector	chainface;
   double	orient;
   long		orientdist[180];

   // segment statistics variables
   long		previous, previous_id;
   long		nlength[MAXNMOLSITES], n_sep[MAXNMOLSITES];
   long		head, tail, seg_on, seg_id, nseg, nxtal, length;
   long		segment[MAXNMOLS][MAXNMOLSITES];	// segments identification on a chain
   long		nsegment[MAXNMOLS];			// # of segments on a chain
   long		seg_stat[MAXNMOLS][MAXNMOLSITES];	// xtal segment stat for chains
   long		loose_stat[MAXNMOLS][MAXNMOLSITES];	// loose segment stat
   long		nloop, nbridge, ntail, nxseg;		// # of total loops, etc
   double	lloop, lbridge, ltail, lxseg;		// average length of loops, etc
   long		nlloop[MAXNMOLSITES], nlbridge[MAXNMOLSITES];	// segment length distribution
   long		nltail[MAXNMOLSITES], nlxseg[MAXNMOLSITES];

   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];		// group beads in the same nucleus

   vector	com[MAXNMOLS], temp;		// average center of mass
   long		ncom[MAXNMOLS];

   double	imagex, imagey, imagez;		// variables for creating pbc images
   long		imagen;

   if (argc<2) {
      printf("extractlammps (c) 2011 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\textractlammps [-option] [nshots= 1] [dn= 2] chainid= 0 lammpsdumpfile\n\n");

      printf("Notes:\n");
      printf("\t* -option = \n");
      printf("\t* -poly: polydisperse system, has molid in dump file\n");
      printf("\t* -conf: configuration file output\n");
      printf("\t* -car: car file output\n");
      printf("\t* -pdb: pdb file output\n");
      printf("\t* nshots=: number of snapshots\n");
      printf("\t* dn=: only analyze every dn frames\n");
      printf("\t* chainid: id of the chain we extract, id starts from 0\n");
      printf("\t* nfiles=: number of input files if more than 1 (must be <=8)\n");
      printf("\t* \"=\" must immediately follow x or y or z or n or dn\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-poly"))	polydisperse	=	1;
      else if (samestr(par, "-conf"))	confflag	=	1;
      else if (samestr(par, "-car"))	carflag		=	1;
      else if (samestr(par, "-pdb"))	pdbflag		=	1;
      else if (samestr(par, "nshots="))	nshots		=	atol(argv[i+1]);
      else if (samestr(par, "dn=")) 	dnframe		=	atol(argv[i+1]);
      else if (samestr(par, "nfiles="))	nfiles		=	atol(argv[i+1]);
      else if (samestr(par, "chainid="))chainid		=	atol(argv[i+1]);
      
   }
   for (i=0; i<nfiles; i++) {
      strcpy(infile[nfiles-1-i], argv[argc-1-i]);	// get input filenames
   }

   // Open output files
 
   if (nfiles==1)	strcpy(filein, infile[0]);
   else			strcpy(filein, "multi");

   strcpy(filename, filein);
   strcat(filename, ".car");
   if (carflag && (fout=fopen(filename, "w"))==NULL )
      Exit("lammps2hst", "main", "open car file failed.");

   strcpy(filename, filein);
   strcat(filename, ".pdb");
   if (pdbflag && ((fpdb=fopen(filename, "w"))==NULL || (fdat=fopen("vmd.dat","w"))==NULL))
      Exit("lammps2hst", "main", "open pdb file failed.");

   if (!DEBUG) {
      strcpy(filename, filein);
      strcat(filename, ".out");
      freopen(filename, "w", stdout);	// redirect standard output stream to a file
   }

   printf("nshots = %d dn = %d chainid = %d infile[0] = %s\n", 
		nshots, dnframe, chainid, infile[0]);

   ////////////////////
   // Initialization //
   ////////////////////

   if (DEBUG)	printf("Initialization: start ... \n");

   moln	=	(molstruct *) calloc(nshots, sizeof(molstruct));
  
   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   InitUnits();				// initialize units for calculation
   InitForcefield();			// initialize Lennard-Jones potential mixing rule

   system	=	0;		// for now 2/14/08, only one system
   					// initialize sampling
					// initialize chain rotation distribution
				   	// initialize segment variables
 
   nframe	=	-1;

for (ifile=0; ifile<nfiles; ifile++){		// multiple input files

   fin	=	fopen(infile[ifile], "r");

   while (!feof(fin)) {

      if (DEBUG)	printf("Read one configuration: start ...\n");

      /* Read in one configuration from dump file */

      if (!fgets(dummy, sizeof(dummy), fin))	break;			// end of file
      fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &xlo, &xhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &ylo, &yhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &zlo, &zhi);	fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);

      BOX[system].lx	=	xhi-xlo;
      BOX[system].ly	=	yhi-ylo;
      BOX[system].lz	=	zhi-zlo;

      LENGTH		=	NSITES/NMOLS;	// monodisperse system for now (4/26/2008)
      accum		=	0;

      for (i=0; i<nsites; i++) {
         fscanf(fin, "%ld", &id);
         fscanf(fin, "%ld", &molid);		// Need to check the lammps.dump file format
						// because some early lammps.dump file 
						// does not have molid output
//         if (polydisperse)	fscanf(fin, "%ld", &molid);
         fscanf(fin, "%ld", &type);
         fscanf(fin, "%lf%lf%lf %lf%lf%lf %ld%ld%ld", &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
         fgets(dummy, sizeof(dummy), fin);

         if (polydisperse) {			// polydisperse
            molid	--;
            /*accum	=	0;
            for (j=0; j<molid; j++)
               accum	+=	(mol+j)->nsites;
            siteid	=	id  - accum;
            */
         }
         else {					// monodisperse
            id	--; 				// -1 because lammps index starts from 1
   	    molid	=	(long) (id/LENGTH);
            siteid	=	id % LENGTH;
         }
         mol[molid].box		=	system;		// for now, only one system
         mol[molid].nsites	=	LENGTH;		// for now, Jan/10/2010
         mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx);
         mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly);
         mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);
         mol[molid].type[siteid]=	type - 1;	// -1 because lammps index starts from 1
      }

      for (system=0; system<NSYSTEMS; system++) {
         NMols[system]	=	0;
         NSites[system]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;			// total # of mols in certain system
            NSites[system]	+=	moli->nsites;	// total # of sites in certain system
         }
      }

      for (moli=mol; moli<mol+NMOLS; moli++) { 		
         for (i=0; i<moli->nsites; i++)  {
            moli->flags[i]	=	1;		// activate all the sites on this processor
            moli->parent[i]	=	i-1;		// initialize parent site
         }
         moli->flip		=	0;		// flip to the original direction
         moli->origin		=	CenterofMass(moli);
      }

      ///////////////////////////////////////////////////
      /* Skip repeated frames in different input files */
      ///////////////////////////////////////////////////
      if (ifile>=1 && timestep==previous_timestep) {
         continue;
      }
      previous_timestep	=	timestep;	// check repeat frames in different infiles

      ///////////////////////////
      /* Skip dnframe-1 frames */
      ///////////////////////////

      nframe	++;

      if (mod(nframe, dnframe))	continue;	// analyze every dnframe frames
      
      //////////////////////////
      // Extract the molecule //
      //////////////////////////

      if (DEBUG)	printf("Extract the molecule ...\n");

      moln[ishots]	=	mol[chainid];
      MolInBox2(moln+ishots);
      ishots	++;

      if (ishots >= nshots)	break;

      ///////////////////////////////////////////////////////////////////////////
      // Start: Calculate the average position of center of mass of each chain //
      ///////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////
      /* Convert coordinates and box size from SI to system unit */
      /////////////////////////////////////////////////////////////
      CoorSI2System();
      //_________________________________________________________//

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
      } 

      /////////////////////
      // build cell list //
      /////////////////////

      //////////////////////
      // Perform analysis //
      //////////////////////

      /////////////////////////////
      // Output analysis results //
      /////////////////////////////

      CoorSystem2SI();		// convert coordinates and box size back to SI units

   } 

   fclose(fin);		// close current input file
  }			// multiple input files

   ////////////////////////////////////////////
   // OUTPUT .car file for VMD visualization //
   ////////////////////////////////////////////
      
   if (carflag) {
      fprintf(fout, "!BIOSYM archive 3\n");
      fprintf(fout, "PBC=ON\n");
      fprintf(fout, "!TIMESTEP %d\n", timestep);
      fprintf(fout, "!DATE %s", asctime(localtime(&t)));
      fprintf(fout, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

      n	=	0;

      for (moli=moln; moli<moln+nshots; moli++) {
         if (system==moli->box) {

            //MolInBox2(moli);
            for (i=0; i<moli->nsites; i++) {

               if (moli->nuclid[i]>0)	// crystal like particle
                  sprintf(s, "N%d", n++);	// N: blue color in VMD
               else
                  sprintf(s, "O%d", n++);	// O: red color in VMD

               fprintf(fout, "%-5.5s ", s);
               sprintf(s, "M%d", moli-moln);
               fprintf(fout, "%14.8g %14.8g %14.8g ", moli->p[i].x, moli->p[i].y, moli->p[i].z);
               strcpy(ff, "O");
               fprintf(fout, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-moln, Element(moli->type[i], s));
            } 
         }   
      }
      fprintf(fout, "end\nend\n");
   }
   fflush(fout);
   //___________OUTPUT .car file_____________//
/*
   ///////////////////////////////////////////
   // OUTPUT .pdb file for further analysis //
   ///////////////////////////////////////////

   if (pdbflag) {
      fprintf(fpdb, "HEADER: file created from %s on %s", argv[1], asctime(localtime(&t)));
      fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
        	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

      m		=	0;	// molecule sequence number
      n		=	0;	// atom sequence number
#define SIZECAP	5
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if (system==moli->box) {
            //MolInBox2(moli);

            drawmol	=	0;
            for (i=0; i<moli->nsites; i++) {
               if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {
                  drawmol	=	1;		// participate the biggest nucleus
                  break;
               }
	    }
//temp=CenterofMass(moli);
//if (temp.z < -0.25 * BOX[system].lz) {
            m	++; 
            for (i=0; i<moli->nsites; i++) {
               if (drawmol) {
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {	// nuclid index starts from 1
                     atomname	=	'N';	// N: blue color in VMD
	             fprintf(fdat, " 10");
                  }
                  else {
                     atomname	=	'O';	// O: red color in VMD
                     fprintf(fdat, " 0");
                  }
               }
               else {
	          atomname	=	'C';		// C: cyan color in VMD
	          fprintf(fdat," -1");
	       }

               n	++;
	       fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
               fprintf(fpdb, "%5d ", n);	// atom number
               fprintf(fpdb, " %c  ", atomname);	// atom name
               fprintf(fpdb, " ");		// alternate location indiator
  	       fprintf(fpdb, " C8");		// residue name
	       fprintf(fpdb, " ");		// column 21
               fprintf(fpdb, " ");		// chain identifier, column 22
	       fprintf(fpdb, "%4d", m);	// residue sequence number, 23-26
	       fprintf(fpdb, " ");		// code for insertion of residues, 27
               fprintf(fpdb, "   ");		// column 28-30
               fprintf(fpdb, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z);
               fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
               fprintf(fpdb, "%5.5s", "");
               fprintf(fpdb, "\n"); 

               if (imagen) {			// for image box
                  n	++;
                  fprintf(fpdb, "ATOM  ");
                  fprintf(fpdb, "%5d ", n);
                  fprintf(fpdb, " %c  ", atomname);
                  fprintf(fpdb, " ");
                  fprintf(fpdb, " C8");
	          fprintf(fpdb, " ");
                  fprintf(fpdb, " ");
	          fprintf(fpdb, "%4d", m);
                  fprintf(fpdb, " ");
                  fprintf(fpdb, "   ");
                  fprintf(fpdb, "%8.3f%8.3f%8.3f", 
			moli->p[i].x + imagex, moli->p[i].y+imagey, moli->p[i].z+imagez);
                  fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                  fprintf(fpdb, "%5.5s", "");
                  fprintf(fpdb, "\n"); 
               }
            } 
//}
         }   
      }
      fprintf(fpdb, "END\n");
      fflush(fpdb);
   }	//pdbflag
   //____________OUTPUT .pdb file___________//
*/

   if (DEBUG)	printf("Closing output files ...\n");

   if (carflag)		fclose(fout);
   if (pdbflag)	     {	fclose(fpdb); fclose(fdat);}

   fflush(stdout);
   fclose(stdout);

   return	0;
}
                                                                                                                                                                                                                                                                                                                                                                           src/forcefield.c                                                                                    0000600 0143352 0000144 00000134550 11302536372 013334  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    forcefield.c
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Calculation of the potential of a system

    Modified:	Sept 20, 2007
		Long range correction no longer calculated for each particle,
		but calculated for the whole system at one time.
		
*/

#define __FORCEFIELD_MODULE
#include "forcefield.h"

/* Some basic operations to vstruct and wstruct */

void vstructNull(vstruct *V)
{
   V->lj	=	0.0;
   V->ljcorr	=	0.0;
   V->hs	=	0.0;
   V->stretch	=	0.0;
   V->bending	=	0.0;
   V->torsion	=	0.0;
   V->corr	=	0.0;
   V->bonded	=	0.0;
   V->nonbonded	=	0.0;
   V->tot	=	0.0;
}

void wstructNull(wstruct *VIR)
{
   VIR->lj	=	0.0;
   VIR->stretch	=	0.0;
   VIR->torsion	=	0.0;
   VIR->tot	=	0.0;
}

void vstructNegate(vstruct *V)
{
   V->lj	*=	-1.0;
   V->ljcorr	*=	-1.0;
   V->hs	*=	-1.0;
   V->stretch	*=	-1.0;
   V->bending	*=	-1.0;
   V->torsion	*=	-1.0;
   V->corr	*=	-1.0;
   V->bonded	*=	-1.0;
   V->nonbonded	*=	-1.0;
   V->tot	*=	-1.0;
}

void wstructNegate(wstruct *VIR)
{
   VIR->lj	*=	-1.0;
   VIR->stretch *=	-1.0;
   VIR->torsion	*=	-1.0;
   VIR->tot	*=	-1.0;
}

vstruct vstructSum(vstruct *V1, vstruct *V2)
{
   vstruct	V;
  
   V.lj		=	V1->lj + V2->lj;
   V.ljcorr	=	V1->ljcorr + V2->ljcorr;
   V.hs		=	V1->hs + V2->hs;
   V.stretch	=	V1->stretch + V2->stretch;
   V.bending	=	V1->bending + V2->bending;
   V.torsion	=	V1->torsion + V2->torsion;

   V.corr	=	V1->corr + V2->corr;
   V.bonded	=	V1->bonded + V2->bonded;
   V.nonbonded	=	V1->nonbonded + V2->nonbonded;
   V.tot	=	V1->tot + V2->tot;
   return	V;
}

wstruct wstructSum(wstruct *W1, wstruct *W2)
{
   wstruct	W;

   W.lj		=	W1->lj + W2->lj;
   W.stretch	=	W1->stretch + W2->stretch;	
   W.torsion	=	W1->torsion + W2->torsion;

   W.tot	=	W1->tot + W2->tot;
   return	W;
}

void Printvstruct(vstruct *V)
{
   printf("v.lj\t=\t%f\n", V->lj);
   printf("v.ljcorr\t=\t%f\n", V->ljcorr);
   printf("v.hs\t=\t%f\n", V->hs);
   printf("v.stretch\t=\t%f\n", V->stretch);
   printf("v.bending\t=\t%f\n", V->bending);
   printf("v.torsion\t=\t%f\n", V->torsion);
   printf("v.corr\t=\t%f\n", V->corr);
   printf("v.bonded\t=\t%f\n", V->bonded);
   printf("v.nonbonded\t=\t%f\n", V->nonbonded);
   printf("v.tot\t=\t%f\n", V->tot);
}


/* Lennard Jones combining rules */

/* SIGMA_eff 	= 0.5 * (SIGMA_1 + SIGMA_2) */
/* EPSILON_eff 	= sqrt( EPSILON_1 * EPSILON_2 ) */

void CalcMixSigma(long typei, long typej)
{
   if (typei == typej)
      type[typei].mix[typei].SIGMA	=	type[typei].SIGMA;
   else {
      type[typei].mix[typej].SIGMA	
		=	0.5 * (type[typei].SIGMA + type[typej].SIGMA);
      type[typej].mix[typei].SIGMA
		=	type[typei].mix[typej].SIGMA;
   }
}		

void CalcMixEpsilon(long typei, long typej)
{
   if (typei == typej)
      type[typei].mix[typej].EPSILON	=	type[typei].EPSILON;
   else {
      type[typei].mix[typej].EPSILON
		=	sqrt( type[typei].EPSILON * type[typej].EPSILON );
      type[typej].mix[typei].EPSILON
		=	type[typei].mix[typej].EPSILON;
   }
}

void InitForcefield()
{
   long		i, j;
   
   for (i=0; i<NTYPES; i++)
      for (j=i; j<NTYPES; j++) {	// symmetric combining rule
         CalcMixSigma(i, j);
         CalcMixEpsilon(i, j);
      }
}

// OPLS() and OPLS2() are two equivalent forms, with same k1, k2, and k3.

double OPLS(double phi, double k1, double k2, double k3)	// torsional energy, OPLS model
{
   return	0.5 * ( k1*(1-cos(phi)) + k2*(1-cos(2*phi)) + k3*(1-cos(3*phi)) );
}

double OPLS2(double cosphi, double k1, double k2, double k3)
{
   return	0.5* ((k1+2*k2+k3) + (3*k3-k1)*cosphi - 2*k2*cosphi*cosphi -4*k3*cosphi*cosphi*cosphi);
}

double POLY(double cosphi, double k0, double k1, double k2, double k3, double k4, double k5)
{
   double	cosphi2	= cosphi * cosphi;
   double	cosphi4	= cosphi2 * cosphi2;

   return	k0 + k1 * cosphi + k2 * cosphi2 + k3 * cosphi * cosphi2 + k4 * cosphi4 
		+ k5 * cosphi * cosphi4;
}

/* Calculate different types of interactions */
/* V*Site() and V*Mol() do NOT initialize the virial */
/* in their parameter list, so we have to make sure  */
/* they do the right thing to the right virial.      */

double VHSSite(molstruct *molm, long site)
{
   static long		ibox, n;
   static molstruct	*moln;
   static double	r2;
   static vector	pm;
   static typestruct	*typema;
#ifdef CELL_LIST
   static long		i, j;
   static cellstruct	*cellm, *celli;
#endif

   if (!V_HS || 0==(molm->flags[site]) ) 
      return	0.0;
   
   ibox		=	molm->box;			// determine which box
   pm		=	molm->p[site];
   typema	=	type + molm->type[site];

#ifdef CELL_LIST
   cellm	=	molm->cell[site];
   for (j=0; j<cellm->nneigh; j++) {
      
      celli	=	cellm->neigh[j];

      for (i=0; i<celli->nsites; i++) {
         if (moln=celli->mol[i]) {

            n	=	celli->molsite[i];
#else
   for (moln=mol; moln < mol + NMOLS; moln++) {				// search through all mols

      if (moln->box	==	ibox) {					// check if two mols are in the same box

         for (n=0; n<moln->nsites; n++) {				// search through all sites
#endif /* CELL_LIST */
            
            if ( (moln->flags[n]>0) && (molm!=moln || n!=site) ) {	// interaction with itself is excluded

               r2	=	DistSQ(pm, moln->p[n], ibox);
               
               if (r2 < (typema->mix+moln->type[n])->SIGMA * (typema->mix+moln->type[n])->SIGMA)
		  return	1.0e4;
               else
		  return	0.0;
	    }
         }   
      }
   }
}


double VHSMol(molstruct *moli)
{
   long		i;
   double	vhs=0.0;

   if (moli->box <0 || !V_HS )
      return	0.0;

   for (i=0; i<moli->nsites; i++)
      if ( vhs=VHSSite(moli, i) > 1.0)
	 return	vhs;	
   
   return	vhs;
}


void CalcVHS()				// calculate total hard sphere interaction
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++)
      v[ibox].hs	=	0.0;

   for (moli=mol; moli<mol+NMOLS; moli++) {
      if ( (ibox=moli->box) >=0 && v[ibox].hs < 1.0) 
         v[ibox].hs	+=	VHSMol(moli);
   }
}


double VStretchSite(molstruct *molm, long site, double *w)	// V = 0.5 * k * (l-l0)^2
{
   long			i, k=0, system=molm->box;
   double		l0, cosa, v=0.0, f;
   vector		dr[2], *p0, *p1;
   typestruct		*t;

   if ( (!V_STRETCH) || (molm->flags[site]==0) || (site<0) || (site>=molm->nsites) )
      return	0.0;
   p0		=	molm->p+site;
   i		=	site;
   while ( (k<1) && ((i=molm->parent[i])>=0) && (molm->flags[i]>0) ) {
      p1	=	molm->p+i;
      dr[k]	=	V_Subtr(p0, p1);
      p0	=	p1;
      k	++;
   }
   if ( !k )
      return	0.0;
   t		=	type + molm->type[site];
   l0		=	sqrt( V_Dot(dr, dr) );
   f		=	l0 - t->LSTRETCH;
   v		=	0.5 * f * f * t->KSTRETCH;
   if (V_VIRIAL)
      *w	+=	l0 * f * t->KSTRETCH;
   return	v;
}


double VStretchMol(molstruct *molm, double *w)
{
   long		i;
   double	vstretch = 0.0;

   if (!V_STRETCH || molm->box<0) 
      return	0.0;
   for (i=0; i<molm->nsites; i++) 
      vstretch	+=	VStretchSite(molm, i, w);
   return 	vstretch;
}


void CalcVStretch()
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++) {
      v[ibox].stretch		=	0.0;
      vir[ibox].stretch		=	0.0;
   }
   if (!V_STRETCH) 	
      return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox=moli->box) >=0 ) 
         v[ibox].stretch	+=	VStretchMol(moli, &(vir[ibox].stretch));
}


double VLJSite(molstruct *molm, long site, double *w)	// calculate the LJ potential energy of one site
{
   static long		ibox, n, flags_m;
   static vector	dp;
   static double	vlj, wlj;			// lj potential, lj virial
   static molstruct	*moln;
   static double	Sigma, Epsilon, r2, r6i, r12i, rc2;
   static vector	pm;
   static double	ljcut;
   typestruct		*typem = type + molm->type[site];
   mixstruct		*typemix;
#ifdef CELL_LIST
   static long		i, j;
   static cellstruct	*cellm, *celli;
#endif
   int	nn=0;

   if (!V_LJ || (flags_m=molm->flags[site]) == 0 || (typem->EPSILON < ZERO) )
      return	0.0;

   vlj	=	0.0;					// MUST, because vlj is static variable
   wlj	=	0.0;
   
   ibox		=	molm->box;			// determine which box
   rc2		=	BOX[ibox].rc * BOX[ibox].rc;
   pm		=	molm->p[site];
#ifdef CELL_LIST
   cellm	=	molm->cell[site];

   for (j=0; j<cellm->nneigh; j++) {
      
      celli	=	cellm->neigh[j];

      for (i=0; i<celli->nsites; i++) {
         if (moln=celli->mol[i]) {

            n	=	celli->molsite[i];
#else
   for (moln=mol; moln < mol + NMOLS; moln++) {		// search through all mols

      if (moln->box	==	ibox) {			// check if two mols are in the same box

         for (n=0; n<moln->nsites; n++) {		// search through all sites
#endif /* CELL_LIST */
            if ( (moln->flags[n]>0) && (moln==molm ? abs(n-site)>=DLJ : 1) ) {	// LJ interaction condition
									// interaction with itself is excluded
									// automatically
               r2	=	DistSQ(pm, moln->p[n], ibox);	
               if (r2 < rc2) {
		  typemix	=	typem->mix + moln->type[n];
                  Sigma		=	typemix->SIGMA;
	          Epsilon	=	typemix->EPSILON;

         	  r6i 		= 	Sigma * Sigma/r2;
	          r6i 		= 	r6i * r6i * r6i;
		  r12i		=	r6i * r6i;
/*
		  if ( moln==molm && DLJ==abs(n-site) ) {	// LJ 1-4 pair scaling
         	     vlj 	+= 	2.0 * Epsilon * (r12i - r6i);
                     if (V_VIRIAL)
                        wlj	+=	-24.0 * Epsilon * (r12i - 0.5 * r6i);
  		  }
		  else { 
*/         	     vlj 	+= 	4.0 * Epsilon * (r12i - r6i);
                     if (V_VIRIAL)
                        wlj	+=	-48.0 * Epsilon * (r12i - 0.5 * r6i);
  //    		  }

/*
		  if (V_LJSHIFT) {                
      		     ljcut	=	Sigma / BOX[ibox].rc;		// sigma/rc
		     ljcut	=	ljcut * ljcut * ljcut;		// (sigma/rc)^3
		     ljcut	=	ljcut * ljcut;			// (sigma/rc)^6
		     ljcut	=	4.0 * Epsilon * (ljcut * ljcut - ljcut);
//		     if ( moln==molm && DLJ==abs(n-site) ) 	// LJ 1-4 pair scaling
//		        vlj	-=	0.5 * ljcut;
//		     else
			vlj	-=	ljcut;
                  }
*/
               }
	    }
         }   
      }
   }
   if (V_VIRIAL)
      *w	+=	wlj;
   return vlj;
}


double VLJSiteinner(molstruct *molm, long site, double *w)	// same as VLJSite() but w/ a smaller cutoff
{
   static long		ibox, n, flags_m;
   static vector	dp, pm;
   static molstruct	*moln;
   static double	Sigma, Epsilon, r2, r6i, r12i, rc2, ljcut, vlj, wlj;
   typestruct		*typem = type + molm->type[site];
   mixstruct		*typemix;
#ifdef CELL_LIST
   static long		i, j;
   static cellstruct	*cellm, *celli;
#endif
   double		rc2low=Rclow*Rclow;			// smaller cutoff

   if (!V_LJ || (flags_m=molm->flags[site]) == 0 || (typem->EPSILON < ZERO) )
      return	0.0;

   vlj	=	0.0;					// MUST, because vlj is static variable
   wlj	=	0.0;
   
   ibox		=	molm->box;			// determine which box
   rc2		=	BOX[ibox].rc * BOX[ibox].rc;
   pm		=	molm->p[site];

#ifdef CELL_LIST
   cellm	=	molm->cell[site];

   for (j=0; j<cellm->nneigh; j++) {
      
      celli	=	cellm->neigh[j];

      for (i=0; i<celli->nsites; i++) {
         if (moln=celli->mol[i]) {

            n	=	celli->molsite[i];
#else
   for (moln=mol; moln < mol + NMOLS; moln++) {		// search through all mols

      if (moln->box	==	ibox) {			// check if two mols are in the same box

         for (n=0; n<moln->nsites; n++) {		// search through all sites
#endif /* CELL_LIST */

            if ( (moln->flags[n]>0) && (moln==molm ? abs(n-site)>=DLJ : 1) ) {	// LJ interaction condition
									// interaction with itself is excluded
									// automatically
               r2	=	DistSQ(pm, moln->p[n], ibox);	

               if (r2 < rc2low) {
		  typemix	=	typem->mix + moln->type[n];
                  Sigma		=	typemix->SIGMA;
	          Epsilon	=	typemix->EPSILON;

         	  r6i 		= 	Sigma * Sigma/r2;
	          r6i 		= 	r6i * r6i * r6i;
		  r12i		=	r6i * r6i;

                  vlj 	+= 	4.0 * Epsilon * (r12i - r6i);
                  if (V_VIRIAL)
                     wlj	+=	-48.0 * Epsilon * (r12i - 0.5 * r6i);
               }
	    }
         }   
      }
   }
   if (V_VIRIAL)
      *w	+=	wlj;
   return vlj;
}


double VLJSiteouter(molstruct *molm, long site, double *w)	// same as VLJSite() but w/ additional inner cutoff
{
   static long		ibox, n, flags_m;
   static vector	dp, pm;
   static molstruct	*moln;
   static double	Sigma, Epsilon, r2, r6i, r12i, rc2, ljcut, vlj, wlj;
   typestruct		*typem = type + molm->type[site];
   mixstruct		*typemix;
#ifdef CELL_LIST
   static long		i, j;
   static cellstruct	*cellm, *celli;
#endif
   double		rc2low=Rclow*Rclow;			// smaller cutoff

   if (!V_LJ || (flags_m=molm->flags[site]) == 0 || (typem->EPSILON < ZERO) )
      return	0.0;

   vlj	=	0.0;					// MUST, because vlj is static variable
   wlj	=	0.0;
   
   ibox		=	molm->box;			// determine which box
   rc2		=	BOX[ibox].rc * BOX[ibox].rc;
   rc2low	=	Rclow * Rclow;			// smaller cutoff
   pm		=	molm->p[site];

#ifdef CELL_LIST
   cellm	=	molm->cell[site];

   for (j=0; j<cellm->nneigh; j++) {
      
      celli	=	cellm->neigh[j];

      for (i=0; i<celli->nsites; i++) {
         if (moln=celli->mol[i]) {

            n	=	celli->molsite[i];
#else
   for (moln=mol; moln < mol + NMOLS; moln++) {		// search through all mols

      if (moln->box	==	ibox) {			// check if two mols are in the same box

         for (n=0; n<moln->nsites; n++) {		// search through all sites
#endif /* CELL_LIST */

            if ( (moln->flags[n]>0) && (moln==molm ? abs(n-site)>=DLJ : 1) ) {	// LJ interaction condition
									// interaction with itself is excluded
									// automatically
               r2	=	DistSQ(pm, moln->p[n], ibox);	

               if (rc2low <= r2 && r2 < rc2) {
		  typemix	=	typem->mix + moln->type[n];
                  Sigma		=	typemix->SIGMA;
	          Epsilon	=	typemix->EPSILON;

         	  r6i 		= 	Sigma * Sigma/r2;
	          r6i 		= 	r6i * r6i * r6i;
		  r12i		=	r6i * r6i;

                  vlj 	+= 	4.0 * Epsilon * (r12i - r6i);
                  if (V_VIRIAL)
                     wlj	+=	-48.0 * Epsilon * (r12i - 0.5 * r6i);

		  if (V_LJSHIFT) {                
      		     ljcut	=	Sigma / BOX[ibox].rc;		// sigma/rc
		     ljcut	=	ljcut * ljcut * ljcut;		// (sigma/rc)^3
		     ljcut	=	ljcut * ljcut;			// (sigma/rc)^6
		     ljcut	=	4.0 * Epsilon * (ljcut * ljcut - ljcut);
	             vlj	-=	ljcut;
                  }
               }
	    }
         }   
      }
   }
   if (V_VIRIAL)
      *w	+=	wlj;
   return vlj;
}


double VLJoutSite(molstruct *molm, long site, double *w)	// calc. LJ potential energy b/w one site
							// and all other sites NOT on the same chain
{
   static long		ibox, n, flags_m;
   static vector	dp;
   static double	vljout, wljout;			// lj potential, lj virial
   static molstruct	*moln;
   static double	Sigma, Epsilon, r2, r6i, r12i, rc2;
   static vector	pm;
   static double	ljcut;
   typestruct		*typem = type + molm->type[site];
   mixstruct		*typemix;
#ifdef CELL_LIST
   static long		i, j;
   static cellstruct	*cellm, *celli;
#endif

   if (!V_LJ || (flags_m=molm->flags[site]) == 0 || (typem->EPSILON < ZERO) )
      return	0.0;

   vljout	=	0.0;					// MUST, because vlj is static variable
   wljout	=	0.0;
   
   ibox		=	molm->box;			// determine which box this site locates
   rc2		=	BOX[ibox].rc * BOX[ibox].rc;	// LJ interaction cutoff
   pm		=	molm->p[site];			// store the position of this site

#ifdef CELL_LIST
   cellm	=	molm->cell[site];

   for (j=0; j<cellm->nneigh; j++) {
      
      celli	=	cellm->neigh[j];

      for (i=0; i<celli->nsites; i++) {
         if (moln=celli->mol[i]) {

            n	=	celli->molsite[i];
#else
   for (moln=mol; moln < mol + NMOLS; moln++) {		// search through all mols

      if (moln->box	==	ibox) {			// check if two mols are in the same box

         for (n=0; n<moln->nsites; n++) {		// search through all sites
#endif /* CELL_LIST */

            if ( (moln->flags[n]>0) && (moln!=molm) ) {		// LJ interaction condition, not on same chain
								// interaction with itself is excluded
								// automatically
               r2	=	DistSQ(pm, moln->p[n], ibox);	

               if (r2 < rc2) {
		  typemix	=	typem->mix + moln->type[n];
                  Sigma		=	typemix->SIGMA;
	          Epsilon	=	typemix->EPSILON;

         	  r6i 		= 	Sigma * Sigma/r2;
	          r6i 		= 	r6i * r6i * r6i;
		  r12i		=	r6i * r6i;
	          vljout	+= 	4.0 * Epsilon * (r12i - r6i);
                  if (V_VIRIAL)
                     wljout	+=	-48.0 * Epsilon * (r12i - 0.5 * r6i);

		  if (V_LJSHIFT) {                
      		     ljcut	=	Sigma / BOX[ibox].rc;		// sigma/rc
		     ljcut	=	ljcut * ljcut * ljcut;		// (sigma/rc)^3
		     ljcut	=	ljcut * ljcut;			// (sigma/rc)^6
		     ljcut	=	4.0 * Epsilon * (ljcut * ljcut - ljcut);
		     vljout	-=	ljcut;
                  }
               }
	    }
         }   
      }
   }
   if (V_VIRIAL)
      *w	+=	wljout;

   return vljout;
}


double VLJinSite(molstruct *molm, long site, double *w)		// LJ energy b/w one site
			// and other sites on the same chain, no need of using cell list
{
   static long		ibox, m, flags_m;
   static vector	dp, pm;
   static double	vljin, wljin;			// lj potential, lj virial
   static molstruct	*moln;
   static double	Sigma, Epsilon, r2, r6i, r12i, rc2, ljcut;
   typestruct		*typem = type + molm->type[site];
   mixstruct		*typemix;

   if (!V_LJ || (flags_m=molm->flags[site]) == 0 || (typem->EPSILON < ZERO) )
      return	0.0;

   vljin	=	0.0;				// MUST, because vlj is static variable
   wljin	=	0.0;
   
   ibox		=	molm->box;			// determine which box this site locates
   rc2		=	BOX[ibox].rc * BOX[ibox].rc;	// LJ interaction cutoff
   pm		=	molm->p[site];			// store the position of this site

   for (m=0; m<molm->nsites; m++) {			// search thru sites on the same chain
      if ( abs(m-site)>=DLJ && molm->flags[m]>0 ) {
         r2	=	DistSQ(pm, molm->p[m], ibox);

         if (r2 < rc2) {
            typemix	=	typem->mix + molm->type[m];
            Sigma	=	typemix->SIGMA;
            Epsilon	=	typemix->EPSILON;

            r6i 	= 	Sigma * Sigma/r2;
	    r6i 	= 	r6i * r6i * r6i;
	    r12i	=	r6i * r6i;
	    vljin	+= 	4.0 * Epsilon * (r12i - r6i);

            if (V_VIRIAL)
               wljin	+=	-48.0 * Epsilon * (r12i - 0.5 * r6i);

	    if (V_LJSHIFT) {                
      	       ljcut	=	Sigma / BOX[ibox].rc;		// sigma/rc
	       ljcut	=	ljcut * ljcut * ljcut;		// (sigma/rc)^3
	       ljcut	=	ljcut * ljcut;			// (sigma/rc)^6
	       ljcut	=	4.0 * Epsilon * (ljcut * ljcut - ljcut);
	       vljin	-=	ljcut;
            }
         }
      }
   }
   if (V_VIRIAL)
      *w	+=	wljin;

   return	vljin;
}


double VLJMol(molstruct *molm, double *w)
{
   long		i;
   double	vlj = 0.0;

   if (!V_LJ || molm->box<0 )
      return	0.0;
   for (i=0; i<molm->nsites; i++)
      vlj	+=	VLJSite(molm, i, w);

   return	vlj; 
}


double VLJoutMol(molstruct *molm, double *w)
{
   long		i;
   double	vljout = 0.0;

   if (!V_LJ || molm->box<0)
      return	0.0;

   for (i=0; i<molm->nsites; i++)
      vljout	+=	VLJoutSite(molm, i, w);

   return	vljout;
}


double VLJinMol(molstruct *molm, double *w)
{
   long		i;
   double	vljin = 0.0;

   if (!V_LJ || molm->box<0)
      return	0.0;

   for (i=0; i<molm->nsites; i++)
      vljin	+=	VLJinSite(molm, i, w);

   return	vljin;
}


void CalcVLJ()
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++) {
      v[ibox].lj	=	0.0;
      vir[ibox].lj	=	0.0;
   }
   if (!V_LJ)	return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox = moli->box) >=0) 
         v[ibox].lj	+=	0.5 * VLJMol(moli, &(vir[ibox].lj));	// fix double-counting

   if (V_VIRIAL) 
      for (ibox=0; ibox<NBOX; ibox++)
         vir[ibox].lj	*=	0.5;			// fix double-counting of virial
   
   return;
}


void CalcVLJout()			// LJ energy b/w sites on the different chains ONLY
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NSYSTEMS; ibox++) {
      v[ibox].ljout	=	0.0;
      vir[ibox].ljout	=	0.0;
   }
   if (!V_LJ)	return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox=moli->box) >= 0)
         v[ibox].ljout	+=	0.5 * VLJoutMol(moli, &(vir[ibox].ljout));	// fix double-counting

   if (V_VIRIAL)
      for (ibox=0; ibox<NSYSTEMS; ibox++)
         vir[ibox].ljout	*=	0.5;		// fix double-counting

   return;
}


void CalcVLJin()			// LJ energy b/w sites on the same chain ONLY
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NSYSTEMS; ibox++) {
      v[ibox].ljin	=	0.0;
      vir[ibox].ljin	=	0.0;
   }
   if (!V_LJ)	return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox=moli->box) >= 0)
         v[ibox].ljin	+=	0.5 * VLJinMol(moli, &(vir[ibox].ljin));	// fix double-counting

   if (V_VIRIAL)
      for (ibox=0; ibox<NSYSTEMS; ibox++)
         vir[ibox].ljin	*=	0.5;		// fix double-counting

   return;
}


double VBendingSite(molstruct *molm, long site)		// Only calculate the bending
{							// energy on its parent side
   long		i, k=0;
   double	l0, l1, cosa, v=0.0, f;
   vector	dr[2], *p0, *p1;
   typestruct	*t;

   if ( (!V_BENDING) || (0==molm->flags[site]) || (site<0) || (site>=molm->nsites) )
      return	0.0;
   p0		=	molm->p+site;
   i		=	site;
   while ( (k<2) && ((i=molm->parent[i])>=0) && (molm->flags[i]>0)) {
      p1	=	molm->p+i;
      dr[k++]	=	V_Subtr(p0, p1);
      p0	=	p1;
   }
   if (k<2)
      return	0.0;
   t		=	type + molm->type[site];
   l0		=	sqrt(V_Dot(dr, dr));
   l1		=	sqrt(V_Dot(dr+1, dr+1));
   cosa		=	V_Dot(dr, dr+1)/(l0*l1);
   f		=	acos(cosa) - t->THETA;		// OPLS model
//   f		=	cosa - cos(t->THETA);	
   v		=	0.5 * t->KBENDING * f * f;
   return	v;
}


double VBendingMol(molstruct *molm)
{
   long		i;
   double	v_bending = 0.0;

   if (!V_BENDING || (molm->box<0))
      return	0.0;
   for (i=0; i<molm->nsites; i++)
      v_bending	+=	VBendingSite(molm, i);
   return	v_bending;
}


void CalcVBending()
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++) 
      v[ibox].bending	=	0.0;

   if (!V_BENDING)	return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox = moli->box) >=0 )
         v[ibox].bending	+=	VBendingMol(moli);	// no double-counting
}


double VTorsionSite(molstruct *molm, long site)		// no virial contribution
{							// only toward parents side !!
   long		i, k=0;
   double	cosb, b, l0, l1;
   vector	dr[3], n0, n1, *p0, *p1;
   typestruct	*t;
double temp1, temp2;

   if ( (!V_TORSION) || (0==molm->flags[site]) || (site<0) || (site>=molm->nsites) ) 
      return	0.0;
   p0	=	molm->p+site;
   i	=	site;
   while ( (k<3) && ((i=molm->parent[i])>=0) && (molm->flags[i]>0)) {
      p1	=	molm->p+i;
      dr[k++]	=	V_Subtr(p0, p1);
      p0	=	p1;
   }
   if (k<3) {
      return	0.0;
   }
   n0	=	V_Cross(dr+1, dr);
   n1	=	V_Cross(dr+2, dr+1);
   l0	=	sqrt(V_Dot(&n0, &n0));
   l1	=	sqrt(V_Dot(&n1, &n1));
   cosb	=	-V_Dot(&n0, &n1)/(l0*l1);
   t	=	type + molm->type[site];
/*   
   if (fabs(cosb-1)<ZERO)
      b	=	0;
   else if (fabs(cosb+1) <ZERO)
      b	=	M_PI;
   else
      b	=	acos(cosb);
*/
   if (samestr(TORTYPE, "opls")) 
      return	OPLS2(cosb, t->TORSION[1], t->TORSION[2], t->TORSION[3]);
   else if (samestr(TORTYPE, "poly"))
      return	POLY(cosb, t->TORSION[0], t->TORSION[1], t->TORSION[2], t->TORSION[3],
			t->TORSION[4], t->TORSION[5]);
   else {
      printf("Torsion type error!\n");
      exit(-1);
   }
}


double VTorsionMol(molstruct *molm)
{
   long		i;
   double	v_tor = 0.0;

   if (!V_TORSION || (molm->box<0))
      return	0.0;
   for (i=0; i<molm->nsites; i++)
      v_tor	+=	VTorsionSite(molm, i);
   return	v_tor;
}


void CalcVTorsion()
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++) 
      v[ibox].torsion	=	0.0;

   if (!V_TORSION)	return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox = moli->box) >=0 )
         v[ibox].torsion	+=	VTorsionMol(moli);	// no double-counting
}


void CalcVLJCorr()		// recalculate LJ long-range energy correction
{
   long		n, m, flags_m, *nsites, system;
   double	Epsilon, Sigma, rc3, rc9, vol;
   molstruct	*molm, *moln;
   typestruct	*typem;
   mixstruct	*typemix;
   long		NActive[MAXNSYSTEMS];

   if (!V_LJLRC)	return;

   if (!(nsites = (long *) calloc(NSYSTEMS, sizeof(long))))
      Exit("forcefield", "CalcVLJCorr", "Out of memory");

   for (system=0; system<NSYSTEMS; system++) {
      v[system].ljcorr	=	0.0;
      nsites[system]	=	0;
      NActive[system]	=	0;
   }

   for (molm=mol; molm<mol+NMOLS; molm++)
      if ( (system=molm->box) >= 0)
         for (m=0; m<molm->nsites; m++)
            if (molm->flags[m]>0)
                NActive[system]	++;

   for (molm=mol; molm<mol+NMOLS; molm++)
      if ( (system=molm->box) >= 0)
         for (m=0; m<molm->nsites; m++) 
            if ( (molm->flags[m]>0) && ((typem=type+molm->type[m])->EPSILON) )
               for (moln=mol; moln<mol+NMOLS; moln++)
                  if (moln->box == system)
	             for (n=0; n<moln->nsites; n++)
			if ( (moln->flags[n]>0) && (moln==molm ? abs(n-m)>=DLJ : 1) 
				&& (Epsilon=(typemix=(typem->mix)+moln->type[n])->EPSILON) ) {

			   Sigma	=	typemix->SIGMA;
	                   rc3		=	BOX[system].rc / Sigma;
		           rc3		=	rc3 * rc3 * rc3;
		           rc9		=	rc3 * rc3 * rc3;
		           vol		=	BOX[system].vol;
			   if (V_LJSHIFT)
                              v[system].ljcorr	+=	Epsilon/vol * Sigma * Sigma * Sigma * (4.0/(3*rc9)-2.0/rc3);
			   else
                              v[system].ljcorr	+=	Epsilon/vol * Sigma * Sigma * Sigma * (1.0/(3*rc9)-1.0/rc3);
			   nsites[system]	++;
                        }
         
   for (system=0; system<NSYSTEMS; system++)
      v[system].ljcorr	*=	nsites[system] ? 8.0*pi/3.0*NSites[system]*NActive[system]/nsites[system] : 0.0;
   free(nsites);
}   


void CalcVCorr()
{
   long		ib;

   if (V_LJLRC)
      CalcVLJCorr();

   for (ib=0; ib<NBOX; ib++) 
      v[ib].corr	=	v[ib].ljcorr;
}

/*
void CalcV_volchange1(long system)	// Calculate system energy for vol. change, 
				// the relative distance b/w sites on the same chain doesn't 
				// change, then only LJ interaction needs to be recalc.
{
   v[system].nonbonded	-=	(v[system].ljout + v[system].corr);     
   v[system].tot	-=	(v[system].ljout + v[system].corr);
   vir[system].tot	-=	vir[system].ljout;

   if (V_LJ)			// only LJ energy needs to be updated
      CalcVLJout();

   CalcVCorr();

   v[system].nonbonded	+=	(v[system].ljout + v[system].corr);     
   v[system].tot	+=	(v[system].ljout + v[system].corr);
   vir[system].tot	+=	vir[system].ljout;

   return;
}

void CalcV_volchange2()		// for vol. change, if relative distance b/w sites on the same chain
				// also scale with box size, then LJ interaction can be scaled, but
				// bond energy needs to be recalculated
{}
*/

void CalcV_mcvol(double volscale)	// recalc. energy for volume change move, volscale=volnew/volold
{
   long		ib;

   for (ib=0; ib<NBOX; ib++) {		// LJ and HS energy need to be recalculated, others not
      v[ib].lj		=	0.0;
      v[ib].hs		=	0.0;
      vir[ib].lj	=	0.0;
   }

   if (V_LJ)		CalcVLJ();
   if (V_HS)		CalcVHS();

   for (ib=0; ib<NBOX; ib++) {				// energy tail correction simple scaling
      v[ib].ljcorr	/=	volscale;		// if cutoff NOT scale with box size
      v[ib].corr	=	v[ib].ljcorr;
   }

   for (ib=0; ib<NBOX; ib++) {
      v[ib].bonded	=	v[ib].stretch + v[ib].bending + v[ib].torsion;
      v[ib].nonbonded	=	v[ib].lj + v[ib].hs + v[ib].corr;
      v[ib].tot		=	v[ib].bonded + v[ib].nonbonded;

      vir[ib].tot	=	vir[ib].lj + vir[ib].stretch + vir[ib].torsion;
   }
   return;
}


void CalcV()
{
   long		ib;

   for (ib=0; ib<NBOX; ib++) {
      vstructNull(v+ib);
      wstructNull(vir+ib);
   }
   if (V_LJ)		CalcVLJ();
   if (V_HS)		CalcVHS();
   if (V_STRETCH)	CalcVStretch();
   if (V_BENDING)	CalcVBending();
   if (V_TORSION)	CalcVTorsion();

   CalcVCorr();

   for (ib=0; ib<NBOX; ib++) {
      v[ib].bonded	=	v[ib].stretch + v[ib].bending + v[ib].torsion;
      v[ib].nonbonded	=	v[ib].lj + v[ib].hs + v[ib].corr;
      v[ib].tot		=	v[ib].bonded + v[ib].nonbonded;

      vir[ib].tot	=	vir[ib].lj + vir[ib].stretch + vir[ib].torsion;
   }
   return;
}


vstruct CalcVSite(molstruct *moli, long site, wstruct *VIRSite)	// energy and virial
{								// contribution of one
   long		i;						// single site to the system
   vstruct	V;

   vstructNull(&V);
   wstructNull(VIRSite);

   if (V_LJ) 
      V.lj		+=	VLJSite(moli, site, &(VIRSite->lj));
   if (V_HS)
      V.hs		+=	VHSSite(moli, site);
   if (V_STRETCH)
      for (i=0; i<2; i++) 
         V.stretch	+=	VStretchSite(moli, site+i, &(VIRSite->stretch));
   if (V_BENDING) 
      for (i=0; i<3; i++)
         V.bending	+=	VBendingSite(moli, site+i);
   if (V_TORSION)
      for (i=0; i<4; i++)
         V.torsion	+=	VTorsionSite(moli, site+i);

   // tail correction is calculated for the whole system, not for individual site or molecule
   V.bonded	=	V.stretch + V.torsion + V.bending;
   V.nonbonded	=	V.lj + V.hs;
   V.tot	=	V.bonded + V.nonbonded;

   if (V_VIRIAL)
      VIRSite->tot	=	VIRSite->lj + VIRSite->stretch + VIRSite->torsion;

   return	V;   
}


vstruct CalcVSiteinner(molstruct *moli, long site, wstruct *VIRSite)	// energy and virial
{								// contribution of one
   long		i;						// single site to the system
   vstruct	V;

   vstructNull(&V);
   wstructNull(VIRSite);

   if (V_LJ) 
      V.lj		+=	VLJSiteinner(moli, site, &(VIRSite->lj));
   if (V_HS)
      V.hs		+=	VHSSite(moli, site);
   if (V_STRETCH)
      for (i=0; i<2; i++) 
         V.stretch	+=	VStretchSite(moli, site+i, &(VIRSite->stretch));
   if (V_BENDING) 
      for (i=0; i<3; i++)
         V.bending	+=	VBendingSite(moli, site+i);
   if (V_TORSION)
      for (i=0; i<4; i++)
         V.torsion	+=	VTorsionSite(moli, site+i);

   // tail correction is calculated for the whole system, not for individual site or molecule
   V.bonded	=	V.stretch + V.torsion + V.bending;
   V.nonbonded	=	V.lj + V.hs;
   V.tot	=	V.bonded + V.nonbonded;

   if (V_VIRIAL)
      VIRSite->tot	=	VIRSite->lj + VIRSite->stretch + VIRSite->torsion;

   return	V;   
}


vstruct CalcVSiteouter(molstruct *moli, long site, wstruct *VIRSite)	// energy and virial
{								// contribution of one
   long		i;						// single site to the system
   vstruct	V;

   vstructNull(&V);
   wstructNull(VIRSite);

   if (V_LJ) 
      V.lj		+=	VLJSiteouter(moli, site, &(VIRSite->lj));

   // tail correction is calculated for the whole system, not for individual site or molecule
   V.bonded	=	V.stretch + V.torsion + V.bending;
   V.nonbonded	=	V.lj + V.hs;
   V.tot	=	V.bonded + V.nonbonded;

   if (V_VIRIAL)
      VIRSite->tot	=	VIRSite->lj + VIRSite->stretch + VIRSite->torsion;

   return	V;   
}


/* VDeleteSites, VAddSites, and grow update ALL energy and virial components */
/* bonded, nonbonded, tot, EXCEPT for long-range energy correction */

double VDeleteSites(molstruct *moli, long i_0, long i_n)
{
   long		i, j;
   vstruct	*v_new		=	v+moli->box;
   wstruct	*vir_new	=	vir+mol->box;
   double	v_old		=	v_new->tot;

   if (V_VIRIAL)
      wstructNegate(vir_new);		// need to be paired up in the end

   for (i=i_0; i<=i_n; i++) {
      
      if (V_STRETCH) 			// energy of other site(s) caused by this site
         for (j=1; j<2; j++)
            v_new->stretch	-=	VStretchSite(moli, i+j, &(vir_new->stretch));
      if (V_BENDING)
         for (j=1; j<3; j++)
            v_new->bending	-=	VBendingSite(moli, i+j);
      if (V_TORSION)
         for (j=1; j<4; j++)
            v_new->torsion	-=	VTorsionSite(moli, i+j);

      if (moli->flags[i] > 0) {		// energy of this site itself
         if (V_LJ) 
            v_new->lj		-=	VLJSite(moli, i, &(vir_new->lj));
         if (V_HS)
            v_new->hs		-=	VHSSite(moli, i);
         if (V_STRETCH) 
            v_new->stretch	-=	VStretchSite(moli, i, &(vir_new->stretch));
         if (V_BENDING)
	    v_new->bending	-=	VBendingSite(moli, i);
         if (V_TORSION)
            v_new->torsion	-=	VTorsionSite(moli, i);

#ifdef CELL_LIST
	 CL_Delete(moli, i);			// Unregister site.
#endif
         moli->flags[i]	=	0;		// Deactivate site. Important!! to avoid double-counting
      }
   }
   v_new->bonded	=	v_new->stretch + v_new->bending + v_new->torsion;
   v_new->nonbonded	=	v_new->lj + v_new->hs + v_new->corr;
   v_new->tot		=	v_new->bonded + v_new->nonbonded;
   if (V_VIRIAL) {				// update virial
      wstructNegate(vir_new);			// VIR = -VIR
      vir_new->tot	=	vir_new->lj + vir_new->stretch + vir_new->torsion;
   }
   return	v_new->tot - v_old;
}


double VAddSites(molstruct *moli, long i_0, long i_n)
{
   long		i, j;
   vstruct	*v_new		=	v+moli->box;
   wstruct	*vir_new	=	vir+moli->box;
   double	v_old		=	v_new->tot;

   for (i=i_0; i<=i_n; i++) {			// must be from parent to children

      if (moli->flags[i] == 0) {
         moli->flags[i]		=	1;		// activate this site
#ifdef CELL_LIST
	 CL_Add(moli, i);
#endif
         if (V_LJ)
            v_new->lj		+=	VLJSite(moli, i, &(vir_new->lj));
         if (V_HS)
            v_new->hs		+=	VHSSite(moli, i);
         if (V_STRETCH)
            v_new->stretch	+=	VStretchSite(moli, i, &(vir_new->stretch));
         if (V_BENDING)
	    v_new->bending	+=	VBendingSite(moli, i);
         if (V_TORSION)
            v_new->torsion	+=	VTorsionSite(moli, i);
      }
      if (V_STRETCH)
         for (j=1; j<2; j++)
            v_new->stretch	+=	VStretchSite(moli, i+j, &(vir_new->stretch));
      if (V_BENDING)
         for (j=1; j<3; j++)
            v_new->bending	+=	VBendingSite(moli, i+j);
      if (V_TORSION)
         for (j=1; j<4; j++)
            v_new->torsion	+=	VTorsionSite(moli, i+j);
   }
   v_new->bonded	=	v_new->stretch + v_new->bending + v_new->torsion;
   v_new->nonbonded	=	v_new->lj + v_new->hs + v_new->corr;
   v_new->tot		=	v_new->bonded + v_new->nonbonded;
   if (V_VIRIAL)
      vir_new->tot	=	vir_new->lj + vir_new->stretch + vir_new->torsion;
   return	v_new->tot - v_old;
}


long Select(double *wt, double sumw, long ntrial)
{
   double	ws, cumw;
   long		n;

   ws		=	ran1(seed) * sumw;
   n		=	0;
   cumw		=	wt[0];

   while (cumw < ws) {
      n		++;
      cumw	+=	wt[n];
   }
   if (n>=ntrial)
      Exit("forcefield", "Select", "n>ntrial.");
   return	n;
}


double grow(char *flag, molstruct *molm, long site)
{
   long		i, j, k, n, ib=molm->box;
   double	W, sumw, lbond;
   vector	dp;
   sphere	s;
   static long		init = 1;
   static double	*wt;
   static vector	*pt, *p0;
   static vstruct	*Vt;
   static wstruct	*VIRt;
   long		ntrial, f, old, new;
   vstruct		vinner, vouter;
   wstruct		winner, wouter;
   double		margin=1e-8;

   k	=	NTRIALCONF;
   f	=	NTRIALFIRSTBEAD;

   if (!strcmp(flag, "old"))		{	old	=	1;	new	=	0; }
   else if (!strcmp(flag, "new"))	{	new	=	1;	old	=	0; }
   else					Exit("forcefield", "grow", "flag not recognized");

   if (init) {
      ntrial	=	MAX(NTRIALCONF, NTRIALFIRSTBEAD);

      if (! (wt=(double*) calloc(ntrial, sizeof(double))) )		// allocate for trial variables
         Exit("forcefield", "grow", "out of memory!");
      if (! (pt=(vector *) calloc(ntrial, sizeof(vector))) )		// trial position
         Exit("forcefield", "grow", "out of memory!");
      if (! (Vt=(vstruct *) calloc(ntrial, sizeof(vstruct))) )		// trial energy
         Exit("forcefield", "grow", "out of memory!");
      if (! (VIRt=(wstruct *) calloc(ntrial, sizeof(wstruct))) )	// trial virial
         Exit("forcefield", "grow", "out of memory!");

      init	=	0;
   }

   if ( old ) {			// for energy and virial update for old and new conf.
      vstructNegate(v+ib);			// need to be paired up in the end
      wstructNegate(vir+ib);
   }

   for (i=site; i<molm->nsites; i++) {
      if (1 == molm->flags[i]) {		// if active site
         molm->flags[i]	=	0;		// deactivate sites
#ifdef CELL_LIST
         CL_Delete(molm, i);  			// remove these sites from cells
#endif
      }
   }

   W	=	0.0;				// Rosenbluth factor

   for (i=site; i<molm->nsites; i++) {		// grow the chain

      if (i==0)		ntrial	=	NTRIALFIRSTBEAD;	// f 
      else		ntrial	=	NTRIALCONF;		// k

      // Generate trial positions for each bead

      for (j=0; j<ntrial; j++) {

         molm->flags[i]	=	1;		// activate this site for energy calc.
						// in trail conf. generation

         if (old && j==0) 			// first trial for old conf. is itself
	    pt[j]	=	molm->p[i];	// in fact, pt[0]
         else {
            if (i==0) {				// first bead trial position
	       pt[j].x	=	(ran1(seed)-0.5) * BOX[ib].lx;
	       pt[j].y	=	(ran1(seed)-0.5) * BOX[ib].ly;
	       pt[j].z	=	(ran1(seed)-0.5) * BOX[ib].lz;
            }
            else {				// non-first bead trial position
/*
	       dp	=	tors_bonda(molm, i);
               lbond	=	bondl_g(type[0].LSTRETCH, BOX[ib].temp);
               dp	=	V_Mult(lbond, &dp);
	       pt[j]	=	V_Add(molm->p+i-1, &dp);	// trial position
*/
	       s.d	=	bondl_g(type[0].LSTRETCH, BOX[ib].temp);
	       if (i==1) {
		  dp	=	ranor();
	          dp	=	V_Mult(s.d, &dp);
		  pt[j]	=	V_Add(molm->p+i-1, &dp);
	       }
	       else if (i==2) {
		  s.alpha	=	bonda_g(type[0].THETA, BOX[ib].temp);
		  s.beta	=	(ran1(seed)-0.5) * 2 * M_PI;
	          pt[j]		=	SiteCartesian(molm, i, s);
	       }
	       else {
                  bonda_tors(type[0].THETA, BOX[ib].temp, &(s.alpha), &(s.beta));
		  //s.alpha	=	bonda_g(type[0].THETA, BOX[ib].temp);
		  //s.beta	=	tors(BOX[ib].temp);
	          pt[j]		=	SiteCartesian(molm, i, s);
	       }
	    }
            molm->p[i]	=	pt[j];		// pass trial conf. to mol. for energy calc.
         }

#ifdef CELL_LIST
         CL_Add(molm, i);			// add trial site into cell, for external
						// energy calc., it in principle should go
						// with activation, but the energy calc. evolved
						// in trial conf. generation doesn't require cell
						// list, so we don't have to call it so often
#endif
         Vt[j]	=	CalcVSite(molm, i, VIRt+j);	// need cell list from previous step
         //Vt[j]	=	CalcVSiteinner(molm, i, VIRt+j);	// smaller LJ cutoff to save time
         wt[j]	=	exp(-Vt[j].nonbonded/BOX[ib].temp);

         molm->flags[i]	=	0;		// deactivate this site
#ifdef CELL_LIST
         CL_Delete(molm, i);			// delete trial site from cell list
#endif
      }
      // Pick one trial conf. for each bead with probability

      sumw	=	0.0;				// sum of wt
      for (j=0; j<ntrial; j++) 
         sumw	+=	wt[j];
         
      W	+=	log(sumw);				// +log(sumw) rather than *sumw is to avoid 
							// extremely big number
      if (old)						// grow old conf.
         n	=	0;				// pick old position
      else						// grow new conf.
         n	=	Select(wt, sumw, ntrial);	// select one trial pos.

      molm->p[i]	=	pt[n];

#ifdef CELL_LIST
      CL_Add(molm, i);				// add trial site into cell
#endif
      molm->flags[i]	=	1;			

      //vouter		=	CalcVSiteouter(molm, i, &wouter);	// outer LJ layer
      v[ib]		=	vstructSum(v+ib, Vt+n);
      //v[ib]		=	vstructSum(v+ib, &vouter);
      if (V_VIRIAL) {
         vir[ib]	=	wstructSum(vir+ib, VIRt+n);
         //vir[ib]	=	wstructSum(vir+ib, &wouter);
      }
      //W	+=	(-vouter.nonbonded)/BOX[ib].temp;			// fix Rosenbluth factor
   }	// chain grow complete
   if (old) {
      vstructNegate(v+ib);
      wstructNegate(vir+ib);
   }
   return	W;
}

/*
double grow1(char *flag, molstruct *molm, long site)
{
   long		i, j, k, n, ib=molm->box;
   double	W, sumw, lbond;
   vector	dp;
   sphere	s;
   static long		init = 1;
   static double	*wt;
   static vector	*pt, *p0;
   static vstruct	*Vt;
   static wstruct	*VIRt;
   long		ntrial, f;

   k	=	NTRIALCONF;
   f	=	NTRIALFIRSTBEAD;

   if (init) {
      ntrial	=	MAX(NTRIALCONF, NTRIALFIRSTBEAD);

      if (! (wt=(double*) calloc(ntrial, sizeof(double))) )		// allocate for trial variables
         Exit("forcefield", "grow", "out of memory!");
      if (! (pt=(vector *) calloc(ntrial, sizeof(vector))) )		// trial position
         Exit("forcefield", "grow", "out of memory!");
      if (! (Vt=(vstruct *) calloc(ntrial, sizeof(vstruct))) )		// trial energy
         Exit("forcefield", "grow", "out of memory!");
      if (! (VIRt=(wstruct *) calloc(ntrial, sizeof(wstruct))) )	// trial virial
         Exit("forcefield", "grow", "out of memory!");

      init	=	0;
   }

   if (!strcmp(flag,"old") ) {				// for energy and virial update for old conf.
      vstructNegate(v+ib);			// need to be paired up in the end
      wstructNegate(vir+ib);
   }

   for (i=site; i<molm->nsites; i++) {
      if (1 == molm->flags[i]) {			// if active site
         molm->flags[i]	=	0;			// deactivate sites
#ifdef CELL_LIST
         CL_Delete(molm, i);  				// remove these sites from cells
#endif
      }
   }

   W	=	0.0;					// Rosenbluth factor

   for (i=site; i<molm->nsites; i++) {
      if (0==i) {				// grow a whole chain, first atom
         if (!strcmp(flag,"new") ) {				// new conf.
	    molm->p[0].x	=	(ran1(seed)-0.5) * BOX[ib].lx;
	    molm->p[0].y	=	(ran1(seed)-0.5) * BOX[ib].ly;
	    molm->p[0].z	=	(ran1(seed)-0.5) * BOX[ib].lz;
         }

#ifdef CELL_LIST
         CL_Add(molm, i);			// add into cell
#endif
         molm->flags[i]	=	1;		// activated this site 

         Vt[0]		=	CalcVSite(molm, i, VIRt);	// energy and virial contribution
         wt[j]	=	exp(-Vt[j].nonbonded/BOX[ib].temp);

         molm->flags[i]	=	0;		// deactivate this site
#ifdef CELL_LIST
         CL_Delete(molm, i);			// delete trial site from cell list
#endif
      }
      // Pick one trial conf. for each bead with probability

      sumw	=	0.0;				// sum of wt
      for (j=0; j<ntrial; j++) 
         sumw	+=	wt[j];
         
      W	+=	log(sumw);				// +log(sumw) rather than *sumw is to avoid 
							// extremely big number
      if (old)						// grow old conf.
         n	=	0;				// pick old position
      else						// grow new conf.
         n	=	Select(wt, sumw, ntrial);	// select one trial pos.

      molm->p[i]	=	pt[n];
      v[ib]		=	vstructSum(v+ib, Vt+n);
      if (V_VIRIAL)
         vir[ib]	=	wstructSum(vir+ib, VIRt+n);

#ifdef CELL_LIST
      CL_Add(molm, i);				// add trial site into cell
#endif
      molm->flags[i]	=	1;			
   }	// chain grow complete
   if (old) {
      vstructNegate(v+ib);
      wstructNegate(vir+ib);
   }
   return	W;
}


double grow1(char *flag, molstruct *molm, long site)
{
   long		i, j, k, n, ib=molm->box;
   double	W, sumw, lbond;
   vector	dp;
   sphere	s;
   static long		init = 1;
   static double	*wt;
   static vector	*pt, *p0;
   static vstruct	*Vt;
   static wstruct	*VIRt;
   long		ntrial, f;

   k	=	NTRIALCONF;
   f	=	NTRIALFIRSTBEAD;

   if (init) {
      ntrial	=	MAX(NTRIALCONF, NTRIALFIRSTBEAD);

      if (! (wt=(double*) calloc(ntrial, sizeof(double))) )		// allocate for trial variables
         Exit("forcefield", "grow", "out of memory!");
      if (! (pt=(vector *) calloc(ntrial, sizeof(vector))) )		// trial position
         Exit("forcefield", "grow", "out of memory!");
      if (! (Vt=(vstruct *) calloc(ntrial, sizeof(vstruct))) )		// trial energy
         Exit("forcefield", "grow", "out of memory!");
      if (! (VIRt=(wstruct *) calloc(ntrial, sizeof(wstruct))) )	// trial virial
         Exit("forcefield", "grow", "out of memory!");

      init	=	0;
   }

   if (!strcmp(flag,"old") ) {				// for energy and virial update for old conf.
      vstructNegate(v+ib);			// need to be paired up in the end
      wstructNegate(vir+ib);
   }

   for (i=site; i<molm->nsites; i++) {
      if (1 == molm->flags[i]) {			// if active site
         molm->flags[i]	=	0;			// deactivate sites
#ifdef CELL_LIST
         CL_Delete(molm, i);  				// remove these sites from cells
#endif
      }
   }

   W	=	0.0;					// Rosenbluth factor

   for (i=site; i<molm->nsites; i++) {
      if (0==i) {				// grow a whole chain, first atom
         if (!strcmp(flag,"new") ) {				// new conf.
	    molm->p[0].x	=	(ran1(seed)-0.5) * BOX[ib].lx;
	    molm->p[0].y	=	(ran1(seed)-0.5) * BOX[ib].ly;
	    molm->p[0].z	=	(ran1(seed)-0.5) * BOX[ib].lz;
         }

#ifdef CELL_LIST
         CL_Add(molm, i);			// add into cell
#endif
         molm->flags[i]	=	1;		// activated this site 

         Vt[0]		=	CalcVSite(molm, i, VIRt);	// energy and virial contribution
         v[ib]		=	vstructSum(v+ib, Vt);		// of this site to the system

         if (V_VIRIAL)
            vir[ib]	=	wstructSum(vir+ib, VIRt);

         W	=	log(k) - Vt[0].nonbonded/BOX[ib].temp;
      }
      else {						// not the first atom

         for (j=0; j<k; j++) {				// search thru trial orient.

            molm->flags[i]	=	1;	// activate this site, before tors_bonda()

            if (!strcmp(flag, "old") && (0==j))
               pt[0]	=	molm->p[i];
            else {
	       // generate trial conf. according to Vbonded

	       dp	=	tors_bonda(molm, i);
               lbond	=	bondl_g(type[0].LSTRETCH, BOX[ib].temp);
               dp	=	V_Mult(lbond, &dp);

	       pt[j]	=	V_Add(molm->p+i-1, &dp);	// trial position
               //MapInBox2(pt+j, PBC, BOX[ib].lbox);
               molm->p[i]	=	pt[j];

	       //s.d	=	bondl_g(type[0].LSTRETCH, BOX[ib].temp);
               //s.alpha	=	(i>1 ? bonda_g(type[0].THETA, BOX[ib].temp) : 0.0);
	       //s.beta	=	(i>2 ? tors(BOX[ib].temp) : 0.0);
		
	       //pt[j]	=	SiteCartesian(molm, i, s);
               //molm->p[i]	=	pt[j];
            }

#ifdef CELL_LIST
	    CL_Add(molm, i);				// add trial site into cell
#endif

            Vt[j]	=	CalcVSite(molm, i, VIRt+j);
            wt[j]	=	exp(-Vt[j].nonbonded/BOX[ib].temp);

            molm->flags[i]	=	0;		// deactivate this site
#ifdef CELL_LIST
	    CL_Delete(molm, i);				// delete trial site from cell list
#endif
            //printf("%s\ttrial %d\t of %d\t, vt[j].tot=%f\tvt[j].nonbonded=%f\n", s, j, i, Vt[j].tot, Vt[j].nonbonded);
         }

         sumw	=	0.0;				// sum of wt
         for (j=0; j<k; j++) 
	    sumw	+=	wt[j];
         
         W	+=	log(sumw);

	 if (!strcmp(flag, "old"))				// grow old conf.
            n	=	0;				// pick old position
         else						// grow new conf.
            n	=	Select(wt, sumw, k);		// select one trial pos.

         molm->p[i]	=	pt[n];
         v[ib]		=	vstructSum(v+ib, Vt+n);
         if (V_VIRIAL)
            vir[ib]	=	wstructSum(vir+ib, VIRt+n);

#ifdef CELL_LIST
         CL_Add(molm, i);				// add trial site into cell
#endif
         molm->flags[i]	=	1;			

	 //printf("Adding %d done, trial #%d.\n", i, n);
      }		// not the atom #0 done
   }		// all atoms done
   if (!strcmp(flag, "old")) {
      vstructNegate(v+ib);
      wstructNegate(vir+ib);
   }
   return	W;
}
*/

void EnergyCheck()		// Check step-by-step energy update correctness
{
   double	*V, *VIR;
   long		i;

   V	=	(double *) calloc(NBOX, sizeof(double));
   VIR	=	(double *) calloc(NBOX, sizeof(double));

   for (i=0; i<NBOX; i++) {
      V[i]	=	v[i].tot;
      VIR[i]	=	vir[i].tot;
   }
   CalcV();

   for (i=0; i<NBOX; i++) {
      if ( fabs(V[i] - v[i].tot) > 1.0e-6)
          fprintf(foutput, "\n#############total energy problem in box %d.###############\n", i);
      if ( fabs(VIR[i] - vir[i].tot) > 1.0e-6)
          fprintf(foutput, "\n#############total virial problem in box %d.###############\n", i);

      fprintf(foutput, "\n\tBOX[%d]\n\n", i);
      fprintf(foutput, "\tTotal energy end of simulation:\t%f\n", V[i]);
      fprintf(foutput, "\t      energy recalculated:\t%f\n", v[i].tot);
      fprintf(foutput, "\tTotal virial end of simulation:\t%f\n", VIR[i]);
      fprintf(foutput, "\t      virial recalculated:\t%f\n", vir[i].tot);
      fprintf(foutput, "\n");
   }

   free(VIR);
   free(V);
   return;
}
                                                                                                                                                        src/histogram.c                                                                                     0000600 0143352 0000144 00000004360 10714730275 013226  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
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


                                                                                                                                                                                                                                                                                src/history.c                                                                                       0000600 0143352 0000144 00000067735 10772306334 012750  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	history.c
    author:	Peng Yi borrowed from Pieter J. in 't Veld
    date:	January 10, 2008
    purpose:	Binary history file module.
*/
#define __HISTORY_MODULE
#include "history.h"

int fputd(double *d, FILE *stream)
{
  fwrite(d, sizeof(double), 1, stream);
  return ferror(stream)!=0;
}


int fgetd(double *d, FILE *stream)
{
  fread(d, sizeof(double), 1, stream);
  return feof(stream)!=0;
}


int fputl(long *l, FILE *stream)
{
  fwrite(l, sizeof(long), 1, stream);
  return ferror(stream)!=0;
}


int fgetl(long *l, FILE *stream)
{
  fread(l, sizeof(long), 1, stream);
  return feof(stream)!=0;
}


int fput_vector(vector *r, FILE *stream)
{
  if (fputd(&(r->x), stream)) return 1;
  if (fputd(&(r->y), stream)) return 1;
  if (fputd(&(r->z), stream)) return 1;
  return 0;
}


int fget_vector(vector *r, FILE *stream)
{
  if (fgetd(&(r->x), stream)) return 1;
  if (fgetd(&(r->y), stream)) return 1;
  if (fgetd(&(r->z), stream)) return 1;
  return 0;
}


int fput_matrix(matrix *m, FILE *stream)
{
  if (fput_vector(&(m->x), stream)) return 1;
  if (fput_vector(&(m->y), stream)) return 1;
  if (fput_vector(&(m->z), stream)) return 1;
  return 0;
}


int fget_matrix(matrix *m, FILE *stream)
{
  if (fget_vector(&(m->x), stream)) return 1;
  if (fget_vector(&(m->y), stream)) return 1;
  if (fget_vector(&(m->z), stream)) return 1;
  return 0;
}


int fput_sphere(sphere *s, FILE *stream)
{
  if (fputd(&(s->d), stream)) return 1;
  if (fputd(&(s->alpha), stream)) return 1;
  if (fputd(&(s->beta), stream)) return 1;
  return 0;
}


int fget_sphere(sphere *s, FILE *stream)
{
  if (fgetd(&(s->d), stream)) return 1;
  if (fgetd(&(s->alpha), stream)) return 1;
  if (fgetd(&(s->beta), stream)) return 1;
  return 0;
}

/* More customized functions */

int fput_mol(molstruct *mol, FILE *stream)
{
  int			i;
  
  if (fputl(&(mol->nsites), stream)) return 1;
  if (fputl(&(mol->box), stream)) return 1;
  if (fputl(&(mol->fix), stream)) return 1;
  if (fputl(&(mol->flip), stream)) return 1;
  for (i=0; i<mol->nsites; ++i)
  {
    if (fputl(mol->type+i, stream)) return 1;
    if (fputl(mol->parent+i, stream)) return 1;
    if (fputl(mol->flags+i, stream)) return 1;
    if (fput_vector(mol->p+i, stream)) return 1;
  }
  return 0;
}
  

int fget_mol(molstruct *mol, FILE *stream, long version)
{
  int			i;
  
  if (fgetl(&(mol->nsites), stream)) return 1;
  if (fgetl(&(mol->box), stream)) return 1;
  if (fgetl(&(mol->fix), stream)) return 1;
  if (fgetl(&(mol->flip), stream)) return 1;
  if ((version<17)&&fget_vector(mol->p, stream)) return 1;
  for (i=0; i<mol->nsites; ++i)
  {
    if (fgetl(mol->type+i, stream)) return 1;
    if (fgetl(mol->parent+i, stream)) return 1;
    if (fgetl(mol->flags+i, stream)) return 1;
    if (fget_vector(mol->p+i, stream)) return 1;
  }
  return 0;
}


int fput_type(typestruct *type, FILE *stream)
{
  if (fputd(&(type->M), stream)) return 1;
  if (fputd(&(type->Q), stream)) return 1;
  if (fputd(&(type->SIGMA), stream)) return 1;
  if (fputd(&(type->EPSILON), stream)) return 1;
  if (fputd(&(type->LSTRETCH),stream)) return 1;
  if (fputd(&(type->KSTRETCH),stream)) return 1;
  if (fputd(&(type->KBENDING),stream)) return 1;
  if (fputd(&(type->THETA),stream)) return 1;
  if (fputd(type->TORSION,stream)) return 1;
  if (fputd(type->TORSION+1,stream)) return 1;
  if (fputd(type->TORSION+2,stream)) return 1;
  if (fputd(type->TORSION+3,stream)) return 1;
  if (fputd(type->TORSION+4,stream)) return 1;
  if (fputd(type->TORSION+5,stream)) return 1;
  return 0;
}


int fget_type(typestruct *type, FILE *stream, long version)
{
  if (fgetd(&(type->M), stream)) return 1;
  if (fgetd(&(type->Q), stream)) return 1;
  if (fgetd(&(type->SIGMA), stream)) return 1;
  if (fgetd(&(type->EPSILON), stream)) return 1;
  if (fgetd(&(type->LSTRETCH),stream)) return 1;
  if (fgetd(&(type->KSTRETCH),stream)) return 1;
  if (fgetd(&(type->KBENDING),stream)) return 1;
  if (fgetd(&(type->THETA),stream)) return 1;
  if (fgetd(type->TORSION,stream)) return 1;
  if (fgetd(type->TORSION+1,stream)) return 1;
  if (fgetd(type->TORSION+2,stream)) return 1;
  if (fgetd(type->TORSION+3,stream)) return 1;
  if (fgetd(type->TORSION+4,stream)) return 1;
  if (fgetd(type->TORSION+5,stream)) return 1;
  return 0;
}


int fput_system(systemstruct *system, FILE *stream)
{
  if (fputl(&(system->n), stream)) return 1;
  if (fputl(system->nmols, stream)) return 1;
  if (fputl(system->nsites, stream)) return 1;
  if (fputd(system->pres, stream)) return 1;
  if (fputd(system->vol, stream)) return 1;
  if (fputd(system->temp, stream)) return 1;
  if (fputd(system->drmax, stream)) return 1;
  if (fputd(system->dlmax, stream)) return 1;
  if (fputd(system->damax, stream)) return 1;
  /*
  if (fput_vector(&(system->a), stream)) return 1;
  if (fput_vector(&(system->b), stream)) return 1;
  if (fput_vector(&(system->c), stream)) return 1;
  */
  if (fputd(system->lx, stream)) return 1;
  if (fputd(system->ly, stream)) return 1;
  if (fputd(system->lz, stream)) return 1;
  return 0;
}


int fget_system(systemstruct *system, FILE *stream, long version)
{
  if (fgetl(&(system->n), stream)) return 1;
  if (fgetl(system->nmols, stream)) return 1;
  if (fgetl(system->nsites, stream)) return 1;
  if (fgetd(system->pres, stream)) return 1;
  if (fgetd(system->vol, stream)) return 1;
  if (fgetd(system->temp, stream)) return 1;
  if (fgetd(system->drmax, stream)) return 1;
  if (fgetd(system->dlmax, stream)) return 1;
  if (fgetd(system->damax, stream)) return 1;
  /*
  if (fget_vector(&(system->a), stream)) return 1;
  if (fget_vector(&(system->b), stream)) return 1;
  if (fget_vector(&(system->c), stream)) return 1;
  */
  if (fgetd(system->lx, stream)) return 1;
  if (fgetd(system->ly, stream)) return 1;
  if (fgetd(system->lz, stream)) return 1;

  return 0;
}


int fput_command(commandstruct *command, FILE *stream)
{
  if (fputl(command->hs, stream)) return 1;
  if (fputl(command->lj, stream)) return 1;
  if (fputl(command->ljshift, stream)) return 1;
  if (fputl(command->ljlrc, stream)) return 1;
  if (fputl(command->stretch, stream)) return 1;
  if (fputl(command->bending, stream)) return 1;
  if (fputl(command->torsion, stream)) return 1;
  if (fputl(command->virial, stream)) return 1;
  if (fputl(command->scalecut, stream)) return 1;
  if (fputl(command->nvt, stream)) return 1;
  if (fputl(command->npt, stream)) return 1;
  if (fputl(command->gibbs, stream)) return 1;
  if (fputl(command->mpi, stream)) return 1;
  if (fputl(command->density, stream)) return 1;
  if (fputl(command->energy, stream)) return 1;
  if (fputl(command->pressure, stream)) return 1;
  if (fputl(command->drift, stream)) return 1;
  if (fputl(command->torsional, stream)) return 1;
  if (fputl(command->bonda, stream)) return 1;
  if (fputl(command->bondl, stream)) return 1;
  if (fputl(command->radial, stream)) return 1;
  if (fputl(command->localp2, stream)) return 1;
  if (fputl(command->xtalsize, stream)) return 1;

/*
  if (fputl(command->hs, stream)) return 1;
  if (fputl(command->lj, stream)) return 1;
  if (fputl(command->coulomb, stream)) return 1;
  if (fputl(command->polymer, stream)) return 1;
  if (fputl(command->mpi, stream)) return 1;
  if (fputl(command->async, stream)) return 1;
  if (fputl(command->bias, stream)) return 1;
  if (fputl(command->nvt, stream)) return 1;
  if (fputl(command->npt, stream)) return 1;
  if (fputl(command->gibbs, stream)) return 1;
  if (fputl(command->insert, stream)) return 1;
  if (fputl(command->widom, stream)) return 1;
  if (fputl(command->canonical, stream)) return 1;
  if (fputl(command->cavity, stream)) return 1;
  if (fputl(command->density, stream)) return 1;
  if (fputl(command->density3d, stream)) return 1;
  if (fputl(command->densfree, stream)) return 1;
  if (fputl(command->denstalobr, stream)) return 1;
  if (fputl(command->orient, stream)) return 1;
  if (fputl(command->orientfree, stream)) return 1;
  if (fputl(command->orienttalobr, stream)) return 1;
  if (fputl(command->orientcorr01, stream)) return 1;
  if (fputl(command->orientcorr02, stream)) return 1;
  if (fputl(command->orientcorr03, stream)) return 1;
  if (fputl(command->orientcorr04, stream)) return 1;
  if (fputl(command->orientcorr05, stream)) return 1;
  if (fputl(command->orientcorrCR, stream)) return 1;
  if (fputl(command->tails_etc, stream)) return 1;
  if (fputl(command->radial, stream)) return 1;
  if (fputl(command->energy, stream)) return 1;
  if (fputl(command->hs_dens, stream)) return 1;
  if (fputl(command->temper, stream)) return 1;
  if (fputl(command->jacob, stream)) return 1;
  if (fputl(command->d_bridge, stream)) return 1;
  if (fputl(command->torsion, stream)) return 1;
  if (fputl(command->re_torsion, stream)) return 1;
  if (fputl(command->e_profile, stream)) return 1;
  if (fputl(command->n_profile, stream)) return 1;
  if (fputl(command->b_length, stream)) return 1;
  if (fputl(command->b_angle, stream)) return 1;
  if (fputl(command->ptorsion, stream)) return 1;
  if (fputl(command->loopreentry, stream)) return 1;
  if (fputl(command->virial, stream)) return 1;
  if (fputl(command->e_n_function, stream)) return 1;
  if (fputl(command->monodisperse, stream)) return 1;
  if (fputl(command->w_profile, stream)) return 1;
  if (fputl(command->stretch, stream)) return 1;
*/
  return 0;
}


int fget_command(commandstruct *command, FILE *stream, long version)
{
  long			extra_l;

  if (fgetl(command->hs, stream)) return 1;
  if (fgetl(command->lj, stream)) return 1;
  if (fgetl(command->ljshift, stream)) return 1;
  if (fgetl(command->ljlrc, stream)) return 1;
  if (fgetl(command->stretch, stream)) return 1;
  if (fgetl(command->bending, stream)) return 1;
  if (fgetl(command->torsion, stream)) return 1;
  if (fgetl(command->virial, stream)) return 1;
  if (fgetl(command->scalecut, stream)) return 1;
  if (fgetl(command->nvt, stream)) return 1;
  if (fgetl(command->npt, stream)) return 1;
  if (fgetl(command->gibbs, stream)) return 1;
  if (fgetl(command->mpi, stream)) return 1;
  if (fgetl(command->density, stream)) return 1;
  if (fgetl(command->energy, stream)) return 1;
  if (fgetl(command->pressure, stream)) return 1;
  if (fgetl(command->drift, stream)) return 1;
  if (fgetl(command->torsional, stream)) return 1;
  if (fgetl(command->bonda, stream)) return 1;
  if (fgetl(command->bondl, stream)) return 1;
  if (fgetl(command->radial, stream)) return 1;
  if (fgetl(command->localp2, stream)) return 1;
  if (fgetl(command->xtalsize, stream)) return 1;
/*  
  if (fgetl(command->hs, stream)) return 1;
  if (fgetl(command->lj, stream)) return 1;
  if (fgetl(command->coulomb, stream)) return 1;
  if (fgetl(command->polymer, stream)) return 1;
  if ((version>19)&&fgetl(command->mpi, stream)) return 1;
  if ((version>19)&&fgetl(command->async, stream)) return 1;
  if (fgetl(command->bias, stream)) return 1;
  if (fgetl(command->nvt, stream)) return 1;
  if (fgetl(command->npt, stream)) return 1;
  if (fgetl(command->gibbs, stream)) return 1;
  if (fgetl(command->insert, stream)) return 1;
  if (fgetl(command->widom, stream)) return 1;
  if (fgetl(command->canonical, stream)) return 1;
  if (fgetl(command->cavity, stream)) return 1;
  if (fgetl(command->density, stream)) return 1;
  if ((version>23)&&fgetl(command->density3d, stream)) return 1;
  if ((version>22)&&fgetl(command->densfree, stream)) return 1;
  if ((version>22)&&fgetl(command->denstalobr, stream)) return 1;
  if ((version>19)&&fgetl(command->orient, stream)) return 1;
  if ((version>22)&&fgetl(command->orientfree, stream)) return 1;
  if ((version>22)&&fgetl(command->orienttalobr, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr01, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr02, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr03, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr04, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr05, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorrCR, stream)) return 1;
  if (fgetl(command->tails_etc, stream)) return 1;
  if (fgetl(command->radial, stream)) return 1;
  if (fgetl(command->energy, stream)) return 1;
  if (fgetl(command->hs_dens, stream)) return 1;
  if (fgetl(command->temper, stream)) return 1;
  if (version<2)
  {
    if (!fgetl(&extra_l, stream)) return 1;
    if (!fgetl(&extra_l, stream)) return 1;
    if (!fgetl(&extra_l, stream)) return 1;
    return fgetl(&extra_l, stream);
  }
  if (version<7)
    return 0;
  if (fgetl(command->jacob, stream)) return 1;
  if (fgetl(command->d_bridge, stream)) return 1;
  if (version<8)
    return 0;
  if (fgetl(command->torsion, stream)) return 1;
  if ((version>24)&&fgetl(command->re_torsion, stream)) return 1;
  if (fgetl(command->e_profile, stream)) return 1;
  if (fgetl(command->n_profile, stream)) return 1;
  if (version<9)
    return 0;
  if (fgetl(command->b_length, stream)) return 1;
  if (fgetl(command->b_angle, stream)) return 1;
  if ((version>9)&&fgetl(command->ptorsion, stream)) return 1;
  if ((version>10)&&fgetl(command->loopreentry, stream)) return 1;
  if ((version>11)&&fgetl(command->virial, stream)) return 1;
  if ((version>13)&&fgetl(command->e_n_function, stream)) return 1;
  if ((version>17)&&fgetl(command->monodisperse, stream)) return 1;
  if ((version>21)&&fgetl(command->w_profile, stream)) return 1;
  if (fgetl(command->stretch, stream)) return 1;
*/
  return 0;
}


int fput_n(nstruct *n, FILE *stream)
{
  if (fputl(n->cycle, stream)) return 1;
  if (fputl(n->tape, stream)) return 1;
  if (fputl(n->systems, stream)) return 1;
  if (fputl(n->mols, stream)) return 1;
  if (fputl(n->sites, stream)) return 1;
  if (fputl(n->types, stream)) return 1;
  if (fputl(n->volchange, stream)) return 1;
  if (fputl(n->swap, stream)) return 1;
  if (fputl(n->displace, stream)) return 1;
  if (fputl(n->cbmc, stream)) return 1;
/*
  if (fputl(n->cycle, stream)) return 1;
  if (fputl(n->systems, stream)) return 1;
  if (fputl(n->mols, stream)) return 1;
  if (fputl(n->sites, stream)) return 1;
  if (fputl(n->types, stream)) return 1;
  if (fputl(n->volumes, stream)) return 1;
  if (fputl(n->swaps, stream)) return 1;
  if (fputl(n->inserts, stream)) return 1;
  if (fputl(n->cavinserts, stream)) return 1;
  if (fputl(n->blocks, stream)) return 1;
  if (fputl(n->cycles, stream)) return 1;
  if (fputl(n->box, stream)) return 1;
  if (fputl(n->rotation, stream)) return 1;
  if (fputl(n->reptation, stream)) return 1;
  if (fputl(n->endbridge, stream)) return 1;
  if (fputl(n->rebridge, stream)) return 1;
  if (fputl(n->bridges, stream)) return 1;
  if (fputl(n->fixed, stream)) return 1;
  if (fputl(n->semifixed, stream)) return 1;
  if (fputl(n->seed, stream)) return 1;
  if (fputl(n->temper, stream)) return 1;
  if (fputl(n->sample, stream)) return 1;
  if (fputl(n->stretch, stream)) return 1;
  if (fputl(n->free, stream)) return 1;
  if (fputl(n->ends, stream)) return 1;
  if (fputl(n->system, stream)) return 1;
*/
  return 0;
}


int fget_n(nstruct *n, FILE *stream, long version)
{
  long			extra_l;

  if (fgetl(n->cycle, stream)) return 1;
  if (fgetl(n->tape, stream)) return 1;
  if (fgetl(n->systems, stream)) return 1;
  if (fgetl(n->mols, stream)) return 1;
  if (fgetl(n->sites, stream)) return 1;
  if (fgetl(n->types, stream)) return 1;
  if (fgetl(n->volchange, stream)) return 1;
  if (fgetl(n->swap, stream)) return 1;
  if (fgetl(n->displace, stream)) return 1;
  if (fgetl(n->cbmc, stream)) return 1;
/*  if (fgetl(n->cycle, stream)) return 1;
  if (fgetl(n->systems, stream)) return 1;
  if (fgetl(n->mols, stream)) return 1;
  if (fgetl(n->sites, stream)) return 1;
  if (fgetl(n->types, stream)) return 1;
  if (fgetl(n->volumes, stream)) return 1;
  if (fgetl(n->swaps, stream)) return 1;
  if (fgetl(n->inserts, stream)) return 1;
  if (fgetl(n->cavinserts, stream)) return 1;
  if (fgetl(n->blocks, stream)) return 1;
  if (fgetl(n->cycles, stream)) return 1;
  if (fgetl(n->box, stream)) return 1;
  if (fgetl(n->rotation, stream)) return 1;
  if (fgetl(n->reptation, stream)) return 1;
  if (fgetl(n->endbridge, stream)) return 1;
  if (fgetl(n->rebridge, stream)) return 1;
  if ((version<=19)&&fgetl(&extra_l, stream)) return 1; // loops
  if ((version<=19)&&fgetl(&extra_l, stream)) return 1; // tails
  if (fgetl(n->bridges, stream)) return 1;
  if (fgetl(n->fixed, stream)) return 1;
  if ((version>20)&&fgetl(n->semifixed, stream)) return 1;
  if (fgetl(n->seed, stream)) return 1;
  if (fgetl(n->temper, stream)) return 1;
  if (version<2)
  {
    if (fgetl(&extra_l, stream)) return 1;
    if (fgetl(&extra_l, stream)) return 1;
    if (fgetl(&extra_l, stream)) return 1;
    return fgetl(&extra_l, stream);
  }
  if (fgetl(n->sample, stream)) return 1;
  if (version<9) 
    return 0;
  if (fgetl(n->stretch, stream)) return 1;
  if (version<15) 
    return 0;
  if (fgetl(n->free, stream)) return 1;
  if (fgetl(n->ends, stream)) return 1;
  if ((version>19)&&fgetl(n->system, stream)) return 1;
*/
  return 0;
}


int fput_av(avstruct *av, FILE *stream)
{
  if (fputl(&(av->move), stream)) return 1;
  if (fputl(&(av->acc_move), stream)) return 1;
  if (fputl(&(av->rot), stream)) return 1;
  if (fputl(&(av->acc_rot), stream)) return 1;
  if (fputl(&(av->vol), stream)) return 1;
  if (fputl(&(av->acc_vol), stream)) return 1;
  if (fputl(&(av->cbmc), stream)) return 1;
  if (fputl(&(av->acc_cbmc), stream)) return 1;
  if (fputl(&(av->rep), stream)) return 1;
  if (fputl(&(av->acc_rep), stream)) return 1;
  if (fputl(&(av->seq), stream)) return 1;
  if (fputl(&(av->acc_seq), stream)) return 1;
  return 0;
}


int fget_av(avstruct *av, FILE *stream, long version)
{
  if (fgetl(&(av->move), stream)) return 1;
  if (fgetl(&(av->acc_move), stream)) return 1;
  if (fgetl(&(av->rot), stream)) return 1;
  if (fgetl(&(av->acc_rot), stream)) return 1;
  if (fgetl(&(av->vol), stream)) return 1;
  if (fgetl(&(av->acc_vol), stream)) return 1;
  if (fgetl(&(av->cbmc), stream)) return 1;
  if (fgetl(&(av->acc_cbmc), stream)) return 1;
  if (fgetl(&(av->rep), stream)) return 1;
  if (fgetl(&(av->acc_rep), stream)) return 1;
  if (fgetl(&(av->seq), stream)) return 1;
  if (fgetl(&(av->acc_seq), stream)) return 1;
  return 0;
}


int fput_v(vstruct *v, FILE *stream)
{
  if (fputd(&(v->tot), stream)) return 1;
  if (fputd(&(v->bonded), stream)) return 1;
  if (fputd(&(v->nonbonded), stream)) return 1;
  if (fputd(&(v->hs), stream)) return 1;
  if (fputd(&(v->lj), stream)) return 1;
  if (fputd(&(v->ljcorr), stream)) return 1;
  if (fputd(&(v->corr), stream)) return 1;
  if (fputd(&(v->stretch), stream)) return 1;
  if (fputd(&(v->bending), stream)) return 1;
  if (fputd(&(v->torsion), stream)) return 1;
  return 0;
}


int fget_v(vstruct *v, FILE *stream, long version)
{
  if (fgetd(&(v->tot), stream)) return 1;
  if (fgetd(&(v->bonded), stream)) return 1;
  if (fgetd(&(v->nonbonded), stream)) return 1;
  if (fgetd(&(v->hs), stream)) return 1;
  if (fgetd(&(v->lj), stream)) return 1;
  if (fgetd(&(v->ljcorr), stream)) return 1;
  if (fgetd(&(v->corr), stream)) return 1;
  if (fgetd(&(v->stretch), stream)) return 1;
  if (fgetd(&(v->bending), stream)) return 1;
  if (fgetd(&(v->torsion), stream)) return 1;
  return 0;
}


int fput_w(wstruct *w, FILE *stream)
{
  if (fputd(&(w->tot), stream)) return 1;
  if (fputd(&(w->lj), stream)) return 1;
  if (fputd(&(w->stretch), stream)) return 1;
  if (fputd(&(w->torsion), stream)) return 1;
  return 0;
}


int fget_w(wstruct *w, FILE *stream, long version)
{
  if (fgetd(&(w->tot), stream)) return 1;
  if (fgetd(&(w->lj), stream)) return 1;
  if (fgetd(&(w->stretch), stream)) return 1;
  if (fgetd(&(w->torsion), stream)) return 1;
  return 0;
}


int fput_dist(diststruct *dist, int top, FILE *stream)
{
  long                  i;

  if (fputl(&(dist->nbins), stream)) return 1;
  if (fputl(&(dist->n), stream)) return 1;
  if (fputl(&(dist->ncount), stream)) return 1;
  if (fputl(&(dist->startbin), stream)) return 1;
  if (fputl(&(dist->level), stream)) return 1;
  if (top)
    for (i=0; i<=dist->level; ++i)
      if (fputd(dist->binsize+i, stream)) return 1;

  if (dist->level)				// distribution fork
  {
    for (i=0; i<dist->nbins; ++i)
      if (fput_dist(dist->dist+i, FALSE, stream)) return 1;
  }
  else if (dist->nbins)				// data fork
  {
    for (i=0; i<D_NAVERAGE; ++i)
      if (fputd(dist->average+i, stream)) return 1;
    for (i=0; i<dist->nbins; ++i)
    {
      if (fputl(dist->bin+i, stream)) return 1;
      if (fputd(dist->data+i, stream)) return 1;
      if (fputd(dist->cweight+i, stream)) return 1;
    }
  }
  return 0;
}


int fget_dist(diststruct *dist, double *binsize, FILE *stream, long version)
{
  long                  i;

  if (fgetl(&(dist->nbins), stream)) return 1;
  if (fgetl(&(dist->n), stream)) return 1;
  if ((version>19)&&fgetl(&(dist->ncount), stream)) return 1;
  if (fgetl(&(dist->startbin), stream)) return 1;

  if (version<24) dist->level = 0;
  else if (fgetl(&(dist->level), stream)) return 1;
  if (!binsize)
  {
    for (i=0; i<=dist->level; ++i)
      if (fgetd(dist->binsize+i, stream)) return 1;
  }
  else dist->binsize	= binsize;

  if (dist->level)				// distribution fork
  {
    if (dist->nbins)
    {
      if (!(dist->dist = calloc(dist->nbins, sizeof(diststruct))))
	Exit("history", "fget_dist", "distribution fork calloc error");
      for (i=0; i<dist->nbins; ++i)
	if (fget_dist(dist->dist+i, dist->binsize+1, stream, version)) return 1;
    }
  }
  else						// data fork
  {
    if ((version<24)||dist->nbins)
    {
      if (!(dist->average = (double *) calloc(D_NAVERAGE, sizeof(double))))
        Exit("history", "fget_dist", "average calloc error");
      for (i=0; i<D_NAVERAGE; ++i)
        if (fgetd(dist->average+i, stream)) return 1;
    }
    if (dist->nbins)
    {
      if (!((dist->bin = (long *) calloc(dist->nbins, sizeof(long)))&&
            (dist->data = (double *) calloc(dist->nbins, sizeof(double)))&&
            (dist->cweight = (double *) calloc(dist->nbins, sizeof(double)))))
        Exit("history", "fget_dist", "data fork calloc error");
      for (i=0; i<dist->nbins; ++i)
      {
        if (fgetl(dist->bin+i, stream)) return 1;
        if (fgetd(dist->data+i, stream)) return 1;
        if (fgetd(dist->cweight+i, stream)) return 1;
      }
    }
  }
  return 0;
}


int fput_all_dist(bridgestruct *c, FILE *stream)
{
  long			i, system;

  for (i=0; i<D_NDIST; ++i) 
    for (system=0; c->dist[i]&&(system<NSYSTEMS); ++system)
      if (fput_dist(c->dist[i]+system, TRUE, stream)) return 1;
  return 0;
}


int fget_all_dist(bridgestruct *c, FILE *stream, long version)
{
  long			i, system;

  if (version<8) return 0;
  ReinitSample(c->dist);
  for (i=0; i<D_NDIST; ++i)
    for (system=0; c->dist[i]&&(system<*(c->n.systems)); ++system)
    {
      D_Reset(c->dist[i]+system);
      if (fget_dist(c->dist[i]+system, NULL, stream, version)) return 1;
    }
  return 0;
}


int fput_all(bridgestruct *c, FILE *stream)
{
  int			flag, i, j;
  
//  if (!(flag = fputl(c->d_type, stream))) {}
  if (!(flag = fput_n(&(c->n), stream)))			// 0
      
    flag		= fput_command(&(c->command), stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 1
    flag		= fput_system(c->system+i, stream);
  for (i=0; (i<*(c->n.types))&&!flag; ++i)			// 2
    flag		= fput_type(c->type+i, stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 3
    for (j=0; (j<*(c->system[i].nmols))&&!flag; ++j)
      flag		= fput_mol(c->mol[i]+j, stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 4
    flag		= fput_av(c->av[i], stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 5
    flag		= fput_v(c->v[i], stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 6
    flag		= fput_w(c->vir[i], stream);
  if (!flag) 							// 7
    flag 		= fput_all_dist(c, stream);
  return flag;
}


int fget_all(bridgestruct *c, FILE *stream, long version)
{
  int			flag = 0, i = 0;
  long			cnt = 0, debug = 0;	// cnt keeps track of where
						// read in fails
  long			data;
  
//  if (version>19) flag = fgetl(c->d_type, stream);
//  if (!flag)
  flag = fget_n(&(c->n), stream, version);			// 0
  if (!flag)
    flag		= fget_command(&(c->command), stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 1
    flag		= fget_system(c->system+i, stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.types))&&!flag; ++i)			// 2
    flag		= fget_type(c->type+i, stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.mols))&&!flag; ++i)			// 3
    flag		= fget_mol(c->mol[0]+i, stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 4
    flag		= fget_av(c->av[i], stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 5
    flag		= fget_v(c->v[i], stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 6
    flag		= fget_w(c->vir[i], stream, version);
  if (!flag) ++cnt;

  if (!flag)							// 7
    flag		= fget_all_dist(c, stream, version);
  if (debug&&flag)
    fprintf(stderr, "history::fget_all: read error at identity #%d, i = %d\n", 
      cnt, i);
 
/* 
  if (!flag) {				// copy some values to corresponding system variable
    for (i=0; i<*(c->n.systems); i++) {
       BOX[i].lx	=	c->system[i].a.x;
       BOX[i].ly	=	c->system[i].b.y;
       BOX[i].lz	=	c->system[i].c.z;
    }
  }
*/
/*
  if (!flag)
  {
    U_SYSTEM		= TRUE;
    D_TYPE		= NTYPES+1;
    SetBoxVectors(c->system->a, c->system->b, c->system->c);
    CalcUnits(-1);
    if (version<16) AllCartesian();
    InitForceField();
    if (version<15)
    {
      for (i=0; i<*(c->n.systems); ++i)
      {
        av[i].v		-= v[i].corr*av[i].total/NActive[i];
        av[i].vnb	-= v[i].corr*av[i].total/NActive[i];
        v[i].total	-= v[i].corr;
        v[i].nonbonded	-= v[i].corr;
      }
      CalcVCorr();
      for (i=0; i<*(c->n.systems); ++i)
      {
        av[i].v		+= v[i].corr*av[i].total/NActive[i];
        av[i].vnb	+= v[i].corr*av[i].total/NActive[i];
        v[i].total	+= v[i].corr;
        v[i].nonbonded	+= v[i].corr;
      }
    }   
  }
*/
  return flag;
}


FILE *fcreate_hist(char *name)
{
  FILE			*fp;
  long			version	= HIST_VERSION, i;
  char			ident[5] = HIST_IDENT;
  
  if (fp = fopen(name, "w"))
  {
    for (i=0; i<4; ++i)
      fputc((int) ident[i], fp);
    fputl(&version, fp);
    return fp;
  }
  return NULL;
}


FILE *fread_hist(char *name, long *version)
{
  char			c, ident[5] = HIST_IDENT;
  long			i, flag = 1;
  FILE			*fp;
  
  if (!name[0])
    fp			= stdin;
  else 
    fp			= fopen(name, "r");
  if (fp)
  {
    for (i=0; (i<4)&&flag; ++i)
    {
      c			= fgetc(fp);
      flag		= (c==ident[i]) ? 1 : 0;
    }
    if (flag)
    {
      fgetl(version, fp);
      if (*version<=HIST_VERSION)
        return fp;
    }
  }
  return NULL;
}


bridgestruct		*VarBridge = NULL;

int StartNewHistoryFile(char *name, long flag_mpi)
{
  FILE			*fp;
  
  if (fp = fcreate_hist(name))
  {
    VarBridge		= BridgeMap(flag_mpi);
    fclose(fp);
    return 1;
  }
  return 0;
}


int H_StoreBridge(char *name, bridgestruct *bridge)
{
  int			flag;
  FILE			*fp;

  if (fp = fopen(name, "a"))
  {
    flag		= fput_all(bridge, fp);
    fclose(fp);
    return flag;
  }
  return 0;
}


int H_StoreCurrent(char *name)
{
  return H_StoreBridge(name, VarBridge);
}


int H_GetHeader(char *name, long *version)
{
  long			flag, i;
  FILE			*fp;

  if (fp = fread_hist(name, version))
  {
    VarBridge		= BridgeMap(TRUE);
    flag		= 1;
    if ((*version>19)&&fgetl(VarBridge->d_type, fp)) flag = 0;
    if (flag&&fget_n(&(VarBridge->n), fp, *version)) flag = 0;
    if (flag&&fget_command(&(VarBridge->command), fp, *version)) flag = 0;
    for (i=0; (i<*(VarBridge->n.systems))&&flag; ++i)
      flag		= !fget_system(VarBridge->system+i, fp, *version);
    fclose(fp);
    return flag ? 0 : 1;
  }
  return 1;
}
  

FILE *H_GetFirst(char *name, long *version, long flag_mpi)
{
  FILE			*fp;
  
  if (fp = fread_hist(name, version))
  {
    VarBridge		= BridgeMap(flag_mpi);
    if (!fget_all(VarBridge, fp, *version))
    {
      VarBridge		= BridgeMap(flag_mpi);
      return fp;
    }
  }
  return NULL;
}


bridgestruct *H_Bridge()		// not used
{
  return VarBridge;
}


int H_GetNext(FILE *fp, long version)
{
  if (VarBridge)
    fget_all(VarBridge, fp, version);
  return feof(fp)!=0;
}

/*
void H_InitSystem(char *argv[])		// not used
{
  InitUnits();
  InitMove();
  //InitEnsembles();
  //InitForceField();
  InitCavity(argv);
  InitSample(argv);
}
*/
                                   src/hstbk.c                                                                                         0000600 0143352 0000144 00000027614 11103062060 012332  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    hst.c
    author:     Peng Yi at MIT
    date:       January 12, 2008
    purpose:    get info from binary history files and do analysis
*/

#define __MAIN_PROGRAM
#include "header.h"

#include "correlation.h"

long	timestep = 0,
	min=MAXNMOLS, 
	max=0,
	nnmax[MAXNMOLS],		// nnmax[i] is the counts of conf. with nmax=i
	n[1000][MAXNMOLS];		// n[i][j] is the counts of nuclei of size j when nmax=i
double	pn[1000][MAXNMOLS],		// pn[i][j] is the prob. of nuclei of size j when nmax=i
	Gn[MAXNMOLS],			// free energy of one nucleus
	G[MAXNMOLS];			// free energy of the system

// some functions declared here
void printout(FILE *);
void hst2conf();
void hst2lammpsdump(FILE *fPtr);

int main(int argc, char * argv[])
{
   FILE		*fhist, *fp, *fp1, *fsize, *fconf, *frdf, *flammps;
   long		addflag=0, splitflag=0, normalflag=0, 
		Nsplit, 
		Nblock = -1;		// # of records in each block, could be set by command line
   long		version, system;
   char		file_hst1[256],
		comments[1024];
   double	dummy;
   long		i, j, k; 
   molstruct	*moli;
   double	var, bias; 
   static long	init = 1;

   if (argc<3) {
      printf("hst (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\thst -normal history_file [N]\n\n");
      printf("\thst -add history_file1 history_file2\n\n");
      printf("\thst -split history_file N\n\n");
      printf("Notes:\n");
      printf("\t* require setup file to get e.g. Temp for analysis\n\n");
      printf("\t* -normal: normal history file analysis\n");
      printf("\t\t if N is given, then do analysis in blocks each with N records, otherwise as one big block\n\n");
      printf("\t* -add: add history_file2 to the end of history_file1\n\n");
      printf("\t* -split: split history_file to binsplit1 and binsplit2 and the former contains N timesteps\n\n");
      exit(1);
   }

   if (!strcmp(argv[1], "-add")) {
      addflag	=	1;
      sprintf(file_hst1, argv[2]);
      sprintf(file_hst, argv[3]);
   }
   else if (!strcmp(argv[1], "-split")) {
      splitflag	=	1;
      sprintf(file_hst, argv[2]);
      Nsplit		=	atol(argv[3]);
   }
   else if (!strcmp(argv[1], "-normal")) {
      normalflag	=	1;
      sprintf(file_hst, argv[2]);
      if (argc == 4)					// optional
	 Nblock		=	atol(argv[3]);
   }

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate max. memory for molecules
   GetSetup(argv);			// most of the variables read from setup 
					// will be updated by the reading from history file
   InitUnits();				// so far the binary file is in system units (4/26/08)	
   InitForcefield();

   if (normalflag) {
      if (!(fsize = fopen("Psize","w"))) 		// open output file for size distribution information
         Exit("hst", "main","Size prob. distribution file failed to open");
   }
   else if (addflag) {
      if (!(fp1=fopen(file_hst1, "a")))
         Exit("hst", "main", "add target history file fails to open");
   }

   if (!(fdump=fopen("dump", "w"))) 			// to convert binary hist file to dump file through Write_Conf
      Exit("hst", "main", "dump file failed to open");
   if (!(flammps=fopen("lammps.dump", "w")))		// to convert to lammps dump file for lammps2pdb
      Exit("hst", "main", "lammps dump file failed to open");


   if (!(fp = H_GetFirst(file_hst, &version, FALSE))) {		// open history file, get some info. stored in bin file
      printf("File format not recognized.\n\n");
      exit(-1);
   }
   InitSample(argv);		// InitSample will need some commands read in from the first record
   fclose(fp);			// Or we can move the pointer to the beginning of the history file

   timestep	=	0;
   if (!(fp = H_GetFirst(file_hst, &version, FALSE))) {	// reopen history file
      printf("File format not recognized.\n\n");
      exit(-1);
   }

   //while (!H_GetNext(fp, version)) {	
   do {					// do analysis and read in next record/configuration

      // some pre-processing of data for analysis
      for (i=0; i<NSYSTEMS; i++) {
         NMols[i]	=	0;
         NSites[i]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++)
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;				// total # of mols in certain system
            NSites[system]	+=	moli->nsites;		// total # of sites in certain system
         }
      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox	=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol	=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].pres 	=	P;
         BOX[i].temp 	=	kT;
         BOX[i].drmax 	=	DRMAX;
         BOX[i].dlmax 	=	DLMAX;
         BOX[i].damax	=	DAMAX;
         BOX[i].rc	=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb	=	Rb;
         BOX[i].rv	=	Rv;
      } 
      for (i=0; i<NMOLS; i++) { 				// note NMOLS might change due to MPI
         for (j=0; j<mol[i].nsites; j++)  {
            mol[i].flags[j]	=	1;		// activate all the sites on this processor
            mol[i].parent[j]	=	j-1;		// initialize parent sites
         }
         mol[i].flip		=	0;		// flip to the original direction
         mol[i].origin		=	CenterofMass(mol+i);
      }
#ifdef CELL_LIST
      if (init) {
         CL_Init();				// Cell list should be ready before Verlet list
         init	=	0;
      }
      else
         CL_Destroy();
      CL_Build();
#endif	/* CELL_LIST */

      CalcV();					// need cell list
      if (dynvar==3)
         SampleP2All();				// calculate global P2 and local p2
      SampleSpherical();
      Find_Nuclei();

printf("%d %f %f %d\n", timestep, v[0].tot/NSITES, NSITES/BOX[0].vol, nmax[0][0]);

      // calculate bias for each record/configuration 
      if (dynvar==1)
         var	=	MAXSIZE[0];
      else if (dynvar==3)
	 var	=	P2M[0];
      bias	=	exp(0.5*kP2*(var-P2middle)*(var-P2middle));

      // collect probability distribution of cluster size 
      nnmax[MAXSIZE[0]]	++;
      for (i=1; i<=MAXSIZE[0]; i++) {
	 pn[0][i]	+=	bias * sizedist[i];		// probability distribution of cluster size
	 n[0][i]	+=	sizedist[i];

         pn[MAXSIZE[0]][i]	+=	bias * sizedist[i];
         n[MAXSIZE[0]][i]	+=	sizedist[i];
      }
      if (MAXSIZE[0] > max)
	 max	=	MAXSIZE[0];
      if (MAXSIZE[0] < min) 
         min	=	MAXSIZE[0];

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
      } 
      if (!((timestep*ITAPE)%IRADIAL))
         radial("sample");		// sample rdf

      correlation();			// calculate correlation function

//      if (normalflag && MAXSIZE[0]==17 && secondNmax[0]<=6) 
//      if (normalflag && MAXSIZE[0] >=15 && MAXSIZE[0] <=17 )
	// hst2conf(fconf);		// convert to configuration file
      
      hst2lammpsdump(flammps);		// convert to lammps dump format to use lammps2pdb in the future
      Write_Conf(timestep);		// convert to ASCII dump file
      timestep	+=	ITAPE;		// because the binary hist file was stored every ITAPE cycles

      if (Nblock!=-1 && timestep%Nblock==0) {		// summarize the results for this block and output
         for (i=1; i<=max; i++) 
            Gn[i]	=	-log(pn[0][i]);
         for (i=min; i<=max; i++) {
            for (j=1; j<=i; j++) {
               G[i]	+=	n[i][j] * Gn[j];
            }
         }
	 printout(fsize);

         // clear up the counters
         for (i=0; i<MAXNMOLS; i++) {
	    nnmax[i]	=	0;
	    Gn[i]	=	0.0;
	    G[i]	=	0.0;
	 }
         for (i=0; i<100; i++) { 
            for (j=0; j<MAXNMOLS; j++) {
               pn[i][j]	=	0.0;
	       n[i][j]	=	0;
            }
         }
         max	=	0;
	 min	=	MAXNMOLS;
      }
   } while (!H_GetNext(fp, version));		// read in next record/configuration

   /* normalize correlation functions and print out */

   corr_norm();
   corr_print();
   radial("print");

   // summarize results for all records and print out size distribution

   if (Nblock==-1) {
      for (i=1; i<=max; i++) 
         Gn[i]	=	-log(pn[0][i]);			// free energy of nucleus
      for (i=min; i<=max; i++) {
         for (j=1; j<=i; j++) {
            G[i]	+=	n[i][j] * Gn[j];	// free energy of the system
         }
      }
      printout(fsize);
   }

   fclose(fdump);
   fclose(fp);

   if (normalflag) {
      fclose(fsize);
   }
   return EXIT_SUCCESS;
}


void printout(FILE *fPtr)
{
   long		i, j;

   fprintf(fPtr, "\nOverall nuclei size distribution P(n):\n\n");
   fprintf(fPtr, "n    N(n)   fixed p(n)   -log(fixed pn[i])\n\n");
   for (i=1; i<=max; i++) 
      fprintf(fPtr, "%4d   %8d   %14.4f   %f\n", i, n[0][i], pn[0][i], Gn[i]);

   fprintf(fPtr, "P(n) for copy purpose...\n");
   for (i=1; i<=max; i++)
      fprintf(fPtr, "%f\n", Gn[i]);

   fprintf(fPtr, "\nDistribution of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, "n_max    N(n_max)\n");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%4d  %8d\n", i, nnmax[i]);

   fprintf(fPtr, "\nSize distribution N(n) for different value of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, " n_max");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%8d  ", i);
   fprintf(fPtr, "\n   n \n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%8d  ", n[j][i]);
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n  Normalized to N(n_max)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", (nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n  -Logrithm of N(n)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", -log(nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
}


void hst2conf()
{
   long		i;
   molstruct	*moli;
   static FILE	*fPtr;
   static long	init =1;

   vector	rA, rB, rAB;
   double	r2, L;
   long		id, print=0;

   if (init) {
      if (!(fPtr=fopen("hst.conf", "w")) )
          Exit("hst", "main", "open conf file failed.");
      init	=	0;
   }
/*
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (sizeofnucl[moli->nuclid[0]]==MAXSIZE[0]) {
         id	=	moli->nuclid[0];
	 rA	=	CenterofMass(moli);
         break;
      }
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (moli->nuclid[0] == id) {
         rB	=	CenterofMass(moli);
	 rAB	=	V_Subtr(&rA, &rB);
         r2	=	V_Dot(&rAB, &rAB);
         L	=	MIN(BOX[0].lx, MIN(BOX[0].ly, BOX[0].lz));
         if (r2 > L * L/2) {
            print	=	1;
	    break;
         }
      }

   if (!print)
      return;
*/
   fprintf(fPtr, "!TIMESTEP %d Nmax = %d 2ndNmax = %d\n", timestep, MAXSIZE[0], secondNmax[0]);
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);

   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx*unit.LENGTH, BOX[i].ly*unit.LENGTH, BOX[i].lz*unit.LENGTH);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
   }
}

void hst2lammpsdump(FILE *fPtr)
{
   long		i;
   molstruct	*moli;

   fprintf(fPtr, "ITEM: TIMESTEP\n");	
   fprintf(fPtr, "%d\n", timestep);		
   fprintf(fPtr, "ITEM: NUMBER OF ATOMS\n");
   fprintf(fPtr, "%d\n", NSITES);
   fprintf(fPtr, "ITEM: BOX BOUNDS\n");
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lx*unit.LENGTH, 0.5*BOX[0].lx*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].ly*unit.LENGTH, 0.5*BOX[0].ly*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lz*unit.LENGTH, 0.5*BOX[0].lz*unit.LENGTH); 
   fprintf(fPtr, "ITEM: ATOMS\n");
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         fprintf(fPtr, "%d %d %f %f %f %f %f %f %d %d %d\n", (moli-mol)*moli->nsites+i, moli->type[i], 
              moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH, 0, 0, 0, 0, 0, 0);
   return;
}

                                                                                                                    src/hst.c                                                                                           0000600 0143352 0000144 00000042063 11450206644 012025  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    hst.c
    author:     Peng Yi at MIT
    date:       January 12, 2008
    purpose:    get info from binary history files and do analysis
    note:	April 15, 2010 add biggest nucleus shape analysis
*/

#define __MAIN_PROGRAM
#include "header.h"

#include "correlation.h"

long	nrecords,			// # of records
	min=MAXNMOLS, 
	max=0,
	nnmax[MAXNMOLS],		// nnmax[i] is the counts of conf. with nmax=i
	n[1000][MAXNMOLS];		// n[i][j] is the counts of nuclei of size j when nmax=i
double	pn[1000][MAXNMOLS],		// pn[i][j] is the prob. of nuclei of size j when nmax=i
	Gn[MAXNMOLS],			// free energy of one nucleus
	G[MAXNMOLS];			// free energy of the system
long	lencnt[MAXNSYSTEMS][MAXNMOLSITES+1];	// number of chains of different length
long	BIN;				// binsize for nucleus size statistics

// some functions declared here
void printout(FILE *);			// print out probability distribution
void hst2conf();			// convert to configuration file
void hst2lammpsdump(FILE *fPtr);	// convert to dump file in lammps format

int main(int argc, char * argv[])
{
   FILE		*fhist, *fp, *fp1, *fsize, *fconf, *flammps;
   char		file_hst1[256],
		comments[1024];
   long		addflag=0, splitflag=0, normalflag=0, 
		Nsplit, 
		Nblock = -1;		// # of records in each block, could be set by command line
   long		version, system;
   long		i, j, k, m, nuclid; 
   double	dummy;
   double	var, bias; 		// bias variable
   static long	init = 1;
   molstruct	*moli;
   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];	// contain the beads in the biggest nucleus

   if (argc<3) {
      printf("hst (c) 2010 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\thst -normal history_file BINSIZE [N]\n\n");
      printf("\thst -add history_file1 history_file2\n\n");
      printf("\thst -split history_file N\n\n");
      printf("Notes:\n");
      printf("\t* require setup file to get information e.g. Temp for analysis\n\n");
      printf("\t* -normal: normal history file analysis\n");
      printf("\t\t BINSIZE is the binsize for nucleus increment\n\n");
      printf("\t\t if N is given, then do analysis in blocks each with N records, otherwise as one big block\n\n");
      printf("\t* -add: add history_file2 to the end of history_file1\n\n");
      printf("\t* -split: split history_file to binsplit1 and binsplit2 and the former contains N records\n\n");
      exit(1);
   }

printf("start\n");
   if (!strcmp(argv[1], "-add")) {
      addflag	=	1;
      sprintf(file_hst1, argv[2]);
      sprintf(file_hst, argv[3]);
   }
   else if (!strcmp(argv[1], "-split")) {
      splitflag	=	1;
      sprintf(file_hst, argv[2]);
      Nsplit		=	atol(argv[3]);
   }
   else if (!strcmp(argv[1], "-normal")) {
      normalflag	=	1;
      sprintf(file_hst, argv[2]);
      BIN	=	atol(argv[3]);
      if (argc == 5)					// optional
	 Nblock		=	atol(argv[4]);
   }

   InitMols(MAXNMOLS, MAXNMOLS);// allocate max. memory for molecules
   GetSetup(argv);		// most of the variables read from setup file
				// will be updated by the reading from history file
   InitUnits();			// calculate the transformation b/w SI units and system units	
   InitForcefield();		// initialize forcefield
   InitSample(argv);		// initialize sampling

   if (normalflag) {		// open output file for distribution output
      if (!(fsize = fopen("Psize","w"))) 
         Exit("hst", "main","Size prob. distribution file failed to open");
   }
   else if (addflag) {
      if (!(fp1=fopen(file_hst1, "a")))
         Exit("hst", "main", "add target history file fails to open");
   }

/*
   if (!(fdump=fopen("dump", "w"))) 		// to convert binary hist file to dump file through Write_Conf
      Exit("hst", "main", "dump file failed to open");
*/

   if (!(flammps=fopen("lammps.dump", "w")))	// to convert to lammps dump file for lammps2pdb
      Exit("hst", "main", "lammps dump file failed to open");

   if (!(fhist=fopen(file_hst, "r")))
      Exit("hst", "main", "file_hst failed to open.");

   Read_MultiConf(fhist);			// read the first conf., no analysis for this one

   nrecords	=	0;

   while (Read_MultiConf(fhist)) {		// read next conf., start analysis

      printf("read in %d configuration...\n", nrecords);

      /////////////////////////////////////////
      /* PRE-PROCESSING of data for analysis */
      /////////////////////////////////////////
 
      CoorSI2System();				// convert from SI to system units

      for (i=0; i<NSYSTEMS; i++) {
         NMols[i]	=	0;
         NSites[i]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++)
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;			// total # of mols in certain system
            NSites[system]	+=	moli->nsites;	// total # of sites in certain system
            lencnt[system][moli->nsites]	++;	// polydispersity
         }
      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox	=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol	=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].pres 	=	P;
         BOX[i].temp 	=	kT;
         BOX[i].drmax 	=	DRMAX;
         BOX[i].dlmax 	=	DLMAX;
         BOX[i].damax	=	DAMAX;
         BOX[i].rc	=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb	=	Rb;
         BOX[i].rv	=	Rv;
      } 
      for (i=0; i<NMOLS; i++) { 			// note NMOLS might change due to MPI
         for (j=0; j<mol[i].nsites; j++)  {
            mol[i].flags[j]	=	1;		// activate all the sites on this processor
            mol[i].parent[j]	=	j-1;		// initialize parent sites
         }
         mol[i].flip		=	0;		// flip to the original direction
         mol[i].origin		=	CenterofMass(mol+i);
      }
#ifdef CELL_LIST				// cell list needed for energy calc.
      if (init) {
         CL_Init();				// Cell list should be ready before Verlet list
         init	=	0;
      }
      else
         CL_Destroy();
      CL_Build();
#endif	/* CELL_LIST */

      ///////////////////////
      /*  Perform Analysis */
      ///////////////////////

      printf("start performing analysis...\n");

      CalcV();					// need cell list
      if (V_VIRIAL)	Sample_Pressure();
      SampleP2All();				// calculate global P2 and local p2
      Dist_p2();				// collect distribution of local p2
      SampleSpherical();			// sample spherical coordinates
      Dist_Spherical();				// collect distribution

      //////////////////////////////////////////////////
      /* Calculate Bias for each record/configuration */
      //////////////////////////////////////////////////

      printf("calculate bias ...\n");

      Find_Nuclei_p2(1);

//      Find_Nuclei(1);				// sample nuclei for bias calculation
//printf("%d %d %d ", nrecords, nmax[0][0], nmax[0][1]);

      if (dynvar==1)
         var	=	MAXSIZE[0];
      else if (dynvar==3)
	 var	=	P2M[0];
      else if (dynvar==6)
	 var	=	MAXSIZE[0];
    
      bias	=	exp(0.5*kP2*(var-P2middle)*(var-P2middle));

      ////////////////////////////////////////////
      /* Find the nuclid of the biggest nucleus */
      ////////////////////////////////////////////

      sizeofnucl	=	sizeofnuclp2;	// Find_Nuclei_p2() uses sizeofnucleip2
      sizedist		=	sizedistp2;

      system	=	0;		// for now only one system, 4/16/2010
      i	= 1;				// nucleus id starts from 1
      while (sizeofnucl[i] != nmax[system][0]) {	// find the biggest nucleus
         i ++;
      }
      nuclid	=	i;		// store the id of the biggest nucleus
    
      ////////////////////////////////////////////////////////////////////////////
      /* Group the segments belong to the biggest nucleus based on sizeofnucl[] */
      ////////////////////////////////////////////////////////////////////////////

      m	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if (system == moli->box) {
            for (i=0; i<moli->nsites; i++) {
               if (moli->nuclid[i] == nuclid) {
                  nucleus[m].moli	=	moli;
	          nucleus[m].site	=	i;
                  m	++;
	       }
            }
         }
      }

      ////////////////////////////////
      /* Nucleus structure analysis */
      ////////////////////////////////

//      cylindershape(nucleus, nmax[system][0], nuclid); 	// calc. length and radius
							// assuming cylindrical shape

      //////////////////////////////////////////////////////
      /* Collect probability distribution of cluster size */
      //////////////////////////////////////////////////////

      printf("collect histogram ...\n");

      //Find_Nuclei_p2(1);			// Needed if the nucleus definition used
						// in biasing is different from the nucleus
						// definition used to calculate the nucleus 
						// size distribution

      //printf("%d %d\n ", nmax[0][0], nmax[0][1]);

      nnmax[MAXSIZE[0]]	++;
      for (i=1; i<=MAXSIZE[0]; i++) {		// index starts from 1
         if (BIN==1) {
	   pn[0][i] += bias * sizedistp2[i];	// probability distribution of cluster size
	   n[0][i]  += sizedistp2[i];
	 }
         else {
           pn[0][i/BIN+1] += bias * sizedistp2[i];
	   n[0][i/BIN+1]  += sizedistp2[i];
         }
         pn[MAXSIZE[0]][i] += bias * sizedistp2[i];
         n[MAXSIZE[0]][i]  += sizedistp2[i];
      }
      if (MAXSIZE[0] > max)
	 max	=	MAXSIZE[0];
      if (MAXSIZE[0] < min) 
         min	=	MAXSIZE[0];

      ////////////////////////////
      /* Calculate correlations */
      ////////////////////////////

      printf("calculate correlation ...\n");

      correlation();			// calculate correlation function
      if (!((nrecords*ITAPE)%IRADIAL))
         radial("sample");		// sample radial distribution function

      ///////////////////////////////////////
      /* Output to other trajactory format */
      ///////////////////////////////////////

//      hst2lammpsdump(flammps);	// convert to lammps dump format to use lammps2pdb 

//      Write_Conf(TIMESTEP);		// convert to ASCII dump file
//      if (normalflag && MAXSIZE[0]==17 && secondNmax[0]<=6) 
//      if (normalflag && MAXSIZE[0] >=15 && MAXSIZE[0] <=17 )
	// hst2conf(fconf);		// convert to configuration file    

      /////////////////////////////////////////////////////
      /* Summarize the results for this block and output */
      /////////////////////////////////////////////////////

      if (Nblock!=-1 && TIMESTEP%Nblock==0) {
         for (i=1; i<=max; i++) 
            Gn[i]	=	-log(pn[0][i]);
         for (i=min; i<=max; i++) {
            for (j=1; j<=i; j++) {
               G[i]	+=	n[i][j] * Gn[j];
            }
         }
	 printout(fsize);

         // clear up the counters

         for (i=0; i<MAXNMOLS; i++) {
	    nnmax[i]	=	0;
	    Gn[i]	=	0.0;
	    G[i]	=	0.0;
	 }
         for (i=0; i<100; i++) { 
            for (j=0; j<MAXNMOLS; j++) {
               pn[i][j]	=	0.0;
	       n[i][j]	=	0;
            }
         }
         max	=	0;
	 min	=	MAXNMOLS;
      }
      nrecords	++;
   }

   /////////////////////////////////////////////////
   /* Print out sampling results and correlations */
   /////////////////////////////////////////////////

   printf("print out sampling results and correlations...\n");

   S_PrintAll();
   corr_norm();
   corr_print();
   radial("print");

   ///////////////////////////////////////////////////////////////////////
   /* Summarize results for all records and print out size distribution */
   ///////////////////////////////////////////////////////////////////////
 
   if (Nblock==-1) {
      for (i=1; i<=max; i++) 
         Gn[i]	=	-log(pn[0][i]);			// free energy of nuclei
      for (i=min; i<=max; i++) {
         for (j=1; j<=i; j++) {
            G[i]	+=	n[i][j] * Gn[j];	// free energy of system
         }
      }
      printout(fsize);
   }
   if (normalflag) {
      fclose(fsize);
   }
   fclose(fhist);
   return EXIT_SUCCESS;
}


void printout(FILE *fPtr)		// output main information and distribution
{
   long		i, j;

   fprintf(fPtr, "ITAPE = %d (dump file written every ITAPE MC cycles)\n", ITAPE); 
   fprintf(fPtr, "Number of records: %d\n\n", nrecords);
/*
   fprintf(fPtr, "bias nucleus definition = ");
   switch (clusdef) {
      case	1:	fprintf(fPtr, "%s\n", );	break;
      case	1:	fprintf(fPtr, "%s\n", );	break;
      case	1:	fprintf(fPtr, "%s\n", );	break;
      default:		break;
   }
   fprintf(fPtr, "kP2 = %f\t P2middle = %f\n",kP2, P2middle);
   fprintf(fPtr, "distribution nucleus definition = ");
   switch (dynvar) {
      case	1:	fprintf(fPtr, "%s\n", );	break;
      case	1:	fprintf(fPtr, "%s\n", );	break;
      case	1:	fprintf(fPtr, "%s\n", );	break;
      default:		break;
   }
*/  

   /* Output Polydispersity */

   fprintf(fPtr, "\n#Polydispersity\n\n");
   for (i=0; i<NSYSTEMS; i++) {
      fprintf(fPtr, "length      number of chains (SYSTEM #%d)\n\n", i);
         for (j=0; j<=MAXNMOLSITES; j++)
            if (lencnt[i][j]!=0)
               fprintf(fPtr, "%d    %d\n", j, lencnt[i][j]/nrecords);
   }

   /* Output nuclei size distribution */

   fprintf(fPtr, "\n#Overall nuclei size distribution P(n):\n\n");
   fprintf(fPtr, "n    N(n)   fixed p(n)   -log(fixed pn[i])\n\n");
   for (i=1; i<=max/BIN; i++) 
      fprintf(fPtr, "%4d   %8d   %14.4f   %f\n", i*BIN, n[0][i], pn[0][i], Gn[i]-Gn[1]);

   fprintf(fPtr, "\n\n#Throw away poorly sampled nucleus size\n\n");
   fprintf(fPtr, "n    N(n)   fixed p(n)   -log(fixed pn[i])\n\n");
   for (i=1; i<=max/BIN; i++) {
      if (n[0][i] <10) {			// poorly sampled nucleus size
         fprintf(fPtr, "%4d   %8d   nan      nan\n", i*BIN, n[0][i]);
      }
      else {
         fprintf(fPtr, "%4d   %8d   %14.4f   %f\n", i*BIN, n[0][i], pn[0][i], Gn[i]-Gn[1]);
      }
   } 

   fprintf(fPtr, "\n#Nucleus Free energy G(n) for copy purpose...\n");
   for (i=1; i<=max/BIN; i++) {
      if (n[0][i]<10)
         fprintf(fPtr, "nan\n");		// poorly sampled nucleus size
      else
         fprintf(fPtr, "%f\n", Gn[i]-Gn[1]);
   }

   fprintf(fPtr, "\n#Distribution of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, "n_max    N(n_max)\n");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%4d  %8d\n", i, nnmax[i]);

   fprintf(fPtr, "\n#Size distribution N(n) for different value of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, " n_max");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%8d  ", i);
   fprintf(fPtr, "\n   n \n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%8d  ", n[j][i]);
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n#Normalized to N(n_max)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", (nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n#-Logrithm of N(n)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", -log(nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
}


void hst2conf()
{
   long		i;
   molstruct	*moli;
   static FILE	*fPtr;
   static long	init =1;

   vector	rA, rB, rAB;
   double	r2, L;
   long		id, print=0;

   if (init) {
      if (!(fPtr=fopen("hst.conf", "w")) )
          Exit("hst", "main", "open conf file failed.");
      init	=	0;
   }
/*
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (sizeofnucl[moli->nuclid[0]]==MAXSIZE[0]) {
         id	=	moli->nuclid[0];
	 rA	=	CenterofMass(moli);
         break;
      }
   for (moli=mol; moli<mol+NMOLS; moli++)
      if (moli->nuclid[0] == id) {
         rB	=	CenterofMass(moli);
	 rAB	=	V_Subtr(&rA, &rB);
         r2	=	V_Dot(&rAB, &rAB);
         L	=	MIN(BOX[0].lx, MIN(BOX[0].ly, BOX[0].lz));
         if (r2 > L * L/2) {
            print	=	1;
	    break;
         }
      }

   if (!print)
      return;
*/
   fprintf(fPtr, "!TIMESTEP %d Nmax = %d 2ndNmax = %d\n", TIMESTEP, MAXSIZE[0], secondNmax[0]);
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);

   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx*unit.LENGTH, BOX[i].ly*unit.LENGTH, BOX[i].lz*unit.LENGTH);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
   }
}

void hst2lammpsdump(FILE *fPtr)
{
   long		i, accum;
   molstruct	*moli;

   fprintf(fPtr, "ITEM: TIMESTEP\n");	
   fprintf(fPtr, "%d\n", TIMESTEP);		
   fprintf(fPtr, "ITEM: NUMBER OF ATOMS\n");
   fprintf(fPtr, "%d\n", NSITES);
   fprintf(fPtr, "ITEM: BOX BOUNDS\n");
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lx*unit.LENGTH, 0.5*BOX[0].lx*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].ly*unit.LENGTH, 0.5*BOX[0].ly*unit.LENGTH); 
   fprintf(fPtr, "%f %f\n", -0.5*BOX[0].lz*unit.LENGTH, 0.5*BOX[0].lz*unit.LENGTH); 
   fprintf(fPtr, "ITEM: ATOMS\n");
   
   accum	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         fprintf(fPtr, "%d %d %d ", accum+i+1, moli-mol+1, moli->type[i]); 
         fprintf(fPtr, "%f %f %f ", moli->p[i].x * unit.LENGTH, moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
         fprintf(fPtr, "%f %f %f %d %d %d\n", 0.0, 0.0, 0.0, 0, 0, 0);
      }
      accum	+=	moli->nsites;
   }
   return;
}

 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             src/hstcomb.c                                                                                       0000600 0143352 0000144 00000017604 11006366521 012667  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    hstcomb.c
    author:     Peng Yi at MIT
    date:       May 1, 2008
    purpose:    combine two binary history files
*/

#define __MAIN_PROGRAM
#include "header.h"

#include "correlation.h"

long	timestep = 0,
	min=MAXNMOLS, max=0,
	nnmax[MAXNMOLS],		// nnmax[i] is the counts of conf. with nmax=i
	n[100][MAXNMOLS];		// n[i][j] is the counts of nuclei of size j when nmax=i
double	pn[100][MAXNMOLS],		// pn[i][j] is the prob. of nuclei of size j when nmax=i
	Gn[MAXNMOLS],			// free energy of one nucleus
	G[MAXNMOLS];			// free energy of the system

// some functions declared here
void printout(FILE *);
void hst2conf(FILE *fPtr);
void hst2dump(FILE *fPtr);

int main(int argc, char * argv[])
{
   FILE		*fp, *fp1, *fsize, *fconf, *fcomb;
   long		addflag=0, splitflag=0, normalflag=0, 
		Nsplit, 
		Nblock = -1;		// # of records in each block, could be set by command line
   long		version;
   char		file_hst1[256], file_hst2[256],
		comments[1024];
   char		data;
   double	dummy;
   long		i, j, k; 
   molstruct	*moli;
   double	var, bias; 

   if (argc<2) {
      printf("hstcombine (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\thstcombine history_file1 history_file2\n\n");
      printf("\t\tcombine two history files and generate a combined one binhstcomb\n\n");
      exit(1);
   }

   sprintf(file_hst1, argv[1]);
   sprintf(file_hst2, argv[2]);
   sprintf(file_hst, "binhstcomb");

/*
   InitMols(MAXNMOLS, MAXNMOLS);	// allocate max. memory for molecules
   GetSetup(argv);			// most of the variables read from setup 
					// will be updated by the reading from history file
   CalcUnits(ConvertUnits);		// so far the binary file is in system units (4/26/08)
*/

   if (!(fp = fread_hist(file_hst1, &version)))		// open the first history file
      Exit("hstcombine", "main", "first history file fails to open.");

   if (!(fcomb=fcreate_hist(file_hst)))			// open the combo file, must be after one history file
      Exit("hstcombine", "main", "history combo file fails to open.");

   data = fgetc(fp);
   while(!feof(fp)) {					// copy 1st history file to combo file
      fputc((int)data, fcomb);
      data = fgetc(fp);
   }
   fclose(fp);

   if (!(fp = fread_hist(file_hst2, &version)))		// open the second history file
      Exit("hstcombine", "main", "second history file fails to open.");
   
   data = fgetc(fp);
   while(!feof(fp)) {					// copy 2nd history file to combo file
      fputc((int)data, fcomb);
      data = fgetc(fp);
   }
   fclose(fp);
   fclose(fcomb);


/*
   if (!(fp = H_GetFirst(file_hst, &version, FALSE))) {
      printf("File format not recognized.\n\n");
   }
   InitSample(argv);		// InitSample will need some commands read in from the first record
   fclose(fp);			// Or we can move the pointer to the beginning of the history file

   if (!(fp = H_GetFirst(file_hst, &version, FALSE))) {	// reopen history file
      printf("File format not recognized.\n\n");
   }
   while (!H_GetNext(fp, version)) {			// read in next record/configuration

      if (dynvar==3)
         SampleP2All();					// calculate global P2 and local p2
      SampleSpherical();
      Find_Nuclei();

      // calculate bias for each record/configuration
      if (dynvar==1)
         var	=	MAXSIZE[0];
      else if (dynvar==3)
	 var	=	P2M[0];

      bias	=	exp(0.5*kP2*(var-P2middle)*(var-P2middle));

      // collect probability distribution of cluster size
 
      nnmax[MAXSIZE[0]]	++;
      for (i=1; i<=MAXSIZE[0]; i++) {
	 pn[0][i]	+=	bias * sizedist[i];		// probability distribution of cluster size
	 n[0][i]	+=	sizedist[i];

         pn[MAXSIZE[0]][i]	+=	bias * sizedist[i];
         n[MAXSIZE[0]][i]	+=	sizedist[i];
      }
      if (MAXSIZE[0] > max)
	 max	=	MAXSIZE[0];
      if (MAXSIZE[0] < min) 
         min	=	MAXSIZE[0];

      correlation();			// calculate correlation function

      // if (normalflag) hst2conf(fconf);		// convert to configuration file
      // hst2dump();			// convert to lammps dump format to use lammps2pdb in the future

      timestep	++;

      if (Nblock!=-1 && timestep%Nblock==0) {		// summarize the results for this block and output
         for (i=1; i<=max; i++) 
            Gn[i]	=	-log(pn[0][i]);
         for (i=min; i<=max; i++) {
            for (j=1; j<=i; j++) {
               G[i]	+=	n[i][j] * Gn[j];
            }
         }

	 printout(fsize);

         // clear up the counters
         for (i=0; i<MAXNMOLS; i++) {
	    nnmax[i]	=	0;
	    Gn[i]	=	0.0;
	    G[i]	=	0.0;
	 }
         for (i=0; i<100; i++) { 
            for (j=0; j<MAXNMOLS; j++) {
               pn[i][j]	=	0.0;
	       n[i][j]	=	0;
            }
         }
         max	=	0;
	 min	=	MAXNMOLS;
      }
   }

   // normalize correlation functions and print out

   corr_norm();
   corr_print();

   // summarize results for all records and print out size distribution

   if (Nblock==-1) {
      for (i=1; i<=max; i++) 
         Gn[i]	=	-log(pn[0][i]);			// free energy of nucleus
      for (i=min; i<=max; i++) {
         for (j=1; j<=i; j++) {
            G[i]	+=	n[i][j] * Gn[j];	// free energy of the system
         }
      }
      printout(fsize);
   }

   fclose(fp);

   if (normalflag) {
      fclose(fsize);
      fclose(fconf);
   }
*/
   return EXIT_SUCCESS;
}

/*
void printout(FILE *fPtr)
{
   long		i, j;

   fprintf(fPtr, "\nOverall nuclei size distribution\n\n");
   fprintf(fPtr, "n    N(n)   fixed p(n)   -log(fixed pn[i])\n\n");
   for (i=1; i<=max; i++) 
      fprintf(fPtr, "%4d   %8d   %14.4f   %f\n", i, n[0][i], pn[0][i], Gn[i]);

   fprintf(fPtr, "\nDistribution of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, "n_max    N(n_max)\n");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%4d  %8d\n", i, nnmax[i]);

   fprintf(fPtr, "\nSize distribution N(n) for different value of n_max in the range of [%-3d, %3d]\n\n", min, max);
   fprintf(fPtr, " n_max");
   for (i=min; i<=max; i++)
      fprintf(fPtr, "%8d  ", i);
   fprintf(fPtr, "\n   n \n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%8d  ", n[j][i]);
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n  Normalized to N(n_max)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", (nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
   fprintf(fPtr, "\n  -Logrithm of N(n)\n\n");
   for (i=1; i<=max; i++) {				// list thru n
      fprintf(fPtr, "%4d  ", i);
      for (j=min; j<=max; j++)				// list thru n_max
         fprintf(fPtr, "%10.3f  ", -log(nnmax[j]>0 ? (double)n[j][i]/nnmax[j] : 0));
      fprintf(fPtr, "\n");
   }
}


void hst2conf(FILE *fPtr)
{
   long		i;
   molstruct	*moli;

   fprintf(fPtr, "!TIMESTEP %d\n", timestep);
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);

   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx*unit.LENGTH, BOX[i].ly*unit.LENGTH, BOX[i].lz*unit.LENGTH);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
			moli->p[i].y * unit.LENGTH, moli->p[i].z * unit.LENGTH);
   }
}


void hst2dump(FILE *fPtr)
{
   long		i;
   molstruct	*moli;

   printf("ITEM: TIMESTEP\n");	
   printf("%d\n", timestep);		
   printf("ITEM: NUMBER OF ATOMS\n");
   printf("%d\n", NSITES);
   printf("ITEM: BOX BOUNDS\n");
   printf("%f %f\n", -BOX[0].lx, BOX[0].lx); 
   printf("%f %f\n", -BOX[0].ly, BOX[0].ly); 
   printf("%f %f\n", -BOX[0].lz, BOX[0].lz); 
   printf("ITEM: ATOMS\n");
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         printf("%d %d %f %f %f %f %f %f %d %d %d\n", (moli-mol)*moli->nsites+i, moli->type[i], 
              moli->p[i].x, moli->p[i].y, moli->p[i].z, 0, 0, 0, 0, 0, 0);
}
*/
                                                                                                                            src/init.c                                                                                          0000600 0143352 0000144 00000015347 11535543255 012205  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    init.c
    author:     Peng Yi at MIT
    date:       October 23, 2006
    purpose:    Suggested initialization sequence
*/

#define __INIT_MODULE
#include "init.h"


void InitMols(long nmols, long noldmols)	// allocate memory for molecules
{
   if (nmols > 0)
      mol		=	(molstruct *) calloc(nmols, sizeof(molstruct));
   if (noldmols > 0)
      oldmol	=	(molstruct *) calloc(noldmols, sizeof(molstruct));
   if ( (nmols>0&&mol==NULL) || (noldmols>0&&oldmol==NULL) )
      Exit("init", "InitMols", "Out of memory.");
}


void GetCoordinates(char * filename)	// read coordinates from initial file
{
   long		i, j, id, system;
   long		nsystems, nmols, nsites;
   molstruct	*moli;
  
   Read_Conf(filename);			// read in WITHOUT unit conversion

   for (i=0; i<NMOLS; i++)		// convert coordinates to system units
      for (j=0; j<mol[i].nsites; j++)
         mol[i].p[j]		=	V_Mult(1.0/unit.LENGTH, mol[i].p+j);

   for (i=0; i<NSYSTEMS; i++) {		// convert box size to system units
      BOX[i].lx	/=	unit.LENGTH;
      BOX[i].ly	/=	unit.LENGTH;
      BOX[i].lz	/=	unit.LENGTH; 
   }

#ifdef MPI						// if parallelized
   system		=	CommRank -1 ;		// system id named by processor id

   nsites		=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {		// move mols in this processor
      if (moli->box == system) {			// up front in the mol array
         nsites		+=	moli->nsites;		// and change NMOLS and NSITES
         *oldmol	=	mol[nmols];		// such that in the following
         mol[nmols++]	=	*moli;			// this processor can only "see"
         *moli		=	*oldmol;		// these mols and sites
      }
   }
   NMOLS		=	nmols;			// # of mols in this processor
   NSITES		=	nsites;			// # of sites in this processor
   NMols[system]	=	nmols;
   NSites[system]	=	nsites;
#else							// not parallelized, but still multiple systems
   for (system=0; system<NSYSTEMS; system++) {
      NMols[system]	=	0;
      NSites[system]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if ( (system=moli->box) >= 0) {
         NMols[system]	++;				// total # of mols in certain system
         NSites[system]	+=	moli->nsites;		// total # of sites in certain system
      }
   }
#endif /* MPI */

   for (i=0; i<NSYSTEMS; i++) {
      BOX[i].lbox	=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
      BOX[i].vol	=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
      BOX[i].pres 	=	P;
      BOX[i].temp 	=	kT;
      BOX[i].drmax 	=	DRMAX;
      BOX[i].dlmax 	=	DLMAX;
      BOX[i].damax	=	DAMAX;
      BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
      BOX[i].rb		=	Rb;
      BOX[i].rv		=	Rv;
   } 

   for (i=0; i<NMOLS; i++) { 				// note NMOLS might change due to MPI
      for (j=0; j<mol[i].nsites; j++)  {
         mol[i].flags[j]	=	1;		// activate all the sites on this processor
         mol[i].parent[j]	=	j-1;		// initialize parent sites
      }
      mol[i].flip		=	0;		// flip to the original direction
      mol[i].origin		=	CenterofMass(mol+i);
   }
}


void InitSystem()
{
   if (!strcmp(INITCONF, "initconf"))			// read in coordinates from initconf
      GetCoordinates("initconf");
   else						
      Exit("init", "InitSystem", "initconf not found");

   AllSpherical();					// calculate spherical coordinates
}


void CheckSetup()					// check setup parameters
{
   if ( NTYPES>MAXNTYPES )	
      Exit("init", "CheckSetup", "NTYPES invalid!");
   if (NMOLS >MAXNMOLS)			Exit("init", "CheckSetup", "NMOLS invalid!");
   if (NSYSTEMS >MAXNSYSTEMS)		Exit("init", "CheckSetup", "NSYSTEMS invalid!");
}


void PrintInit()					// print out initial system info.
{
   long		i, j;

   PrintUnits();

   fprintf(foutput, "\nRandom number test:\n\n");
   for (i=0; i<5; i++)
      fprintf(foutput, "\trandom number #%d:\t%f\n", i+1, ran1(seed));

   fprintf(foutput, "\nInitial box setup (system units):\n");
   for (i=0; i<NSYSTEMS; i++) {
      fprintf(foutput, "\n\tBOX[%d]\n\n", i);
      fprintf(foutput, "\tNMols:\t%d\n", NMols[i]);
      fprintf(foutput, "\tNSites:\t%d\n", NSites[i]);
      fprintf(foutput, "\tLBOX:\t%f\n", BOX[i].lbox);
      fprintf(foutput, "\tLBOXX:\t%f\n", BOX[i].lx);
      fprintf(foutput, "\tLBOXY:\t%f\n", BOX[i].ly);
      fprintf(foutput, "\tLBOXZ:\t%f\n", BOX[i].lz);
      fprintf(foutput, "\tRc:\t%f\n", BOX[i].rc);
      fprintf(foutput, "\tkT:\t%f\n", BOX[i].temp);
#ifdef CELL_LIST
      fprintf(foutput, "\tCell[%d]:\t%3d %3d %3d\n", i, M[i].x, M[i].y, M[i].z);
#endif
   }
   fprintf(foutput, "\nInitial potential energy and virial (system units):\n");
   for (i=0; i<NSYSTEMS; i++) {
      fprintf(foutput, "\n\tBOX[%d]\n\n", i);
      fprintf(foutput, "\tlj energy (without long range):\t%f\n", v[i].lj);
      fprintf(foutput, "\tlj energy (long range part):\t%f\n", v[i].ljcorr);
      fprintf(foutput, "\tstretch energy:\t%f\n", v[i].stretch);
      fprintf(foutput, "\tbending energy:\t%f\n", v[i].bending);
      fprintf(foutput, "\ttorsion energy:\t%f\n", v[i].torsion);
      fprintf(foutput, "\ttotal energy:\t%f\n", v[i].tot);
      fprintf(foutput, "\tlj virial:\t%f\n", vir[i].lj);
      fprintf(foutput, "\tstretch virial:\t%f\n", vir[i].stretch);
      fprintf(foutput, "\ttotal virial:\t%f\n", vir[i].tot);
   }
   fflush(foutput);
}


void InitAll(char * argv[])
{
   InitMols(MAXNMOLS, MAXNMOLS);	// allocate max. memory for molecules
   InitFile();				// initialize the i/o files
   GetSetup(argv);			// read setup parameters
   CheckSetup();			// check setup parameters
 
   InitUnits();				// convert units, read in coord. will need units
   PrintSetup();			// print out setup file
   InitSystem();			// prepare coordinates and systems (boxes), and convert their units
  
#ifdef CELL_LIST
   CL_Init();				// Cell list should be ready before Verlet list
   CL_Build();
#endif	/* CELL_LIST */

   InitForcefield();			// initialize mixing rule
   CalcV();				// calculate total energy and virial
   if (UMBRELLA) {
      New_Qlm(l_of_Ylm);		// calculate initial qlm and Qlm
      New_Clist();			// set up connected neighbor list
      Find_Nuclei(dynvar);
      Calc_Qlm(6);			// calc. Q6
   }
   InitSample(argv);			// initialize distributions
   PrintInit();				// print initial system info.
}


// below is some old stuff...

void InitConfig()   //initialize the whole system, old stuff
{
  long		i, j, k;
  FILE 		*fconf;
  molstruct	*moli;

  S_STOP	=	0;  

  if (V_SCALECUTOFF) {
     Rc	=	MIN(0.5 * LBOX, Rc / LBOX0 * LBOX);	// initialize cutoff radii, with Rc0, Rb0 and Rv0 ...
     Rb	=	Rb / LBOX0 * LBOX;			// ... the desired value at equilibrium, and LBOX0
     Rv	=	Rv / LBOX0 * LBOX;			// ... is the estimated LBOX at equilibrium
  }
#ifdef CELL_LIST
  CL_Init();
  CL_Build();
#endif

#ifdef VERLET_LIST
  vlistplus		=	calloc(2*MAXVERLETNEIGH+1, sizeof(long));	// allocate for move-affected Verlet list
  New_Vlist();				// get Verlet list for each particle, using array
  maxnverlet = 0;
#endif	/* VERLET_LIST */
}
                                                                                                                                                                                                                                                                                         src/io.c                                                                                            0000600 0143352 0000144 00000141317 11312062016 011626  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
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
  fprintf(foutput, "[%d] %s::%s: %s\n", CommRank, module, procedure, error);
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
	
      rewind(fPtr);
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
         else
            fgets(comments, sizeof(comments), fPtr);
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

      fprintf(foutput, "%s", title1);
      fprintf(foutput, "%s", title2);
      fprintf(foutput, "Program Version: \t%s\n", VERSION);
      for (i=0; i<80; i++)	
	 fprintf(foutput, "*");
      fprintf(foutput, "\n");

      fprintf(foutput, "\n/* Basic system property. */\n\n");

      fprintf(foutput, "Molecule Type\t=\t%s\n", moltype);
      fprintf(foutput, "ConvertUnits\t=\t%s\n", (ConvertUnits ? "Yes":"No"));
      fprintf(foutput, "Rho 	\t=\t%f\n", Rho);
      fprintf(foutput, "Dpoly   \t=\t%f\n", Dpoly);
      fprintf(foutput, "kT	\t=\t%f\n", kT);
      fprintf(foutput, "P 	\t=\t%f\n", P);
      fprintf(foutput, "NBOX	\t=\t%d\n", NBOX);
      fprintf(foutput, "NMOLS    \t=\t%d\n", NMOLS);
      fprintf(foutput, "NSITES    \t=\t%d\n", NSITES);
      fprintf(foutput, "NPARTS    \t=\t%d\n", NPARTS);

      fprintf(foutput, "\n/* Basic simulation setup */\n\n");
      
      fprintf(foutput, "INITCONF  \t=\t%s\n", INITCONF);
      fprintf(foutput, "Stage 	\t=\t%d\n", Stage);
      fprintf(foutput, "PBC 	\t=\t%d\n", PBC);
#ifdef CELL_LIST
      fprintf(foutput, "CELL_LIST \t=\t%s\n", "Yes");
#else
      fprintf(foutput, "CELL_LIST \t=\t%s\n", "No");
#endif
      fprintf(foutput, "Nequil 	\t=\t%d\n", Nequil);
      fprintf(foutput, "Nprod 	\t=\t%d\n", Nprod);
      fprintf(foutput, "SEQUENCE  \t=\t%d\n", SEQUENCE);
      fprintf(foutput, "TRIALRUN  \t=\t%d\n", TRIALRUN);
      fprintf(foutput, "PROD 	\t=\t%d\n", PROD);
      fprintf(foutput, "ITAPE 	\t=\t%d\n", ITAPE);
      fprintf(foutput, "ICONF 	\t=\t%d\n", ICONF);
      fprintf(foutput, "NBLOCKS   \t=\t%d\n", NBLOCKS);
      fprintf(foutput, "NGSAMPLE  \t=\t%d\n", NGSAMPLE);
      fprintf(foutput, "IRADIAL   \t=\t%d\n", IRADIAL);
      fprintf(foutput, "NGRBINS   \t=\t%d\n", NGRBINS);
      fprintf(foutput, "CRIT 	\t=\t%f\n", CRIT);
      fprintf(foutput, "Alpha 	\t=\t%f\n", Alpha);

      fprintf(foutput, "\n/* Forcefield model */\n\n");

      fprintf(foutput, "NTYPES	\t=\t%d\n", NTYPES);
      fprintf(foutput, "DLJ	\t=\t%d\n", DLJ);
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
      fprintf(foutput, "critconn\t=\t%d\n", critconnect);
      fprintf(foutput, "SCALECUTOFF \t=\t%s\n", (V_SCALECUTOFF ? "Yes" : "No"));

      fprintf(foutput, "\n/* Ensemble setup */\n\n");
      
      fprintf(foutput, "E_NVT	\t=\t%s\n", (E_NVT ? "Yes" : "No"));
      fprintf(foutput, "E_NPT	\t=\t%s\n", (E_NPT ? "Yes" : "No"));
      fprintf(foutput, "E_GIBBS \t=\t%s\n", (E_GIBBS ? "Yes" : "No"));
      fprintf(foutput, "NDISPLACE \t=\t%d\n", NDISPLACE);
      fprintf(foutput, "NSWAP	\t=\t%d\n", NSWAP);
      fprintf(foutput, "NREPTATION\t=\t%d\n", NREPTATION);
      fprintf(foutput, "NENDROT   \t=\t%d\n", NENDROT);
      fprintf(foutput, "NCBMC	\t=\t%d\n", NCBMC);
      fprintf(foutput, "NENDBR	\t=\t%d\n", NENDBR);
      fprintf(foutput, "NREBR	\t=\t%d\n", NREBR);
      fprintf(foutput, "NFLIP	\t=\t%d\n", NFLIP);
      fprintf(foutput, "NDB  	\t=\t%d\n", NDB);
      fprintf(foutput, "NIDR 	\t=\t%d\n", NIDR);
      fprintf(foutput, "NVOLCHANGE\t=\t%d\n", NVOLCHANGE);
      fprintf(foutput, "NGIBBSVOL \t=\t%d\n", NGIBBSVOL);
      fprintf(foutput, "NTRIALCONF\t=\t%d\n", NTRIALCONF);
      fprintf(foutput, "NTRIALFIRSTBEAD\t=\t%d\n", NTRIALFIRSTBEAD);
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
      fprintf(foutput, "NMAXmiddle\t=\t%d\n", NMAXmiddle);
      fprintf(foutput, "NMAXbinsize\t=\t%d\n", NMAXbinsize);
      fprintf(foutput, "NMAXbins  \t=\t%d\n", NMAXbins);
      fprintf(foutput, "kN	\t=\t%f\n", kN);
      fprintf(foutput, "Qlmiddle  \t=\t%f\n", Qlmiddle);
      fprintf(foutput, "Qlbinsize \t=\t%f\n", Qlbinsize);
      fprintf(foutput, "Qlbins	\t=\t%d\n", Qlbins);
      fprintf(foutput, "kQ 	\t=\t%f\n", kQ);
      fprintf(foutput, "l_of_Ylm \t=\t%d\n", l_of_Ylm);
      fprintf(foutput, "critqlproduct \t=\t%f\n", critqlproduct);
      fprintf(foutput, "critconnect   \t=\t%d\n", critconnect);
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


void InitFile()					//initialize files
{
   long	i;

   if ((fhst=fopen("histogram.dat","w"))==NULL)
      Exit("io", "InitFile", "Histogram file cannot open!");
   if ((fdump=fopen("dump","w"))==NULL)
      Exit("io", "InitFile", "Dump file cannot open!");

   sprintf(file_hst, "binhist");
   frame = 0;

   if ((foutput=fopen("output","a"))==NULL)
      Exit("io", "InitFile", "Output file cannot open!");

   curtime	=	time(NULL);
   loctime	=	localtime(&curtime);
   for (i=0; i<80; i++)		fprintf(foutput, "*");
   fprintf(foutput, "\nOutput opened at:\t%s\n", asctime(loctime));
}


void CloseFile()				//close files
{
   long		i;
 
   for (i=0; i<NBOX; i++) {
      fprintf(foutput, "\n\tBOX[%d]:\n", i);
      fprintf(foutput, "\t# of mols:\t%d\n", NMols[i]);
      fprintf(foutput, "\t# of sites:\t%d\n", NSites[i]);
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
      fprintf(foutput, "\n\tBOX[%d]:\n", i);
      fprintf(foutput, "\tmove:\t%d\t acc:\t%f\n", av[i].move, ((double)av[i].acc_move)/av[i].move);
      fprintf(foutput, "\tvol: \t%d\t acc:\t%f\n", av[i].vol, ((double)av[i].acc_vol)/av[i].vol);
      fprintf(foutput, "\tcbmc: \t%d\t acc:\t%f\n", av[i].cbmc, ((double)av[i].acc_cbmc)/av[i].cbmc);
      fprintf(foutput, "\trep: \t%d\t acc:\t%f\n", av[i].rep, ((double)av[i].acc_rep)/av[i].rep);
      fprintf(foutput, "\tendrot:\t%d\t acc:\t%f\n", av[i].erot, ((double)av[i].acc_erot)/av[i].erot);
      fprintf(foutput, "\tendbr: \t%d\t acc:\t%f\n", av[i].eb, ((double)av[i].acc_eb)/av[i].eb);
      fprintf(foutput, "\trebr: \t%d\t acc:\t%f\n", av[i].re, ((double)av[i].acc_re)/av[i].re);
      fprintf(foutput, "\tflip: \t%d\t acc:\t%f\n", av[i].flip, ((double)av[i].acc_flip)/av[i].flip);
      fprintf(foutput, "\tdb: \t%d\t acc:\t%f\n", av[i].db, ((double)av[i].acc_db)/av[i].db);
      fprintf(foutput, "\tidr: \t%d\t acc:\t%f\n", av[i].idr, ((double)av[i].acc_idr)/av[i].idr);
      fprintf(foutput, "\tswap: \t%d\t acc:\t%f\n", av[i].swap, ((double)av[i].acc_swap)/av[i].swap);
      fprintf(foutput, "\trot: \t%d\t acc:\t%f\n", av[i].rot, ((double)av[i].acc_rot)/av[i].rot);
      fprintf(foutput, "\tmovemol: \t%d\t acc:\t%f\n", av[i].movemol, ((double)av[i].acc_movemol)/av[i].movemol);
      fprintf(foutput, "\tseq: \t%d\t acc:\t%f\n", av[i].seq, ((double)av[i].acc_seq)/av[i].seq);
   }
   if (NCBMC) {
      fprintf(foutput, "\tcbmc statistics:\n");
      for (i=0; i<NSITES/NMOLS; i++)
         fprintf(foutput, "\tcut point=\t%d\t accepted:\t%d\n", i, cbmcsucc[i]);
   }
   curtime	=	time(NULL);
   loctime	=	localtime(&curtime);

   fprintf(foutput, "\nOutput closed at:\t%s\n", asctime(loctime));
   for (i=0; i<80; i++) {
      fprintf(foutput, "*");
   }
   fprintf(foutput, "\n");

   fclose(foutput);

   fclose(fhst);
   fclose(fdump);
   fflush(stdin);
}


/****************************************************************************************/
/*	Print_Nuclei()									*/
/*											*/
/*	Print out nuclei information.							*/
/****************************************************************************************/

void Print_Nuclei()		//print out nuclei information
{
  long	size, id;
  long 	i;
  fprintf(foutput,"Printout nuclei size information......\n");
  /*
  for (id=1; id<NPARTS; id++) {			//nuclei index starts from 1
    if (sizeofnucl[id] != 0) {
      fprintf(foutput,"Size of nucleus #%d\t=%d\n", id, sizeofnucl[id]);
    }
  }
  */
  /*
  for (i=0; i<NPARTS; i++) {		//print out nuclei info. of each particle
    fprintf(foutput, "particle #%d\tbelongs to nuclei #%d\n", i+1, part[i].nuclid);
  }
  */
  for (i=1; i<NPARTS+1; i++) {		//printout size distribution
    if (sizedist[i]!=0) {
      fprintf(foutput, "Pring_Nuclei: Number of nuclei of size\t%d\t=%d\n", i, sizedist[i]);
    }
  }
  fprintf(foutput, "Maxsize:\t%d\n", MAXSIZE);
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
    fprintf(foutput,"Particle #%d has %d Verlet neighbors:\n",i+1, part[i].nverlet);  
    fprintf(foutput,"They are:\t");
    for (jj=0; jj<part[i].nverlet; jj++) {
      fprintf(foutput,"%d\t", part[i].vlist[jj]+1);
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
   fprintf(foutput,"%d\t",vlistplus[i]+1); 
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
      fprintf(foutput,"Particle #%d is connected to %d neighbors. They are:\t", i+1, part[i].nconnect);
      for (jj=0; jj<part[i].nconnect; jj++) {
        fprintf(foutput,"%d\t", part[i].clist[jj]+1);
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
      fprintf(foutput,"q(l=%d,m=%d,%d)=\t%+f+\t%+f*i\n", l, m, i+1, creal(part[i].qlm[m+l]), cimag(part[i].qlm[m+l]));
    }
    fprintf(foutput,"q(l=%d,%d)=\t%f\n\n", l, i+1, ql(l, i));
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
    fprintf(foutput, "aveQ(l=%d, m=%d)=\t%+f+\t%+f*i\n", l, m, creal(aveQlm(l, m)), cimag(aveQlm(l, m)));
  } 
  fprintf(foutput, "Q(l=%d)=\t%f\n", l, CalcQl(l));
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
        fprintf(foutput, "%d -> %d\t%+f+\t%+f*i\n", i+1, k+1, creal(term), cimag(term));
      }
    }
  }
  fprintf(foutput, "Printout ql product finished.\n\n");
  return;
}
*/
/****************************************************************************************/
/* 	Visualize(int type)								*/
/*											*/
/* 	Print out crystal nuclei configuration in pdb format for VMD visualization.	*/
/*	type = 0 (all particles) or 1 (biggest crystals only)				*/
/*	Color: O (red), N (blue), F (Green), H (White), C (blue-green)			*/
/****************************************************************************************/

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
   sprintf(tempname,"%d",frame);
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
   sprintf(tempname,"%d",frame);
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
      fprintf(fhere, "#MAXSIZE=%d\t# of biggest nuclei=%d\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);

      for (i=0; i<NPARTS; i++) {
         if (sizeofnucl[part[i].nuclid] == MAXSIZE) {           //print out only the biggest nucleus (nuclei)   

            pvisual     =       part[i].p;

            fprintf(fhere,"ATOM%7d  O   ACE%6d    %8.3f%8.3f%8.3f%6.2f%6.2f\n", part[i].nuclid, i+1, pvisual.x, pvisual.y, pvisual.z, 1.0, 0.0);
         }      //print out only the biggest nucleus/nuclei
      }
   }
   else if (vtype == 0) {                                       //all particles in different colors
/*
      fprintf(fhere, "#MAXSIZE=%d\t# of biggest nuclei=%d\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);
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
      fprintf(fhere, "#MAXSIZE=%d\t# of biggest nuclei=%d\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);
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
      fprintf(fhere, "#MAXSIZE=%d\t# of biggest nuclei=%d\tLBOX=%f\n", MAXSIZE, sizedist[MAXSIZE], LBOX);
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


int Read_Conf(char * filename)		// read in conf. file in SI units, NO unit conversion
{
   long		i, j, id;
   FILE *	fconf;
   char		dummy[80];

   if ((fconf=fopen(filename, "r"))==NULL)
      Exit("io", "Read_Conf", "Input conf. file not found or failed to open!");
   else {
      fscanf(fconf, "%s%ld", dummy, &TIMESTEP);			// read in timestep
      fscanf(fconf, "%ld%ld%ld", &NSYSTEMS, &NMOLS, &NSITES);	// # of systems, mols, and sites
      NBOX	=	NSYSTEMS;				// temporary
      for (i=0; i<NSYSTEMS; i++)
         fscanf(fconf, "%lf%lf%lf", &(BOX[i].lx), &(BOX[i].ly), &(BOX[i].lz));	// box dimension
      for (i=0; i<NMOLS; i++) {
         fscanf(fconf, "%ld%ld%ld", &id, &(mol[i].box), &(mol[i].nsites));	// molecule info.
         //fscanf(fconf, "%ld%ld%ld%ld%ld", &id, &(mol[i].box), &(mol[i].nsites), &(mol[i].fix), &(mol[i].flip));
         if (id!=i)	Exit("io", "Read_Conf", "molecule id mismatch");
         for (j=0; j<mol[i].nsites; j++)					// site info.
            fscanf(fconf, "%ld%lf%lf%lf", &(mol[i].type[j]), 
		&(mol[i].p[j].x), &(mol[i].p[j].y), &(mol[i].p[j].z));
      }
      //CoorSI2System();
   }
   fclose(fconf);
   return	1;
}

int Read_MultiConf(FILE *fPtr)
{
   long		i, j, id;
   char		dummy[80];

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
int Read_Conf(char * filename)			// read in LBOX and coord. of NMOLS molecules
{
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

int Write_Conf(long timestep)			// output conf. in SI units
{
   long		i, j;
   FILE		*fPtr;
   molstruct	*moli;
   double	dL = unit.LENGTH;

   if (timestep==-1) {				// open a file for a single conf.
      if ((fPtr=fopen("finalconf", "w"))==NULL)		
         Exit("io", "Write_Conf", "finalconf file cannot open");
   }
   else {
      fPtr	=	fdump;			// goes to dump file for continuous dumping
   }

   fprintf(fPtr, "TIMESTEP	%d\n", timestep);
   fprintf(fPtr, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
   for (i=0; i<NSYSTEMS; i++)
      fprintf(fPtr, "%f\t%f\t%f\n", BOX[i].lx * dL, BOX[i].ly * dL, BOX[i].lz * dL);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%d\t%f\t%f\t%f\n", moli->type[i], 
		moli->p[i].x * dL, moli->p[i].y * dL, moli->p[i].z * dL);
   }

   if (timestep==-1) {
      fclose(fPtr);
   }
   return	1;
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
//   sprintf(tempname,"%d",frame);

   if (timestep!=-1) {
      sprintf(tempname, "%d", timestep);

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
      fprintf(fconf, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);

      for (i=0; i<NSYSTEMS; i++)
         fprintf(fconf, "%f\t%f\t%f\n", BOX[i].lx*unit.LENGTH, BOX[i].ly*unit.LENGTH, BOX[i].lz*unit.LENGTH);

      for (moli=mol; moli<mol+NMOLS; moli++) {
         fprintf(fconf, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
         //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
         //MolInBox(moli);
         for (i=0; i<moli->nsites; i++) 
            fprintf(fconf, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x * unit.LENGTH, 
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
   sprintf(tempname,"%d",frame);
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
   sprintf(tempname,"%d",frame);
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
         fprintf(fhere, "ATOM%7d  C%12d    %8.3f%8.3f%8.3f\n", i*(moli->nsites)+j+1, i+1, l.x, l.y, l.z);
      }
      moli	++;
   }
   fclose(fhere);
   return	1;
}


void Printout()			// output all
{
  long 		i, j;
  vector 	dp;
  double 	r2, costheta, sintheta, phi, alpha;
 
  //fprintf(foutput, "Box length=\t%f\n", LBOX);

  /*	
  fprintf(foutput, "Mol#\tx\t\ty\t\tz\t\tnbond\tql\t\tnconnect\n");
  for (i=0; i<NPARTS; i++) 		//print out particles coordinates 
    fprintf(foutput, "%d\t%f\t%f\t%f\t%d\t%f\t%d\n", i+1, part[i].p.x, part[i].p.y, part[i].p.z, part[i].nbond, part[i].ql,part[i].nconnect);
  for (i=0; i<NPARTS-1; i++) {
    for (j=i+1; j<NPARTS; j++) {
      dp.x	=	part[j].p.x	-	part[i].p.x;
      dp.y	=	part[j].p.y	-	part[i].p.y;
      dp.z	=	part[j].p.z	-	part[i].p.z;
      MapInBox(&dp);
      r2	=	dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;
      fprintf(foutput, "r(%d,%d)=%f\n",i+1,j+1,sqrt(r2));
      costheta	=	dp.z/sqrt(r2);
      sintheta	=	sqrt(1-costheta*costheta);
        phi	=	atan2( dp.y, dp.x );
      alpha	=	pow((sqrt(r2)-BONDCUT)/(sigma-BONDCUT),2);
      fprintf(foutput, "%d->%d:\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",i+1, j+1, dp.x, dp.y, dp.z, costheta, sintheta, phi, alpha);
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
                                                                                                                                                                                                                                                                                                                 src/lammps2car.c                                                                                    0000600 0143352 0000144 00000014162 10755232145 013271  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
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
                                                                                                                                                                                                                                                                                                                                                                                                              src/lammps2hst.c                                                                                    0000600 0143352 0000144 00000112774 11623746455 013343  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	lammps2hst.c
	author:		Peng Yi at MIT
	date:		Feb. 14, 2008
	purpose:	read lammps dump file, do some sampling 
			and output histogram file
	note:		
			require setup file
			July 28, 2009	Add average center of mass calculation
 			March 7, 2011	Add segments statistics
			March 9, 2011	Allow skipping frames
			April 28, 2011	Allow multiple input data files 
			June 15, 2011	Allow triclinic box, though energy and pressure calculate
					not fixed yet
*/

#define __MAIN_PROGRAM
#include "header.h"

#define CAP	5		// largest nuclei considered still the melt
#define p2threshold	0.4
#define DEBUG	0

#include "correlation.h"

long		timestep;
double		rshift2;	// shift of the biggest nucleus
long		nsitesp2;	// # of sites with p2 greater than a threshold
long		nmolsp2;
long		nmaxp2_1[10];	// nmax using p2 nucleus definition
long		nmaxp2_2[10];
long		nmaxp2_3[10];
long		nmaxp2_4[10];

char *Element(long t, char *s)
{
   double	mass = type[t].M;

   s[0]		=	0;
/*
   if (14==mass || 15==mass)	strcpy(s, "C");
   else if (1.01==mass)		strcpy(s, "H");
   else if (28.086==mass)	strcpy(s, "Si");
   else if (26.982==mass)	strcpy(s, "Al");
   else if (16==mass)		strcpy(s, "O");
*/
   strcpy(s, "C");
   return	s;
}


void Print_hst(FILE *fPtr)
{
   long		i;
   static long	init = 1;

   if (init) {
      // print title

      fprintf(fPtr, "Timestep ");
      fprintf(fPtr, "E_pot ");
      fprintf(fPtr, "Volume ");
      fprintf(fPtr, "Pressure ");

      fprintf(fPtr, "Lx ");	// box dimension
      fprintf(fPtr, "Ly ");
      fprintf(fPtr, "Lz ");
      fprintf(fPtr, "xy ");	// triclinic box parameter
      fprintf(fPtr, "xz ");
      fprintf(fPtr, "yz ");

      fprintf(fPtr, "P2 ");
      fprintf(fPtr, "P2m ");
      fprintf(fPtr, "P2z ");
      fprintf(fPtr, "Transfrac ");
      fprintf(fPtr, "Xtal ");
      fprintf(fPtr, "RealXtal ");
      fprintf(fPtr, "Nnucl ");

fprintf(fPtr, "Nmaxp2 ");
fprintf(fPtr, "2ndNmaxp2 ");
/*
fprintf(fPtr, "Nmaxp2_1 ");
fprintf(fPtr, "2ndNmaxp2_2 ");
fprintf(fPtr, "Nmaxp2_2 ");
fprintf(fPtr, "2ndNmaxp2_2 ");
fprintf(fPtr, "Nmaxp2_3 ");
fprintf(fPtr, "2ndNmaxp2_3 ");
fprintf(fPtr, "Nmaxp2_4 ");
fprintf(fPtr, "2ndNmaxp2_4 ");
*/
      fprintf(fPtr, "Q6 ");
      fprintf(fPtr, "rshift2 ");
      fprintf(fPtr, "nsitesp2 ");

fprintf(fPtr, "E_lj ");
fprintf(fPtr, "E_ljcorr ");
fprintf(fPtr, "E_bond ");
fprintf(fPtr, "E_angle ");
fprintf(fPtr, "E_tors");

      fprintf(fPtr, "\n");
      init	=	0;
   }

   // output data

   fprintf(fPtr, "%-6d ", timestep);
   fprintf(fPtr, "%8.4f ", v[0].tot);
   fprintf(fPtr, "%8.4f ", BOX[0].vol);
   fprintf(fPtr, "%8.4f ", BOX[0].pres);

   fprintf(fPtr, "%5.3f ", BOX[0].lx);
   fprintf(fPtr, "%5.3f ", BOX[0].ly);
   fprintf(fPtr, "%5.3f ", BOX[0].lz);
   fprintf(fPtr, "%5.3f ", BOX[0].xy);
   fprintf(fPtr, "%5.3f ", BOX[0].xz);
   fprintf(fPtr, "%5.3f ", BOX[0].yz);

   fprintf(fPtr, "%6.4f ", P2[0]);
   fprintf(fPtr, "%6.4f ", P2M[0]);
   fprintf(fPtr, "%6.4f ", P2z[0]);
   fprintf(fPtr, "%6.4f ", transfrac[0]);
   fprintf(fPtr, "%4d ", Xtal[0]);
   fprintf(fPtr, "%4d ", realXtal[0]);
   fprintf(fPtr, "%4d ", Nnucl[0]);

for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_3[i]);
/*
for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_1[i]);
for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_2[i]);
for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_3[i]);
for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_4[i]);
*/

   fprintf(fPtr, "%6.4f ", Q6[0]);
   fprintf(fPtr, "%6.4f ", rshift2);
   fprintf(fPtr, "%5d", nsitesp2);

fprintf(fPtr, " %8.4f", v[0].lj);
fprintf(fPtr, " %8.4f", v[0].ljcorr);
fprintf(fPtr, " %8.4f", v[0].stretch);
fprintf(fPtr, " %8.4f", v[0].bending);
fprintf(fPtr, " %8.4f", v[0].torsion);

   fprintf(fPtr, "\n");
   fflush(fPtr);
}

void shiftbox(long system, beadstruct *nucleus, long nsites) //copied from conf2car.c (07/13/09)
{
   molstruct	*moli;
   long		i, n, site;
   vector	rA, rO, rBA, rOA;

   // step 1: find one segment that belongs to the biggest nucleus
   
   rA	=	nucleus[0].moli->p[nucleus[0].site];
  
   // step 2: calc. the shift of com of nucleus to this segment
   
   V_Null(&rBA);
   V_Null(&rOA);

   for (n=0; n<nsites; n++) { 
      moli	=	nucleus[n].moli;
      site	=	nucleus[n].site;
      rBA	=	V_Subtr(moli->p+site, &rA);
      rBA	=	MapInBox2(&rBA, PBC, system);
      rOA	=	V_Add(&rOA, &rBA);
   }
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);

   // step 3: every segment shift and move to central box
  
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (site=0; site<moli->nsites; site++) {
         moli->p[site]	=	V_Subtr(moli->p+site, &rO);
         moli->p[site]	=	MapInBox2(moli->p+site, PBC, system);
      }
   }
} 

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   molstruct	*moli;
   long		i, j, k, system, m, n, nuclid, nsites;
   long		previous_timestep=-1;
   long		nx, ny, nz, id, siteid, molid, type;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo,			// box upper boundary
		xy=0, xz=0, yz=0;		// triclinic parameter
   double	temp1, temp2, temp3;		// dummy variables

   char		infile[8][255], filein[80], filename[80];
   char		s[80], ff[80], par[80], dummy[255];
   FILE		*fin, *fhst, *fout, *fconf, *fpdb, *fdat;
   long		LENGTH, accum, confflag=0, carflag=0, pdbflag=0, polydisperse=0, drawmol;
   long		nframe, dnframe=1;		// analyze only every dnframe
   long		nfiles=1, ifile;		// number of input files 
   char		atomname;
   static long	init=1;	

   vector		con;	// center of nucleus
   static vector	rO;	// original position of nucleus

   // chain rotational angle distribution
   vector	chainface;
   double	orient;
   long		orientdist[180];

   // segment statistics variables
   long		previous, previous_id;
   long		nlength[MAXNMOLSITES], n_sep[MAXNMOLSITES];
   long		head, tail, seg_on, seg_id, nseg, nxtal, length;
   long		segment[MAXNMOLS][MAXNMOLSITES];	// segments identification on a chain
   long		loose_stat[MAXNMOLS][MAXNMOLSITES];	// loose segment stat
   long		nloop, nbridge, ntail, nxseg;		// # of total loops, etc
   double	lloop, lbridge, ltail, lxseg;		// average length of loops, etc
   long		nlloop[MAXNMOLSITES], nlbridge[MAXNMOLSITES];	// segment length distribution
   long		nltail[MAXNMOLSITES], nlxseg[MAXNMOLSITES];

   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];		// group beads in the same nucleus

   vector	com[MAXNMOLS], temp;		// average center of mass
   long		ncom[MAXNMOLS];

   long		imagen;				// variables for creating pbc images
   double	imagex, imagey, imagez;	

   if (argc<2) {
      printf("lammps2hst (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tlammps2hst [-option] [x= 0.1 y= 0.1 z= 0.1 n= 1] [dn= 2 nfiles= 2] lammpsdumpfile\n\n");
      printf("Notes:\n");
      printf("\t* -option = \n");
      printf("\t* -poly: polydisperse system, has molid in dump file\n");
      printf("\t* -conf: configuration file output\n");
      printf("\t* -car: car file output\n");
      printf("\t* -pdb: pdb file output\n");
      printf("\t* x= y= z=: duplicate the system and shift unit vector\n");
      printf("\t* n=: multiple of shift vector\n");
      printf("\t* dn=: only analyze every dn frames\n");
      printf("\t* nfiles=: number of input files if more than 1 (must be <=8)\n");
      printf("\t* \"=\" must immediately follow x or y or z or n or dn\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-poly"))	polydisperse	=	1;
      else if (samestr(par, "-conf"))	confflag	=	1;
      else if (samestr(par, "-car"))	carflag		=	1;
      else if (samestr(par, "-pdb"))	pdbflag		=	1;
      else if (samestr(par, "x="))	imagex		=	atof(argv[i+1]);
      else if (samestr(par, "y=")) 	imagey		=	atof(argv[i+1]);
      else if (samestr(par, "z=")) 	imagez		=	atof(argv[i+1]);
      else if (samestr(par, "n=")) 	imagen		=	atol(argv[i+1]);
      else if (samestr(par, "dn=")) 	dnframe		=	atol(argv[i+1]);
      else if (samestr(par, "nfiles="))	nfiles		=	atol(argv[i+1]);
   }
   for (i=0; i<nfiles; i++) {
      strcpy(infile[nfiles-1-i], argv[argc-1-i]);	// get input filenames
   }

   // Open output files
 
   if (nfiles==1)	strcpy(filein, infile[0]);
   else			strcpy(filein, "multi");

   strcpy(filename, filein);
   strcat(filename, ".hst");
   if ((fhst=fopen(filename, "w"))==NULL )
      Exit("lammps2hst", "main", "open hst file failed.");

   strcpy(filename, filein);
   strcat(filename, ".conf");
   if (confflag && (fconf=fopen(filename, "w"))==NULL )
      Exit("lammps2hst", "main", "open conf file failed.");   

   strcpy(filename, filein);
   strcat(filename, ".car");
   if (carflag && (fout=fopen(filename, "w"))==NULL )
      Exit("lammps2hst", "main", "open car file failed.");

   strcpy(filename, filein);
   strcat(filename, ".pdb");
   if (pdbflag && ((fpdb=fopen(filename, "w"))==NULL || (fdat=fopen("vmd.dat","w"))==NULL))
      Exit("lammps2hst", "main", "open pdb file failed.");

   if (!DEBUG) {
      strcpy(filename, filein);
      strcat(filename, ".out");
      freopen(filename, "w", stdout);	// redirect standard output stream to a file
   }

   ////////////////////
   // Initialization //
   ////////////////////

   if (DEBUG)	printf("Initialization: start ... \n");

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   InitUnits();				// initialize units for calculation
   InitForcefield();			// initialize Lennard-Jones potential mixing rule

   system	=	0;		// for now 2/14/08, only one system

   InitSample();			// initialize sampling

   for (i=0; i<180; i++)		// initialize chain rotation distribution
      orientdist[i]	=	0;
   for (i=0; i<MAXNMOLSITES; i++) {	// initialize segment variables
      n_sep[i]		=	0;	
      nlength[i]	=	0;
      nlloop[i]		=	0;
      nltail[i]		=	0;
      nlbridge[i]	=	0;
      nlxseg[i]		=	0;
   }
 
   nframe	=	-1;

   ///////////////////////////
   // Start data processing //
   ///////////////////////////

for (ifile=0; ifile<nfiles; ifile++){		// allow multiple input files
   fin	=	fopen(infile[ifile], "r");

   while (!feof(fin)) {

      if (DEBUG)	printf("Read one configuration: start ...\n");

      // Read in one configuration from dump file //

      if (!fgets(dummy, sizeof(dummy), fin)) {	// line 1
         break;					// end of file
      }
      fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);	// line 2
      fgets(dummy, sizeof(dummy), fin);							// line 3
      fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);	// line 4

      fgets(dummy, sizeof(dummy), fin);							// line 5
      if (strstr(dummy, "xy"))	PBC = 3;	// PBC=3 for triclinic box

      fscanf(fin, "%lf%lf", &xlo, &xhi);	if (PBC==3) fscanf(fin, "%lf", &xy); 	// line 6
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &ylo, &yhi);	if (PBC==3) fscanf(fin, "%lf", &xz); 	// line 7
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &zlo, &zhi);	if (PBC==3) fscanf(fin, "%lf", &yz);	// line 8
      fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);							// line 9

      xlo	-=	MIN(0.0, MIN(xy, MIN(xz, xy+xz)));	// triclinic box in general
      xhi	-=	MAX(0.0, MAX(xy, MAX(xz, xy+xz)));
      ylo	-=	MIN(0.0, yz);
      yhi	-=	MAX(0.0, yz);

      BOX[system].xy	=	xy;
      BOX[system].xz	=	xz;
      BOX[system].yz	=	yz;
      
      BOX[system].lx	=	xhi-xlo;
      BOX[system].ly	=	yhi-ylo;
      BOX[system].lz	=	zhi-zlo;

      LENGTH		=	NSITES/NMOLS;	// monodisperse system for now (4/26/2008)
      accum		=	0;

      for (i=0; i<nsites; i++) {
         fscanf(fin, "%ld", &id);
         fscanf(fin, "%ld", &molid);		// Need to check the lammps.dump file format
						// because some early lammps.dump file 
						// does not have molid output
//         if (polydisperse)	fscanf(fin, "%ld", &molid);
         fscanf(fin, "%ld", &type);
         fscanf(fin, "%lf%lf%lf %lf%lf%lf %ld%ld%ld", &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
         fgets(dummy, sizeof(dummy), fin);

         if (polydisperse) {			// polydisperse
            molid	--;
            /*accum	=	0;
            for (j=0; j<molid; j++)
               accum	+=	(mol+j)->nsites;
            siteid	=	id  - accum;
            */
         }
         else {					// monodisperse
            id	--; 				// -1 because lammps index starts from 1
   	    molid	=	(long) (id/LENGTH);
            siteid	=	id % LENGTH;
         }
         mol[molid].box		=	system;		// for now, only one system
         mol[molid].nsites	=	LENGTH;		// for now, Jan/10/2010

	 // triclinic box in general
         mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx) + ny*(BOX[system].xy) + nz*(BOX[system].xz);
         mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly) + nz*(BOX[system].yz);
         mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);
         
         mol[molid].type[siteid]=	type - 1;	// -1 because lammps index starts from 1
      }

      for (system=0; system<NSYSTEMS; system++) {
         NMols[system]	=	0;
         NSites[system]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;			// total # of mols in certain system
            NSites[system]	+=	moli->nsites;	// total # of sites in certain system
         }
      }

      for (moli=mol; moli<mol+NMOLS; moli++) { 		
         for (i=0; i<moli->nsites; i++)  {
            moli->flags[i]	=	1;		// activate all the sites on this processor
            moli->parent[i]	=	i-1;		// initialize parent site
         }
         moli->flip		=	0;		// flip to the original direction
         moli->origin		=	CenterofMass(moli);
      }

      ///////////////////////////////////////////////////
      // Skip repeated frames in different input files //
      ///////////////////////////////////////////////////

      if (ifile>=1 && timestep==previous_timestep) {
         continue;
      }
      previous_timestep	=	timestep;	// check repeat frames in different infiles

      ///////////////////////////
      // Skip dnframe-1 frames //
      ///////////////////////////

      nframe	++;
      if (mod(nframe, dnframe))	continue;	// analyze every dnframe frames
      
      ///////////////////////////////////////////////////////////////////////////
      // Start: Calculate the average position of center of mass of each chain //
      ///////////////////////////////////////////////////////////////////////////

      for (moli=mol; moli<mol+NMOLS; moli++) {
      	  temp	=	CenterofMass(moli);
          if (fabs(temp.x) < 0.45 * BOX[moli->box].lx && 
		fabs(temp.y) < 0.45 * BOX[moli->box].ly &&
		fabs(temp.z) < 0.45 * BOX[moli->box].lz) {
	     com[moli-mol]	=	V_Add(&temp, com+(moli-mol));
             ncom[moli-mol]	++;
          }
      }

      /////////////////////////////////////////////////////////////
      // Convert coordinates and box size from SI to system unit //
      /////////////////////////////////////////////////////////////
      CoorSI2System();
      //_________________________________________________________//

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
      } 
      /////////////////////
      // Build cell list //
      /////////////////////

#ifdef CELL_LIST	
/*
      if (init) {
         CL_Init();		
         init	=	0;
      } 
      else 
         CL_Destroy();

      CL_Build();
*/
 
      if (!init) {
	 CL_Destroy();
         init=0;
      }
      CL_Init();
      CL_Build();

/*
printf("Number of cells = %d (%d %d %d)\n", NCELLS, M[0].x, M[0].y, M[0].z);
for (i=0; i<NCELLS; i++) {
   printf("%f %f %f %f %f %f %d\n", Cell[i].p_min.x, Cell[i].p_max.x, Cell[i].p_min.y, Cell[i].p_max.y, Cell[i].p_min.z, Cell[i].p_max.z, Cell[i].nsites);
}
*/
#endif	/* CELL_LIST */


      //////////////////////
      // Perform Analysis //
      //////////////////////

      if (DEBUG)	printf("Analysis: start ...\n");

      CalcV();				// calc energy and virial
      if (V_VIRIAL)	Sample_Pressure();


      if (DEBUG)	printf("...sampling p2\n");

      SampleP2All();			// sample P2 and P2m and local p2
      Dist_p2();			// put local p2 into distribution

//SampleConnection();
//continue;
/*
printf("perform sampleM_Q\n");
SampleM_Q();
   printf("nconnect vs p2\n");
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++) 
         printf("%d %f\n", moli->np2[i], moli->p2[i]);
*/

      if (DEBUG)	printf("...sampling spherical coordinates\n");

      SampleSpherical();		// sample spherical coordinates
      Dist_Spherical();			// put spherical coord. into distribution

      if (DEBUG)	printf("...find nuclei\n");

      Find_Nuclei_p2(1);
//sizeofnucl	=	sizeofnuclp2;
//sizedist	=	sizedistp2;

//      Calc_Qlm(6);			// calc. Qlm for LJ system, require CELL_LIST

      ////////////////////////////////////////////
      /* Calculate the drift of biggest nucleus */
      ////////////////////////////////////////////
      /*
      if (timestep==0) {
         rO		=	CoM_MaxNucleus(0);
         rshift2	=	0;
      }
      else {
         con	=	CoM_MaxNucleus(0);
         con	=	V_Subtr(&con, &rO);
         con	=	MapInBox2(&con, PBC, system);
         rshift2	=	V_Dot(&con, &con);
      }*/

      /////////////////////////////////////////////////////////////////////
      // Compute chain orientation distribution, test the rotator phase. //
      /////////////////////////////////////////////////////////////////////

      if (DEBUG)	printf("...computer chain orientation\n");

      for (moli=mol; moli<mol+NMOLS; moli++) {
         V_Null(&chainface);

         for (i=0; i<moli->nsites; i++) {
            if (mod(i,2)) {
	       chainface	=	V_Add(&chainface, moli->p+i);
            }
	    else {
	       chainface	=	V_Subtr(&chainface, moli->p+i);
	    }
	 }
	 orient	=	atan2(chainface.y, chainface.x)	+ M_PI;
         //temp=CenterofMass(moli);
         //if (temp.z < -0.25 * BOX[system].lz && temp.z <0) 
	 orientdist[(int)(orient/(2*M_PI)*36)]	++;
      }

      /////////////////////////////////////////////////////////////////////
      // calculate the # of sites with local p2 greater than a threshold //
      /////////////////////////////////////////////////////////////////////

      nsitesp2	=	0;
      nmolsp2	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
//       nsitesp2	=	0;
         for (i=0; i<moli->nsites; i++) {
            if (moli->p2[i] > 0.4)
               nsitesp2	++;
         }
//            if (nsitesp2>1 && nsitesp2<6)
//	       nmolsp2	++;
      }

      ////////////////////////////////////////////////
      // Find the nucleus id of the biggest nucleus //
      ////////////////////////////////////////////////

      sizeofnucl	=	sizeofnuclp2;
      sizedist		=	sizedistp2;

      system	=	0;		// for now only one system, 4/16/2010
      i	= 1;
      while (sizeofnucl[i] != nmax[system][0]) {
         i ++;
      }
      nuclid	=	i;
   
      ////////////////////////////////////////////////////////////////////////////
      // Group the segments belong to the biggest nucleus based on sizeofnucl[] //////
      ////////////////////////////////////////////////////////////////////////////

      m	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if (system == moli->box) {
            for (i=0; i<moli->nsites; i++) {
               if (moli->nuclid[i] == nuclid) {
                  nucleus[m].moli	=	moli;
	          nucleus[m].site	=	i;
                  m	++;
	       }
            }
         }
      }
//      cylindershape(nucleus, nmax[system][0], nuclid); 
      //_______Group segments in the biggest nucleus______//

      ////////////////////////////////////////////////////////
      //  Segments statistics: counting loops, bridges, etc //
      ////////////////////////////////////////////////////////

      for (i=0; i<MAXNMOLS; i++) {
         for (j=0; j<MAXNMOLSITES; j++){
            segment[i][j]	=	-1;	// xtal seg identification
            loose_stat[i][j]	=	-1;	// loose seg stat, each seg needs 3 elem
         }
      }

      // calculate the distance between two nearby xtal-like beads
/*
      for (moli=mol; moli<mol+NMOLS; moli++) {
         previous_id	=	-1;			// previous xtal bead id
	 previous	=	0;
         for (i=0; i<moli->nsites; i++) {
            if (moli->nuclid[i] > 0) {			// xtal like bead
               if (previous_id > 0) {
                  n_sep[i - previous]	++;
               }
               previous_id	=	moli->nuclid[i];
	       previous	=	i;
            }
         }
      }
*/
      // coarse graining the xtal segment
/*
      for (moli=mol; moli<mol+NMOLS; moli++) {
         for (i=0; i<moli->nsites-8; i++) {
            head	=	-1;
	    tail	=	-1;
	    nxtal	=	0;
            for (j=i; j<i+8; j++) {			// look a segment of 8 beads
               if (moli->nuclid[j] >0) {
                  nxtal	++;
		  tail	=	j;			// final xtal bead in this seg
		  if (head == -1) {
	             head	=	j;		// first xtal bead in this seg
                  }
               }
            }
            if (nxtal >= 5) {				// 5 out of 8 are xtal beads
               seg_id	=	moli->nuclid[head];	// seg_id = nucl_id
               for (j=head; j<=tail; j++) {		// belong to same xtal seg
                  segment[moli-mol][j]	=	seg_id;

                  if (moli->nuclid[j]>0 && moli->nuclid[j] != seg_id) {	// sanity check
                     printf("Error: more than two nuclei in one xtal segment!\n");
                     printf("mol [%d] head %d tail %d\n", moli-mol, head, tail);
                  }
               }
            }
         }
      }
*/
      // Segment statistics: xtal segments and loose segments
/*
      // Based on segment[mol][site]
      for (moli=mol; moli<mol+NMOLS; moli++) {
         j		=	moli-mol;
         previous	=	segment[j][0];

         nseg		=	0;

         seg_stat[j][0]	=	0;			// head position of 1st segment
         seg_stat[j][1]	=	previous;		// segment id

         for (i=0; i<moli->nsites; i++) {
            if (segment[j][i] != previous) {
               seg_stat[j][nseg*3+2]	=	i-1;		// tail position
               seg_stat[j][nseg*3+1]	=	previous;	// seg id
               nseg	++;
               seg_stat[j][nseg*3]	=	i;	// head position of next seg

               if (i==moli->nsites-1) {
                  seg_stat[j][nseg*3+2]	=	i;
                  seg_stat[j][nseg*3+1]	=	segment[j][i];
		  nseg	++;
               }
	    }
            else if (i==moli->nsites-1) {
               seg_stat[j][nseg*3+2]	=	i;
               seg_stat[j][nseg*3+1]	=	segment[j][i];
	       nseg	++;
            } 
	    previous	=	segment[j][i];
	 }
         nsegment[j]	=	nseg;
      }
*/
      // Based on nuclid

//      Find_segments();			// find all segments, do fluctuation analysis if needed
//      Seg_smooth();

      for (i=0; i<NMOLS; i++) {
         k	=	0;
         for (j=0; j<nsegment[i]; j++) {
            length	=	seg_stat[i][j*3+2]-seg_stat[i][j*3]+1;
            k		+=	length;	// for sanity check

            if (seg_stat[i][j*3+1] == -1) 	// loose segment
               n_sep[length]	++;
            else 				// xtal segment
               nlength[length]	++;
         }
         if (k!=(mol+i)->nsites) 			// sanity check
            printf("Error; k!=(mol+i)->nsites, k=%d, moli->nsites=%d\n",k, (mol+i)->nsites);         
      }

      // Segment analysis for the BIGGEST nucleus

      nloop	=	0;			// # of loop for the biggest nucleus
      nbridge	=	0;			// # of bridge
      ntail	=	0;			// # of tail
      nxseg	=	0;			// # of xtal segment
      lloop	=	0.0;			// average loop length
      ltail	=	0.0;			// average tail length
      lbridge	=	0.0;			// average bridge length
      lxseg	=	0.0;			// average xseg length

      for (moli=mol; moli<mol+NMOLS; moli++) {
         k	=	moli-mol;

         if (nsegment[k] >1) {				// containing xtal segment
            for (j=0; j<nsegment[k]; j++) {
               length=	seg_stat[k][j*3+2]-seg_stat[k][j*3]+1;

               if (seg_stat[k][j*3+1]==-1) {		// loose segment
		  head	=	(j==0 ? -1 : seg_stat[k][(j-1)*3+1]);	
                  tail	=	(j==nsegment[k]-1 ? -1 : seg_stat[k][(j+1)*3+1]);

                  if ((head==-1 && tail==nuclid) || (head==nuclid && tail==-1)) {
                     ntail		++;
                     ltail		+=	length;
                     nltail[length]	++;
                  }
                  else if (head==nuclid && tail==nuclid) {
                     if (length >= 4) { 
                        nloop		++;
                        lloop		+=	length;
                     }
                     nlloop[length]	++;
                  }
                  else if (head==nuclid || tail==nuclid) {
                     if (length >= 4) {
                        nbridge	++;
                        lbridge		+=	length;
                     }
                     nlbridge[length]	++;
                  }   
               }
               else if (seg_stat[k][j*3+1]==nuclid) {	// nseg xtal segment
                  if (length >= 4) {
                     nxseg	++;
                     lxseg	+=	length;
                  }
                  nlxseg[length]	++;
               }
            }
         }
      }
      printf("%5d ", nmax[0][0]);
      printf("%5d %5d %5d %5d ", ntail, nloop, nbridge, nxseg);
      printf("%8.3f %8.3f %8.3f %8.3f\n", ltail/ntail, lloop/nloop, lbridge/nbridge, lxseg/nxseg);

      //////////////////////////////////////////////////
      // use another nuclei definition for comparison //
      //////////////////////////////////////////////////
      if (DEBUG)	printf("...use another nuclei definition\n");

      nmaxp2_3[0]	=	nmax[0][0];	// Rp=2.5, Rconn=1.3, critp2=0.4
      nmaxp2_3[1]	=	nmax[0][1];
/*
      temp1	=	Rp;		// store original criterion
      temp2	=	Rconn;
      temp3	=	critp2;

      Rp	=	1.5;		// new definition
      Rconn	=	1.5;
      critp2	=	0.6;
      SampleP2All();			// sample P2 and P2m and local p2
      Find_Nuclei_p2(1);
      nmaxp2_1[0]	=	nmax[0][0];
      nmaxp2_1[1]	=	nmax[0][1];

      Rp	=	2.0;		// another definition
      Rconn	=	1.5;
      critp2	=	0.5;
      SampleP2All();			// sample P2 and P2m and local p2
      Find_Nuclei_p2(1);
      nmaxp2_2[0]	=	nmax[0][0];
      nmaxp2_2[1]	=	nmax[0][1];

      Rp	=	3.0;		// another definition
      Rconn	=	1.3;
      critp2	=	0.4;
      SampleP2All();			// sample P2 and P2m and local p2
      Find_Nuclei_p2(1);
      nmaxp2_4[0]	=	nmax[0][0];
      nmaxp2_4[1]	=	nmax[0][1];

      Rp	=	temp1;		// restore the original criterion
      Rconn	=	temp2;
      critp2	=	temp3;
*/
      //_______Another nuclei definition______________//

      /////////////////////////////////////////
      // Correlation and Radial distribution //
      /////////////////////////////////////////

      if (DEBUG)	printf("Calculate correlation: start...\n");

      correlation();			// calculate correlation, in system units
    //  if (!(timestep%IRADIAL)) {
//         radial("sample");		// sample radial distribution function
//         sq(stdout, "sample");
    //  }

//      shiftbox(system, nucleus, n);

      /////////////////////////////
      // Output analysis results //
      /////////////////////////////
      
      if (DEBUG)	printf("Output: start ...\n");

      Print_hst(fhst);		// print out histgram
      CoorSystem2SI();		// convert coordinates and box size back to SI units

      ////////////////////////////////////////////
      // OUTPUT .car file for VMD visualization //
      ////////////////////////////////////////////
      
      if (carflag) {
         fprintf(fout, "!BIOSYM archive 3\n");
         fprintf(fout, "PBC=ON\n");
         fprintf(fout, "!TIMESTEP %d\n", timestep);
         fprintf(fout, "!DATE %s", asctime(localtime(&t)));
         fprintf(fout, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

	 n	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {

               MolInBox2(moli);
               for (i=0; i<moli->nsites; i++) {
                  //moli->p[i] = MapInBox2(moli->p+i, PBC, system); //temp

                  if (moli->nuclid[i]>0)	// crystal like particle
                     sprintf(s, "N%d", n++);	// N: blue color in VMD
                  else
	             sprintf(s, "O%d", n++);	// O: red color in VMD

/*
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])		// note nuclid index starts from 1
                     sprintf(s, "M%d", n++);
                  else if (sizeofnucl[moli->nuclid[i]] -MAXSIZE[system] > -3 && MAXSIZE[system]>=10)
                     sprintf(s, "C%d", n++);
	          else if (moli->nuclid[i] >= 1)
                     sprintf(s, "O%d", n++);
                  else
                     sprintf(s, "H%d", n++);
*/
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
      fflush(fout);
      //___________OUTPUT .car file_____________//


      ///////////////////////////////////////////
      // OUTPUT .pdb file for further analysis */
      ///////////////////////////////////////////

      if (pdbflag) {
         fprintf(fpdb, "HEADER: file created from %s on %s", argv[1], asctime(localtime(&t)));
         fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
        	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

	 m		=	0;	// molecule sequence number
         n		=	0;	// atom sequence number
#define SIZECAP	5
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {
               //MolInBox2(moli);

               drawmol	=	0;
               for (i=0; i<moli->nsites; i++) {
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {
	             drawmol	=	1;		// participate the biggest nucleus
		     break;
                  }
	       }
//temp=CenterofMass(moli);
//if (temp.z < -0.25 * BOX[system].lz) {
               m	++; 
               for (i=0; i<moli->nsites; i++) {
                  if (drawmol) {
                     if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {	// nuclid index starts from 1
                        atomname	=	'N';	// N: blue color in VMD
		        fprintf(fdat, " 10");
                     }
                     else {
                        atomname	=	'O';	// O: red color in VMD
                        fprintf(fdat, " 0");
                     }
                  }
                  else {
		     atomname	=	'C';		// C: cyan color in VMD
		     fprintf(fdat," -1");
	          }

                  n	++;
	          fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
                  fprintf(fpdb, "%5d ", n);	// atom number
                  fprintf(fpdb, " %c  ", atomname);	// atom name
                  fprintf(fpdb, " ");		// alternate location indiator
  	          fprintf(fpdb, " C8");		// residue name
	          fprintf(fpdb, " ");		// column 21
                  fprintf(fpdb, " ");		// chain identifier, column 22
	          fprintf(fpdb, "%4d", m);	// residue sequence number, 23-26
	          fprintf(fpdb, " ");		// code for insertion of residues, 27
                  fprintf(fpdb, "   ");		// column 28-30
                  fprintf(fpdb, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z);
                  fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                  fprintf(fpdb, "%5.5s", "");
                  fprintf(fpdb, "\n"); 

                  if (imagen) {			// for image box
                     n	++;
                     fprintf(fpdb, "ATOM  ");
                     fprintf(fpdb, "%5d ", n);
                     fprintf(fpdb, " %c  ", atomname);
                     fprintf(fpdb, " ");
                     fprintf(fpdb, " C8");
	             fprintf(fpdb, " ");
                     fprintf(fpdb, " ");
	             fprintf(fpdb, "%4d", m);
                     fprintf(fpdb, " ");
                     fprintf(fpdb, "   ");
                     fprintf(fpdb, "%8.3f%8.3f%8.3f", 
			moli->p[i].x + imagex, moli->p[i].y+imagey, moli->p[i].z+imagez);
                     fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                     fprintf(fpdb, "%5.5s", "");
                     fprintf(fpdb, "\n"); 
                  }
               } 
//}
            }   
         }
         fprintf(fpdb, "END\n");
         fflush(fpdb);
      }	//pdbflag
      //____________OUTPUT .pdb file___________//

      ////////////////////////////////////////////////////
      // OUTPUT configuration file for further analysis //
      ////////////////////////////////////////////////////

      if (confflag) { 
         fprintf(fconf, "!TIMESTEP %d\n", timestep);
         fprintf(fconf, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
         for (i=0; i<NSYSTEMS; i++)
            fprintf(fconf, "%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

         for (moli=mol; moli<mol+NMOLS; moli++) {
            fprintf(fconf, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
            //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
            //MolInBox(moli);
            for (i=0; i<moli->nsites; i++) 
               fprintf(fconf, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
         }
         fflush(fconf);
      }
      //____________OUTPUT configuration file___________//
   } 

   fclose(fin);		// close current input file
  }			// multiple input files

   /////////////////////////////////////////////
   // output average center of mass of chains //
   /////////////////////////////////////////////

   if (pdbflag) {
      fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
      fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
 	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

      n		=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         i	=	moli-mol;
         if (ncom[i]>0) {
//if( fabs(com[i].x/ncom[i])<3 && fabs(com[i].y/ncom[i])<3 && fabs(com[i].z/ncom[i])<7 ) {
            n	++;
            atomname	=	'N';
            fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
            fprintf(fpdb, "%5d ", n);		// atom number
            fprintf(fpdb, " %c  ", atomname);	// atom name
            fprintf(fpdb, " ");		// alternate location indiator
            fprintf(fpdb, "   ");	// residue name
            fprintf(fpdb, " ");		// column 21
            fprintf(fpdb, " ");		// chain identifier, column 22
            fprintf(fpdb, "    ");	// residue sequence number, 23-26
            fprintf(fpdb, " ");		// code for insertion of residues, 27
            fprintf(fpdb, "   ");	// column 28-30
            fprintf(fpdb, "%8.3f%8.3f%8.3f", com[i].x/ncom[i], com[i].y/ncom[i], com[i].z/ncom[i]);
            fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
            fprintf(fpdb, "\n"); 
         }
//}
      }
      fprintf(fpdb, "END\n");
      fflush(fpdb);
   }
   //______OUTPUT average center of mass of chains to .pdb file__//

   S_PrintAll();		// print out distributions
   corr_norm();			// normalize correlation function
   corr_print();		// print out correlation function
   radial("print");		// print out radial distribution function
   sq(stdout, "print");

/*
   for (i=0; i<180; i++) {
      printf("%d %d\n", i, orientdist[i]);
   }
*/
   /////////////////////////////////////
   // Output segment analysis results //
   /////////////////////////////////////

   printf("**********Segment length statistics**********\n");
   printf("length\tn_sep\tnlength\tnltail\tnlloop\tnlbridge\tnlxseg\n");
   for (i=0; i<NSITES/NMOLS; i++) {
      printf("%d\t %d\t %d\t %d\t %d\t %d\t %d\n", 
	i, n_sep[i], nlength[i], nltail[i], nlloop[i], nlbridge[i], nlxseg[i]);
      if (n_sep[i]-nltail[i]-nlloop[i]-nlbridge[i] < 0)	printf("Error 1\n");
      if (nlength[i] < nlxseg[i])			printf("Error 2\n");
   }

   printf("**********Chain segment information**********\n");
   for (i=0; i<NMOLS; i++) {
      if (nsegment[i] == 1) {			// condition for output
         printf("chain %d has %d segments\n", i, nsegment[i]);
         for (j=0; j<(mol+i)->nsites; j++) {
            if (segment[i][j]>0 )
               printf("%d_%d  ", j, segment[i][j]);
         }
         printf("\n");
         for (j=0; j<nsegment[i]; j++) {
            printf("head = %d\t id = %d\t tail = %d\t length = %d\n", seg_stat[i][j*3], seg_stat[i][j*3+1], seg_stat[i][j*3+2], seg_stat[i][j*3+2]-seg_stat[i][j*3]+1);
         }
         printf("\n");
      }
   }

   printf("**********Biggest nucleus segment analysis**********\n");
   for (i=0; i<NMOLS; i++) {
      for (j=0; j<nsegment[i]; j++) {
         if (seg_stat[i][j*3+1] == nuclid) {		// related to BIGGEST nucleus
            printf("chain #%d\n", i);
            printf("head\t id\t tail\t length\n");
            for (k=0; k<nsegment[i]; k++) {
               printf("%d\t%d\t%d\t%d\n", 
			seg_stat[i][k*3], seg_stat[i][k*3+1], seg_stat[i][k*3+2], seg_stat[i][k*3+2]-seg_stat[i][k*3]+1);
            }
            printf("\n");
            break;
         }
      }
   }

   /////////////
   // Closing //
   /////////////

   if (DEBUG)	printf("Closing output files ...\n");

   fclose(fhst);	// close output files
   if (carflag)		fclose(fout);
   if (confflag)    	fclose(fconf);
   if (pdbflag)	     {	fclose(fpdb); fclose(fdat);}

   fflush(stdout);
   fclose(stdout);

   return	0;
}
    src/lammpsdump.c                                                                                    0000600 0143352 0000144 00000004505 11034433270 013400  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
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
            fprintf(fout, "ITEM: ATOMS\n");
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
                                                                                                                                                                                           src/lammpstest.c                                                                                    0000600 0143352 0000144 00000050270 11453124066 013417  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	lammpstest.c
	author:		Peng Yi at MIT
	date:		Oct. 05, 2010
	purpose:	read lammps dump file, do some simple analysis
	note:		
			require setup file
*/

#define __MAIN_PROGRAM
#include "header.h"

#define CAP	5		// largest nuclei considered still the melt
#define p2threshold	0.4

#include "correlation.h"

long		timestep;
double		rshift2;		// shift of the biggest nucleus
long		nsitesp2;	// # of sites with p2 greater than a threshold
long		nmolsp2;
long		nmaxp2_1[10];		// nmax using p2 nucleus definition
long		nmaxp2_2[10];
long		nmaxp2_3[10];
long		nmaxp2_4[10];

char *Element(long t, char *s)
{
   double	mass = type[t].M;

   s[0]		=	0;
/*
   if (14==mass || 15==mass)	strcpy(s, "C");
   else if (1.01==mass)		strcpy(s, "H");
   else if (28.086==mass)	strcpy(s, "Si");
   else if (26.982==mass)	strcpy(s, "Al");
   else if (16==mass)		strcpy(s, "O");
*/
   strcpy(s, "C");
   return	s;
}


void Print_hst(FILE *fPtr)
{
   long		i;
   static long	init = 1;

   if (init) {
      fprintf(fPtr, "timestep ");
      fprintf(fPtr, "poteng ");
      fprintf(fPtr, "volume ");
      fprintf(fPtr, "pressure ");
      fprintf(fPtr, "P2 ");
      fprintf(fPtr, "P2m ");
      fprintf(fPtr, "P2z ");
      fprintf(fPtr, "transfrac ");
      fprintf(fPtr, "Xtal ");
      fprintf(fPtr, "realXtal ");
      fprintf(fPtr, "Nnucl ");
      fprintf(fPtr, "Nmaxp2_1 ");
      fprintf(fPtr, "2ndNmaxp2_2 ");

fprintf(fPtr, "Nmaxp2_2 ");
fprintf(fPtr, "2ndNmaxp2_2 ");

fprintf(fPtr, "Nmaxp2_3 ");
fprintf(fPtr, "2ndNmaxp2_3 ");

fprintf(fPtr, "Nmaxp2_4 ");
fprintf(fPtr, "2ndNmaxp2_4 ");

      fprintf(fPtr, "Q6 ");
      fprintf(fPtr, "rshift2 ");
      fprintf(fPtr, "nsitesp2");
      fprintf(fPtr, "\n");
      init	=	0;
   }
   fprintf(fPtr, "%-6d ", timestep);
   fprintf(fPtr, "%8.4f ", v[0].tot);
   fprintf(fPtr, "%8.4f ", BOX[0].vol);
   fprintf(fPtr, "%8.4f ", BOX[0].pres);
   fprintf(fPtr, "%6.4f ", P2[0]);
   fprintf(fPtr, "%6.4f ", P2M[0]);
   fprintf(fPtr, "%6.4f ", P2z[0]);
   fprintf(fPtr, "%6.4f ", transfrac[0]);
   fprintf(fPtr, "%4d ", Xtal[0]);
   fprintf(fPtr, "%4d ", realXtal[0]);
   fprintf(fPtr, "%4d ", Nnucl[0]);
   for (i=0; i<2; i++)
      fprintf(fPtr, "%4d ", nmaxp2_1[i]);
   for (i=0; i<2; i++)
      fprintf(fPtr, "%4d ", nmaxp2_2[i]);
   for (i=0; i<2; i++)
      fprintf(fPtr, "%4d ", nmaxp2_3[i]);
   for (i=0; i<2; i++)
      fprintf(fPtr, "%4d ", nmaxp2_4[i]);

   fprintf(fPtr, "%6.4f ", Q6[0]);
   fprintf(fPtr, "%6.4f ", rshift2);
   fprintf(fPtr, "%5d", nsitesp2);
fprintf(fPtr, " %8.4f", v[0].lj);
fprintf(fPtr, " %8.4f", v[0].ljcorr);
fprintf(fPtr, " %8.4f", v[0].stretch);
fprintf(fPtr, " %8.4f", v[0].bending);
fprintf(fPtr, " %8.4f", v[0].torsion);
   fprintf(fPtr, "\n");
   fflush(fPtr);
}

void shiftbox(long system, beadstruct *nucleus, long nsites) //copied from conf2car.c (07/13/09)
{
   molstruct	*moli;
   long		i, n, site;
   vector	rA, rO, rBA, rOA;

   // step 1: find one segment that belongs to the biggest nucleus
   
   rA	=	nucleus[0].moli->p[nucleus[0].site];
  
   // step 2: calc. the shift of com of nucleus to this segment
   
   V_Null(&rBA);
   V_Null(&rOA);

   for (n=0; n<nsites; n++) { 
      moli	=	nucleus[n].moli;
      site	=	nucleus[n].site;
      rBA	=	V_Subtr(moli->p+site, &rA);
      rBA	=	MapInBox2(&rBA, PBC, system);
      rOA	=	V_Add(&rOA, &rBA);
   }
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);

   // step 3: every segment shift and move to central box
  
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (site=0; site<moli->nsites; site++) {
         moli->p[site]	=	V_Subtr(moli->p+site, &rO);
         moli->p[site]	=	MapInBox2(moli->p+site, PBC, system);
      }
   }
} 

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   long		i, j, k, system, m, n, nuclid;
   long		nsites, record;
   long		nx, ny, nz, id, siteid, molid, type;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo;			// box upper boundary
   molstruct	*moli;
   double	temp1, temp2, temp3;		// temperory variables

   char		filein[255], s[80], ff[80], par[80], dummy[255];
   FILE		*fin, *fhst, *fout, *fconf, *fpdb, *fdat;
   long		LENGTH, accum, confflag=0, carflag=0, pdbflag=0, polydisperse=0, drawmol;
   char		atomname;
   static long	init=1;

   vector		con;	// center of nucleus
   static vector	rO;	// original position of nucleus

   vector		chainface;
   double		orient;
   long			orientdist[180];	// chain orientation distribution

   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];

   vector	com[MAXNMOLS], temp;		// average center of mass
   long		ncom[MAXNMOLS];

   double	Rg2,				// radius of gyration square
		R02;				// end-to-end distance square
   double	imagex, imagey, imagez;		// variables for creating pbc images
   long		imagen;

   if (argc<2) {
      printf("lammpstest (c) 2010 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tlammpstest [-option] [x= 0.1 y= 0.1 z= 0.1 n= 1] lammpsdumpfile\n\n");
      printf("Notes:\n");
      printf("\t* -option = \n");
      printf("\t* -poly: polydisperse system, has molid in dump file\n");
      printf("\t* -conf: configuration file output\n");
      printf("\t* -car: car file output\n");
      printf("\t* -pdb: pdb file output\n");
      printf("\t* x= y= z=: duplicate the system and shift unit vector\n");
      printf("\t* n=: multiple of shift vector\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-poly"))
         polydisperse	=	1;
      else if (samestr(par, "-conf")) 
         confflag	=	1;
      else if (samestr(par, "-car"))
         carflag	=	1;
      else if (samestr(par, "-pdb"))
         pdbflag	=	1;
      else if (samestr(par, "x=")) 
         imagex	=	atof(argv[i+1]);
      else if (samestr(par, "y=")) 
         imagey	=	atof(argv[i+1]);
      else if (samestr(par, "z=")) 
         imagez	=	atof(argv[i+1]);
      else if (samestr(par, "n=")) 
         imagen	=	atol(argv[i+1]);
   }
   strcpy(filein, argv[argc-1]);

   if ( (fin=fopen(filein, "r"))==NULL )
      Exit("lammpstest", "main", "open input file failed.");
   if ((fhst=fopen("lammps.hst", "w"))==NULL )
      Exit("lammpstest", "main", "open hst file failed.");
   if (confflag && (fconf=fopen("lammps.conf", "w"))==NULL )
      Exit("lammpstest", "main", "open conf file failed.");   
   if (carflag && (fout=fopen("lammps.car", "w"))==NULL )
      Exit("lammpstest", "main", "open car file failed.");
   if (pdbflag && ((fpdb=fopen("lammps.pdb", "w"))==NULL || (fdat=fopen("vmd.dat","w"))==NULL))
      Exit("lammpstest", "main", "open pdb file failed.");

   printf("Initialization: start ... \n");

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   InitUnits();				// initialize units for calculation
   InitForcefield();			// initialize Lennard-Jones potential mixing rule

   system	=	0;		// for now 2/14/08, only one system

   InitSample();			// initialize sampling

   for (i=0; i<180; i++)		// initialize chain rotation distribution
      orientdist[i]	=	0;

   while (!feof(fin)) {

      printf("Read one configuration: start ...\n");

      /* Read in one configuration from dump file */

      if (!fgets(dummy, sizeof(dummy), fin))	break;			// end of file
      fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &xlo, &xhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &ylo, &yhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &zlo, &zhi);	fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);

      BOX[system].lx	=	xhi-xlo;
      BOX[system].ly	=	yhi-ylo;
      BOX[system].lz	=	zhi-zlo;

      LENGTH		=	NSITES/NMOLS;	// monodisperse system for now (4/26/2008)
      accum		=	0;

      for (i=0; i<nsites; i++) {
         fscanf(fin, "%ld", &id);
         fscanf(fin, "%ld", &molid);		// Need to check the lammps.dump file format
						// because some early lammps.dump file 
						// does not have molid output
//         if (polydisperse)	fscanf(fin, "%ld", &molid);
         fscanf(fin, "%ld", &type);
         fscanf(fin, "%lf%lf%lf %lf%lf%lf %ld%ld%ld", &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
         fgets(dummy, sizeof(dummy), fin);

         if (polydisperse) {			// polydisperse
            molid	--;
            /*accum	=	0;
            for (j=0; j<molid; j++)
               accum	+=	(mol+j)->nsites;
            siteid	=	id  - accum;
            */
         }
         else {					// monodisperse
            id	--; 				// -1 because lammps index starts from 1
   	    molid	=	(long) (id/LENGTH);
            siteid	=	id % LENGTH;
         }
         mol[molid].box		=	system;		// for now, only one system
         mol[molid].nsites	=	LENGTH;		// for now, Jan/10/2010
         mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx);
         mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly);
         mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);
         mol[molid].type[siteid]=	type - 1;	// -1 because lammps index starts from 1
      }

      for (system=0; system<NSYSTEMS; system++) {
         NMols[system]	=	0;
         NSites[system]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;			// total # of mols in certain system
            NSites[system]	+=	moli->nsites;	// total # of sites in certain system
         }
      }

      for (moli=mol; moli<mol+NMOLS; moli++) { 		
         for (i=0; i<moli->nsites; i++)  {
            moli->flags[i]	=	1;		// activate all the sites on this processor
            moli->parent[i]	=	i-1;		// initialize parent site
         }
         moli->flip		=	0;		// flip to the original direction
         moli->origin		=	CenterofMass(moli);
      }

      /////////////////////////////////////////////////////////////
      /* Convert coordinates and box size from SI to system unit */
      /////////////////////////////////////////////////////////////
      CoorSI2System();
      //_________________________________________________________//

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
      } 

      /////////////////////
      /* Build cell list */
      /////////////////////
#ifdef CELL_LIST	
      if (init) {
         CL_Init();		
         init	=	0;
      } 
      else 
         CL_Destroy();

      CL_Build();
#endif	/* CELL_LIST */


      //////////////////////
      /* Perform Analysis */
      //////////////////////
      printf("Analysis: start ...\n");
/*
      CalcV();				// calc energy and virial
      if (V_VIRIAL)	Sample_Pressure();

printf("...sampling p2\n");
      SampleP2All();			// sample P2 and P2m and local p2
      Dist_p2();			// put local p2 into distribution

printf("...sampling spherical coordinates\n");
      SampleSpherical();		// sample spherical coordinates
      Dist_Spherical();			// put spherical coord. into distribution
printf("...find nuclei\n");
      Find_Nuclei_p2(1);
*/

//sizeofnucl	=	sizeofnuclp2;
//sizedist	=	sizedistp2;


//      shiftbox(system, nucleus, n);

      /////////////////////////////
      /* Conformation properties */
      /////////////////////////////

      Rg2	=	0;
      R02	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         Rg2	+=	R2_gyration(moli);
	 R02	+=	R2_n2n(moli);
      }
printf("Rg2 = %f\tR02 = %f\n", Rg2/NMOLS, R02/NMOLS);

      /////////////////////////////
      /* Output analysis results */
      /////////////////////////////
      
      printf("Output: start ...\n");
      Print_hst(fhst);		// print out histgram
      CoorSystem2SI();		// convert coordinates and box size back to SI units

      ////////////////////////////////////////////
      /* OUTPUT .car file for VMD visualization */
      ////////////////////////////////////////////
      
      if (carflag) {
         fprintf(fout, "!BIOSYM archive 3\n");
         fprintf(fout, "PBC=ON\n");
         fprintf(fout, "!TIMESTEP %d\n", timestep);
         fprintf(fout, "!DATE %s", asctime(localtime(&t)));
         fprintf(fout, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

	 n	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {

               MolInBox2(moli);
               for (i=0; i<moli->nsites; i++) {

                  if (moli->nuclid[i]>0)	// crystal like particle
                     sprintf(s, "N%d", n++);	// N: blue color in VMD
                  else
	             sprintf(s, "O%d", n++);	// O: red color in VMD

/*
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])		// note nuclid index starts from 1
                     sprintf(s, "M%d", n++);
                  else if (sizeofnucl[moli->nuclid[i]] -MAXSIZE[system] > -3 && MAXSIZE[system]>=10)
                     sprintf(s, "C%d", n++);
	          else if (moli->nuclid[i] >= 1)
                     sprintf(s, "O%d", n++);
                  else
                     sprintf(s, "H%d", n++);
*/
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
      fflush(fout);
      //___________OUTPUT .car file_____________//


      ///////////////////////////////////////////
      // OUTPUT .pdb file for further analysis */
      ///////////////////////////////////////////

      if (pdbflag) {
         fprintf(fpdb, "HEADER: file created from %s on %s", argv[1], asctime(localtime(&t)));
         fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
        	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

	 m		=	0;	// molecule sequence number
         n		=	0;	// atom sequence number
#define SIZECAP	5
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {
               //MolInBox2(moli);

               drawmol	=	0;
               for (i=0; i<moli->nsites; i++) {
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {
	             drawmol	=	1;		// participate the biggest nucleus
		     break;
                  }
	       }
//temp=CenterofMass(moli);
//if (temp.z < -0.25 * BOX[system].lz) {
               m	++; 
               for (i=0; i<moli->nsites; i++) {
                  if (drawmol) {
                     if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {	// nuclid index starts from 1
                        atomname	=	'N';	// N: blue color in VMD
		        fprintf(fdat, " 10");
                     }
                     else {
                        atomname	=	'O';	// O: red color in VMD
                        fprintf(fdat, " 0");
                     }
                  }
                  else {
		     atomname	=	'C';		// C: cyan color in VMD
		     fprintf(fdat," -1");
	          }

                  n	++;
	          fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
                  fprintf(fpdb, "%5d ", n);	// atom number
                  fprintf(fpdb, " %c  ", atomname);	// atom name
                  fprintf(fpdb, " ");		// alternate location indiator
  	          fprintf(fpdb, " C8");		// residue name
	          fprintf(fpdb, " ");		// column 21
                  fprintf(fpdb, " ");		// chain identifier, column 22
	          fprintf(fpdb, "%4d", m);	// residue sequence number, 23-26
	          fprintf(fpdb, " ");		// code for insertion of residues, 27
                  fprintf(fpdb, "   ");		// column 28-30
                  fprintf(fpdb, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z);
                  fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                  fprintf(fpdb, "%5.5s", "");
                  fprintf(fpdb, "\n"); 

                  if (imagen) {			// for image box
                     n	++;
                     fprintf(fpdb, "ATOM  ");
                     fprintf(fpdb, "%5d ", n);
                     fprintf(fpdb, " %c  ", atomname);
                     fprintf(fpdb, " ");
                     fprintf(fpdb, " C8");
	             fprintf(fpdb, " ");
                     fprintf(fpdb, " ");
	             fprintf(fpdb, "%4d", m);
                     fprintf(fpdb, " ");
                     fprintf(fpdb, "   ");
                     fprintf(fpdb, "%8.3f%8.3f%8.3f", 
			moli->p[i].x + imagex, moli->p[i].y+imagey, moli->p[i].z+imagez);
                     fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                     fprintf(fpdb, "%5.5s", "");
                     fprintf(fpdb, "\n"); 
                  }
               } 
//}
            }   
         }
         fprintf(fpdb, "END\n");
         fflush(fpdb);
      }	//pdbflag
      //____________OUTPUT .pdb file___________//

      ////////////////////////////////////////////////////
      /* OUTPUT configuration file for further analysis */
      ////////////////////////////////////////////////////

      if (confflag) { 
         fprintf(fconf, "!TIMESTEP %d\n", timestep);
         fprintf(fconf, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
         for (i=0; i<NSYSTEMS; i++)
            fprintf(fconf, "%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

         for (moli=mol; moli<mol+NMOLS; moli++) {
            fprintf(fconf, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
            //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
            //MolInBox(moli);
            for (i=0; i<moli->nsites; i++) 
               fprintf(fconf, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
         }
         fflush(fconf);
      }
      //____________OUTPUT configuration file___________//
   } 

   /////////////////////////////////////////////
   /* output average center of mass of chains */
   /////////////////////////////////////////////

   if (pdbflag) {
      fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
      fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
 	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

      n		=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         i	=	moli-mol;
         if (ncom[i]>0) {
//if( fabs(com[i].x/ncom[i])<3 && fabs(com[i].y/ncom[i])<3 && fabs(com[i].z/ncom[i])<7 ) {
            n	++;
            atomname	=	'N';
            fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
            fprintf(fpdb, "%5d ", n);		// atom number
            fprintf(fpdb, " %c  ", atomname);	// atom name
            fprintf(fpdb, " ");		// alternate location indiator
            fprintf(fpdb, "   ");	// residue name
            fprintf(fpdb, " ");		// column 21
            fprintf(fpdb, " ");		// chain identifier, column 22
            fprintf(fpdb, "    ");	// residue sequence number, 23-26
            fprintf(fpdb, " ");		// code for insertion of residues, 27
            fprintf(fpdb, "   ");	// column 28-30
            fprintf(fpdb, "%8.3f%8.3f%8.3f", com[i].x/ncom[i], com[i].y/ncom[i], com[i].z/ncom[i]);
            fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
            fprintf(fpdb, "\n"); 
         }
//}
      }
      fprintf(fpdb, "END\n");
      fflush(fpdb);
   }
   //______OUTPUT com of chains to .pdb file__//

   S_PrintAll();		// print out distributions
   corr_norm();			// normalize correlation function
   corr_print();		// print out correlation function
   //radial("print");		// print out radial distribution function

   for (i=0; i<180; i++) {
      printf("%d %d\n", i, orientdist[i]);
   }

   printf("Closing files ...\n");

   fclose(fin);		// close files
   fclose(fhst);
   if (carflag)		fclose(fout);
   if (confflag)    	fclose(fconf);
   if (pdbflag)	     {	fclose(fpdb); fclose(fdat);}

   return	0;
}
                                                                                                                                                                                                                                                                                                                                        src/lists.c                                                                                         0000600 0143352 0000144 00000104416 11565350547 012377  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    lists.c
    author:     Peng Yi at MIT
    date:       April 10, 2007
    purpose:    build verlet list using link-list operations
*/

#define __LISTS_MODULE
#include "lists.h"

/*
   lvector	M[MAXNBOX];		// for debug
   long		NCELLS;			// debug
   cellstruct	*Cell;			// debug
*/

long elem_index(long * array, long arraylength, long element)	//find long type element in array
{
   long		i;
   for (i=0; i<arraylength; i++) {
      if (array[i] == element)
	 return	i;			//if found, return the index
   }
   return	-1;			//if not found, return -1
}
		
long List_Length(liststruct *currentPtr)	//the length of the linked list
{
   long		length;
   length	=	0;

   if (currentPtr == NULL) {
      return	0;
   }
   else {
      while (currentPtr!=NULL) {
	 length		++;
         currentPtr	=	currentPtr->nextPtr;
      }
      return	length;
   }
}

int List_Insert(liststruct ** sPtr, long value)	//insert one element into the list, from small to big
{
   liststruct *	newPtr, * currentPtr, * previousPtr;
   newPtr	=	malloc( sizeof(liststruct));

   if (newPtr != NULL) {			//if memory available
      newPtr->neighbor	=	value;
      newPtr->nextPtr	=	NULL;

      previousPtr	=	NULL;
      currentPtr	=	*sPtr;

      while (currentPtr != NULL && value>currentPtr->neighbor) {
	 previousPtr	=	currentPtr;		//move to ...
	 currentPtr	=	currentPtr->nextPtr;	// ... next list node
      }

      if (currentPtr != NULL && value == currentPtr->neighbor) {	//if same element found in the list
	 printf("Element already in the linked list!\n");
	 return	0;
      }
      else {
         if (previousPtr == NULL) {			//if value is the smallest in the list...
            newPtr->nextPtr	=	*sPtr;		// ... speciall care must be taken to update
	    *sPtr			=	newPtr;		// ... sPtr, because sPtr is supposed to be
         }							// ... pointing to the FIRST element in list
         else {
	    previousPtr->nextPtr	=	newPtr;
	    newPtr->nextPtr	=	currentPtr;
         }
         return	1;
      }
   }
   else {
      printf("Add in list failed.  No memory available!\n");
      return	0;
   }
}

int List_Remove(liststruct ** sPtr, long value)
{
   liststruct * previousPtr, * currentPtr, * tempPtr;

   if (value == (*sPtr)->neighbor) {
      tempPtr	=	*sPtr;
      *sPtr	=	(*sPtr)->nextPtr;
      free(tempPtr);
      return	1;
   }
   else {
      previousPtr	=	*sPtr;
      currentPtr	=	(*sPtr)->nextPtr;

      while (currentPtr != NULL && currentPtr->neighbor != value) {
	 previousPtr	=	currentPtr;
	 currentPtr	=	currentPtr->nextPtr;
      }

      if (currentPtr != NULL) {
	 tempPtr		=	currentPtr;
	 previousPtr->nextPtr	=	currentPtr->nextPtr;
	 free(tempPtr);
	 return	1;
      }
      else {
	 printf("Remove from list failed.  Not found in this list!\n");
	 return 0;
      }
   }
}

void Free_List(liststruct ** sPtr)
{
   liststruct * tempPtr;

   if ( (*sPtr) != NULL) {
      Free_List( &((*sPtr)->nextPtr) );
      tempPtr	=	*sPtr;
      (*sPtr)	=	NULL;
      free(tempPtr);
   }
}

void Printlist(liststruct * currentPtr)		//print a link-list
{
   if (currentPtr == NULL) {
      printf("List is empty!\n");
   }
   else {
      printf("The list is:\n");
      while (currentPtr!=NULL) {
         printf("%d-->", currentPtr->neighbor);
         currentPtr	=	currentPtr->nextPtr;
      }
      printf("NULL \n\n");
   }
}

int List_is_Empty(liststruct * sPtr)
{
   return	sPtr==NULL;
}

#ifdef CELL_LIST

long CL_Neighbor(long i, long j, long n)		// determine neighboring relationship
{
   if (i>j)
      return	j ? (i-j<2) : (i-j<2) || (n-i<2);
   else
      return	i ? (j-i<2) : (j-i<2) || (n-j<2);
} 


void CL_Init()				// Initialize cell lists for all boxes
{
   long		celli, cellj, cellstart, i, ib;
   vector	cellsize;		// cellsize in x, y, z-direction
   long		ix, iy, iz, jx, jy, jz;
   long		nx, ny, nz;
   static long	init=1;

   // Determine the total # of cells in all boxes

   NCELLS	=	0;

   for (ib=0; ib<NBOX; ib++) {

#ifdef VERLET_LIST
      cellsize.x	=	BOX[ib].rv;
#else
      cellsize.x	=	BOX[ib].rc;
#endif
      cellsize.y	=	cellsize.x;
      cellsize.z	=	cellsize.x;

      M[ib].x		=	(int) (BOX[ib].lx/cellsize.x);
      M[ib].y		=	(int) (BOX[ib].ly/cellsize.y);
      M[ib].z		=	(int) (BOX[ib].lz/cellsize.z);

      if (mod(M[ib].x, 2))	M[ib].x	--;		// make M even number
      if (M[ib].x ==2 )		M[ib].x	= 1;		// too few cells

      if (mod(M[ib].y, 2))	M[ib].y	--;
      if (M[ib].y ==2 )		M[ib].y	= 1;

      if (mod(M[ib].z, 2))	M[ib].z	--;
      if (M[ib].z ==2 )		M[ib].z	= 1;

      if (PBC==1)	NCELLS	+=	M[ib].x * M[ib].y * M[ib].z;
   }



   if (!init)
      free(Cell);

   Cell		=	(cellstruct *) calloc (NCELLS, sizeof(cellstruct));
   init		=	0;

   // determine neighboring cells

   celli	=	0;

   for (ib=0; ib<NBOX; ib++) {
      if (PBC==1) {
         nx		=	M[ib].x;
         ny		=	M[ib].y;
         nz		=	M[ib].z;

         cellsize.x	=	BOX[ib].lx / nx;			// recalculate cell size
         cellsize.y	=	BOX[ib].ly / ny;
         cellsize.z	=	BOX[ib].lz / nz;
   
         cellstart	=	celli;
   
         for (iz=0; iz<nz; iz++) {
            for (iy=0; iy<ny; iy++) {
               for (ix=0; ix<nx; ix++) {

                  i		=	celli;
	 	  Cell[i].box	=	ib;

		  Cell[i].p_min.x	=	((double) ix/nx - 0.5) * BOX[ib].lx;
		  Cell[i].p_min.y	=	((double) iy/ny - 0.5) * BOX[ib].ly;
		  Cell[i].p_min.z	=	((double) iz/nz - 0.5) * BOX[ib].lz;

		  Cell[i].center.x	=	((double) (ix+0.5)/nx - 0.5) * BOX[ib].lx;
		  Cell[i].center.y	=	((double) (iy+0.5)/ny - 0.5) * BOX[ib].ly;
		  Cell[i].center.z	=	((double) (iz+0.5)/nz - 0.5) * BOX[ib].lz;

		  Cell[i].p_max.x	=	((double) (ix+1)/nx - 0.5) * BOX[ib].lx;
		  Cell[i].p_max.y	=	((double) (iy+1)/ny - 0.5) * BOX[ib].ly;
		  Cell[i].p_max.z	=	((double) (iz+1)/nz - 0.5) * BOX[ib].lz;

		  // determine neighboring cells

		  Cell[i].nneigh	=	1;
		  Cell[i].neigh[0]	=	Cell + i;

		  cellj		=	cellstart;

		  for (jz=0; jz<nz; jz++)
		     for (jy=0; jy<ny; jy++)
			for (jx=0; jx<nx; jx++) {
			   if ( (ix!=jx || iy!=jy || iz!=jz) && CL_Neighbor(ix, jx, nx) 
				&& CL_Neighbor(iy, jy, ny) && CL_Neighbor(iz, jz, nz) ) { 
			      Cell[i].neigh[Cell[i].nneigh]	=	Cell + cellj;
	                      Cell[i].nneigh	++;
                           }
                           cellj		++;
                        }

		  celli	++;
               }
	    }
         }
      }	// PBC == 1
   }	// for all boxes
   if (celli != NCELLS)
      Exit("list", "CL_Init", "celli!=NCELLS");
}


inline long CL_InsideCell(vector *p, vector *p_min, vector *p_max)
{
   return ( (p->x >= p_min->x) && (p->y >= p_min->y) && (p->z >= p_min->z) &&
	    (p->x < p_max->x) && (p->y < p_max->y) && (p->z < p_max->z) );
}


long CL_Findcell(molstruct *moli, long site, long ib, long PBC)	// find a cell for certain particle
{
   long		i;
   cellstruct	*celli; 
   double	x, y, z;
   vector	pimg;

   celli	=	Cell;
   while ( (celli->box != ib) && celli < Cell+NCELLS )
      celli	++;

   x	=	moli->p[site].x;
   y	=	moli->p[site].y;
   z	=	moli->p[site].z;

   if (PBC==1) {
      celli	+=	((int) (x / BOX[ib].lx * M[ib].x + M[ib].x * 0.5))
		+	((int) (y / BOX[ib].ly * M[ib].y + M[ib].y * 0.5)) * M[ib].x
		+	((int) (z / BOX[ib].lz * M[ib].z + M[ib].z * 0.5)) * M[ib].x * M[ib].y;
   }
/*			// haven't figure out the truncated oct. with lx!=ly!=lz
   if (PBC==2) {
      if ( p.z<0 ) {
         pimg.z	=	p.z	+	0.5 * BOX[ib].lbox;
	 pimg.x	=	p.x	+	(p.x >=0 ? -0.5 : 0.5) * BOX[ib].lbox;
	 pimg.y	=	p.y	+	(p.y >=0 ? -0.5 : 0.5) * BOX[ib].lbox;
      }
      else {
	 pimg	=	p;
      }
      icell	=	((int) (pimg.x * M[ib]/BOX[ib].lbox + M[ib]/2))
			+ ((int) (pimg.y * M[ib]/BOX[ib].lbox + M[ib]/2)) * M[ib]
			+ ((int) (pimg.z * M[ib]/BOX[ib].lbox)) * M[ib] * M[ib]
			+ NCELLS;
   }
*/
   if (celli >= Cell+NCELLS) {
      Exit("list", "CL_Findcell", "celli > NCELLS.");
   }
   return	celli - Cell;
}


void CL_Add(molstruct *moli, long site)
{
   static long		ib, n;
   static cellstruct	*celli;
   static vector	l;
   static double	x, y, z;

   if ( (ib = moli->box) < 0)
      Exit("list", "CL_Add", "ib < 0.");

   l	=	moli->p[site];
   l	=	MapInBox2(&l, PBC, ib);

   celli	=	Cell;
   while ( (celli->box != ib) && celli < Cell+NCELLS )
      celli	++;

   celli	+=	((int) (l.x / BOX[ib].lx * M[ib].x + M[ib].x * 0.5))
		+	((int) (l.y / BOX[ib].ly * M[ib].y + M[ib].y * 0.5)) * M[ib].x
		+	((int) (l.z / BOX[ib].lz * M[ib].z + M[ib].z * 0.5)) * M[ib].x * M[ib].y;

   if (celli-Cell >= NCELLS) {
      printf("site coordinates: %f\t%f\t%f\n", l.x, l.y, l.z);
      printf("box dimension: %f\t%f\t%f\n", BOX[ib].lx, BOX[ib].ly, BOX[ib].lz);
      printf("NCELLS=%d\t, celli-Cell=%d\n", NCELLS, celli-Cell);
      Exit("list", "CL_Add", "celli >= NCELLS!");
   }

   if (celli->nempty) {
      celli->nempty	--;
      n	=	celli->empty[celli->nempty];
   }
   else if (celli->nsites < MAXNCELLSITES) {
      n	=	celli->nsites;
      celli->nsites	++;
   }
   else {
      //fprintf(foutput, "Cell[0].nsites=%d\n", Cell[0].nsites);
      Exit("list", "CL_Add", "MAXNCELLSITES exceeded.");
   }

   celli->mol[n]	=	moli;
   celli->molsite[n]	=	site;

   moli->cell[site]	=	celli;
   moli->cellsite[site]	=	n;
}


void CL_Delete(molstruct *moli, long site)
{
   static cellstruct	*celli;
   static long		i;

   if ( (celli = moli->cell[site]) && celli->nsites ) {
      // celli != NULL and celli has at least one site

      i	=	moli->cellsite[site];				// position of this site in its cell

      celli->empty[celli->nempty]	=	i;		// put atom in stack
      celli->nempty	++;

      celli->mol[i]	=	NULL;
      moli->cell[site]	=	NULL;
      // because of the stack, we here don't update celli->nsites
   }
}


// Relink molecule sites to cell list

void CL_Relink(molstruct *moli)
{
   static cellstruct	*celli;
   static long		i;
   long			j;

   for (i=0; i<moli->nsites; i++)
      if ( (moli->flags[i]>0) && (celli = moli->cell[i]) ) {
         j		=	moli->cellsite[i];
         celli->mol[j]	=	moli;
         celli->molsite[j]	=	i;
      }
}  


void CL_Build()				// place molecule atoms in cells
{
   long		i;
   molstruct	*moli;

   for (i=0; i<NCELLS; i++) {
      Cell[i].nsites	=	0;
      Cell[i].nempty	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         CL_Add(moli, i);
      }
   }
}


void CL_Destroy()			// reverse of CL_Build()
{
   long		i;
   molstruct	*moli;
   cellstruct	*celli;

   for (moli=mol; moli<mol+NMOLS; moli++) {		// detach mols from cell
      for (i=0; i<moli->nsites; i++) {
         moli->cell[i]		=	NULL;
         moli->cellsite[i]	=	-1;
      }
   }
   for (celli=Cell; celli<Cell+NCELLS; celli++) {	// empty cells
      celli->nsites	=	0;
      celli->nempty	=	0;
   }    
}

#endif	/* CELL_LIST */

long NeighborList(molstruct *molm, neighborlist *list)	
{						// build neighbor list for bridging moves
  long			i, k, n;
  long			site, flag, system, reverse;
  double		d2_min, d2_max, d2, alpha_min, alpha_max;
  vector		dr, dr_old, p[7];
  sphere		s[7];
  molstruct		*moln;
  register vector	*r, *q;
#ifdef CELL_LIST
  long			j;
  cellstruct		*cellm, *celli;
#endif

  if ((!(i = molm->nsites-1)))
    return list->n = 0;

  system                = molm->box;
  alpha_min		= (alpha_max = M_PI-type[0].THETA);
  d2_min		= (d2_max = type[0].LSTRETCH);

  for (k=1; (k<NTYPES); ++k)			// Determine bond angle and
  {						// bond length extremes
    if (alpha_min>M_PI-type[k].THETA) alpha_min = M_PI-type[k].THETA;
    if (alpha_max<M_PI-type[k].THETA) alpha_max = M_PI-type[k].THETA;
    if (d2_min>type[k].LSTRETCH) d2_min = type[k].LSTRETCH;
    if (d2_max<type[k].LSTRETCH) d2_max = type[k].LSTRETCH;
  }
//printf("%f %f %f %f\n", alpha_min, alpha_max, d2_min, d2_max);

//  flag                  = molm->fix&7;		// Calculate min and max span
  d2_max                = 4.0*sin(0.5*alpha_min)*d2_max;
  d2_min                = 0.0; //4.0*cos(0.5*alpha_max)*cos(alpha_max)*d2_min;
						// 10% fudge factor (1.21=1.1^2)
  d2_min                *= d2_min/1.21;		// Equal bond angles and
  d2_max                *= 1.21*d2_max;		// lengths assumed
  r                     = molm->p+molm->nsites-1;	// position of this end
  n                     = 0;
#ifdef CELL_LIST_1
  cellm			= molm->cell[i];
  for (i=0; i<cellm->nneigh; ++i)
  {
    celli		= cellm->neigh[i];
    for (j=0; j<celli->nsites; ++j)
    {
      if (moln = celli->mol[j])
      {
        site		= celli->molsite[j];
//	if ((moln!=molm)&&
//	    (flag ? (((moln->fix&7)==0)||((moln->fix&7)==3)):(moln->fix&7)<3)&&
//	    (site>=SITE_GAP)&&(site<moln->nsites-SITE_GAP)&&
//	    (!moln->flags[site]))
        if (moln!=molm)
#else
  for (moln=mol; moln<mol+NMOLS; ++moln)
//    if ((moln->box==molm->box)&&(moln!=molm)&&
//	(flag ? (((moln->fix&7)==0)||((moln->fix&7)==3)) : (moln->fix&7)<3))
    if ((moln->box == molm->box) && moln!=molm)
      for (site=0; site<moln->nsites; ++site)
//        if (!moln->flags[site])
#endif
	{
          dr.x		= (q = moln->p+site)->x - r->x;
	  dr.y		= q->y - r->y;
	  dr.z		= q->z - r->z;
	  dr_old        = dr;
          MapInBox2(&dr, PBC, BOX[system].lbox);
          d2            = dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;	// distance 
          if ((d2>=d2_min)&&(d2<=d2_max))
          {
	    dr.x	-= dr_old.x;		// map to the same box image
	    dr.y	-= dr_old.y;		// as the end site
	    dr.z	-= dr_old.z;
	    for (reverse=0; reverse<2; ++reverse)	// try both directions
	    {					// Hard-coded minimum
//	      if (((reverse ? moln->nsites-site-1 : site)-3>=
//		   (moln->fix&7 ? NMINSITES : E_NMINFREE))&&
//		  !RebridgeSetup(moln, site+1+4*reverse, reverse, p, s))
/*
	      if (((reverse ? moln->nsites-1-site : site)-3 >= E_NMINFREE) &&
		  !RebridgeSetup(moln, site+1+4*reverse, reverse, p, s))
	      {
  		p[0]	= *(r-1);		
  		p[1]	= *r;			// end site
		for (k=2; k<7; ++k)		// Translate to molm frame
		{
		  p[k].x += dr.x;
		  p[k].y += dr.y;
		  p[k].z += dr.z;
		}
		if (Feasible(p, s))		// Check for solutions
		{
*/
	          if (n>=MAXNNEIGHBORS)
        	    Exit("lists", "NeighborList", "MAXNNEIGHBORS exceeded");
		  list->dr[n] = dr;
                  list->mol[n] = moln;
	          list->reverse[n] = reverse;
                  list->site[n++] = site;
/*		}
	      }
*/
	    }
	  }
        }
#ifdef CELL_LIST_1
      }
    }
  }
#endif
  return list->n = n;
}


long DB_NeighborList(molstruct *molm, long sitem, neighborlist *list)
{				// build neighbor list for double bridging moves
  long			i, k, n;
  long			site, flag, system, reverse, sitem2, siten, siten2;
  double		d2_min, d2_max, d2, alpha_min, alpha_max;
  vector		dr, dr_old, pm[7], pn[7];
  sphere		sm[7], sn[7];
  molstruct		*moln;
  register vector	*r, *q;
#ifdef CELL_LIST
  long			j;
  cellstruct		*cellm, *celli;
#endif

  if (sitem<2 || sitem>molm->nsites-3)	// at least two bonds from the end
    return list->n = 0;

  system        = molm->box;

  alpha_max     = M_PI-type[0].THETA;
  alpha_min	= alpha_max;
  d2_max	= type[0].LSTRETCH;
  d2_min	= d2_max;

  for (k=1; (k<NTYPES); ++k)	// Determine bond angle and bond length extremes
  {
    if (alpha_min>M_PI-type[k].THETA) alpha_min = M_PI-type[k].THETA;
    if (alpha_max<M_PI-type[k].THETA) alpha_max = M_PI-type[k].THETA;
    if (d2_min>type[k].LSTRETCH) d2_min = type[k].LSTRETCH;
    if (d2_max<type[k].LSTRETCH) d2_max = type[k].LSTRETCH;
  }

  d2_max                = 4.0*sin(0.5*alpha_min)*d2_max;
  d2_min                = 0.0; 		//4.0*cos(0.5*alpha_max)*cos(alpha_max)*d2_min;
						// 10% fudge factor (1.21=1.1^2)
  d2_min                *= d2_min/1.21;		// Equal bond angles and
  d2_max                *= 1.21*d2_max;		// lengths assumed

  r                     = molm->p+sitem;	// position of this site
  n                     = 0;			// # of neighbors

#ifdef CELL_LIST
  cellm			= molm->cell[sitem];
  for (i=0; i<cellm->nneigh; ++i) {
    celli		= cellm->neigh[i];
    for (j=0; j<celli->nsites; ++j) {
      if (moln = celli->mol[j]) {
        siten		= celli->molsite[j];
        if (moln!=molm)
#else
  for (moln=mol; moln<mol+NMOLS; ++moln)
    if ((moln->box == molm->box) && moln!=molm)
      for (siten=0; siten<moln->nsites; ++siten)
#endif
	{
          if (siten < 2 || siten > moln->nsites-3)	continue;

          q		= moln->p+siten;
          dr.x		= q->x - r->x;
	  dr.y		= q->y - r->y;
	  dr.z		= q->z - r->z;
	  dr_old        = dr;
          MapInBox2(&dr, PBC, BOX[system].lbox);	// minimum distance image
          d2            = dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;	// distance 

          if ((d2>=d2_min)&&(d2<=d2_max))
          {
	    dr.x	-= dr_old.x;		// now dr is a vector pointing from
	    dr.y	-= dr_old.y;		// siten's box to sitem's box
	    dr.z	-= dr_old.z;

            for (reverse=0; reverse<4; reverse++) {
              switch (reverse) {
                case	0:  sitem2 = sitem + 4;	siten2 = siten + 4;	break; 
                case	1:  sitem2 = sitem - 4;	siten2 = siten + 4; 	break; 
                case	2:  sitem2 = sitem + 4;	siten2 = siten - 4; 	break; 
                case	3:  sitem2 = sitem - 4;	siten2 = siten - 4; 	break; 
                default:  break;
              }
              if (sitem2<2 || sitem2 > molm->nsites-3 || siten2<2 || siten2>moln->nsites-3)
                continue;
   
              d2	=	DistSQ(molm->p[sitem2], moln->p[siten2], system);
              if (d2 >= d2_min && d2 <= d2_max ) {
/*
              	RebridgeSetup(molm, MAX(sitem, sitem2)+1, 0, pm, sm);
		RebridgeSetup(moln, MAX(siten, siten2)+1, 0, pn, sn);
                pm[0]	=	V_Add(moln->p+MIN(siten, siten2)-1, &dr);
                pm[1]	=	V_Add(moln->p+MIN(siten, siten2)-1, &dr);
                pn[0]	=	V_Subtr(moln->p+MIN(sitem, sitem2)-1, &dr);
                pn[1]	=	V_Subtr(moln->p+MIN(sitem, sitem2)-1, &dr);

                if (Feasible(pm, sm) && Feasible(pn, sn)) {
*/
                  if (n>=MAXNNEIGHBORS)
                    Exit("lists", "DB_NeighborList", "MAXNNEIGHBORS exceeded");
		  list->dr[n] 		= dr;
                  list->mol[n] 		= moln;
	          list->reverse[n] 	= reverse;
                  list->site[n++] 	= site;
//                }
	      }
	    }
	  }
        }
#ifdef CELL_LIST
      }
    }
  }
#endif
  return list->n = n;
}


long IDR_NeighborList(molstruct *molm, long sitem, neighborlist *list)
{			// build neighbor list for intramolecular double rebridging
  long			i, k, n;
  long			site, flag, system, reverse, sitem2, siten, siten2;
  double		d2_min, d2_max, d2, alpha_min, alpha_max;
  vector		dr, dr_old, pm[7], pn[7];
  sphere		sm[7], sn[7];
  molstruct		*moln;
  register vector	*r, *q;
#ifdef CELL_LIST
  long			j;
  cellstruct		*cellm, *celli;
#endif

  if (sitem<2 || sitem>molm->nsites-3)	// at least two bonds from the end
    return list->n = 0;

  system        = molm->box;

  alpha_max     = M_PI-type[0].THETA;
  alpha_min	= alpha_max;
  d2_max	= type[0].LSTRETCH;
  d2_min	= d2_max;

  for (k=1; (k<NTYPES); ++k)	// Determine bond angle and bond length extremes
  {
    if (alpha_min>M_PI-type[k].THETA) alpha_min = M_PI-type[k].THETA;
    if (alpha_max<M_PI-type[k].THETA) alpha_max = M_PI-type[k].THETA;
    if (d2_min>type[k].LSTRETCH) d2_min = type[k].LSTRETCH;
    if (d2_max<type[k].LSTRETCH) d2_max = type[k].LSTRETCH;
  }

  d2_max                = 4.0*sin(0.5*alpha_min)*d2_max;
  d2_min                = 0.0; //4.0*cos(0.5*alpha_max)*cos(alpha_max)*d2_min;
						// 10% fudge factor (1.21=1.1^2)
  d2_min                *= d2_min/1.21;		// Equal bond angles and
  d2_max                *= 1.21*d2_max;		// lengths assumed

  r                     = molm->p+sitem;	// position of this site
  n                     = 0;			// # of neighbors

  for (siten=0; siten<molm->nsites; ++siten) {
    if (siten < 2 || siten > molm->nsites-3)	continue;
    if (abs(siten-sitem) <=4)			continue;

    q		= molm->p+siten;
    dr.x	= q->x - r->x;
    dr.y	= q->y - r->y;
    dr.z	= q->z - r->z;
    dr_old      = dr;
    MapInBox2(&dr, PBC, BOX[system].lbox);	// minimum distance image
    d2            = dr.x*dr.x+dr.y*dr.y+dr.z*dr.z;	// distance 

    if ((d2>=d2_min)&&(d2<=d2_max)) {
      dr.x	-= dr_old.x;		// now dr is a vector pointing from
      dr.y	-= dr_old.y;		// siten's box to sitem's box
      dr.z	-= dr_old.z;

      for (reverse=0; reverse<2; reverse++) {
        sitem2	=	sitem + 4 - reverse * 8;
        siten2	=	siten + 4 - reverse * 8;
 
        if (sitem2<2 || sitem2 > molm->nsites-3 || siten2<2 || siten2>molm->nsites-3)
          continue;
   
        d2	=	DistSQ(molm->p[sitem2], molm->p[siten2], system);
        if (d2 >= d2_min && d2 <= d2_max ) {
/*
          RebridgeSetup(molm, MAX(sitem, sitem2)+1, 0, pm, sm);
	  RebridgeSetup(molm, MAX(siten, siten2)+1, 0, pn, sn);
          pm[0]	=	V_Add(pn+6, &dr);
          pm[1]	=	V_Add(pn+5, &dr);
          pn[0]	=	V_Subtr(pm+6, &dr);
          pn[1]	=	V_Subtr(pm+5, &dr);

          if (Feasible(pm, sm) && Feasible(pn, sn)) {
*/            if (n>=MAXNNEIGHBORS)
              Exit("lists", "DB_NeighborList", "MAXNNEIGHBORS exceeded");
	      list->dr[n] 		= dr;
              list->mol[n] 		= molm;
	      list->reverse[n] 	= reverse;
              list->site[n++] 	= site;
//            }
          }
        }
      }
  }
  return list->n = n;
}

#ifdef VERLET_LIST

void New_Vlist()		//build a new Verlet list for all particles using array
{
  long		i, j, jj, k;
  double 	r2;
  long		icell, jcell;

  for (i=0; i<NPARTS; i++) {
    part[i].nverlet	=	0;
    part[i].pv		=	part[i].p;	//save particle coordinates
  }
#ifdef CELL_LIST				//given a cell list
  for (i=0; i<NPARTS-1; i++) {
     icell	=	part[i].icell;		//determine cell number

     for (k=0; k<Cell[icell].nneigh; k++) {	//loop over the neighbor cells, including itself
        jcell	=	Cell[icell].neigh[k];

        for (jj=0; jj<Cell[jcell].sites; jj++) {
	   j	=	Cell[jcell].list[jj];

           if (j > i && DistSQ(part[j].p, part[i].p) < Rv2) {
#else
  for (i=0; i<NPARTS-1; i++) {
    for (j=i+1; j<NPARTS; j++) {
      r2	=	DistSQ(part[j].pv, part[i].pv);	//it's pv matters, though now pv=p
      if (r2 < Rv*Rv) {
      {
#endif	/* CELL_LIST */
	part[i].vlist[part[i].nverlet]	=	j;
	part[j].vlist[part[j].nverlet]	=	i;
	part[i].nverlet			++;
	part[j].nverlet			++;
      }
    }
  }}

  for (i=0; i<NPARTS; i++) {
     if (part[i].nverlet > MAXVERLETNEIGH) 
	printf("Verlet neighbor # exceeds limit.\n");
  }
  return;
}

void New_Vlist_LL()		//build a new Verlet list for all particles, using linked list
{
   long		i, j;
   vector 	dp;
   double 	r2;

   for (i=0; i<NPARTS; i++) {
      part[i].pv	=	part[i].p;	//save particle coordinates
      Free_List(&(part[i].startPtr));
      if (part[i].startPtr != NULL)
	 printf("Verlet linked list initialization failed!\n");
   }
   for (i=0; i<NPARTS-1; i++) {
      for (j=i+1; j<NPARTS; j++) {
         r2	=	DistSQ(part[j].p, part[i].p);
         if (r2 < Rv*Rv) {
            List_Insert(&part[i].startPtr, j);	//add j to i's verlet list, using linked list
            List_Insert(&part[j].startPtr, i);
         }
      }
   }
}

void Update_Vlist(long n, vector pv_old, vector pv_new)
{
   long		i, j, k, kcell, index;
   double	r2, r2old;

   part[n].pv	=	pv_new;

#ifdef CELL_LIST					//given a cell list
   for (k=0; k<nneighcellplus; k++) {			//loop over expanded cell neighbors
      kcell	=	neighcellplus[k];

      for (j=0; j<Cell[kcell].sites; j++) {
         i	=	Cell[kcell].list[j];
#else
   for (i=0; i<NPARTS; i++) {
      {					//match the number of loops

#endif	/* CELL_LIST */

      if (i!=n) {
         r2	=	DistSQ(part[i].pv, pv_new);	// It's pv matters in vlist updating
         r2old	=	DistSQ(part[i].pv, pv_old);	// It's pv matters in vlist updating
         if (r2 < Rv*Rv && r2old >= Rv*Rv) {
	    index	=	elem_index(part[n].vlist, part[n].nverlet, i);
	    if (index != -1) {
		   printf("Update_Vlist error, Was already in list, oldr2=%f\tr2=%f\tRv=%f\n", r2old, r2, Rv);
            }
	    else {
	       part[n].vlist[part[n].nverlet]	=	i;
	       part[n].nverlet			++;
	       part[i].vlist[part[i].nverlet]	=	n;
	       part[i].nverlet			++;
	       if (part[n].nverlet > MAXVERLETNEIGH || part[i].nverlet > MAXVERLETNEIGH)
		  printf("Verlet neighbor number exceed limits.\n");
	    } 
	 }
	 else if (r2 >= Rv*Rv && r2old < Rv*Rv) {
	    index	=	elem_index(part[n].vlist, part[n].nverlet, i);
	    if (index == -1) {
	       printf("Update_Vlist error, Was not in list, oldr2=%f\tr2=%f\tRv=%f\n", r2old, r2, Rv);
	    }
	    else {
	       part[n].vlist[index]	=	part[n].vlist[part[n].nverlet-1];
	       part[n].nverlet		--;  
            }		    
	    index	=	elem_index(part[i].vlist, part[i].nverlet, n);
	    if (index == -1) {
	       printf("Update_Vlist error, Was not in list, oldr2=%f\tr2=%f\tRv=%f\n", r2old, r2, Rv);
            }
            else {
	       part[i].vlist[index]	=	part[i].vlist[part[i].nverlet-1];
	       part[i].nverlet		--;  
            }		    
	 }
      }
   }}		//need to match the number of loops
   return;
}

void Update_Vlist2(long n, vector pv_old, vector pv_new)	//update vlist due to displacement of n
{				
   long		i;
   double	r2, r2old;

   part[n].pv	=	pv_new;

   for (i=0; i<NPARTS; i++) {
      if (i != n) {
         r2	=	DistSQ(part[i].pv, pv_new);	// It's pv matters in vlist updating
         r2old	=	DistSQ(part[i].pv, pv_old);	// It's pv matters in vlist updating
         if (r2 < Rv*Rv && r2old >= Rv*Rv) {		//even if it was already in the vlist
            if (!List_Insert( &(part[n].startPtr), i ))
	       printf("%f\t%f\n", r2old, r2);
	    if (!List_Insert( &(part[i].startPtr), n ))
	       printf("%f\t%f\n", r2old, r2);
         }
         else if (r2 >= Rv*Rv && r2old < Rv*Rv) {		//even if it wasn't there in the vlist
            if (!List_Remove( &(part[n].startPtr), i ))
	       printf("%f\t%f\n", r2old, r2);
            if (!List_Remove( &(part[i].startPtr), n ))
	       printf("%f\t%f\n", r2old, r2);
         }
      }
   }
}

#endif	/* VERLET_LIST */

/*******************************************************************************************/

void New_ConnectList()		// build a new connection neighbor list for all particles
{
   molstruct	*moli, *molj;
   long		i, j, k, n, system;

#ifdef CELL_LIST
   cellstruct	*celli, *celln;
#endif

   // Initialization
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         moli->nconn[i]	=	0;

   // Building the connection list
   for (moli=mol; moli<mol+NMOLS; moli++) {
      system	=	moli->box;
      for (i=0; i<moli->nsites; i++) {

#ifdef CELL_LIST
	 celli	=	moli->cell[i];			// cell id of moli->i
	 for (n=0; n<celli->nneigh; n++) {
	    celln	=	celli->neigh[n];	// neighboring cell id

	    for (k=0; k<celln->nsites; k++) {		// sites in neighboring cell
	       if (molj=celln->mol[k]) {
		  j	=	celln->molsite[k];	// pick up one
                  if (moli!=molj || i!=j) {
#else
	 for (molj=mol; molj<mol+NMOLS; molj++) {
	    if (molj->box == system) {
               for (j=0; j<molj->nsites; j++) {
                  if (molj>moli || j>i) {
#endif
		     if (DistSQ(moli->p[i], molj->p[j], system) < Rb*Rb && qlproductSQ(l_of_Ylm, moli, i, molj, j) >= critqlproductSQ ) {

			moli->connmol[i][moli->nconn[i]]	=	molj;
			moli->connsite[i][moli->nconn[i]]	=	j;
			moli->nconn[i]	++;
#ifdef CELL_LIST
#else
			molj->connmol[j][molj->nconn[j]]	=	moli;
			molj->connsite[j][molj->nconn[j]]	=	i;
			molj->nconn[j]	++;
#endif
                        if (moli->nconn[i] > MAXCONNEIGH || molj->nconn[j] > MAXCONNEIGH)
			   Exit("list.c", "ConnectList", "Connected neighbor # exceeds limit.");
		     }
         }  }  }  }
      }
   }
}

#ifdef TEST
void New_Clist()	//build a new connection neighbor list for all particles
{
   long	i, jj, k, kk, jcell, icell;
   for (i=0; i<NPARTS; i++) {			//initialize connected neighbor number
      part[i].nconnect	=	0;
   }

#ifdef VERLET_LIST
   for (i=0; i<NPARTS-1; i++) {
      for (jj=0; jj<part[i].nverlet; jj++) {	//search from within the Verlet list
         k		=	part[i].vlist[jj];
	 {
#elif CELL_LIST
   for (i=0; i<NPARTS-1; i++) {
      icell	=	part[i].icell;

      for (jj=0; jj<Cell[icell].nneigh; jj++) {
         jcell	=	Cell[icell].neigh[jj];

         for (kk=0; kk<Cell[jcell].sites; kk++) {
	    k	=	Cell[jcell].list[kk]; 
#else
   for (i=0; i<NPARTS-1; i++) {
      for (k=i+1; k<NPARTS; k++) {
	 {
#endif
         if ( k>i && DistSQ(part[i].p, part[k].p, part[k].box)<Rb*Rb && qlproductSQ(l_of_Ylm, i, k) >= critqlproductSQ ) {
	    part[i].clist[part[i].nconnect]	=	k;	//being connected neighbor to each other
            part[i].nconnect		++;
	    part[k].clist[part[k].nconnect]	=	i;
	    part[k].nconnect		++;
         }						//k > i to avoid double counting
      }
   }
   }

   for (i=0; i<25; i++) {		//sample distribution prob. of neighbor connections
      cnndist[i]	=	0;
   }
   for (i=0; i<NPARTS; i++) {
      if (part[i].nconnect > MAXCONNEIGH)
	 printf("Connected neighbor # exceeds limit.\n");
      cnndist[part[i].nconnect]	++;	//number distribution of connected neighbors
   }
   return;
} 
#endif /* TEST */

/****************************************************************************************/

#ifdef TEST

void  Update_Clist(long n, vector p_new)	//require Vlist2, update the Clist of all particles in Vlist2 ...
{							// ... because these are the particles that have q6 changes
   long		jj, kk, ll;
   long		i, j, k, index, index2;
   long		kcell, jcell;

   //remember that before we do Update_Clist, Vlist and Qlm have already been updated.

#ifdef VERLET_LIST
   for (kk=0; kk<nverletplus-1; kk++) {		//update particle n first.  We update particle n and other particles ...
						//...in vlistplus separately because when we update n, we need to check 
						//...vlistplus while when we update other particles in vlistplus, we 
						//...only have to check their verlet list.
      j	=	vlistplus[kk];
      index	=	elem_index(part[n].clist, part[n].nconnect, j);

      if (index != -1 && (DistSQ(p_new, part[j].p) >= Rb2 || qlproductSQ(l_of_Ylm, n, j) < critqlproductSQ)) {
						//j WAS in n's clist but no longer qualified
         if (index != part[n].nconnect-1)	//check if j was the last one or the only one in n's clist, optional
            part[n].clist[index]	=	part[n].clist[part[n].nconnect-1];
	 part[n].nconnect	--;

	 index2	=	elem_index(part[j].clist, part[j].nconnect, n);
	 if (index2 != part[j].nconnect-1)
	    part[j].clist[index2]	=	part[j].clist[part[j].nconnect-1];
	 part[j].nconnect	--;
      }
      if (index == -1 && DistSQ(p_new, part[j].p) < Rb*Rb && qlproductSQ(l_of_Ylm, n, j) >= critqlproductSQ) {
         part[n].clist[part[n].nconnect]	=	j;
         part[j].clist[part[j].nconnect]	=	n;
         part[n].nconnect		++;
         part[j].nconnect		++;
      }	
   }

   for (kk=0; kk<nverletplus-1; kk++) {		//update other particles in vlistplus
      i	=	vlistplus[kk];

      for (ll=0; ll<part[i].nverlet; ll++) {
         j	=	part[i].vlist[ll];

	 if (j!=n) {				//n has been dealt with, no need to repeat here
            index	=	elem_index(part[i].clist, part[i].nconnect, j);		//look for j in i's clist

	    if (index != -1 && (DistSQ(part[i].p, part[j].p) >= Rb*Rb || qlproductSQ(l_of_Ylm, i, j) < critqlproductSQ)) {
               if (index != part[i].nconnect-1)
	          part[i].clist[index]	=	part[i].clist[part[i].nconnect-1];
	        part[i].nconnect	--;

	       index2	=	elem_index(part[j].clist, part[j].nconnect, i);
	       if (index2 != part[j].nconnect-1)
	          part[j].clist[index2]	=	part[j].clist[part[j].nconnect-1];
	       part[j].nconnect	--;
	    }
            if (index == -1 && DistSQ(part[i].p, part[j].p) < Rb*Rb && qlproductSQ(l_of_Ylm, i, j) >= critqlproductSQ) {
	       part[i].clist[part[i].nconnect]	=	j;
	       part[j].clist[part[j].nconnect]	=	i;
	       part[i].nconnect		++;
	       part[j].nconnect		++;
	    }	
         }
      }
   }
#elif CELL_LIST
   for (k=0; k<nneighcellplus; k++) {		//loop over n's expanded cell neighbors
      kcell	=	neighcellplus[k];

      for (kk=0; kk<Cell[kcell].sites; kk++) {	//loop over all particles in one neighbor cell
         j	=	Cell[kcell].list[kk];

	 if (j!=n) {
            index	=	elem_index(part[n].clist, part[n].nconnect, j);

            if (index != -1 && (DistSQ(p_new, part[j].p) >= Rb2 
			|| qlproductSQ(l_of_Ylm, n, j) < critqlproductSQ)) {
						//j WAS in n's clist but no longer qualified
               part[n].clist[index]	=	part[n].clist[part[n].nconnect-1];
               part[n].nconnect	--;

               index2	=	elem_index(part[j].clist, part[j].nconnect, n);
	       part[j].clist[index2]	=	part[j].clist[part[j].nconnect-1];
	       part[j].nconnect	--;
            }
            if (index == -1 && DistSQ(p_new, part[j].p) < Rb*Rb 
			&& qlproductSQ(l_of_Ylm, n, j) >= critqlproductSQ) {
               part[n].clist[part[n].nconnect]	=	j;
               part[j].clist[part[j].nconnect]	=	n;
               part[n].nconnect		++;
               part[j].nconnect		++;
            }	
         }
      }
   }
   for (k=0; k<nneighcellplus; k++) {		//update other particles in expanded neighbor cells
      kcell	=	neighcellplus[k]; 

      for (kk=0; kk<Cell[kcell].sites; kk++) {
         i	=	Cell[kcell].list[kk];
         if (i!=n) {				//pick up one particle other than particle n

            for (jj=0; jj<Cell[kcell].nneigh; jj++) {
               jcell	=	Cell[kcell].neigh[jj];

	       for (ll=0; ll<Cell[jcell].sites; ll++) {
		  j	=	Cell[jcell].list[ll];
		  if (j!=i && j!=n) {

                     index	=	elem_index(part[i].clist, part[i].nconnect, j);		//look for j in i's clist

                     if (index != -1 && (DistSQ(part[i].p, part[j].p) >= Rb*Rb 
				|| qlproductSQ(l_of_Ylm, i, j) < critqlproductSQ)) {
	          
                        part[i].clist[index]	=	part[i].clist[part[i].nconnect-1];
            	        part[i].nconnect	--;

	                index2	=	elem_index(part[j].clist, part[j].nconnect, i);
	                part[j].clist[index2]	=	part[j].clist[part[j].nconnect-1];
                        part[j].nconnect	--;
                     }
                     if (index == -1 && DistSQ(part[i].p, part[j].p) < Rb*Rb 
				&& qlproductSQ(l_of_Ylm, i, j) >= critqlproductSQ) {
                        part[i].clist[part[i].nconnect]	=	j;
                        part[j].clist[part[j].nconnect]	=	i;
                        part[i].nconnect		++;
                        part[j].nconnect		++;
                     }	
                  }
               }
            }
         }
      }
   }
#endif

   for (i=0; i<25; i++) {
      cnndist[i]	=	0;
   }
   for (i=0; i<NPARTS; i++) {
      if (part[i].nconnect > MAXCONNEIGH)
         printf("Connected neighbor # exceeds limit.\n");
      cnndist[part[i].nconnect]	++;	//number distribution of connected neighbors
   }
   return;
}

#endif /* TEST */
/**********************************************************************************************************/
                                                                                                                                                                                                                                                  src/lmp_conf.c                                                                                      0000600 0143352 0000144 00000050055 11577465122 013034  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	lmp_conf.c
	author:		Peng Yi at MIT
	date:		March 2, 2011
	purpose:	Read lammps dump file, do quick analysis related to chain conformation
			or box shape, etc.  It is meant to do quick analysis, so complicated 
			analysis will still be done using lammps2hst.
	note:		
			require setup file
*/

#define __MAIN_PROGRAM
#include "header.h"

#define DEBUG		0
#define NEWVERSION	1
#define Rcut2		9	// 3sigma squared

#include "correlation.h"

long		timestep;
double		rshift2;		// shift of the biggest nucleus
long		nsitesp2;	// # of sites with p2 greater than a threshold
long		nmolsp2;
double		Rg2, Ro2, DRg2, DRo2;		// radius of gyration, n2n distance 
long		nRg2, nRo2;
double		nchains;

void Print_hst(FILE *fPtr)
{
   long		i;
   static long	init = 1;

   if (init) {
      fprintf(fPtr, "timestep ");
      fprintf(fPtr, "volume ");
      fprintf(fPtr, "lx ");
      fprintf(fPtr, "ly ");
      fprintf(fPtr, "lz ");
      fprintf(fPtr, "Rg2 ");		// average radius of gyration squared
      fprintf(fPtr, "Ro2 ");		// average n2n distance squared
      fprintf(fPtr, "nchains");
      fprintf(fPtr, "\n");
      init	=	0;
   }
   fprintf(fPtr, "%-6d ", timestep);
   fprintf(fPtr, "%8.4f ", BOX[0].vol);
   fprintf(fPtr, "%8.4f ", BOX[0].lx);
   fprintf(fPtr, "%8.4f ", BOX[0].ly);
   fprintf(fPtr, "%8.4f ", BOX[0].lz);
   fprintf(fPtr, "%6.4f ", Rg2);
   fprintf(fPtr, "%6.4f ", Ro2);
   fprintf(fPtr, "%6.4f", nchains);
   fprintf(fPtr, "\n");
   fflush(fPtr);
}

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   long		i, j, k, system, m, n, nuclid;
   long		nsites, record;
   long		nx, ny, nz, id, siteid, molid, type;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo;			// box upper boundary
   molstruct	*moli, *molj;
   double	temp1, temp2, temp3;		// temperory variables
   double	R2;

   char		filein[255], s[80], ff[80], par[80], dummy[255];
   FILE		*fin, *fhst, *fout, *fconf, *fpdb, *fdat;
   long		LENGTH, accum, confflag=0, carflag=0, pdbflag=0, polydisperse=0, drawmol;
   char		atomname;
   static long	init=1;

   vector		con;	// center of nucleus
   static vector	rO;	// original position of nucleus

   vector	chainface;
   double	orient;
   long		orientdist[180];	// chain orientation distribution

   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];

   vector	com[MAXNMOLS], temp;		// average center of mass
   long		ncom[MAXNMOLS];

   double	imagex, imagey, imagez;		// variables for creating pbc images
   long		imagen;

   if (argc<2) {
      printf("lmp_conf (c) 2011 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tlmps2hst [-option] [x= 0.1 y= 0.1 z= 0.1 n= 1] lammpsdumpfile\n\n");
      printf("Notes:\n");
      printf("\t* -option = \n");
      printf("\t* -poly: polydisperse system, has molid in dump file\n");
      printf("\t* -conf: configuration file output\n");
      printf("\t* x= y= z=: duplicate the system and shift unit vector\n");
      printf("\t* n=: multiple of shift vector\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-poly"))
         polydisperse	=	1;
      else if (samestr(par, "-conf")) 
         confflag	=	1;
      else if (samestr(par, "x=")) 
         imagex	=	atof(argv[i+1]);
      else if (samestr(par, "y=")) 
         imagey	=	atof(argv[i+1]);
      else if (samestr(par, "z=")) 
         imagez	=	atof(argv[i+1]);
      else if (samestr(par, "n=")) 
         imagen	=	atol(argv[i+1]);
   }
   strcpy(filein, argv[argc-1]);		// input file

   fin		=	fopen(filein, "r");
   fhst		=	fopen("lmp_conf.hst", "w");
   if (confflag)
      fconf	=	fopen("lammps.conf", "w");

   freopen("lammps.out", "w", stdout);	// redirect standard output stream to a file

   if (DEBUG)	printf("Initialization: start ... \n");

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   InitUnits();				// initialize units for calculation
   InitForcefield();			// initialize Lennard-Jones potential mixing rule

   system	=	0;		// for now 2/14/08, only one system

   InitSample();			// initialize sampling

   for (i=0; i<180; i++)		// initialize chain rotation distribution
      orientdist[i]	=	0;

   while (!feof(fin)) {

      if (DEBUG)	printf("Read one configuration: start ...\n");

      /* Read in one configuration from dump file */

      if (!fgets(dummy, sizeof(dummy), fin))	break;			// end of file
      fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &xlo, &xhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &ylo, &yhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &zlo, &zhi);	fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);

      BOX[system].lx	=	xhi-xlo;
      BOX[system].ly	=	yhi-ylo;
      BOX[system].lz	=	zhi-zlo;

      LENGTH		=	NSITES/NMOLS;	// monodisperse system for now (4/26/2008)
      accum		=	0;

      for (i=0; i<nsites; i++) {
         fscanf(fin, "%ld", &id);
         fscanf(fin, "%ld", &molid);
         fscanf(fin, "%ld", &type);
         fscanf(fin, "%lf%lf%lf %lf%lf%lf %ld%ld%ld", &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
         fgets(dummy, sizeof(dummy), fin);

         molid	--;				// note that lammps id starts from 1, not 0
         siteid	=	(id-1) % LENGTH;

         mol[molid].box		=	system;		// for now, only one system
         mol[molid].nsites	=	LENGTH;		// for now, Jan/10/2010
         mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx);
         mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly);
         mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);
         mol[molid].type[siteid]=	type - 1;	// -1 because lammps index starts from 1
      }

      for (system=0; system<NSYSTEMS; system++) {
         NMols[system]	=	0;
         NSites[system]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;			// total # of mols in certain system
            NSites[system]	+=	moli->nsites;	// total # of sites in certain system
         }
      }

      for (moli=mol; moli<mol+NMOLS; moli++) { 		
         for (i=0; i<moli->nsites; i++)  {
            moli->flags[i]	=	1;		// activate all the sites on this processor
            moli->parent[i]	=	i-1;		// initialize parent site
         }
         moli->flip		=	0;		// flip to the original direction
         moli->origin		=	CenterofMass(moli);
      }

      ///////////////////////////////////////////////////////////////////////////
      /* Start: Calculate the average position of center of mass of each chain */
      ///////////////////////////////////////////////////////////////////////////
      for (moli=mol; moli<mol+NMOLS; moli++) {
      	  temp	=	CenterofMass(moli);
          if (fabs(temp.x) < 0.45 * BOX[moli->box].lx && 
		fabs(temp.y) < 0.45 * BOX[moli->box].ly &&
		fabs(temp.z) < 0.45 * BOX[moli->box].lz) {
	     com[moli-mol]	=	V_Add(&temp, com+(moli-mol));
             ncom[moli-mol]	++;
          }
      }

      /////////////////////////////////////////////////////////////
      /* Convert coordinates and box size from SI to system unit */
      /////////////////////////////////////////////////////////////
      CoorSI2System();
      //_________________________________________________________//

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
      } 

      /////////////////////
      /* Build cell list */
      /////////////////////
#ifdef CELL_LIST	
      if (init) {
         CL_Init();		
         init	=	0;
      } 
      else 
         CL_Destroy();

      CL_Build();
#endif	/* CELL_LIST */


      //////////////////////
      /* Perform Analysis */
      //////////////////////
      if (DEBUG)	printf("Analysis: start ...\n");
/*
      CalcV();				// calc energy and virial
      if (V_VIRIAL)	Sample_Pressure();

      if (DEBUG)	printf("...sampling p2\n");

      SampleP2All();			// sample P2 and P2m and local p2
      Dist_p2();			// put local p2 into distribution

      if (DEBUG)	printf("...sampling spherical coordinates\n");

      SampleSpherical();		// sample spherical coordinates
      Dist_Spherical();			// put spherical coord. into distribution
*/
      /////////////////////////////////////////////////////////////////////
      /* Compute chain orientation distribution, test the rotator phase. */
      /////////////////////////////////////////////////////////////////////

      if (DEBUG)	printf("...computer chain orientation\n");

      for (moli=mol; moli<mol+NMOLS; moli++) {
         V_Null(&chainface);

         for (i=0; i<moli->nsites; i++) {
            if (mod(i,2)) {
	       chainface	=	V_Add(&chainface, moli->p+i);
            }
	    else {
	       chainface	=	V_Subtr(&chainface, moli->p+i);
	    }
	 }
	 orient	=	atan2(chainface.y, chainface.x)	+ M_PI;
         //temp=CenterofMass(moli);
         //if (temp.z < -0.25 * BOX[system].lz && temp.z <0) 
	 orientdist[(int)(orient/(2*M_PI)*36)]	++;
      }

      //////////////////////////////////////////////
      /* Chain conformation statistics (3/2/2011) */
      //////////////////////////////////////////////

      Rg2 	= 	0.0;	DRg2	=	0.0;      nRg2	=	0.0;
      Ro2	=	0.0;    DRo2	=	0.0;      nRo2	=	0.0;

      for (moli=mol; moli<mol+NMOLS; moli++) {
         Rg2	+=	R2_gyration(moli);
	 Ro2	+=	R2_n2n(moli);
	 nRg2	++;
	 nRo2	++;
      }
      Rg2	/=	nRg2;
      Ro2	/=	nRo2;

      /////////////////////////////////////////////////////////////////////////////////
      /* calculate on average how many chains appear in the neighborhood of one bead */
      /////////////////////////////////////////////////////////////////////////////////
 
      nchains	=	0.0;
      system	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         for (i=0; i<moli->nsites; i++) {

             for (molj=mol; molj<mol+NMOLS; molj++) {
                if (moli==molj)	continue;

                for (j=0; j<molj->nsites; j++) {

                   R2	=	DistSQ(moli->p[i], molj->p[j], system);
                   if (R2 < Rcut2) {
                      nchains	+=	1;
                      break;
                   }
                }
             }
         }
      }
      nchains	/=	NSITES;

      /////////////////////////////////////////////////////////////////////
      /* calculate the # of sites with local p2 greater than a threshold */
      /////////////////////////////////////////////////////////////////////

      nsitesp2	=	0;
      nmolsp2	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
//       nsitesp2	=	0;
         for (i=0; i<moli->nsites; i++) {
            if (moli->p2[i] > 0.4)
               nsitesp2	++;
         }
//            if (nsitesp2>1 && nsitesp2<6)
//	       nmolsp2	++;
      }

      /////////////////////////////////////////
      /* Correlation and Radial distribution */
      /////////////////////////////////////////

      if (DEBUG)	printf("Calculate correlation: start...\n");

      correlation();			// calculate correlation, in system units
      if (!(timestep%IRADIAL)) {
         radial("sample");		// sample radial distribution function
      }

//      shiftbox(system, nucleus, n);

      /////////////////////////////
      /* Output analysis results */
      /////////////////////////////
      
      if (DEBUG)	printf("Output: start ...\n");

      Print_hst(fhst);		// print out histgram
      CoorSystem2SI();		// convert coordinates and box size back to SI units

      ////////////////////////////////////////////
      /* OUTPUT .car file for VMD visualization */
      ////////////////////////////////////////////
      
      if (carflag) {
         fprintf(fout, "!BIOSYM archive 3\n");
         fprintf(fout, "PBC=ON\n");
         fprintf(fout, "!TIMESTEP %d\n", timestep);
         fprintf(fout, "!DATE %s", asctime(localtime(&t)));
         fprintf(fout, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

	 n	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {

               MolInBox2(moli);
               for (i=0; i<moli->nsites; i++) {

                  if (moli->nuclid[i]>0)	// crystal like particle
                     sprintf(s, "N%d", n++);	// N: blue color in VMD
                  else
	             sprintf(s, "O%d", n++);	// O: red color in VMD

/*
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])		// note nuclid index starts from 1
                     sprintf(s, "M%d", n++);
                  else if (sizeofnucl[moli->nuclid[i]] -MAXSIZE[system] > -3 && MAXSIZE[system]>=10)
                     sprintf(s, "C%d", n++);
	          else if (moli->nuclid[i] >= 1)
                     sprintf(s, "O%d", n++);
                  else
                     sprintf(s, "H%d", n++);
*/
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
      fflush(fout);
      //___________OUTPUT .car file_____________//


      ///////////////////////////////////////////
      // OUTPUT .pdb file for further analysis */
      ///////////////////////////////////////////

      if (pdbflag) {
         fprintf(fpdb, "HEADER: file created from %s on %s", argv[1], asctime(localtime(&t)));
         fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
        	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

	 m		=	0;	// molecule sequence number
         n		=	0;	// atom sequence number
#define SIZECAP	5
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {
               //MolInBox2(moli);

               drawmol	=	0;
               for (i=0; i<moli->nsites; i++) {
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {
	             drawmol	=	1;		// participate the biggest nucleus
		     break;
                  }
	       }
//temp=CenterofMass(moli);
//if (temp.z < -0.25 * BOX[system].lz) {
               m	++; 
               for (i=0; i<moli->nsites; i++) {
                  if (drawmol) {
                     if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {	// nuclid index starts from 1
                        atomname	=	'N';	// N: blue color in VMD
		        fprintf(fdat, " 10");
                     }
                     else {
                        atomname	=	'O';	// O: red color in VMD
                        fprintf(fdat, " 0");
                     }
                  }
                  else {
		     atomname	=	'C';		// C: cyan color in VMD
		     fprintf(fdat," -1");
	          }

                  n	++;
	          fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
                  fprintf(fpdb, "%5d ", n);	// atom number
                  fprintf(fpdb, " %c  ", atomname);	// atom name
                  fprintf(fpdb, " ");		// alternate location indiator
  	          fprintf(fpdb, " C8");		// residue name
	          fprintf(fpdb, " ");		// column 21
                  fprintf(fpdb, " ");		// chain identifier, column 22
	          fprintf(fpdb, "%4d", m);	// residue sequence number, 23-26
	          fprintf(fpdb, " ");		// code for insertion of residues, 27
                  fprintf(fpdb, "   ");		// column 28-30
                  fprintf(fpdb, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z);
                  fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                  fprintf(fpdb, "%5.5s", "");
                  fprintf(fpdb, "\n"); 

                  if (imagen) {			// for image box
                     n	++;
                     fprintf(fpdb, "ATOM  ");
                     fprintf(fpdb, "%5d ", n);
                     fprintf(fpdb, " %c  ", atomname);
                     fprintf(fpdb, " ");
                     fprintf(fpdb, " C8");
	             fprintf(fpdb, " ");
                     fprintf(fpdb, " ");
	             fprintf(fpdb, "%4d", m);
                     fprintf(fpdb, " ");
                     fprintf(fpdb, "   ");
                     fprintf(fpdb, "%8.3f%8.3f%8.3f", 
			moli->p[i].x + imagex, moli->p[i].y+imagey, moli->p[i].z+imagez);
                     fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                     fprintf(fpdb, "%5.5s", "");
                     fprintf(fpdb, "\n"); 
                  }
               } 
//}
            }   
         }
         fprintf(fpdb, "END\n");
         fflush(fpdb);
      }	//pdbflag
      //____________OUTPUT .pdb file___________//

      ////////////////////////////////////////////////////
      /* OUTPUT configuration file for further analysis */
      ////////////////////////////////////////////////////

      if (confflag) { 
         fprintf(fconf, "!TIMESTEP %d\n", timestep);
         fprintf(fconf, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
         for (i=0; i<NSYSTEMS; i++)
            fprintf(fconf, "%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

         for (moli=mol; moli<mol+NMOLS; moli++) {
            fprintf(fconf, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
            //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
            //MolInBox(moli);
            for (i=0; i<moli->nsites; i++) 
               fprintf(fconf, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
         }
         fflush(fconf);
      }
      //____________OUTPUT configuration file___________//
   } 

   /////////////////////////////////////////////
   /* output average center of mass of chains */
   /////////////////////////////////////////////

   if (pdbflag) {
      fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
      fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
 	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

      n		=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         i	=	moli-mol;
         if (ncom[i]>0) {
//if( fabs(com[i].x/ncom[i])<3 && fabs(com[i].y/ncom[i])<3 && fabs(com[i].z/ncom[i])<7 ) {
            n	++;
            atomname	=	'N';
            fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
            fprintf(fpdb, "%5d ", n);		// atom number
            fprintf(fpdb, " %c  ", atomname);	// atom name
            fprintf(fpdb, " ");		// alternate location indiator
            fprintf(fpdb, "   ");	// residue name
            fprintf(fpdb, " ");		// column 21
            fprintf(fpdb, " ");		// chain identifier, column 22
            fprintf(fpdb, "    ");	// residue sequence number, 23-26
            fprintf(fpdb, " ");		// code for insertion of residues, 27
            fprintf(fpdb, "   ");	// column 28-30
            fprintf(fpdb, "%8.3f%8.3f%8.3f", com[i].x/ncom[i], com[i].y/ncom[i], com[i].z/ncom[i]);
            fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
            fprintf(fpdb, "\n"); 
         }
//}
      }
      fprintf(fpdb, "END\n");
      fflush(fpdb);
   }
   //______OUTPUT com of chains to .pdb file__//

   S_PrintAll();		// print out distributions
   corr_norm();			// normalize correlation function
   corr_print();		// print out correlation function
   //radial("print");		// print out radial distribution function

   for (i=0; i<180; i++) {
      printf("%d %d\n", i, orientdist[i]);
   }

   if (DEBUG)	printf("Closing files ...\n");

   fclose(fin);		// close files
   fclose(fhst);
   if (carflag)		fclose(fout);
   if (confflag)    	fclose(fconf);
   if (pdbflag)	     {	fclose(fpdb); fclose(fdat);}

   fflush(stdout);
   fclose(stdout);

   return	0;
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   src/main1.c                                                                                         0000600 0143352 0000144 00000024674 10714730275 012250  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    main.c
    author:     Peng Yi at MIT
    date:       October 23, 2006
    purpose:    Single processor NPT ensemble simulation program
*/

#define __MAIN_PROGRAM
#include "header.h"

#define VERSION		"1.1"

int main(int argc, char * argv[])
{
   long		i;
   int		k;
   long		REALNCYCLE;

   long		jj;	//debug
   liststruct * currentPtr;
   long		j, max;

   tim=(int *)malloc(sizeof(int));     	//random number generator
   seed=(long *)malloc(sizeof(long));
   *tim=(int)time(NULL);
   *seed= -1*(*tim);           		//seed must start out as a negative long

   InitAll(argv);

   radial(0);				//measure radial distribution function
   radial(1);
   radial(2);
   Print_gr();

   Visualize(0);

   if (Stage == 1 || Stage ==	2) {			//1: melting process
							//2: quench the melt, bring system to desired dynvar window
      if (Stage==1)	kT	=	kT_high;
      else		kT	=	kT_low;

      Init_Sample();
      ResetAcceptance();
      radial(0);

      for (counter=0; counter<NCYCLE; counter++) {
         Cycle();					//carry out MC moves
         Adjust_Stepsize();				//adjust MC move step sizes
         Sample_Energy();

	 fprintf(fhst, "%d\t%f\t%6.4f\t", counter, VSYSTEM.tot, Ql);
	 fprintf(fhst, "%d\t%d\t%d\t", MAXSIZE, Xtal, Nnucl);
	 fprintf(fhst, "%4.2f\t%4.2f\t%5.3f\t%5.3f\t", LBOX, Rv, DRMAX, DLMAX);
#ifdef VERLET_LIST
	 fprintf(fhst, "%d\t", maxnverlet);
#endif

         if (Stage==2) {				//more output for quenching process, optional
            if (dynvar==2 && counter>=NCYCLE/2) {
               for (i=0; i<Qlbins; i++)					//sample prob. distribution p(Ql)
                  pQ[i]	*=	(double) (counter-NCYCLE/2) / (counter-NCYCLE/2+1);
               pQ[(int) ((Ql - Qlmiddle + Qlbins/2.0 * Qlbinsize)/Qlbinsize)] += 1.0/(counter-NCYCLE/2+1);
            }
            else if (dynvar==1 && counter>=NCYCLE/2) {
               for (i=0; i<NMAXbins; i++)				//sample prob. distribution p(NMAX)
                  p[i]	*=	(double) (counter-NCYCLE/2) / (counter-NCYCLE/2+1);
               p[(int) ((MAXSIZE - NMAXmiddle + NMAXbins/2.0 * NMAXbinsize)/NMAXbinsize)] += 1.0/(counter-NCYCLE/2+1);
            }

	    if (dynvar==2) {
               for (i=0; i<Qlbins; i++)
                  fprintf(fhst, "%5.4f\t", pQ[i]);
               for (i=0; i<Qlbins; i++)
                  fprintf(fhst, "%6.3f\t", etaQ[i]);
	       for (i=0; i<Qlbins; i++)
	          fprintf(fhst, "%6.3f\t", (pQ[i]>ZERO ? -log(pQ[i]) : 10000));	//Stage 2 MC moves don't depend on eta
            }
	    else if (dynvar==1) {
               for (i=0; i<NMAXbins; i++)
                  fprintf(fhst, "%5.4f\t", p[i]);
               for (i=0; i<NMAXbins; i++)
                  fprintf(fhst, "%6.3f\t", eta[i]);
	       for (i=0; i<NMAXbins; i++)
	          fprintf(fhst, "%6.3f\t", (p[i]>ZERO) ? -log(p[i]) : 10000);
	    }
	 }
         fprintf(fhst, "\n");

         if (counter >= NCYCLE*3/4) {			//after equilibrium is achieved	
            if (mod(counter+1, NCYCLE/4/NRADIAL)==0) {	//total NRADIAL samplings of rdf and take average
               radial(1);
	    }
         }
	 if (mod(counter+1, ITAPE)==0) {		//write out configuration file every ITAPE cycles
	    Write_Conf();
	 }
      }   
      radial(2);
      Print_gr();
      Sample_Done();
   }

   if (Stage ==	3) {			//system in desired dynvar region, do trial runs and production run
      kT	=	kT_low;
//      Update_Eta(0);			//initialize p and eta	;
      radial(0);			//initialize rdf sampling
      trial		=	0;	//in case we need to reduce damping, we use trial
      REALNCYCLE	=	NCYCLE;

      for (k=0; k<=TRIALRUN; k++) {	//make trial runs TRIALRUN times and production run once
         trial	++;
         Init_Sample();
	 
	 if (k==TRIALRUN) {		//calculate production run length
	    REALNCYCLE	=	NCYCLE * PROD;
         }
         for (counter=0; counter<REALNCYCLE; counter++) {
            Cycle();
            Adjust_Stepsize();
//            Sample_Energy();

	    if (dynvar == 2) {
	       if (k<TRIALRUN && counter >= REALNCYCLE/2) {	//sample prob. dist. after equil. achieved
                  for (i=0; i<Qlbins; i++)				//sample prob. distribution p(Ql)
                     pQ[i]	*=	(double)(counter-REALNCYCLE/2) / (counter-REALNCYCLE/2+1);
                  pQ[(int) ((Ql - Qlmiddle + Qlbins/2.0 * Qlbinsize)/Qlbinsize)] += 1.0/(counter-REALNCYCLE/2+1);
	       }
	       else if (k==TRIALRUN) {		//no eta update after last trial run, so prod. run starts from equil.
                  for (i=0; i<Qlbins; i++)
                     pQ[i]	*=	(double)(counter) / (counter+1);
                  pQ[(int) ((Ql - Qlmiddle + Qlbins/2.0 * Qlbinsize)/Qlbinsize)] += 1.0/(counter+1);	
	       }
	    }
	    else if (dynvar == 1) {
               if (k<TRIALRUN && counter >= REALNCYCLE/2) {
                  for (i=0; i<NMAXbins; i++)				//sample prob. distribution p(NMAX)
                     p[i]	*=	(double) (counter-REALNCYCLE/2) / (counter-REALNCYCLE/2+1);
                  p[(int) ((MAXSIZE - NMAXmiddle + NMAXbins/2.0 * NMAXbinsize)/NMAXbinsize)] += 1.0/(counter-REALNCYCLE/2+1);
	       }
	       else if (k==TRIALRUN) {
                  for (i=0; i<NMAXbins; i++)				//sample prob. distribution p(NMAX)
                     p[i]	*=	(double)counter / (counter+1);
                  p[(int) ((MAXSIZE - NMAXmiddle + NMAXbins/2.0 * NMAXbinsize)/NMAXbinsize)] += 1.0/(counter+1);
	       }
            }


	    fprintf(fhst, "%d\t%f\t%6.4f\t", counter + k*NCYCLE, VSYSTEM.tot, Ql);
	    fprintf(fhst, "%d\t%d\t%d\t", MAXSIZE, Xtal, Nnucl);
	    fprintf(fhst, "%4.2f\t%4.2f\t%5.3f\t%5.3f\t", LBOX, Rv, DRMAX, DLMAX);
#ifdef VERLET_LIST
	    fprintf(fhst, "%d\t", maxnverlet);
#endif
	    if (dynvar == 2) {
               for (i=0; i<Qlbins; i++)
	          fprintf(fhst, "%5.4f\t", pQ[i]);
               for (i=0; i<Qlbins; i++)
	          fprintf(fhst, "%6.3f\t", etaQ[i]);
	       for (i=0; i<Qlbins; i++)
		  fprintf(fhst, "%6.3f\t", (pQ[i]>ZERO ? -log(pQ[i])+etaQ[i] : 10000));
            }
            else if (dynvar == 1) {
	       for (i=0; i<NMAXbins; i++)
                  fprintf(fhst, "%5.4f\t", p[i]);
               for (i=0; i<NMAXbins; i++)
	          fprintf(fhst, "%6.3f\t", eta[i]);
	       for (i=0; i<NMAXbins; i++)
		  fprintf(fhst, "%6.3f\t", (p[i]>ZERO) ? -log(p[i])+eta[i] : 10000);
            }
            fprintf(fhst, "\n");


            if (k==TRIALRUN && counter >= REALNCYCLE*3/4) {			//after equilibrium is reached	
	       if (mod(counter, REALNCYCLE/4/NRADIAL)==0) {	//total NRADIAL samplings of rdf and take average
                  radial(1);
	       }
	       if (mod(counter, REALNCYCLE/4/NGSAMPLE)==0) {	//sample free energy and do statistics
		  //Sample_G( (counter-REALNCYCLE*3/4)/(REALNCYCLE/4/NGSAMPLE) );
	       }
            }
	    if (k==TRIALRUN && mod(counter+1, ITAPE)==0) {	//write out configuration file every ITAPE cycles
	       Write_Conf();
	    }
         }
         if (k<TRIALRUN-1) {			//no update eta between last trial run and production run
            Update_Eta(1);			//update eta
         }
	 Sample_Done();
      }
      radial(2);
      Print_gr();
   }
   
/*
  // 1.0 Melt the system
  kT	=	kT_high;	//set temperature high, start melting
  Update_Eta(DYNVAR, 0);

  counter	=	0;
  while(PhaseNotDone(DYNVAR,0)) {
    counter	++;
    if (ran1(seed) < 0.08 )	// possibility of mcvol = 8%
      mcvol();
    else
      mcmove();

    if (DYNVAR == "MAXSIZE") {		//update the histogram
      p[MAXSIZE]	++;
    }
    else if (DYNVAR == "Ql") {
      pQ[Qstatefinder(Ql)]	++;
    }

    if (mod(counter, 200) == 0) {
      Adjust_Stepsize();	//adjust MC move step size
      //printf("DLMAX=%f\tDRMAX=%f\n",DLMAX, DRMAX);
      //printf("Ql=%f\tQstate=%d\tMAXSIZE=%d\n", Ql, Qstatefinder(Ql), MAXSIZE);
    }
    if (mod(counter, 1000) == 0) {
      //printf("counter=%d\tMAXSIZE=%d\tVol=%f\n",counter, MAXSIZE, VOL);
      //Print_Histogram();
    }
  }
  Visualize(0);
  radial(0);
  radial(1);
  Print_gr();
  
  // 2.0 Quench the system
  Update_Eta(DYNVAR, 1);

  for (trial=0; trial<TRIALRUN; trial++) {	//make trial runs to determine weighting factors

    //printf("Maxsize=%d\tsizelower=%d\tsizeupper=%d\n", MAXSIZE, sizelower, sizeupper);
    kT		=	kT_low;			//set temperature low, start quenching
    cycle	=	0;
    for (counter=1; counter<=MCS; counter++) {
      if (ran1(seed) < 0.08) {
        mcvol(); 
      }
      else {
        mcmove();
      }

      if (DYNVAR == "MAXSIZE") {		//update the histogram
        p[MAXSIZE]	++;
        if (MAXSIZE == sizeright && mod(cycle, 2)==0)		cycle ++;
        else if (MAXSIZE == sizeleft && mod(cycle, 2)==1)	cycle ++;
      }
      else if (DYNVAR == "Ql") {
        pQ[Qstatefinder(Ql)]	++;
        if (Qstatefinder(Ql) == Qrightid && mod(cycle, 2)==0)		cycle ++;
        else if (Qstatefinder(Ql) == Qleftid && mod(cycle, 2)==1)	cycle ++;
      }

      if (mod(counter, 200) == 0) {
	Adjust_Stepsize();      
      }
      if (mod(counter*5, MCS) == 0) {
        Print_Histogram();
      }
    }
    Update_Eta(DYNVAR, 2);
  }

  //  3.0 Production run
  cycle	=	0;

  for (counter=1; counter<=MCS*PRODRUN; counter++) {
    if (ran1(seed) < 0.08) {
      mcvol(); 
    }
    else {
      mcmove();
    }
      
    if (DYNVAR == "MAXSIZE") {			//update the histogram
      p[MAXSIZE]	++;
      if (MAXSIZE == sizeright && mod(cycle, 2)==0)		cycle ++;
      else if (MAXSIZE == sizeleft && mod(cycle, 2)==1)		cycle ++;
    }
    else if (DYNVAR == "Ql") {
      pQ[Qstatefinder(Ql)]	++;
      if (Qstatefinder(Ql) == Qrightid && mod(cycle, 2)==0)	cycle ++;
      else if (Qstatefinder(Ql) == Qleftid && mod(cycle, 2)==1)	cycle ++;
    }

    if (mod(counter, 200)==0) {
      Adjust_Stepsize();
    }
    if (mod(counter*1000, MCS*PRODRUN) == 0) {
      Update_Eta(DYNVAR, 3);			//calculate Gibbs free energy 
      Print_Histogram();
    }
    if (mod(counter*10, MCS*PRODRUN) == 0) {
      Visualize(0);
    }
  }
*/

  CloseFile();
  return;
}


/****************************************************************************************/
/*	PhaseNotDone(char *, int)							*/
/*											*/
/*	Check whether each phase is done						*/
/*	phase = 0 (melting), 1 (quenching trial run), 2 (quenching production run)	*/
/****************************************************************************************/
/*
int PhaseNotDone(char * var, int phase)
{
  if (phase == 0) {
    if (var == "MAXSIZE") {
      if (MAXSIZE <= sizeleft)		//melting done
	return 0;
      else
	return 1;
    }
    else if (var == "Ql") {
      if (Qstatefinder(Ql) <= Qleftid) {	//melting done
	return 0;
      }
      else {
	return 1;
      }
    }
  }
}
*/
                                                                    src/main.c                                                                                          0000600 0143352 0000144 00000013124 11274661755 012163  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    main.c
    author:     Peng Yi at MIT
    date:       October 23, 2006
    purpose:    Single processor NPT ensemble simulation program
    note:	Nov. 5, 2009, prepare to output more variables, and use another
			program to view a selected list, like g_energy in 
			Gromacs does.
*/

#define __MAIN_PROGRAM
#include "header.h"

void Print_hist(FILE *fPtr)
{
   long		i;
   static long	init=1;

   if (init) {						// Print out item name line
      fprintf(fPtr, "counter  ");
      for (i=0; i<NBOX; i++) {
         fprintf(fPtr, "vtot[%d]  ", i);		// 1. total potential energy
         fprintf(fPtr, "vbnd[%d]  ", i);		// 2. bond energy
         fprintf(fPtr, "vang[%d]  ", i);		// 3. angle energy
         fprintf(fPtr, "vtors[%d]  ", i);		// 4. torsion energy
         fprintf(fPtr, "vnbnd[%d]  ", i);		// 5. nonbonded energy
         fprintf(fPtr, "vol[%d]  ", i); 		// 6. volume
         fprintf(fPtr, "prestot[%d]  ", i);		// 7. pressure
         fprintf(fPtr, "x/z[%d]  ", i);
         fprintf(fPtr, "y/z[%d]  ", i);
         fprintf(fPtr, "drift2[%d]  ", i);
         fprintf(fPtr, "P2M[%d]  ", i);
         fprintf(fPtr, "P2[%d]  ", i);
         fprintf(fPtr, "P2z[%d]  ", i);
         fprintf(fPtr, "transfrac[%d]  ", i);
         fprintf(fPtr, "MAXSIZE[%d]  ", i);
         fprintf(fPtr, "Xtal[%d]  ", i);
         fprintf(fPtr, "realXtal[%d]  ", i);
         fprintf(fPtr, "2ndNmax[%d]  ", i);
         fprintf(fPtr, "Nnucl[%d] ", i);
         fprintf(fPtr, "Q6[%d]", i);
      }
      fprintf(fPtr, "\n");
      init	=	0;
   }

   fprintf(fPtr, "%-6d  ", counter);				// print out value
   for (i=0; i<NBOX; i++) {
      fprintf(fPtr, "%8.3f  ", v[i].tot);	
      fprintf(fPtr, "%8.3f  ", v[i].stretch);
      fprintf(fPtr, "%8.3f  ", v[i].bending);
      fprintf(fPtr, "%8.3f  ", v[i].torsion);
      fprintf(fPtr, "%8.3f  ", v[i].nonbonded);
      fprintf(fPtr, "%8.3f  ", BOX[i].vol);
      fprintf(fPtr, "%8.3f  ", BOX[i].pres);
      fprintf(fPtr, "%6.4f  ", BOX[i].lx/BOX[i].lz);
      fprintf(fPtr, "%6.4f  ", BOX[i].ly/BOX[i].lz);
      fprintf(fPtr, "%8.4f  ", drift2[i]);		// m.s.a. of c.o.m. displacement
      fprintf(fPtr, "%6.4f  ", P2M[i]);
      fprintf(fPtr, "%6.4f  ", P2[i]);
      fprintf(fPtr, "%6.4f  ", P2z[i]);
      fprintf(fPtr, "%6.4f  ", transfrac[i]);
      fprintf(fPtr, "%4d  ", MAXSIZE[i]);
      fprintf(fPtr, "%4d  ", Xtal[i]);
      fprintf(fPtr, "%4d  ", realXtal[i]);
      fprintf(fPtr, "%4d  ", secondNmax[i]);
      fprintf(fPtr, "%4d  ", Nnucl[i]);
      fprintf(fPtr, "%6.4f", Q6[i]);
   }
   fprintf(fPtr, "\n");
   fflush(fPtr);
}


int main(int argc, char * argv[])
{
   long		i, ii, k;
   long		REALNCYCLE;

   long		jj;	//debug
   long		j, max;

   molstruct	*moli;

   tim=(int *)malloc(sizeof(int));     	//random number generator
   seed=(long *)malloc(sizeof(long));
   *tim=(int)time(NULL);
   *seed= -1*(*tim);           		//seed must start out as a negative long

   InitAll(argv);

#ifdef MPI				
   // to be added ...
#else
   //StartNewHistoryFile(file_hst, FALSE);	// open binary history file
#endif

   /* do the first sampling */

   for (i=0; i<NMOLS; i++)
      mol[i].origin	=	CenterofMass(mol+i);
         
   SampleP2All();				// calculate global P2 and local p2
   SampleSpherical();
   Find_Nuclei(dynvar);
   Dist_p2();					// because p2 is used in biasing
   Dist_Spherical();				// because Ptrans is used in biasing
   //S_PrintAll();
   //H_StoreCurrent(file_hst);			// write to binary history file
   Write_Conf(0);				// dump conf.

   Sample_All();				// initial sampling, energy, pressure, etc
   Print_hist(fhst);				// output variables to histogram file

   /* do the cycles */

   for (ii=0; ii<2; ii++) {
      if (ii==0) {
	 NCYCLE		=	Nequil;			// Equilibration run
         fprintf(foutput, "\nStart equilibration ......\n");
      }
      else if (ii==1) {
	 NCYCLE		=	Nprod;			// Production run
         fprintf(foutput, "\nStart production .......\n");
         InitSample();					// clear the distributions to zero
      }
      InitEnsembles();					// reset acceptance, etc

      for (counter=0; counter<NCYCLE; counter++) {	// do the MC cycles

         if (counter==10)				// in case initial configuration results big number
            CalcV(); 
/*
         if (!mod(counter+1, NCYCLE/8)) {
            type[2].EPSILON	=	type[0].EPSILON * (2-1/8-((double) (counter+1))/NCYCLE);
            type[3].EPSILON	=	type[1].EPSILON * (2-1/8-((double) (counter+1))/NCYCLE);
            fprintf(foutput, "\ttype[2].EPSILON = %f\ttype[3].EPSILON = %f\n", type[2].EPSILON,
		type[3].EPSILON);
            InitForcefield();
         }
*/
         Cycle();					// do MC moves
	 Sample_All();					// do sampling
         Print_hist(fhst);				// output variables

         //SampleN2NSQ();				// sample distributions
         if (fabs(kP2) > ZERO) {
            if (dynvar==3)
               Dist_p2();					// because p2 is used in biasing

            Dist_Spherical();				// because Ptrans is used in biasing
            //Find_Nuclei();				// must after p2 calculation
         }
         if (mod(counter+1, ITAPE)==0) {		// Output binary history file every ITAPE sequences
            //H_StoreCurrent(file_hst);
	    Write_Conf(counter+1);				// dump conf.
	 }
	 if (mod(counter+1, ICONF)==0 && ii==1)	{	// Output configuration file for restart
	    //Write_Conf(counter+1);
	 }
      }	// cycles done

      if (counter!=0) {					// energy check
         EnergyCheck();      
      }

      if (ii<=1) {
         S_PrintAll();					// output all distributions
      }

      if (ii==1) {
         Sample_Done();
         Write_Conf(-1);
         //Visualize(1);
      }
   }

   CloseFile();
   return;
}
                                                                                                                                                                                                                                                                                                                                                                                                                                            src/motion.c                                                                                        0000600 0143352 0000144 00000020474 11175611652 012541  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	motion.c
	author:		Peng Yi at MIT
	date:		Oct 26, 2007
    	purpose:    move molecule, sites, functions borrowed and modified from Pieter's
*/
#define __MOTION_MODULE
#include "motion.h"

// Calculate site based spherical coordinates
// Pieter's code didn't consider the effect of periodic boundary condition

sphere SiteSpherical(molstruct *molm, long i0)
{
   long                  i1, i2, i3;
   vector                n1, d1, d2, d3 = {1.0, 0.0, 0.0};
   sphere                s;
   matrix                A;
   //long		ib	=	molm->box;		// added by PY

  if ((i1 = molm->parent[i0])>=0)
  {
    d1                  = V_Subtr(molm->p+i0, molm->p+i1);
    //MapInBox2(&d1, PBC, BOX[ib].lbox);			// added by PY

    s.d                 = sqrt(V_Dot(&d1, &d1));
    if ((i2 = molm->parent[i1])>=0)
    {
      d2                = V_Subtr(molm->p+i1, molm->p+i2);
      //MapInBox2(&d2, PBC, BOX[ib].lbox);		// added by PY

      s.alpha           = acos((d1.x*d2.x+d1.y*d2.y+d1.z*d2.z)/
                                s.d/sqrt(d2.x*d2.x+d2.y*d2.y+d2.z*d2.z));
      n1                = V_Cross(&d2, &d1);
      n1                = V_Mult(1.0/sqrt(V_Dot(&n1, &n1)), &n1);
      if ((i3 = molm->parent[i2])>=0) {
        d3              = V_Subtr(molm->p+i2, molm->p+i3);
        //MapInBox2(&d3, PBC, BOX[ib].lbox);		// added by PY
      }
      A                 = M_Orientation(&d3, &d2);
      d1.y              = A.y.x*n1.x+A.y.y*n1.y+A.y.z*n1.z;
      d1.z              = A.z.x*n1.x+A.z.y*n1.y+A.z.z*n1.z;
      s.beta            = atan2(-d1.y, d1.z);
      return s;
    }
    s.alpha             = acos(d1.x/s.d);
    s.beta              = atan2(d1.z, d1.y);
    return s;
  }
  s.d                   = sqrt(V_Dot(molm->p+i0, molm->p+i0));
  s.alpha               = acos(molm->p[i0].x/s.d);
  s.beta                = atan2(molm->p[i0].z, molm->p[i0].y);
  return s;
}


void MolSpherical(molstruct *molm)
{
  long                  i;

  for (i=1; i<molm->nsites; ++i)		// why not start from i=0?
    molm->s[i]          = SiteSpherical(molm, i);
}


void AllSpherical()
{
  molstruct             *moli;

  for (moli=mol; moli<mol+NMOLS; ++moli)
    MolSpherical(moli);
}


// Calculate site based cartesian coordinates given a spherical coordinate
// Assumes presence of cartesian coordinates of the parent sites

vector SiteCartesian(molstruct *molm, long i0, sphere s)
{
  long                  i1, i2, i3;
  vector                p, d1, d2 = {1.0, 0.0, 0.0};
  matrix                A, B;
  //long		ib=molm->box;			// added by PY

  if ((i1 = molm->parent[i0])<0) return molm->p[i0];
  A                     = M_Rotation(s.alpha, s.beta);
  if ((i2 = molm->parent[i1])>=0)
  {
    d1                  = V_Subtr(molm->p+i1, molm->p+i2);
    //MapInBox2(&d1, PBC, BOX[ib].lbox);			// added by PY
    if ((i3 = molm->parent[i2])>=0) {
      d2                = V_Subtr(molm->p+i2, molm->p+i3);
      //MapInBox2(&d2, PBC, BOX[ib].lbox);		// added by PY
    }
    B                   = M_Orientation(&d2, &d1);
    A                   = M_Dot(&B, &A);
  }
  p                     = molm->p[i1];
  p.x                   += A.x.x * s.d;
  p.y                   += A.x.y * s.d;
  p.z                   += A.x.z * s.d;

  //MapInBox2(&p, PBC, BOX[ib].lbox);			// added by PY
  return p;
}


void MolCartesian(molstruct *molm)
{
  long                  i;

  for (i=0; i<molm->nsites; ++i)
    molm->p[i]          = SiteCartesian(molm, i, molm->s[i]);
}


void AllCartesian()
{
  molstruct             *moli;

  for (moli=mol; moli<mol+NMOLS; ++moli)
    MolCartesian(moli);
}


/* Molecule manipulator */

void SiteCopy(molstruct *molm, long m, molstruct *moln, long n, long p)
{						// copy moln.n to molm.m
   molm->s[m]		=	moln->s[n];
   molm->p[m]		=	moln->p[n];
   molm->flags[m]	=	moln->flags[n];
   molm->type[m]	=	moln->type[n];
   molm->parent[m]	=	moln->parent[n]+p;
#ifdef CELL_LIST
   molm->cell[m]	=	moln->cell[n];
   molm->cellsite[m]	=	moln->cellsite[n];
   // cell list update through CL_Relink() !!
#endif
}


void SiteSwap(molstruct *molm, long m, molstruct *moln, long n)
{
   static sphere	s;
   static vector	p;
   long			i;
#ifdef CELL_LIST
   cellstruct		*cell;
#endif
  
   s			=	molm->s[m];
   molm->s[m]		=	moln->s[n];
   moln->s[n]		=	s;
   p			=	molm->p[m];
   molm->p[m]		=	moln->p[n];
   moln->p[n]		=	p;
   i			=	molm->flags[m];
   molm->flags[m]	=	moln->flags[n];
   moln->flags[n]	=	i;
   i			=	molm->type[m];
   molm->type[m]	=	moln->type[n];
   moln->type[n]	=	i;
#ifdef CELL_LIST
   cell			=	molm->cell[m];
   molm->cell[m]	=	moln->cell[n];
   moln->cell[n]	=	cell;

   i			=	molm->cellsite[m];
   molm->cellsite[m]	=	moln->cellsite[n];
   moln->cellsite[n]	=	i;
   // cell list update through CL_Relink() !!
#endif
}


void MolReverse(molstruct *mol1)
{
   static long		i;
   static molstruct	mol2;

   mol2.flip		=	1 - mol1->flip;		// track the flip
   mol2.nsites		=	mol1->nsites;
   mol2.box		=	mol1->box;

   for (i=0; i<mol1->nsites; i++) {
      SiteCopy(&mol2, i, mol1, mol1->nsites-1-i, 0);	// s, p, flags, parent, cell, type
      mol2.parent[i]	=	i-1;			// update parent site
   }

   // set correct bond angles and lengths
/*
   if (mol2.nsites > 6) {			// long chain
      for (i=0; i<3; i++)
         mol2.s[i]	=	SiteSpherical(&mol2, i);
      for (i=3; i<mol2.nsites-3; i++) {
         mol2.s[i].d		=	mol1->s[mol2.nsites-i].d;
         mol2.s[i].alpha	=	mol1->s[mol2.nsites+1-i].alpha;
         mol2.s[i].beta		=	mol1->s[mol2.nsites+2-i].beta;
      }
      for (i=mol2.nsites-3; i<mol2.nsites; i++)
         mol2.s[i]	=	SiteSpherical(&mol2, i);
   }
   else {					// short chain
      MolSpherical(&mol2);
   }
*/
   *mol1	=	mol2;
}


void MolFlip(molstruct *moli)
{
   MolReverse(moli);
#ifdef CELL_LIST
   CL_Relink(moli);
#endif
}

// Update cell list through CL_Relink after calling MolAdd
// Functions mainly used by ensembles.c
// CL_Relink is called after mol reassignment to update the cell tables

molstruct MolAdd(molstruct *mol1, molstruct *mol2)
{
   long			i;
   static molstruct	molm;

   molm		=	*mol1;
   for (i=0; i<mol2->nsites; i++) {
      SiteCopy(&molm, i+molm.nsites, mol2, i, molm.nsites);
   }
   molm.nsites	+=	mol2->nsites;
//   molm.fix	=	mol1->fix | mol2->fix;
//   Pieter didn't change type here, so need to do it somewhere else
   return	molm;
}


void ChangeAxis(long system, vector scale)
{
   long		j;
   molstruct	*moli;
   vector	pcm, dp, d;

   BOX[system].lx	*=	scale.x;
   BOX[system].ly	*=	scale.y;
   BOX[system].lz	*=	scale.z;
   d.x			=	scale.x-1.0;
   d.y			=	scale.y-1.0;
   d.z			=	scale.z-1.0;

   BOX[system].lbox	=	MIN(MIN(BOX[system].lx, BOX[system].ly), BOX[system].lz);
   BOX[system].vol	*=	scale.x * scale.y * scale.z;
/*
   if (SCALECUTOFF) {
      BOX[system].rc	*=	scale;	
      BOX[system].rv	*=	scale;
      BOX[system].rb	*=	scale;
   }
*/
   for (moli=mol; moli<mol+NMOLS; moli++) {	// search thru mols in this processor
      if (moli->box == system) {
         pcm	=	CenterofMass(moli);	// center of mass of each chain
         moli->origin.x	*=	scale.x;
         moli->origin.y	*=	scale.y;
         moli->origin.z	*=	scale.z;
         dp.x	=	d.x * pcm.x;		// shift of com of each chain due to volume change
	 dp.y	=	d.y * pcm.y;
         dp.z	=	d.z * pcm.z;
         for (j=0; j<moli->nsites; j++) {	// bond length unchanged
	    moli->p[j].x	+=	dp.x;
	    moli->p[j].y	+=	dp.y;
	    moli->p[j].z	+=	dp.z;
            //MapInBox2(moli->p+j, PBC, BOX[system].lbox);
         }
      }
   }
#ifdef CELL_LIST
   CL_Destroy();
   CL_Build();
#endif
}



void ChangeVolume(long system, double scale)	// box dimension, molecule coord. cell list
{
   long		j;
   molstruct	*moli;
   vector	pcm, dp;

   BOX[system].lbox	*=	scale;
   BOX[system].lx	*=	scale;
   BOX[system].ly	*=	scale;
   BOX[system].lz	*=	scale;
   BOX[system].vol	*=	scale * scale * scale;

   if (V_SCALECUTOFF) {
      BOX[system].rc	*=	scale;
      BOX[system].rv	*=	scale;
      BOX[system].rb	*=	scale;
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {	// search thru mols in this processor
      if (moli->box == system) {
         pcm	=	CenterofMass(moli);
         dp	=	V_Mult(scale-1.0, &pcm);	// CoM of chain shift
         for (j=0; j<moli->nsites; j++) {	// bond length unchanged
	    moli->p[j].x	+=	dp.x;
	    moli->p[j].y	+=	dp.y;
	    moli->p[j].z	+=	dp.z;
            //MapInBox2(moli->p+j, PBC, BOX[system].lbox);
         }
      }
   }
#ifdef CELL_LIST
   CL_Destroy();
   CL_Build();
#endif
}


                                                                                                                                                                                                    src/mymath.c                                                                                        0000600 0143352 0000144 00000000322 10714730275 012522  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    mymath.c
    author:     Peng Yi at MIT
    date:       April 10, 2007
    purpose:    my math functions which are heavily used in the program
*/

#define __MYMATH_MODULE
#include "mymath.h"
                                                                                                                                                                                                                                                                                                              src/period.c                                                                                        0000600 0143352 0000144 00000037576 11247056641 012532  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
	program:	period.c
	author:		Peng Yi at MIT
	date:		Aug 29, 2009
	purpose:	create periodic image 
	note:		
*/

#define __MAIN_PROGRAM
#include "header.h"

#define CAP	5		// largest nuclei considered still the melt
#define p2threshold	0.4

#include "correlation.h"

long		timestep;
double		rshift2;		// shift of the biggest nucleus
long		nsitesp2;	// # of sites with p2 greater than a threshold
long		nmolsp2;
long		nmaxp2[10];		// nmax using p2 nucleus definition

char *Element(long t, char *s)
{
   double	mass = type[t].M;

   s[0]		=	0;
/*
   if (14==mass || 15==mass)	strcpy(s, "C");
   else if (1.01==mass)		strcpy(s, "H");
   else if (28.086==mass)	strcpy(s, "Si");
   else if (26.982==mass)	strcpy(s, "Al");
   else if (16==mass)		strcpy(s, "O");
*/
   strcpy(s, "C");
   return	s;
}


int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   long		i, j, k, system, n;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo;			// box upper boundary
   long		nsites, record;
   long		nx, ny, nz, id, siteid, molid, type;
   molstruct	*moli;

   char		filein[255], s[80], ff[80], par[80];
   FILE		*fin, *fhst, *fout, *fconf, *fpdb, *fdat;
   char		a[255], b[255], c[255], d[255], dummy[255];
   long		LENGTH, confflag=1, carflag=1, pdbflag=1, drawmol;
   char		atomname;
   static long	init=1;

   vector		con;	// center of nucleus
   static vector	rO;	// original position of nucleus

   vector		chainface;
   double		orient;
   long			orientdist[180];	// chain orientation distribution

   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];

   vector	com[MAXNMOLS], temp;		// average center of mass
   long		ncom[MAXNMOLS];

   for (i=0; i<180; i++)
	orientdist[i]	=	0;

   if (argc<6) {
      printf("period (c) 2009 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tperiod x y z n dumpfile\n\n");
      printf("Notes:\n");
      printf("\t* x, y, z: the direction to create image\n");
      printf("\t* n: the number of images to create\n");
      exit(1);
   }

   shiftx	=	atof(argv[1]);
   shifty	=	atof(argv[2]);
   shiftz	=	atof(argv[3]);
   n		=	atol(argv[4]);

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-noconf")) 
         confflag	=	0;
      else if (samestr(par, "-nocar"))
         carflag	=	0;
      else if (samestr(par, "-nopdb"))
         pdbflag	=	0;
   }
   strcpy(filein, argv[argc-1]);

   if ( (fin=fopen(filein, "r"))==NULL )
      Exit("lammps2hst", "main", "open input file failed.");
   if ((fhst=fopen("lammps.hst", "w"))==NULL )
      Exit("lammps2hst", "main", "open hst file failed.");
   if (confflag && (fconf=fopen("lammps.conf", "w"))==NULL )
      Exit("lammps2hst", "main", "open conf file failed.");   
   if (carflag && (fout=fopen("lammps.car", "w"))==NULL )
      Exit("lammps2hst", "main", "open car file failed.");
   if (pdbflag && ((fpdb=fopen("lammps.pdb", "w"))==NULL || (fdat=fopen("vmd.dat","w"))==NULL))
      Exit("lammps2hst", "main", "open pdb file failed.");

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   InitUnits();				// initialize units for calculation
   InitForcefield();			// initialize Lennard-Jones potential mixing rule

   system	=	0;		// for now 2/14/08, only one system

   InitSample();			// initialize sampling

   while (!feof(fin)) {

      /* Read in one configuration from dump file */

      if (!fgets(a, sizeof(a), fin))		break;				// end of file
      fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);
      fgets(b, sizeof(b), fin);
      fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);
      fgets(c, sizeof(c), fin);
      fscanf(fin, "%lf%lf", &xlo, &xhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &ylo, &yhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &zlo, &zhi);	fgets(dummy, sizeof(dummy), fin);
      fgets(d, sizeof(d), fin);

      BOX[system].lx	=	xhi-xlo;
      BOX[system].ly	=	yhi-ylo;
      BOX[system].lz	=	zhi-zlo;

      LENGTH		=	NSITES/NMOLS;	// only for monodisperse system (4/26/2008)

      for (i=0; i<nsites; i++) {
         fscanf(fin, "%ld%ld %lf%lf%lf %lf%lf%lf %ld%ld%ld",
			&id, &type, &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
         fgets(dummy, sizeof(dummy), fin);

         id	--;				// because lammps starts from 1 not 0
         type	--;
         molid	=	(long) (id/LENGTH);
         siteid		=	id % LENGTH;
         mol[molid].box		=	system;		// for now, only one system
         mol[molid].nsites	=	LENGTH;
         mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx);
         mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly);
         mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);
         mol[molid].type[siteid]=	type;
      }

      // Start: Calculate the average position of center of mass
      for (moli=mol; moli<mol+NMOLS; moli++) {
      	  temp	=	CenterofMass(moli);
          if (fabs(temp.x) < 0.45 * BOX[moli->box].lx && 
		fabs(temp.y) < 0.45 * BOX[moli->box].ly &&
		fabs(temp.z) < 0.45 * BOX[moli->box].lz) {
	     com[moli-mol]	=	V_Add(&temp, com+(moli-mol));
             ncom[moli-mol]	++;
          }
      }
      // End: Calculate the average position of center of mass

      CoorSI2System();		// convert coordinates and box size from SI to system unit

      for (system=0; system<NSYSTEMS; system++) {
         NMols[system]	=	0;
         NSites[system]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;				// total # of mols in certain system
            NSites[system]	+=	moli->nsites;		// total # of sites in certain system
         }
      }

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
      } 

      for (moli=mol; moli<mol+NMOLS; moli++) { 		
         for (i=0; i<moli->nsites; i++)  {
            moli->flags[i]	=	1;		// activate all the sites on this processor
            moli->parent[i]	=	i-1;		// initialize parent site
         }
         moli->flip		=	0;		// flip to the original direction
         moli->origin	=	CenterofMass(moli);
      }

#ifdef CELL_LIST					// build cell list
      if (init) {
         CL_Init();		
         init	=	0;
      } 
      else 
         CL_Destroy();

      CL_Build();
#endif	/* CELL_LIST */

      // Perform Analysis
      CalcV();				// calc energy and virial
      if (V_VIRIAL)
         Sample_Pressure();

      SampleP2All();			// sample P2 and P2m and local p2
      Dist_p2();			// put local p2 into distribution
SampleM_Q();
      SampleSpherical();		// sample spherical coordinates
      Dist_Spherical();			// put spherical coord. into distribution
      Find_Nuclei_p2(1);
      nmaxp2[0]	=	nmax[0][0];
      nmaxp2[1]	=	nmax[0][1];
      nmaxp2[2]	=	nmax[0][2];
sizeofnucl	=	sizeofnuclp2;
sizedist	=	sizedistp2;
//      Find_Nuclei(dynvar);		// find nuclei
//      Calc_Qlm(6);			// calc. Qlm for LJ system, require CELL_LIST

/*
	 for (moli=mol; moli<mol+NMOLS; moli++)		// calculate the drift of biggest nucleus
            if (sizeofnucl[moli->nuclid[0]] == MAXSIZE[0]) {
               con	=	CenterofNucleus(moli->nuclid[0], moli);
	       break;
	    }
         if (timestep==0) {
	    rO	=	con;
	    rshift2	=	0;
         }
         else {
	    con	=	V_Subtr(&con, &rO);
	    con	=	MapInBox2(&con, PBC, system);
	    rshift2	=	V_Dot(&con, &con);
         }
*/
      if (timestep==0) {
         rO		=	CoM_MaxNucleus(0);
         rshift2	=	0;
      }
      else {
         con	=	CoM_MaxNucleus(0);
         con	=	V_Subtr(&con, &rO);
         con	=	MapInBox2(&con, PBC, system);
         rshift2	=	V_Dot(&con, &con);
      }

      // chain orientation distribution, test the rotator phase

      for (moli=mol; moli<mol+NMOLS; moli++) {
         V_Null(&chainface);

         for (i=0; i<moli->nsites; i++) {
            if (mod(i,2)) {
	       chainface	=	V_Add(&chainface, moli->p+i);
            }
	    else {
	       chainface	=	V_Subtr(&chainface, moli->p+i);
	    }
	 }
	 orient	=	atan2(chainface.y, chainface.x)	+ M_PI;
	 orientdist[(int)(orient/M_PI*90)]	++;
//	 printf("%f %f %f %d\n", chainface.x, chainface.y, orient, orientdist[(int)(orient/M_PI*90)]);
      }

      // calculate the # of sites with local p2 greater than a threshold

      nsitesp2	=	0;
      nmolsp2	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
//       nsitesp2	=	0;
         for (i=0; i<moli->nsites; i++) {
            if (moli->p2[i] > 0.4)
               nsitesp2	++;
         }
//            if (nsitesp2>1 && nsitesp2<6)
//	       nmolsp2	++;
      }

      Print_hst(fhst);		// print out hisgram
/*
Find_Nuclei_p2(-1);
printf("Xtal = %d\n Nmax = %d\n Nnucl = %d\n", Xtal[0], nmax[0][0], Nnucl[0]);
fprintf(fhst, " %4d\n", nmax[0][0]);
*/
      correlation();			// calculate correlation, in system units
      if (!(timestep%IRADIAL)) {
         radial("sample");		// sample radial distribution function
      }

      /* Group the segments belong to the biggest nucleus */
      n	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if (system == moli->box) {
            for (i=0; i<moli->nsites; i++) {
               if (sizeofnucl[moli->nuclid[i]] == nmax[system][0]) {
                  nucleus[n].moli	=	moli;
	          nucleus[n].site	=	i;
                  n	++;
	       }
            }
         }
      }

//      shiftbox(system, nucleus, n);

      /* Output */

      CoorSystem2SI();		// convert coordinates back to SI units

      // OUTPUT .car file for VMD visualization

      if (carflag) {
         fprintf(fout, "!BIOSYM archive 3\n");
         fprintf(fout, "PBC=ON\n");
         fprintf(fout, "!TIMESTEP %d\n", timestep);
         fprintf(fout, "!DATE %s", asctime(localtime(&t)));
         fprintf(fout, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

	 n	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {

               MolInBox2(moli);
               for (i=0; i<moli->nsites; i++) {

                  if (moli->p2[i]>critp2)
                     sprintf(s, "M%d", n++);
                  else
	             sprintf(s, "C%d", n++);

/*
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])		// note nuclid index starts from 1
                     sprintf(s, "M%d", n++);
                  else if (sizeofnucl[moli->nuclid[i]] -MAXSIZE[system] > -3 && MAXSIZE[system]>=10)
                     sprintf(s, "C%d", n++);
	          else if (moli->nuclid[i] >= 1)
                     sprintf(s, "O%d", n++);
                  else
                     sprintf(s, "H%d", n++);
*/
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

      // OUTPUT .pdb file for further analysis
      if (pdbflag) {
         fprintf(fpdb, "HEADER: file created from %s on %s", argv[1], asctime(localtime(&t)));
         fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
        	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

         n		=	0;
#define SIZECAP	5
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {
               MolInBox2(moli);

               drawmol	=	0;
               for (i=0; i<moli->nsites; i++) {
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {
	             drawmol	=	1;		// participate the biggest nucleus
		     break;
                  }
	       }

               for (i=0; i<moli->nsites; i++) {

                  if (drawmol) {
                     if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {	// nuclid index starts from 1
                        atomname	=	'N';
		        fprintf(fdat, " 10");
                     }
                     else {
                        atomname	=	'O';
                        fprintf(fdat, " 0");
                     }
                  }
                  else {
		     atomname	=	'C';
		     fprintf(fdat," -1");
	          }

                  n	++;

	          fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
                  fprintf(fpdb, "%5d ", n);	// atom number
                  fprintf(fpdb, " %c  ", atomname);	// atom name
                  fprintf(fpdb, " ");		// alternate location indiator
  	          fprintf(fpdb, "   ");		// residue name
	          fprintf(fpdb, " ");		// column 21
                  fprintf(fpdb, " ");		// chain identifier, column 22
	          fprintf(fpdb, "    ");		// residue sequence number, 23-26
	          fprintf(fpdb, " ");		// code for insertion of residues, 27
                  fprintf(fpdb, "   ");		// column 28-30
                  fprintf(fpdb, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z);
                  fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                  fprintf(fpdb, "%5.5s", "");
                  fprintf(fpdb, "\n"); 
               } 
            }   
         }
         fprintf(fpdb, "END\n");
      }	//pdbflag

      // OUTPUT configuration file for further analysis

      if (confflag) { 
         fprintf(fconf, "!TIMESTEP %d\n", timestep);
         fprintf(fconf, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
         for (i=0; i<NSYSTEMS; i++)
            fprintf(fconf, "%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

         for (moli=mol; moli<mol+NMOLS; moli++) {
            fprintf(fconf, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
            //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
            //MolInBox(moli);
            for (i=0; i<moli->nsites; i++) 
               fprintf(fconf, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
         }
      }
   } 

   // output average center of mass of chains
   if (pdbflag) {
      fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
      fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
 	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

      n		=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         i	=	moli-mol;
         if (ncom[i]>0) {
            n	++;
            atomname	=	'N';
            fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
            fprintf(fpdb, "%5d ", n);		// atom number
            fprintf(fpdb, " %c  ", atomname);	// atom name
            fprintf(fpdb, " ");		// alternate location indiator
            fprintf(fpdb, "   ");	// residue name
            fprintf(fpdb, " ");		// column 21
            fprintf(fpdb, " ");		// chain identifier, column 22
            fprintf(fpdb, "    ");	// residue sequence number, 23-26
            fprintf(fpdb, " ");		// code for insertion of residues, 27
            fprintf(fpdb, "   ");	// column 28-30
            fprintf(fpdb, "%8.3f%8.3f%8.3f", com[i].x/ncom[i], com[i].y/ncom[i], com[i].z/ncom[i]);
            fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
            fprintf(fpdb, "\n"); 
         }
      }
      fprintf(fpdb, "END\n");
   }

   S_PrintAll();		// print out distributions
   corr_norm();		// normalize correlation function
   corr_print();		// print out correlation function
   radial("print");		// print out radial distribution function

   fclose(fin);		// close files
   fclose(fhst);
   if (carflag)		fclose(fout);
   if (confflag)    	fclose(fconf);
   if (pdbflag)	     {	fclose(fpdb); fclose(fdat);}

   for (i=0; i<180; i++) {
      printf("%d %d\n", i, orientdist[i]);
   }
   return	0;
}


                                                                                                                                  src/position.c                                                                                      0000600 0143352 0000144 00000176504 11624300133 013072  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
  program:      position.c
  author:       Peng Yi at MIT
  date:         October 22, 2006
  purpose:      Position calculation, Xtal nuclei identification
*/

#define __POSITION_MODULE
#include "position.h"

double AdjustAngle(double x)
{
  while (x>=M_PI) x -= 2.0*M_PI;
  while (x<-M_PI) x += 2.0*M_PI;
  return x;
}

matrix  abc;

/*
void GetBoxVectors(vector *a, vector *b, vector *c)
{
  *a                    = abc.x;
  *b                    = abc.y;
  *c                    = abc.z;
}
*/

void MapInBox(vector *p)
{
   double	L = LBOX;

   (p->x)	/=	L;
   (p->y)	/=	L;
   (p->z)	/=	L;

   if (PBC==1) {		//cubic periodic boundary condition

      if ((p->x<-0.5) || (p->x>=0.5))
         p->x	-=	floor(p->x + 0.5);
      if ((p->y<-0.5) || (p->y>=0.5))
         p->y	-=	floor(p->y + 0.5);
      if ((p->z<-0.5) || (p->z>=0.5))
         p->z	-=	floor(p->z + 0.5);

   }
   if (PBC==2) {		// truncated octahedron periodic boundary conditioni, ref: F.23 of CCP5 library
				// LBOX is the length of the cubic which truncates the octahedrona
				// haven't done the very far map back as PBC==1 case (10/31/07)
      if ((p->x) >= 0.5) 	(p->x) -= 1.0;
      else if ((p->x) < -0.5)	(p->x) += 1.0;
      if ((p->y) >= 0.5)   	(p->y) -= 1.0;
      else if ((p->y) < -0.5)	(p->y) += 1.0;
      if ((p->z) >= 0.5)   	(p->z) -= 1.0;
      else if ((p->z) < -0.5)	(p->z) += 1.0;

      if ( (fabs(p->x)+fabs(p->y)+fabs(p->z)) > 0.75) {
	 if ((p->x) >= 0) 	(p->x) -= 0.5;
	 else			(p->x) += 0.5;         	 	
	 if ((p->y) >= 0) 	(p->y) -= 0.5;
	 else			(p->y) += 0.5;         	 	
	 if ((p->z) >= 0) 	(p->z) -= 0.5;
	 else			(p->z) += 0.5;         	 	
      }
      else if ( (p->z) >= 0 && (fabs(p->x)+fabs(p->y)+fabs(p->z)) == 0.75) {
	 (p->x)	+=	((p->x) >=0 ? -0.5 : 0.5); 
	 (p->y)	+=	((p->y) >=0 ? -0.5 : 0.5);
	 (p->z)	-=	0.5; 
      }
   }   
   (p->x)	*=	L;
   (p->y)	*=	L;
   (p->z)	*=	L;
}


vector MapInBox2(vector *p, long PBC, long system)	// If particle escapes the box, map it back, and pick up 
						// the image due to pbc so that the distance is the minimum
{
   static vector	l;
   double		zlo, zhi, ylo, yhi, xlo, xhi, yoff, xoff;
   vector		A, B, C;

   if (PBC==1) {		//cubic periodic boundary condition
      l.x		=	p->x / BOX[system].lx;	// here do not change *p, unlike MapInBox()
      l.y		=	p->y / BOX[system].ly;
      l.z		=	p->z / BOX[system].lz;

      if ((l.x<-0.5) || (l.x>=0.5))         l.x	-=	floor(l.x + 0.5);
      if ((l.y<-0.5) || (l.y>=0.5))         l.y	-=	floor(l.y + 0.5);
      if ((l.z<-0.5) || (l.z>=0.5))         l.z	-=	floor(l.z + 0.5);

      l.x		*=	BOX[system].lx;
      l.y		*=	BOX[system].ly;
      l.z		*=	BOX[system].lz;
   }
   else if (PBC==3) {		// triclinic periodic box
      A.x = BOX[system].lx;	A.y = 0;		A.z = 0;
      B.x = BOX[system].xy;     B.y = BOX[system].ly;	B.z = 0;
      C.x = BOX[system].xz;	C.y = BOX[system].yz;	C.z = BOX[system].lz;

      xhi	=	0.5 * BOX[system].lx;      xlo	=	-xhi;
      yhi	=	0.5 * BOX[system].ly;      ylo	=	-yhi;
      zhi	=	0.5 * BOX[system].lz;      zlo	=	-zhi;

      l		=	*p;

      // bring z coordinate to be b/w zlo and zhi
      while (l.z < zlo) { 
         l.z += C.z;	l.y += C.y;	l.x += C.x;
      }
      while (l.z >= zhi) {
         l.z -= C.z;	l.y -= C.y;	l.x -= C.x;
      }
      yoff	=	BOX[system].yz * (l.z-zlo) / (zhi-zlo);
      xoff	=	BOX[system].xz * (l.z-zlo) / (zhi-zlo);

      // bring y coordinate to be b/w ylo+yoff and yhi+yoff
      while (l.y < ylo+yoff) {
         l.y += B.y;	l.x += B.x;
      }
      while (l.y >= yhi+yoff) {
         l.y -= B.y;	l.x -= B.x;
      }
      xoff	+=	BOX[system].xy * (l.y - (ylo+yoff))/(yhi-ylo);

      // bring x coordinate to the center box
      while (l.x < xlo+xoff) {
         l.x += A.x;
      }
      while (l.y >= xhi+xoff) {
         l.x -= A.x;
      }
   }
   /*
   if (PBC==2) {		// truncated octahedron periodic boundary conditioni, ref: F.23 of CCP5 library
				// LBOX is the length of the cubic which truncates the octahedron
				// haven't done the very far away map back like PBC==1 case (10/31/07)
      if ((p->x) >= 0.5) 	(p->x) -= 1.0;
      else if ((p->x) < -0.5)	(p->x) += 1.0;
      if ((p->y) >= 0.5)   	(p->y) -= 1.0;
      else if ((p->y) < -0.5)	(p->y) += 1.0;
      if ((p->z) >= 0.5)   	(p->z) -= 1.0;
      else if ((p->z) < -0.5)	(p->z) += 1.0;

      if ( (fabs(p->x)+fabs(p->y)+fabs(p->z)) > 0.75) {
	 if ((p->x) >= 0) 	(p->x) -= 0.5;
	 else			(p->x) += 0.5;         	 	
	 if ((p->y) >= 0) 	(p->y) -= 0.5;
	 else			(p->y) += 0.5;         	 	
	 if ((p->z) >= 0) 	(p->z) -= 0.5;
	 else			(p->z) += 0.5;         	 	
      }
      else if ( (p->z) >= 0 && (fabs(p->x)+fabs(p->y)+fabs(p->z)) == 0.75) {
	 (p->x)	+=	((p->x) >=0 ? -0.5 : 0.5); 
	 (p->y)	+=	((p->y) >=0 ? -0.5 : 0.5);
	 (p->z)	-=	0.5; 
      }
   } 
   */  
   return l;
}


/***************************************************************************************/

double DistSQ(vector p, vector q, long system)	//distance square between vector p and q
{
   long		i, j, k;
   double	r2img[3][3][3], r2min;
   double	xij, yij, zij;
   double	r2, r2new;
   vector	dp, pimg, l, A, B, C, r;
/*
   dp.x	=	p.x - q.x;
   dp.y	=	p.y - q.y;
   dp.z	=	p.z - q.z;
   dp	=	MapInBox2(&dp, PBC, system);
   return	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
*/
   l.x	=	p.x - q.x;
   l.y	=	p.y - q.y;
   l.z	=	p.z - q.z;

   if (PBC==1) {		//cubic periodic boundary condition
      l.x	/=	BOX[system].lx;	
      l.y	/=	BOX[system].ly;
      l.z	/=	BOX[system].lz;

      if ((l.x<-0.5) || (l.x>=0.5))
         l.x	-=	floor(l.x + 0.5);
      if ((l.y<-0.5) || (l.y>=0.5))
         l.y	-=	floor(l.y + 0.5);
      if ((l.z<-0.5) || (l.z>=0.5))
         l.z	-=	floor(l.z + 0.5);

      l.x	*=	BOX[system].lx;
      l.y	*=	BOX[system].ly;
      l.z	*=	BOX[system].lz;
      r2	=	l.x * l.x + l.y * l.y + l.z * l.z;
   }
   else if (PBC==3) {				// added 2011.6.15
      A.x = BOX[system].lx;	A.y = 0;		A.z = 0;
      B.x = BOX[system].xy;     B.y = BOX[system].ly;	B.z = 0;
      C.x = BOX[system].xz;	C.y = BOX[system].yz;	C.z = BOX[system].lz;

      r2	=	l.x * l.x + l.y * l.y + l.z * l.z;
    
      for (i=-1; i<=1; i++) {
         for (j=-1; j<=1; j++) {
            for (k=-1; k<=1; k++) {
	        r.x	=	l.x + i * A.x + j * B.x + k * C.x;
	        r.y	=	l.y + i * A.y + j * B.y + k * C.y;
	        r.z	=	l.z + i * A.z + j * B.z + k * C.z;
                r2new	=	r.x * r.x + r.y * r.y + r.z * r.z;
                if (r2new < r2)		
		   r2	=	r2new;
      }  }  }
   }

   return	r2;

/*
   if (PBC==3) {			// triclinic box
     for (i=-1; i<=1; i++) {
       for (j=-1; j<=1; j++) {
         for (k=-1; k<=1; k++) {
           pimg.x		=	p.x + i * a.x + j * b.x + k * c.x;
           pimg.y		=	p.y + i * a.y + j * b.y + k * c.y;
           pimg.z		=	p.z + i * a.z + j * b.z + k * c.z;
           r2[i+1][j+1][k+1]	=	(pimg.x - q.x)*(pimg.x - q.x) + 
					(pimg.y - q.y)*(pimg.y - q.y) +
					(pimg.z - q.z)*(pimg.z - q.z);
         }
       }
     }
   }
*/
/*
   if (PBC==1) {
     
      xij	=	MIN(fabs(p.x-q.x), L-fabs(p.x-q.x));
      yij	=	MIN(fabs(p.y-q.y), L-fabs(p.y-q.y));
      zij	=	MIN(fabs(p.z-q.z), L-fabs(p.z-q.z));
      r2	=	xij * xij + yij * yij + zij * zij;
   }
   if (PBC==2) {			//truncated octahedron box PBC
      xij	=	MIN(fabs(p.x-q.x), L-fabs(p.x-q.x));		//min distance among the cubic images
      yij	=	MIN(fabs(p.y-q.y), L-fabs(p.y-q.y));
      zij	=	MIN(fabs(p.z-q.z), L-fabs(p.z-q.z));
      r2	=	xij * xij + yij * yij + zij * zij;

      pimg.x	=	p.x + ((q.x-p.x) > 0 ? 0.5 : -0.5)*L;	//look around the octahedron images
      pimg.y	=	p.y + ((q.y-p.y) > 0 ? 0.5 : -0.5)*L;
      pimg.z	=	p.z + ((q.z-p.z) > 0 ? 0.5 : -0.5)*L; 
      r2new	=	(pimg.x-q.x) * (pimg.x-q.x) + (pimg.y-q.y) * (pimg.y-q.y) + (pimg.z-q.z)*(pimg.z-q.z); 

      if (r2new < r2)
	 r2	=	r2new;     
   }

   return	r2;
*/
}

/***************************************************************************************/

void InitLattice(long Nmols, long Nsites, double L, long PBC)	// initialize chain molecules on an fcc lattice
{								// in a cubic box with PBC
   long		i, j, k, nmols;
   long		NC;
   double	bondlength, cell, cell2;
   double	x0, y0, z0;

   if (PBC!=1)
      Exit("position", "InitLattice", "PBC!=1.");

   bondlength	=	1.0;
   cell		=	bondlength * sqrt(2.0);
   cell2	=	0.5 * cell;

   nmols	=	0;
   
   // chains start on x-y plane and grow in the z direction
   // chains are on x-z plane

   for (j=0; j<Nmols; j++) {				// place chains along y-axis to the box wall
      if ( (y0 = j * cell) < L-cell2 ) {		// y0 increases by cell

         for (i=0; i<Nmols; i++) {			// place chains along x-axis to the box wall
            if ( (x0 = i * cell2) < L-cell ) {		// x0 increses by cell2

               if (nmols < Nmols) {			// if still more chains to place

                  for (k=0; k< Nsites; k++) {		// grow chains along z-axis

                     z0		=	k * cell2;
                     mol[nmols].p[k].x	=	x0 + mod(k, 2) * cell2;
                     mol[nmols].p[k].y	=	y0 + mod(i, 2) * cell2;
                     mol[nmols].p[k].z	=	z0; 
                  }
		  nmols	++;
               }
            }
         }
      }
   }
   return;
}


vector ranor()						// Random vector on a unit sphere
{							// F&S algorithm 42
   double	rand1, rand2, ranh, ransq;
   vector	unit;

   ransq	=	2.0;

   do {
      rand1	=	1.0 - 2.0 * ran1(seed);		// -1 < ran1 < 1, cos(Phi)sin(theta)
      rand2	=	1.0 - 2.0 * ran1(seed);		// -1 < ran2 < 1, sin(Phi)sin(theta)
      ransq	=	rand1 * rand1 + rand2 * rand2;
   } while (ransq >= 1.0);

   ranh 	=	2.0 * sqrt(1.0-ransq);
   unit.x	=	rand1 * ranh;
   unit.y	=	rand2 * ranh;
   unit.z	=	(1.0 - 2 * ransq);
   
   return	unit;
}


void Amorph(long nmols, long nmolsites, double LX, double LY, double LZ, long PBC)	
					// generate random conf. of chain molecules
{
   molstruct	*moli;
   long		j, ibox;
   double	bondlength;
   vector	dp;

   bondlength	=	type[0].LSTRETCH;

   for (moli=mol; moli<mol+nmols; moli++) {

      if (1==PBC) {
         moli->p[0].x	=	LX * (ran1(seed)-0.5);		// place the first site
         moli->p[0].y	=	LY * (ran1(seed)-0.5);
         moli->p[0].z	=	LZ * (ran1(seed)-0.5);
      }

      for (j=1; j<nmolsites; j++) {

	 dp	=	ranor();
         //dp	=	tors_bonda(moli, j);
	 dp	=	V_Mult(bondlength, &dp);
	    
	 moli->p[j]	=	V_Add(moli->p+j-1, &dp);
	 //MapInBox2(moli->p+j, PBC, L);
      }
      moli->nsites	=	nmolsites;
   }
   return;
}

/********************************************************************************************/
#ifdef TEST
void sc_lattice(long N, double L, long PBC)		// total particle number N = NC * NC * NC
{							// box dimension L
   long		i, NC;
   double	cell;

   if (PBC==1)
      NC	=	(int) rint(pow(N, 1.0/3));

   if (PBC==1 && N!=NC*NC*NC)
      printf("Error, number of particle incomp. with sc lattice.\n");

   if (PBC==1) {
      cell	=	L/NC;

      for (i=0; i<N; i++) {
         part[i].p.x	=	(mod(i, NC) - NC/2) * cell;			//or (int)mod(i,M)/1
         part[i].p.y	=	((int)(mod(i, NC*NC)/NC) - NC/2) * cell;	//or (int)mod(i,M*M)/M
	 part[i].p.z	=	((int)(i/(NC*NC)) - NC/2) * cell;		//or (int)mod(i,M*M*M)/(M*M)
      }
   }
}


void bcc_lattice(long N, double L, long PBC)		// total particle number N = 2 * NC^3
{							// box dimension L
   long		i, j, NC;
   double	cell;

   if (PBC==1)
      NC	=	(int) rint(pow(N/2, 1.0/3));

   if (PBC==1 && N!=2*NC*NC*NC)
      printf("Error, number of particle incomp. with bcc lattice.\n");

   if (PBC==1) {
      cell	=	L/NC;
    
      for (i=0; i<N/2; i++) {
         part[i].p.x	=	(mod(i, NC) - NC/2) * cell;		//or (int)mod(i,M)/1
         part[i].p.y	=	((int)(mod(i, NC*NC)/NC) - NC/2) * cell;	//or (int)mod(i,M*M)/M
	 part[i].p.z	=	((int)(i/(NC*NC)) - NC/2) * cell;		//or (int)mod(i,M*M*M)/(M*M)

         part[i+N/2].p.x	=	part[i].p.x + 0.5 * cell;
         part[i+N/2].p.y	=	part[i].p.y + 0.5 * cell;
	 part[i+N/2].p.z	=	part[i].p.z + 0.5 * cell;
      }
   }
}


void fcc_lattice(long N, double L, long PBC)			//reference: F.23 of Allen and Tildesley
{
   long		i, j, k, m, ref;
   double	cell, cell2;
   long		NC;

   if (PBC==1)
      NC	=	(int) rint(pow(N/4, 1.0/3));		// rint: nearest integer value
   else if (PBC==2)
      NC	=	(int) rint(pow(N/16, 1.0/3));

   if ( (PBC==1 && N!=4*NC*NC*NC) || (PBC==2 && N!=16*NC*NC*NC))
      printf("Error, number of particles incomp. with fcc lattice.\n");

   if (PBC==1) {
      cell	=	(double)L/NC;			// the reason to use plus sign here is to adapt with our
      cell2	=	0.5*cell;			// p.b.c., that is if X<-LBOX/2, then X+=LBOX
							// if our pbc is if X<=-LBOX/2, then X+=LBOX, then minus sign
							// in another words, our box is [-L, +L)

      part[0].p.x	=	0;	part[0].p.y	=	0;	part[0].p.z	=	0;
      part[1].p.x	=	cell2;	part[1].p.y	=	cell2;	part[1].p.z	=	0;
      part[2].p.x	=	0;	part[2].p.y	=	cell2;	part[2].p.z	=	cell2;
      part[3].p.x	=	cell2;	part[3].p.y	=	0;	part[3].p.z	=	cell2;

      m	=	0;

      for (i=0; i<NC; i++) {
         for (j=0; j<NC; j++) {
	    for (k=0; k<NC; k++) {
	       for (ref=0; ref<4; ref++) {
	          part[ref+m].p.x	=	part[ref].p.x	+	cell*k;
		  part[ref+m].p.y	=	part[ref].p.y	+	cell*j;
		  part[ref+m].p.z	=	part[ref].p.z	+	cell*i;
	       }
	       m	+=	4;
	    }
	 }
      }
      for (i=0; i<N; i++) {
	  part[i].p.x	-=	(double)L/2;		//same reason to use += instead of -= as mentioned above
	  part[i].p.y	-=	(double)L/2;
	  part[i].p.z	-=	(double)L/2;
      }
   }
   if (PBC==2) {					//truncated octahedron periodic boundary condition
      printf("Set up fcc lattice in truncated octahedron box.\n");

      cell	=	L / (2*NC);
      cell2	=	0.5 * cell;

      part[0].p.x	=	0;		part[0].p.y	=	0;	part[0].p.z	=	0;
      part[1].p.x	=	cell2;		part[1].p.y	=	cell2;	part[1].p.z	=	0;
      part[2].p.x	=	0;		part[2].p.y	=	cell2;	part[2].p.z	=	cell2;
      part[3].p.x	=	cell2;		part[3].p.y	=	0;	part[3].p.z	=	cell2;

      m	=	0;

      for (i=0; i<NC; i++) {			// (2*NC, 2*NC, NC)
         for (j=0; j<2*NC; j++) {
	    for (k=0; k<2*NC; k++) {
               for (ref=0; ref<4; ref++) {
	          part[ref+m].p.x	=	part[ref].p.x	+	cell*k;
		  part[ref+m].p.y	=	part[ref].p.y	+	cell*j;
		  part[ref+m].p.z	=	part[ref].p.z	+	cell*i;
	       }
	       m	+=	4;
	    }
	 }
      }
     
      for (i=0; i<N; i++) {
	 part[i].p.x	-=	0.5 * L;	
	 part[i].p.y	-=	0.5 * L;
	 if ( (fabs(part[i].p.x) + fabs(part[i].p.y) + fabs(part[i].p.z)) > 0.75*L) {
	    part[i].p.x	+=	((part[i].p.x >= 0)	? -0.5 : 0.5) * L;
	    part[i].p.y	+=	((part[i].p.y >= 0)	? -0.5 : 0.5) * L;
	    part[i].p.z	+=	((part[i].p.z >= 0)	? -0.5 : 0.5) * L;
	 }
	 else if ( (fabs(part[i].p.x) + fabs(part[i].p.y) + fabs(part[i].p.z) == 0.75 * L) && part[i].p.z>=0) {
	    part[i].p.x	+=	((part[i].p.x >= 0)	? -0.5 : 0.5) * L;
	    part[i].p.y	+=	((part[i].p.y >= 0)	? -0.5 : 0.5) * L;
	    part[i].p.z	-=	0.5 * L;
	 }
      }
   }
}


void randomconf(long N, double L, long PBC)		//assign particle coordinates in random
{
   long		i;

   if (PBC==1) {
      for (i=0; i<N; i++) {
         part[i].p.x	=	(ran1(seed)-0.5) * L;
         part[i].p.y	=	(ran1(seed)-0.5) * L;
         part[i].p.z	=	(ran1(seed)-0.5) * L;
      }
   }
   if (PBC==2) {
      printf("Set up random initial configuration in truncated octahedron box.\n");
      for (i=0; i<N; i++) {
         do {
            part[i].p.x	=	(ran1(seed)-0.5) * L;
            part[i].p.y	=	(ran1(seed)-0.5) * L;
            part[i].p.z	=	(ran1(seed)-0.5) * L;
         } while (fabs(part[i].p.x) + fabs(part[i].p.y) + fabs(part[i].p.z) >= L*0.75);
      }
   }
}
#endif

/*********************************************************************************************/

#ifdef TEST
long getnuclsize(long i, long id)	//Get the size of crystal nucleus #id which contains particle #i.
{
   long	icell, jcell, jj, kk, k;
   long	nuclsize	=	1;	//this nucleus at least has one particle, namely, particle #i

   part[i].nuclid	=	id;	//index of crystal nuclei

#ifdef VERLET_LIST
   for (jj=0; jj<part[i].nverlet; jj++) {	//search its verlet neighbors

      k		=	part[i].vlist[jj];
      {{
#elif CELL_LIST
   icell	=	part[i].icell;
   for (jj=0; jj<Cell[icell].nneigh; jj++) {
      jcell	=	Cell[icell].neigh[jj];

      for (kk=0; kk<Cell[jcell].sites; kk++) {
         k	=	Cell[jcell].list[kk];

         if (k!=i) {
#else
   for (k=0; k<NPARTS; k++) {
      if (k!=i) {
      {
#endif

      if (part[k].nuclid == -1 && part[k].nconnect>=critconnect && DistSQ(part[i].p, part[k].p, part[k].box)<Rb2) {	//if neighbor is Xtal-like and not scanned yet
         nuclsize	+=	getnuclsize(k, id);
      }
   }}}
   return	nuclsize;
}

/********************************************************************************************/

void Find_Nuclei1()		// find and label all crystal nuclei, get size distribution, find max size
{				// the first Find_Nuclei function I wrote, using recursion, accurate
  long	i;
  long	id	=	1;	//the nuclei index starts with 1

  for (i=0; i<NPARTS; i++) {
    part[i].nuclid	=	-1;	//none of the particles has been scanned
  }
  for (i=0; i<NPARTS+1; i++) {		
    sizeofnucl[i]	=	0;	//initialize nuclei sizes
    sizedist[i]	=	0;		//initialize nuclei size distribution
  }
  
  for (i=0; i<NPARTS; i++) {
    if (part[i].nuclid	== -1 && part[i].nconnect >= critconnect) {	//crystal-like particle hasn't been scanned
      sizeofnucl[id]	=	getnuclsize(i, id);		//get the size of nucleus #id which contains particle #i
      id		++;
    }
  }

  MAXSIZE	=	0;
  Nnucl		=	0;		//number of Xtal nuclei
  Xtal		=	0;		//number of Xtal-like particles

  for (id=1; id<NPARTS+1; id++) {		//determine the max nucleus size
    if (sizeofnucl[id] != 0) {

      sizedist[sizeofnucl[id]]	++;
      Nnucl			++;
      Xtal	+=	sizeofnucl[id];

      if (MAXSIZE < sizeofnucl[id]){
	MAXSIZE	=	sizeofnucl[id];
      }
    }
  }
  /*
  for (i=0; i<NPARTS; i++) {
    fprintf(foutput,"%d\t%d\n", part[i].nconnect, part[i].nuclid);
  }
  fprintf(foutput,"\n");
   
  for (i=0; i<NPARTS+1; i++) {
    if (sizedist[i]!=0) {
	fprintf(foutput, "Find_Nuclei, Number of nuclei of size %d\t=\t%d\n", i, sizedist[i]);
    }
  }
  fprintf(foutput, "Maxsize=\t%d\n",MAXSIZE); 
  */
  return;
}
#endif

//////////////////////////////////////////////
/* Decide whether one site is in xtal phase */
//////////////////////////////////////////////
long crystal(molstruct *moli, long i)		// determine a xtal site
{
   if (samestr(moltype, "LJ"))			// Lennard Jones system
      return	moli->nconn[i]	> critconnect;
   else if (samestr(moltype, "monochain"))	// monodisperse chain molecules
      return	moli->p2[i] > critp2;
}
/////____________________________________////


//////////////////////////////
/* Find nuclei in LJ system */
//////////////////////////////
void Find_Nuclei_LJ()				// deal with multiple system
{
   long			i, j, k, jj, n, NIT, LIT, Lk, nuclid, system,
			nsites[MAXNSYSTEMS];	// # of sites in each system
   molstruct		*moli;
   static long		**L;
   static sitestruct	**site;
   vector		pj;
   double		r2;
#ifdef CELL_LIST
   cellstruct		*cellj, *cellk;
#endif
   static long		init=1;

   // allocate memory
   if (init) {
      L			=	(long **) calloc(NSYSTEMS, sizeof(long *));
      site		=	(sitestruct **) calloc(NSYSTEMS, sizeof(sitestruct *));
      if (L==NULL || site==NULL)
         Exit("position.c", "Find_Nuclei", "out of memory");

      for (i=0; i<NSYSTEMS; i++) {
	 L[i]		=	(long *) calloc(NSITES, sizeof(long));		// or NSites[i] to save memory, later
         site[i]	=	(sitestruct *) calloc(NSITES, sizeof(sitestruct));
         if (L[i]==NULL || site[i]==NULL)
            Exit("position.c", "Find_Nuclei", "out of memory");
      }

      sizedist		=	(long *) calloc(NSITES+1, sizeof(long));	// not for all systems, later
      sizeofnucl	=	(long *) calloc(NSITES+1, sizeof(long));
      if (sizedist==NULL || sizeofnucl==NULL)
         Exit("position.c", "Find_Nuclei", "out of memory");

      init	=	0;
   } 

   // identify crystal-like sites and register them
   for (system=0; system < NSYSTEMS; system++)
      nsites[system]	=	0;

   for (moli=mol; moli<mol+NMOLS; moli++) {
      system	=	moli->box;

      for (i=0; i<moli->nsites; i++) {
         n	=	nsites[system];
         site[system][n].mol 		=	moli;
     	 site[system][n].mol_site	=	i;
	 site[system][n].p		=	moli->p[i];
	 site[system][n].cell		=	moli->cell[i];
	 site[system][n].cellsite	=	moli->cellsite[i];
        
         if (crystal(moli, i))
            L[system][n]	=	n;		// label xtal sites
         else
	    L[system][n]	=	-1;		// melt sites

         nsites[system]	++;
         moli->nuclid[i]	=	-1;		// initialize nuclei ida
      }
   }

   // scan to determine clusters in all systems 
   for (system=0; system<NSYSTEMS; system++) {
      for (i=0; i<nsites[system]-1; i++) {
         if (L[system][i] == i) {			// xtal site and not scanned yet

            j	=	i;
            pj	=	site[system][i].p;

	    /* could improve by adding cell list implementation here */

	    for (k=i+1; k<nsites[system]; k++) {
	       Lk	=	L[system][k];
	       if (Lk == k && DistSQ(pj, site[system][k].p, system) < Rconn2) {		// xtal site and not scanned yet and neighbor
		  L[system][k]	=	L[system][j]; 		// exchange label
		  L[system][j]	=	Lk;
               }
	    }
            j	=	L[system][j];
            pj	=	site[system][j].p;

	    while (j!=i) {				// this cluster not complete

               /* could improve by adding cell list implementation here */

               for (k=i+1; k<nsites[system]; k++) {
	          Lk	=	L[system][k];
	          if (Lk == k && DistSQ(pj, site[system][k].p, system) < Rconn2) {
		     L[system][k]	=	L[system][j];
		     L[system][j]	=	Lk;
                  }
	       }
               j	=	L[system][j];
               pj	=	site[system][j].p;
            }
         }
      }
   }	// done for all systems

   // analyze nuclei size distribution
   for (i=0; i<NSYSTEMS; i++) {
      MAXSIZE[i]	=	0;
      Nnucl[i]		=	0;
      Xtal[i]		=	0;
      secondNmax[i]	=	0;
   }
   for (i=0; i<NSITES+1; i++) {
      sizedist[i]	=	0;
      sizeofnucl[i]	=	0;
   }

   for (system=0; system<NSYSTEMS; system++) {	// so far only for system=0 due to sizedist and sizeofnucl, 5/28/08
      nuclid		=	1;		// xtal nuclei index starts from 1, not 0.
      for (i=0; i<nsites[system]; i++) {
         if (L[system][i] >= 0) {
            NIT		=	1;
	    LIT		=	L[system][i];
            L[system][i]	=	-1 * (L[system][i]+2);	// make L[system][i] < 0, clear this site

	    site[system][LIT].mol->nuclid[site[system][LIT].mol_site]	=	nuclid;

            while (LIT != i ) {
	       NIT	+=	1;
	       j	=	LIT;
	       LIT	=	L[system][LIT];
	       L[system][j]	=	-1 * (L[system][j]+2);
	       site[system][LIT].mol->nuclid[site[system][LIT].mol_site]	=	nuclid;
            }

            if (D_XTALSIZE)
	       PutInDistribution(D_Xtalsize+system, NIT, 1.0, 1.0);

            Xtal[system]	+=	NIT;
            Nnucl[system]	++;

	    sizeofnucl[nuclid]	=	NIT;
	    sizedist[NIT]	++;

            if (MAXSIZE[system] < NIT)
	       MAXSIZE[system]	=	NIT;
            
	    nuclid	++;
         }
      }
   }	// done for all systems
}
///______Find_Nuclei_LJ_______//////

///////////////////////////////////////////////////
/* Decide whether chains j and k are connected   */
/* based on Esselink's definition (JCPv101p9033) */
///////////////////////////////////////////////////
long connect(long j, long k, long clusdef)
{				
   double	r2;
   vector	pj, pk;
   matrix	Mj, Mk;
   vector	eigj, vj, eigk, vk;
   long		i;
   double	r2cut, anglecut;

   if ((mol+j)->box != (mol+k)->box)
      Exit("position", "connect", "not in the same box");

   r2cut	=	Rconn2;			// Rconn=1.5sigma
   anglecut	=	0.9848;			// cos 10
/*
   if (clusdef==1) {
      r2cut	=	Rconn2;			// Rconn=1.5sigma
      anglecut	=	0.9848;			// cos 10
   }
   else if (clusdef==2) {
      r2cut	=	Rconn2*1.44;		// 1.8^2 : 1.5^2
      anglecut	=	0.9659;			// 15 degrees
   }       
   else if (clusdef==3) {
      r2cut	=	Rconn2*0.75;		// 1.3^2 : 1.5^2
      anglecut	=	0.9962;			// 5 degrees
   } 
*/
   pj	=	CenterofMass(mol+j);		// chain j center of mass
   pk	=	CenterofMass(mol+k);
   r2	=	DistSQ(pj, pk, (mol+j)->box);

   if (r2 < r2cut) {
      Mj	=	InertiaTensor(mol+j);	// chain j moment of inertia tensor, w.r.t. center of mass
      Mk	=	InertiaTensor(mol+k);
      eigj	=	M_eig(Mj);
      eigk	=	M_eig(Mk);

      vj	=	V_eig(Mj, MIN(MIN(eigj.x, eigj.y), eigj.z));
      vk	=	V_eig(Mk, MIN(MIN(eigk.x, eigk.y), eigk.z));

      if ( fabs( V_Dot(&vj, &vk)/sqrt(V_Dot(&vj, &vj)*V_Dot(&vk, &vk)) ) > anglecut) {	// <= 10 degrees
/*
         printf("mol %d and %d are connected\n", j, k);
         printf("mol %d's coordinates:\n", j);
         for (i=0; i<(mol+j)->nsites; i++)
            printf("\t%f\t%f\t%f\t%f\n", (mol+j)->p[i].x, (mol+j)->p[i].y, (mol+j)->p[i].z, type[(mol+j)->type[i]].M);
         printf("mol %d's coordinates:\n", k);
         for (i=0; i<(mol+k)->nsites; i++)
            printf("\t%f\t%f\t%f\t%f\n", (mol+k)->p[i].x, (mol+k)->p[i].y, (mol+k)->p[i].z, type[(mol+k)->type[i]].M);

         printf("mol %d's center of mass\n", j);
	 V_Print(pj); 
         printf("mol %d's center of mass\n", k);
	 V_Print(pk); 

         printf("mol %d's coordinates:\n", j);
         for (i=0; i<(mol+j)->nsites; i++)
            printf("\t%f\t%f\t%f\t%f\n", (mol+j)->p[i].x-pj.x, (mol+j)->p[i].y-pj.y, (mol+j)->p[i].z-pj.z, type[(mol+j)->type[i]].M);
         printf("mol %d's coordinates:\n", k);
         for (i=0; i<(mol+k)->nsites; i++)
            printf("\t%f\t%f\t%f\t%f\n", (mol+k)->p[i].x-pk.x, (mol+k)->p[i].y-pk.y, (mol+k)->p[i].z-pk.z, type[(mol+k)->type[i]].M);

         printf("mol %d's gyration matrix\n", j);
	 M_Print(Mj); 
         printf("mol %d's gyration matrix\n", k);
	 M_Print(Mk); 
   
         printf("mol %d's eig\n", j);
	 V_Print(eigj);
         printf("mol %d's eig\n", k);
	 V_Print(eigk);

         printf("mol %d's main axis\n", j);
         V_Print(vj);
         printf("mol %d's main axis\n", k);
         V_Print(vk);
*/
	 return	1;
      }
      else {
	 return 0;
      }
   }
   else
      return	0;
}
//////____connect()____________________//////

#define SIZECAP	3

//////////////////////////////////////////////////
/* Added on 12/18/07, used by Esselink et al. 	*/
/* only for one system for now.			*/	
//////////////////////////////////////////////////
void Find_Nuclei(long clusdef)			
{					
   vector	pj;		
   long		i, j, k, n, jj, kk, NIT, LIT, Lk, nuclid;
   static long	*L;
   molstruct	*moli;
   double	r2;
   static long	init = 1;
   long		system = 0, nsites = NMOLS;

   if (init) {
      L			=	(long *) calloc(nsites, sizeof(long));	// only for one system now
      sizeofnucl	=	(long *) calloc(nsites+1, sizeof(long));
      sizedist		=	(long *) calloc(nsites+1, sizeof(long));

      if (L==NULL || sizeofnucl==NULL || sizedist==NULL)
         Exit("position", "Find_Nuclei", "out of memory");
      init	=	0;
   }
   for (i=0; i<nsites; i++) {
      L[i]		=	i;				// all crystal phase
      for (n=0; n<mol[i].nsites; n++)
         mol[i].nuclid[n]	=	-1;
   }

   for (i=0; i<nsites-1; i++) {
      if (L[i] == i) {					// Not scanned crystal particle
         j	=	i;

         for (k=i+1; k<nsites; k++) {			// search through the whole list

            Lk	=	L[k];

	    if (Lk == k) {				// also not scanned crystal particle

               if (connect(j, k, clusdef)) {		// if j and k are connected
		  L[k]	=	L[j];			// save the last particle in L[j]
		  L[j]	=	Lk;
	       }
	    }
	 }
 
         j	=	L[j];		

         while (j!=i) {

            for (k=i+1; k<nsites; k++) {		// search through the whole list again

	       Lk	=	L[k];

               if (Lk == k) {
		  if (connect(j, k, clusdef)) {
		     L[k]	=	L[j];
		     L[j]	=	Lk;
		  }
	       }
	    }

	    j	=	L[j];				// try another particle in the same cluster
	 } 
      }
   }

   /* Collect nuclei size distribution */

   for (i=0; i<NSYSTEMS; i++) {
      MAXSIZE[i]	=	0;
      Nnucl[i]		=	0;			// # of Xtal nuclei
      Xtal[i]		=	0;			// # of Xtal-like particles
   }
   for (i=0; i<nsites+1; i++) {		
      sizeofnucl[i]	=	0;			// initialize nuclei sizes
      sizedist[i]	=	0;			// initialize nuclei size distribution
   }
   
   nuclid	=	1;				// nuclei index starts with 1

   for (i=0; i<nsites; i++) {
      if (L[i] >= 0) {
	 NIT	=	1;
         LIT	=	L[i];
	 L[i]	=	-1 * (L[i]+2);			// clear this particle, but they can be easily recovered

         for (n=0; n<mol[LIT].nsites; n++)
            mol[LIT].nuclid[n]	=	nuclid;

         while (LIT != i) {
            NIT		+=	1;
            j		=	LIT;
            LIT		=	L[LIT];
            L[j]	=	-1 * (L[j]+2);
            
            for (n=0; n<mol[LIT].nsites; n++)
               mol[LIT].nuclid[n]	=	nuclid;
         }

         if (D_XTALSIZE)    
            PutInDistribution(D_Xtalsize+system, NIT, 1.0, 1.0);

	 Xtal[system]		+=	NIT;
	 if (NIT>SIZECAP)
	    Nnucl[system]		++;

	 sizeofnucl[nuclid]	=	NIT;
	 sizedist[NIT]		++;

	 if (MAXSIZE[system] < NIT) {
	    MAXSIZE[system]	=	NIT;
         }
	 nuclid		++;
      }
   }

   /* calculate more nuclei variables */

   for (i=0; i<10; i++)				// calc. the size of 10 biggest nuclei
      nmax[system][i]	=	0;
   k	=	0;			
   for (i=MAXSIZE[system]; i>0; i--) {
      for (j=1; j<=sizedist[i]; j++) {
         nmax[system][k]	=	i;
	 k	++;
         if (k>=10)
	    break;
      }
      if (k>=10)
         break;
   }

   realXtal[system]	=	Xtal[system];
   for (i=1; i<=SIZECAP; i++)
      realXtal[system]	-=	sizedist[i] * i;

   secondNmax[system]	=	0;
   if (sizedist[MAXSIZE[system]] > 1)
      secondNmax[system]	=	MAXSIZE[system];
   else {
      for (i=MAXSIZE[system]-1; i>SIZECAP; i--) {
         if (sizedist[i] > 0) {
	          secondNmax[system]	=	i;
	          break;
         }
      }
   }
}
/////_____Find_Nuclei()______/////


//////////////////////////////////////////////////////////////////
/* Added 2/23/09, Decide whether two sites are connected	*/
/* based on their local p2 value				*/
//////////////////////////////////////////////////////////////////
long connect_p2(long j, long k, long clusdef)
{
   long		sitej, sitek, sitespermol = NSITES/NMOLS;
   double	r2;
   vector	pj, pk;
   molstruct	*molj, *molk;

   long			i, n;
   double		cos2, angle;
   static long		dist[18];
   static double	mincos2=1;
   vector		vj, vk;

   if (mol[j/(NSITES/NMOLS)].box != mol[k/(NSITES/NMOLS)].box)
      Exit("position", "connect", "not in the same box");	// only one box for now

   molj		=	mol + j/sitespermol;
   molk		=	mol + k/sitespermol;
   sitej	=	mod(j, sitespermol);
   sitek	=	mod(k, sitespermol);
   pj	=	molj->p[sitej];
   pk	=	molk->p[sitek];
   //pj	=	mol[j/sitespermol].p[mod(j, sitespermol)];
   //pk	=	mol[k/sitespermol].p[mod(k, sitespermol)];

   r2	=	DistSQ(pj, pk, molj->box);
   if (r2 < Rconn2) {
/*
      if (sitej>0 && sitej<molj->nsites-1 && sitek>0 && sitek<molk->nsites-1) {
         vj		=	V_Subtr(molj->p+sitej+1, molj->p+sitej-1);			
         vk		=	V_Subtr(molk->p+sitek+1, molk->p+sitek-1);			
         cos2		=	V_Dot(&vj, &vk);
         cos2		*=	cos2;
         cos2		/=	V_Dot(&vj, &vj) * V_Dot(&vk, &vk);
         angle		=	acos(sqrt(cos2));
         dist[(int)(angle*180/M_PI/5)]	++;
         for (i=0; i<18; i++)
	    printf("%5d", dist[i]);
         printf("\n");
      } 
*/
      return	1;
   }
   else
      return	0;
}
/////////////////_____________________________________////////////

//////////////////////////////////////////////////////////////////////
/* Add on 8/18/2011, find out all the segments, xtal and amorphous  */
//////////////////////////////////////////////////////////////////////

void Find_segments()
{
   molstruct	*moli;
   long		i, j, k, nseg, previous;
   long		seg_start, seg_O, seg_A, seg_B; 
   long		length, len_A, len_B, typ_O, typ_A, typ_B;
   long		temp, max, maxid, xtal;
   long		done, round;
   static long	Lmin = 4;

   // Initialization

   for (i=0; i<MAXNMOLS; i++) {
      nsegment[i]	=	0;	// # of total segs on each chain
      for (j=0; j<MAXNMOLSITES; j++){
         seg_stat[i][j]	=	-1;	// seg stat, each seg needs 3 elem
      }
   }

   // Find segments for all chains, xtal seg and amorphous seg.

   for (moli=mol; moli<mol+NMOLS; moli++) {
      j		=	moli-mol;
      previous	=	moli->nuclid[0];
      nseg	=	0;
      seg_stat[j][nseg*3]	=	0;

      for (i=0; i<moli->nsites; i++) {
         if (moli->nuclid[i]  != previous) {
            seg_stat[j][nseg*3+2]	=	i-1;	// id of last bead of segment
            seg_stat[j][nseg*3+1]	=	previous;	// segment id
            nseg	++;

            seg_stat[j][nseg*3]	=	i;	// id of first bead of segment

            if (i==moli->nsites-1) {		// last site is a one bead segment
               seg_stat[j][nseg*3+2]	=	i;
               seg_stat[j][nseg*3+1]	=	moli->nuclid[i];
               nseg	++;
            }
         }
         else if (i==moli->nsites-1) {
            seg_stat[j][nseg*3+2]	=	i;
	    seg_stat[j][nseg*3+1]	=	moli->nuclid[i];
            nseg	++;
         }
         previous	=	moli->nuclid[i];
      }
      nsegment[j]	=	nseg;		// # of segments in this chain
   }
   //Seg_smooth();
   return;
}

///////////////////////////
/* Deal with fluctuation */
///////////////////////////

vector Seg_smooth(long nmaxid) {		// reduce segment fluctuation

   molstruct	*moli;
   long		i, j, k, nseg, previous;
   long		seg_start, seg_O, seg_A, seg_B; 
   long		length, len_A, len_B, typ_O, typ_A, typ_B;
   long		temp, max, maxid, xtal, xtalnew;
   long		done, round;
   static long	Lmin = 4;
   vector	result;

   max		=	0;			// look for chain with most segments
   maxid	=	0;
   for (j=0; j<NMOLS; j++) {
      if (nsegment[j] > max) {
         max	=	nsegment[j];
         maxid	=	j;
      }
   }

//   printf("max = %d maxid = %d\n", max, maxid);
   xtal	=	0;
   for (i=0; i<nsegment[maxid]; i++) {
      if (seg_stat[maxid][i*3+1] > 0)
         xtal	+=	seg_stat[maxid][i*3+2]-seg_stat[maxid][i*3]+1;
//      printf("%6d%6d%6d\n", seg_stat[maxid][i*3], seg_stat[maxid][i*3+1], seg_stat[maxid][i*3+2]);
   }
//   printf("total xtal segment = %d\n", xtal);

   for (j=0; j<NMOLS; j++) {			// for all chains
    done=0;
    round=0;
    while (done==0 && round<3) {
      for (i=0; i<nsegment[j]; i++) {		// Find 1st seg of Lmin long or longer
         if (seg_stat[j][i*3+2] - seg_stat[j][i*3]+1 >= Lmin) {
            seg_start	=	i;
            break;
         }
      }  	// There should be at least one such segment on any chain
      
      seg_O	=	seg_start;
      while (seg_O < nsegment[j] - 1) {		// proceed downward
         seg_A	=	seg_O + 1;
         seg_B	=	seg_O + 2;
         typ_O	=	seg_stat[j][seg_O*3+1];
         typ_A	=	seg_stat[j][seg_A*3+1];
         len_A	=	seg_stat[j][seg_A*3+2] - seg_stat[j][seg_A*3]+1;

         if (typ_A != typ_O && len_A < Lmin) { 

            if (seg_A == nsegment[j]-1 && len_A < Lmin-1) {	// A is the last seg
                seg_stat[j][seg_A*3+1]	=	typ_O;
            }
            else if (seg_A < nsegment[j]-1 && seg_stat[j][seg_B*3+1]==typ_O &&
                     seg_stat[j][seg_B*3+2]-seg_stat[j][seg_B*3]+1 >= Lmin) {	// convert
                seg_stat[j][seg_A*3+1]	=	typ_O;
            }
            else if (seg_A < nsegment[j]-1 && seg_stat[j][seg_B*3+1]==typ_O &&
                     seg_stat[j][seg_B*3+2]-seg_stat[j][seg_B*3]+1 < Lmin) {	// swap
                len_B	=	seg_stat[j][seg_B*3+2]-seg_stat[j][seg_B*3]+1;
                seg_stat[j][seg_A*3+2]	=	seg_stat[j][seg_A*3] + len_B-1;
                seg_stat[j][seg_B*3]	=	seg_stat[j][seg_B*3+2] - len_A+1;

                seg_stat[j][seg_A*3+1]	=	seg_stat[j][seg_B*3+1];
                seg_stat[j][seg_B*3+1]	=	typ_A;
            }
         } // if O, A different and len_A < Lmin
         seg_O	++;
      }	// proceed downward finished

      seg_O	=	seg_start;
      while (seg_O > 0) {	// proceed upward
         seg_A	=	seg_O - 1;
         seg_B	=	seg_O - 2;
         typ_O	=	seg_stat[j][seg_O*3+1];
	 typ_A	=	seg_stat[j][seg_A*3+1];
         len_A	=	seg_stat[j][seg_A*3+2] - seg_stat[j][seg_A*3]+1;

         if (typ_A != typ_O && len_A < Lmin) { 

            if (seg_A == 0 && len_A < Lmin-1) {
                seg_stat[j][seg_A*3+1]	=	typ_O;
            }
            else if (seg_A > 0 && seg_stat[j][seg_B*3+1]==typ_O &&
                     seg_stat[j][seg_B*3+2]-seg_stat[j][seg_B*3]+1 >= Lmin) {
                seg_stat[j][seg_A*3+1]	=	typ_O;
            }
            else if (seg_A > 0 && seg_stat[j][seg_B*3+1]==typ_O &&
                     seg_stat[j][seg_B*3+2]-seg_stat[j][seg_B*3]+1 < Lmin) {
                len_B	=	seg_stat[j][seg_B*3+2]-seg_stat[j][seg_B*3]+1;
                seg_stat[j][seg_A*3]	=	seg_stat[j][seg_A*3+2] - len_B+1;                
                seg_stat[j][seg_B*3+2]	=	seg_stat[j][seg_B*3] + len_A-1;

                seg_stat[j][seg_A*3+1]	=	seg_stat[j][seg_B*3+1];
                seg_stat[j][seg_B*3+1]	=	typ_A;
            }
         }
         seg_O	--;
      } // proceed upward finished

      // Resorting, basically combining same type of segments
      nseg	=	0;
      for (i=1; i<nsegment[j]; i++) {
         if (seg_stat[j][i*3+1] == seg_stat[j][nseg*3+1] )
            seg_stat[j][nseg*3+2] = seg_stat[j][i*3+2];
         else {
            nseg	++;
            seg_stat[j][nseg*3]		=	seg_stat[j][i*3];
            seg_stat[j][nseg*3+1]	=	seg_stat[j][i*3+1];
            seg_stat[j][nseg*3+2]	=	seg_stat[j][i*3+2];
         }
      }
      nsegment[j]	=	nseg+1;

      done	=	1;
      for (i=0; i<nsegment[j]; i++) {
         length	=	seg_stat[j][i*3+2] - seg_stat[j][i*3] + 1;

         if ( (i==0 || i==nsegment[j]-1) && length < Lmin-1) {
            done	=	0;
            break;
         }
         else if ( i>0 && i<nsegment[j]-1 && length < Lmin) {
            done	=	0;
            break;
         }
      }
      //printf("chain %d done = %d\n", j, done);
      round	++;
    }
    //printf("\n"); 
   }	// for all chains

//   printf("after resorting\n");
//   printf("max = %d maxid = %d\n", max, maxid);
   xtalnew	=	0;
   for (i=0; i<nsegment[maxid]; i++) {
      if (seg_stat[maxid][i*3+1] > 0)
         xtalnew	+=	seg_stat[maxid][i*3+2]-seg_stat[maxid][i*3]+1;
//      printf("%6d%6d%6d\n", seg_stat[maxid][i*3], seg_stat[maxid][i*3+1], seg_stat[maxid][i*3+2]);
   }
//   printf("total xtal segment = %d\n", xtalnew);

   result.x	=	xtalnew - xtal;

   xtalnew	=	0;
   for (j=0; j<NMOLS; j++) {
      for (i=0; i<nsegment[j]; i++) {
         if (seg_stat[j][i*3+1] == nmaxid ) {
            xtalnew	+=	seg_stat[j][i*3+2]-seg_stat[j][i*3]+1;
   }  }  }      
  
   result.y	=	xtalnew - nmax[0][0];

   return result;
}

///////////////////////////////////////////////////////////////////
/* 8/20/2011. Calculate cylinder shape based on segment analysis */
///////////////////////////////////////////////////////////////////
vector cylinder(beadstruct *nucleus, long size, long nuclid)
{
   static molstruct	**molx;
   molstruct		*moli;
   long			i, j, k, m, index, exist, max;
   long			first, firstid, nseg, sizenew;
   vector		rAB, rpre, rtemp, shape, director;
   static vector	*rtot;
   double		LAB, Ltot, thickness, radius, rho;
   static long		init=1;

   // Initialization

   if (init) {
      molx	=	(molstruct **) calloc(NMOLS, sizeof(molstruct *));
      rtot	=	(vector *) calloc(NMOLS, sizeof(vector));
      init	=	0;
   }
   nseg	=	0;
   Ltot	=	0;
   sizenew	=	0;
   V_Null(&director);
   for (i=0; i<NMOLS; i++) {
      V_Null(rtot+i);
   }

   // List all the chains participating this nucleus, stored in molx[]

   index	=	0;		// index is the # of chains in this nucleus
   for (i=0; i<size; i++) {
      exist	=	0;
      for (j=0; j<index; j++) {
         if (molx[j] ==  nucleus[i].moli) {	// site i and j might have same nucleus.mol
            exist	=	1;		// this mol has been recorded before
            break;
         }
      }
      if (!exist) {
         molx[index]	=	nucleus[i].moli;	// this mol hasn't been recorded
         index	++;
      }      
   }

   // Calculate thickness
/*      
   for (k=0; k<index; k++) {			// for all participating chains
      moli	=	molx[k];
      j		=	moli - mol;
      first	=	1;			// first segment of one chain in the nucleus
      for (i=0; i<nsegment[j]; i++) {
         if (seg_stat[j][i*3+1]==nuclid) {
            rAB		=	V_Subtr(moli->p+seg_stat[j][i*3+2], moli->p+seg_stat[j][i*3]);
            rtemp	=	rAB;
            LAB		=	sqrt(V_Dot(&rAB, &rAB));
            if (first) {
               first		=	0;
            }
            else {
               if (V_Dot(&rAB, &rpre) < 0) {
                  rtemp	=	V_Mult(-1.0, &rAB);
                  if (size<=10) {				// small nuclei
                     printf("Cylinder error: folding for size <=10\n");
                  }     
               }
               else if (size > 10) {				// for big nuclei
                  printf("cylinder error: two nearby segments pointing same direction. chainid=%d segment=%d\n", j, i);
                  for (m=0; m<nsegment[j]; m++) {
                      printf("%d %d %d\n", seg_stat[j][m*3], seg_stat[j][m*3+1], seg_stat[j][m*3+2]);
                  }
               }
            }
            rpre	=	rAB;
            rtot[k]	=	V_Add(rtot+k, &rtemp);
            Ltot	+=	LAB;
            sizenew	+=	seg_stat[j][i*3+2] - seg_stat[j][i*3] + 1;
            nseg	++;
         }
      }
   }
   for (k=0; k<index; k++) {
      director	=	V_Add(&director, rtot+k);
   }
   thickness	=	Ltot/nseg;
*/
   max		=	0;
   for (k=0; k<index; k++) {
      moli	=	molx[k];
      j		=	moli - mol;
      for (i=0; i<nsegment[j]; i++) {
         if (seg_stat[j][i*3+1]==nuclid) {
            if (seg_stat[j][i*3+2] - seg_stat[j][i*3] + 1 > max) {
               max	=	seg_stat[j][i*3+2] - seg_stat[j][i*3] + 1;
               rAB	=	V_Subtr(moli->p+seg_stat[j][i*3+2], moli->p+seg_stat[j][i*3]);
            }
         }
      }
   }
   thickness	=	sqrt(V_Dot(&rAB, &rAB));
   rho		=	((double)NSITES) / BOX[0].vol;
   sizenew	=	size;
   radius	=	((double)sizenew)/ (M_PI*thickness*rho);//size??
   radius	=	sqrt(radius);

   shape.x	=	(double)index;
   shape.y	=	thickness;
   shape.z	=	radius;

   return	shape;
}

////////////////////////////////////////////////////////////////////////////////
/* Add on 4/20/2010, calculate the length and radius of a cylindrical nucleus */
////////////////////////////////////////////////////////////////////////////////
#define	AVELMIN	1	// min. length of stem to be taken into thickness average

vector cylindershape(beadstruct *nucleus, long size, long nuclid)
{
   // only one system for now 4/17/2010
   static molstruct 	**molx;		// molecules that belong to this nucleus
   molstruct		*moli;
   static long		*A, *B;		// head/tail bead id of Xtal segment on one chain
   long			i, j, n, index, maxid, exist, temp;
   static double	*l;
   double		lave, lmax, radius, rho;
   static vector	*rAB;
   vector		rABave;
   static long		init=1;
   vector		shape;

   if (init) {
      molx	=	(molstruct **) calloc(NMOLS, sizeof(molstruct *));
      l		=	(double *) calloc(NMOLS, sizeof(double));
      A		=	(long *) calloc(NMOLS, sizeof(long));
      B		=	(long *) calloc(NMOLS, sizeof(long));
      rAB	=	(vector *) calloc(NMOLS, sizeof(vector));

      init	=	0;
   }

   /* list all the chains participating this nucleus, stored in molx[] */

   index	=	0;		// index is the # of chains in this nucleus
   for (i=0; i<size; i++) {
      exist	=	0;
      for (j=0; j<index; j++) {
         if (molx[j] ==  nucleus[i].moli) {	// site i and j might have same nucleus.mol
            exist	=	1;		// this mol has been recorded before
            break;
         }
      }
      if (!exist) {
         molx[index]	=	nucleus[i].moli;	// this mol hasn't been recorded
         index	++;
      }      
   }

   // calculate how many beads in each chain that are xtal
   lmax	=	0.0;
   for (i=0; i<index; i++) {
      l[i]	=	0.0;
      moli	=	molx[i];
      for (j=0; j<moli->nsites; j++) {
         if (moli->nuclid[j] == nuclid) {
            l[i]	+=	1.0;
         }
      }
      if (l[i]>lmax) {			// find the longest stem
         lmax	=	l[i];
         maxid	=	i;		// maxid is the index in molx
      }
   }

   // find the end beads in each chain that are in this nucleus
   for (i=0; i<index; i++) {
      moli	=	molx[i];
      for (j=0; j<moli->nsites; j++) {
         if (moli->nuclid[j] == nuclid) {
            A[i]	=	j;	
            break;
         }
      }
      for (j=moli->nsites-1; j>=0; j--) {
         if (moli->nuclid[j] == nuclid) {
            B[i]	=	j;
	    break;
         }
      }
      rAB[i]	=	V_Subtr(moli->p+A[i], moli->p+B[i]);
   }

   // calculate the average length and average radius of a cylinder model
   n	=	0;
   V_Null(&rABave);
   for (i=0; i<index; i++) {
      moli	=	molx[i];
      if (V_Dot(rAB+i, rAB+maxid) <0) {
         temp	=	A[i];
         A[i]	=	B[i];
         B[i]	=	temp;
         rAB[i]	=	V_Subtr(moli->p+A[i], moli->p+B[i]);
      }
      if (l[i] >= AVELMIN) {
         rABave	=	V_Add(&rABave, rAB+i);
         n	++;
      }
   }

   rABave	=	V_Mult(1.0/n, &rABave);	// thickness of this nucleus
   lave		=	sqrt(V_Dot(&rABave, &rABave));

   rho		=	((double)NSITES) / BOX[0].vol;
   radius	=	((double)size)/ (M_PI*lave*rho);
   radius	=	sqrt(radius);
/*
  printf("nmax=%d\n", nmax);
  printf("index=%d\n", index);
  for (i=0; i<index; i++) { printf("l[%d]=%f  A=%d  B=%d\n",i, l[i], A[i], B[i]); }
  printf("maxid=%d\n", maxid);
  printf("rho=%f\n", rho);
  printf("lave=%f\n", lave);
  printf("radius=%f\n", radius); 
*/
   shape.x	=	(double)index;
   shape.y	=	lave;
   shape.z	=	radius;

   return	shape;
}

//////////////////////////////////////////////////////////////////
/* Added on 2/23/09, find nuclei based on local p2 of each bead	*/
//////////////////////////////////////////////////////////////////
void Find_Nuclei_p2(long clusdef)
{				
   vector	pj;
   long		i, j, k, m, n, jj, kk, NIT, LIT, Lk, nuclid;
   static long	*L;
   molstruct	*moli;
   double	r2, temp;
   static long	init=1;
   long		system = 0, nsites, sitespermol = NSITES/NMOLS, size, sum;
   static long	*id;
   vector	nreduce;

   // variable for shape analysis
   static FILE		*fPtr;
   vector		shape;
   vector		evalue;			// eigenvalue of gytensor
   matrix		gytensor;		// tensor of gyration of one nucleus
   static long		*number;		// count for average
   static double	*nchain;		// # of chains participate in one nucleus
   static double	*thickness, *radius;	// cylindrical model parameters
   static double	*asphericity;		// asphericity, calc. from evalue

   nsites	=	NSITES;				// definition based on beads

   if (init) {
      init	=	0;

      L		=	(long *) calloc(nsites, sizeof(long));	// only for one system now
      sizeofnuclp2	=	(long *) calloc(nsites+1, sizeof(long));
      sizedistp2	=	(long *) calloc(nsites+1, sizeof(long));

      if (L==NULL || sizeofnuclp2==NULL || sizedistp2==NULL)
         Exit("position", "Find_Nuclei", "out of memory");

      fPtr	=	fopen("shape.out", "w");
      fprintf(fPtr, "##### Nucleus shape analysis output from Find_Nuclei_p2 #####\n");

      // Allocate cylindrical model variables
      number		=	(long *) calloc(NSITES, sizeof(long));
      nchain		=	(double *) calloc(NSITES, sizeof(double));
      thickness		=	(double *) calloc(NSITES, sizeof(double));
      radius		=	(double *) calloc(NSITES, sizeof(double));
      asphericity	=	(double *) calloc(NSITES, sizeof(double));

      id	=	(long *) calloc(NSITES, sizeof(long));
   }
   for (i=0; i<nsites; i++) {
      if (mol[i/sitespermol].p2[mod(i, sitespermol)] > critp2)
         L[i]		=	i;				// crystal phase
      else
	 L[i]		=	-1;				// liquid phase

//      for (n=0; n<mol[i].nsites; n++)
//         mol[i].nuclid[n]	=	-1;
      mol[i/sitespermol].nuclid[mod(i,sitespermol)]	=	-1;
   }

   for (i=0; i<nsites-1; i++) {
      if (L[i] == i) {					// Not scanned crystal particle
         j	=	i;

         for (k=i+1; k<nsites; k++) {			// search through the whole list

            Lk	=	L[k];

	    if (Lk == k) {				// also not scanned crystal particle

               if (connect_p2(j, k, clusdef)) {		// if j and k are connected
		  L[k]	=	L[j];			// save the last particle in L[j]
		  L[j]	=	Lk;
	       }
	    }
	 }
 
         j	=	L[j];		

         while (j!=i) {

            for (k=i+1; k<nsites; k++) {		// search through the whole list again

	       Lk	=	L[k];

               if (Lk == k) {
		  if (connect_p2(j, k, clusdef)) {
		     L[k]	=	L[j];
		     L[j]	=	Lk;
		  }
	       }
	    }

	    j	=	L[j];				// try another particle in the same cluster
	 } 
      }
   }

   /* Collect nuclei size distribution */

   for (i=0; i<NSYSTEMS; i++) {
      MAXSIZE[i]	=	0;
      Nnucl[i]		=	0;			// # of Xtal nuclei
      Xtal[i]		=	0;			// # of Xtal-like particles
   }
   for (i=0; i<nsites+1; i++) {		
      sizeofnuclp2[i]	=	0;			// initialize nuclei sizes
      sizedistp2[i]	=	0;			// initialize nuclei size distribution
   }

   nuclid	=	1;				// nuclei index starts with 1

   for (i=0; i<nsites; i++) {
      if (L[i] >= 0) {
	 NIT	=	1;
         LIT	=	L[i];
	 L[i]	=	-1 * (L[i]+2);	// clear this particle, but they can be easily recovered

//         for (n=0; n<mol[LIT].nsites; n++)
//            mol[LIT].nuclid[n]	=	nuclid;
         mol[LIT/(NSITES/NMOLS)].nuclid[mod(LIT,NSITES/NMOLS)]	=	nuclid;	// should be more general

         while (LIT != i) {
            NIT		+=	1;
            j		=	LIT;
            LIT		=	L[LIT];
            L[j]	=	-1 * (L[j]+2);
            
//            for (n=0; n<mol[LIT].nsites; n++)
//               mol[LIT].nuclid[n]	=	nuclid;
            mol[LIT/(NSITES/NMOLS)].nuclid[mod(LIT, NSITES/NMOLS)]	=	nuclid;
         }

         if (D_XTALSIZE)    
            PutInDistribution(D_Xtalsize+system, NIT, 1.0, 1.0);

	 Xtal[system]		+=	NIT;
//	 if (NIT>SIZECAP)
	    Nnucl[system]		++;

	 sizeofnuclp2[nuclid]	=	NIT;
	 sizedistp2[NIT]		++;

	 if (MAXSIZE[system] < NIT) {
	    MAXSIZE[system]	=	NIT;
         }
	 nuclid		++;
      }
   }

   /* Collect more nuclei variables */

   for (i=0; i<10; i++)				// store the size of 10 biggest nuclei
      nmax[system][i]	=	0;
   k	=	0;			
   for (i=MAXSIZE[system]; i>0; i--) {
      for (j=1; j<=sizedistp2[i]; j++) {
         nmax[system][k]	=	i;
	 k	++;
         if (k>=10)
	    break;
      }
      if (k>=10)
         break;
   }

   realXtal[system]	=	Xtal[system];
   for (i=1; i<=SIZECAP; i++)
      realXtal[system]	-=	sizedistp2[i] * i;

   secondNmax[system]	=	0;
   if (sizedistp2[MAXSIZE[system]] > 1)
      secondNmax[system]	=	MAXSIZE[system];
   else {
      for (i=MAXSIZE[system]-1; i>SIZECAP; i--) {
         if (sizedistp2[i] > 0) {
	          secondNmax[system]	=	i;
	          break;
         }
      }
   }

   // Sort nucleus id based on their size, from small to big

   for (i=0; i<NSITES; i++)
      id[i]	=	0;
   for (i=1; i<=Nnucl[system]; i++) {		// search for every nucleus
      size	=	sizeofnuclp2[i];	// find its size

      sum	=	0;
      for (j=0; j<size; j++) {
         sum	+=	sizedistp2[j];
      }
      while (id[sum] !=0) {
         sum	++;
      }
      id[sum]	=	i;    
   }
/*
   printf("Nuclid size (in the order of increasing nucleus size)\n");
   for (i=0; i<Nnucl[system]; i++) {
      printf("%d  %d\n", id[i], sizeofnuclp2[id[i]]);
   }
*/
   /* Store all beads that belong to the biggest nucleus in a list */
   /* If there are two nuclei both of nmax, then only store one */

   Find_segments();
   nreduce = Seg_smooth(id[Nnucl[system]-1]);

   k	=	0;
   //for (nuclid=1; nuclid<=Nnucl[system]; nuclid++) {		// for all nuclei
   for (j=0; j<Nnucl[system]; j++) {				// for all nuclei
      nuclid	=	id[j];      

      // Group the sites for each nucleus

      size	=	0;			
      for (moli=mol; moli<mol+NMOLS; moli++) {
         for (i=0; i<moli->nsites; i++) {
            if (moli->nuclid[i] == nuclid) {
               nucleus[size].moli	=	moli;
               nucleus[size].site	=	i;
               size	++;
      }  }  }

      // Sanity check

      if (size != sizeofnuclp2[nuclid] ) {
         printf("Error, size %3d != sizeofnuclp2[nuclid%d] %3d\n",size, nuclid, sizeofnuclp2[nuclid]);
         exit(0);
      }

      // Shape analysis for every nucleus
/*
      if (size>10 && k==0) {		// could change the real size of a nucleus, read cylinder().
         Seg_smooth();			// Segment smoothing for cylinder shape characterization.
         k	=	1;		// Do it only once for each snapshot.
      }
*/
      //shape	=	cylindershape(nucleus, size, nuclid);
      shape	=	cylinder(nucleus, size, nuclid);
      gytensor	=	groupGyraTensor(nucleus, size);
      evalue	=	M_eig(gytensor);
      evalue	=	fshape(&evalue);		// calculate asphericity

      // Average

      number[size]	++;
      nchain[size]	+=	shape.x;
      thickness[size]	+=	shape.y;
      radius[size]	+=	shape.z;
      asphericity[size]	+=	evalue.z;
   }

   fprintf(fPtr, "\n********** Cylinder Model Analysis (%d < nmax < %d)**********\n", 
		(nmax[0][0]/100)*100, (nmax[0][0]/100+1)*100);
   fprintf(fPtr, "Rp = %f\tRconn = %f\tcritp2 = %f\n", Rp, Rconn, critp2);
   fprintf(fPtr, "nmax = %d nreduce = %d nmaxreduce = %d\n", nmax[0][0], (long)nreduce.x, (long)nreduce.y);
   fprintf(fPtr, "size\t count\t P(n)\t nchains\t thickness\t radius\t asphericity\n");

   for (size=1; size<=nmax[0][0]; size++) {	// increase by one
      if (number[size]==0) {
         fprintf(fPtr, "%d\t %d\t %f\t nan\t nan\t nan\t nan\n", size, number[size], (double)number[size]/number[1]);
      }
      else {
         temp	=	1.0/number[size];
         fprintf(fPtr,"%d\t %d\t %f\t %f\t ", size, number[size], (double)number[size]/number[1], nchain[size]*temp);
         fprintf(fPtr,"%f\t %f\t %f\n", thickness[size]*temp, radius[size]*temp, asphericity[size]*temp);
      }
   }
/*
   fprintf(fPtr, "----------increase by three-----------\n");
   for (size=3; size<=nmax[0][0]; size+=3) {	// increase by three
      if (number[size]==0) {
         fprintf(fPtr, "%d\t %d\t nan\t nan\t nan\t nan\n", size, number[size]);
      }
      else {
         temp	=	1.0/number[size];
         fprintf(fPtr,"%d\t %d\t %f\t %f\t %f\t %f\n", size, number[size], nchain[size]*temp, thickness[size]*temp, radius[size]*temp, asphericity[size]*temp);
      }
   }
*/
   fflush(fPtr);
/*
   nuclid	=	-1;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         if (sizeofnuclp2[moli->nuclid[i]] == nmax[0][0]) {	// find one bead in the biggest
            nuclid	=	moli->nuclid[i];		// nucleus and record the nuclid
	    break;
	 }
         if (nuclid!=-1)	break;
      }
   }

   size		=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         if (moli->nuclid[i] == nuclid) {	// find one bead in the biggest
            nucleus[size].moli	=	moli;
	    nucleus[size].site	=	i;
	    size	++;
         }
      }
   }

   shape	=	cylindershape(nucleus, nmax[0][0], nuclid);
   printf("%8.3f %8.3f %8.3f\n", shape.x, shape.y, shape.z); 
*/

   return;
}
////////__________Find_Nuclei_p2()___________________/////////////

//////////////////////////////////////////////////////////////////
/* Added on 10/7/09, find connections based on chord vectors	*/
//////////////////////////////////////////////////////////////////
void Find_Nuclei_p2_new(long clusdef)
{				
   vector	pj;
   long		i, j, k, n, jj, kk, NIT, LIT, Lk, nuclid;
   static long	*L;
   molstruct	*moli;
   double	r2;
   static long	init=1;
   long		system = 0, nsites, sitespermol = NSITES/NMOLS;

   nsites	=	NSITES;				// definition based on beads

   if (init) {
      L		=	(long *) calloc(nsites, sizeof(long));	// only for one system now
      sizeofnuclp2	=	(long *) calloc(nsites+1, sizeof(long));
      sizedistp2	=	(long *) calloc(nsites+1, sizeof(long));

      if (L==NULL || sizeofnuclp2==NULL || sizedistp2==NULL)
         Exit("position", "Find_Nuclei", "out of memory");
      init	=	0;
   }
   for (i=0; i<nsites; i++) {
      if (mol[i/sitespermol].p2[mod(i, sitespermol)] > -1)	// all crystal
         L[i]		=	i;				// crystal phase
      else
	 L[i]		=	-1;				// liquid phase

//      for (n=0; n<mol[i].nsites; n++)
//         mol[i].nuclid[n]	=	-1;
      mol[i/sitespermol].nuclid[mod(i,sitespermol)]	=	-1;
   }

   for (i=0; i<nsites-1; i++) {
      if (L[i] == i) {					// Not scanned crystal particle
         j	=	i;

         for (k=i+1; k<nsites; k++) {			// search through the whole list

            Lk	=	L[k];

	    if (Lk == k) {				// also not scanned crystal particle

               if (connect_p2_new(j, k, clusdef)) {	// if j and k are connected
		  L[k]	=	L[j];			// save the last particle in L[j]
		  L[j]	=	Lk;
	       }
	    }
	 }
 
         j	=	L[j];		

         while (j!=i) {

            for (k=i+1; k<nsites; k++) {		// search through the whole list again

	       Lk	=	L[k];

               if (Lk == k) {
		  if (connect_p2_new(j, k, clusdef)) {
		     L[k]	=	L[j];
		     L[j]	=	Lk;
		  }
	       }
	    }

	    j	=	L[j];				// try another particle in the same cluster
	 } 
      }
   }
}

long	connect_p2_new(long j, long k, long clusdef)
{
   vector	pj, pk;
   double	r2, costheta;
   long		sitespermol = NSITES/NMOLS;

   if (mol[j/(NSITES/NMOLS)].box != mol[k/(NSITES/NMOLS)].box)
      Exit("position", "connect", "not in the same box");	// only one box for now

   pj	=	mol[j/sitespermol].p[mod(j, sitespermol)];
   pk	=	mol[k/sitespermol].p[mod(k, sitespermol)];

   r2	=	DistSQ(pj, pk, mol[j/sitespermol].box);
    
   if (r2 < Rconn2)
      return	1;
   else
      return	0;
}

//////////////////////////////////////////////////////////
/* CoM_MaxNucleus(): return the center of mass of one	*/
/* (if there are more than one biggest nuclei in the	*/
/* central box).  11/15/08				*/
//////////////////////////////////////////////////////////

vector CoM_MaxNucleus(long system)
{				
   vector	com, 
		rA,			// center of one chain in nucleus
		rO,			// center of nucleus 
		rBA, 			// rB-rA, B and A belong to same nucleus
		rOA;			// rO-rA
   molstruct	*moli;
   long		id, n;

   // Step 1: find one chain A that belongs to the largest nucleus as a reference

   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->box == system && sizeofnucl[moli->nuclid[0]] == MAXSIZE[system]) {
							// find one molecule in the biggest nucleus
	  id	=	moli->nuclid[0];		// nucleus id
	  rA	=	CenterofMass(moli); 		// take one chain as reference point
	  rA	=	MapInBox2(&rA, PBC, system);	// map back to the central box
          break;
      }
   }

   // Step 2: calc. the shift of center of nucleus to this chain

   n	=	0;					// # of chains in this nucleus
   V_Null(&rBA);
   V_Null(&rOA);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->box==system && moli->nuclid[0] == id) {
	 n	++;
	 com	=	CenterofMass(moli);
	 rBA	=	V_Subtr(&com, &rA);
	 rBA	=	MapInBox2(&rBA, PBC, system);
	 rOA	=	V_Add(&rOA, &rBA);
      }
   }
   rOA	=	V_Mult(1.0/n, &rOA);

   // Step 3: calc. the center of nucleus of this nucleus

   rO	=	V_Add(&rA, &rOA);			// center of nucleus
   return	rO;
}	
/////_____CoM_MaxNucleus()___________________________/////	 

#ifdef TEST
void Find_Nuclei()					// Allen and Tildesley, F.34
{							// ref. Stoddard J Comp Phys, 27, 291, 1977
   vector	pj;					// done 9/11/2007
   long		i, j, k,
		icell, jcell, jj, kk,
		NIT, LIT,
   		L[NPARTS],				// linked list
		Lk,
		nuclid;
   double	r2;

   for (i=0; i<NPARTS; i++) {				// sorting linked list initialization
      if (part[i].nconnect < critconnect) {
	 L[i]	=	-1;				// liquid-like particles
      }
      else {
	 L[i]	=	i;
      }
   }

   for (i=0; i<NPARTS-1; i++) {

      if (L[i] == i) {					// Not scanned crystal particle
	 j	=	i;
	 pj	=	part[j].p;

#ifdef CELL_LIST
         icell	=	part[i].icell;
         for (jj=0; jj<Cell[icell].nneigh; jj++) {
	    jcell	=	Cell[icell].neigh[jj];

	    for (kk=0; kk<Cell[jcell].sites; kk++) {
	       k	=	Cell[jcell].list[kk];

	       if (k!=i) {
#else

         for (k=i+1; k<NPARTS; k++) {			// search through the whole list
	 {{

#endif	/* CELL_LIST */

            Lk	=	L[k];

	    if (Lk == k) {				// also not scanned crystal particle

	       r2	=	DistSQ(pj, part[k].p, part[k].box);
               if (r2 <= Rb2) {
		  L[k]	=	L[j];			// save the last particle in L[j]
		  L[j]	=	Lk;
	       }
	    }
	 }
	 }}
 
         j	=	L[j];		
         pj	=	part[j].p;

         while (j!=i) {

#ifdef CELL_LIST

            icell	=	part[j].icell;
            for (jj=0; jj<Cell[icell].nneigh; jj++) {
               jcell	=	Cell[icell].neigh[jj];

	       for (kk=0; kk<Cell[jcell].sites; kk++) {
	          k	=	Cell[jcell].list[kk];

	          if (k!=j) {
#else

            for (k=i+1; k<NPARTS; k++) {			// search through the whole list again
            {{

#endif	/* CELL_LIST */

	       Lk	=	L[k];

               if (Lk == k) {

		  r2	=	DistSQ(pj, part[k].p, part[k].box);
		  if (r2 <= Rb2) {
		     L[k]	=	L[j];
		     L[j]	=	Lk;
		  }
	       }
	    }
	    }}

	    j	=	L[j];				// try another particle in the same cluster
	    pj	=	part[j].p;
	 } 
      }
   }

   /* Analyze nuclei size distribution */

   MAXSIZE	=	0;
   Nnucl	=	0;				//number of Xtal nuclei
   Xtal		=	0;				//number of Xtal-like particles
   for (i=0; i<NPARTS+1; i++) {		
      sizeofnucl[i]	=	0;			// initialize nuclei sizes
      sizedist[i]	=	0;			// initialize nuclei size distribution
   }
   
   nuclid	=	1;				// nuclei index starts with 1

   for (i=0; i<NPARTS; i++) {
      if (L[i] >= 0) {
	 NIT	=	1;
         LIT	=	L[i];
	 L[i]	=	-1 * (L[i]+2);			// clear this particle, but they can be easily recovered

         while (LIT != i) {
            NIT		+=	1;
            j		=	LIT;
            LIT		=	L[LIT];
            L[j]	=	-1 * (L[j]+2);
         }
         sizeofnucl[nuclid]	=	NIT;
	 Xtal		+=	NIT;
	 Nnucl		++;
         sizedist[NIT]	++;
	 if (MAXSIZE < NIT) {
	    MAXSIZE	=	NIT;
         }
	 nuclid		++;
      }
   }
}
#endif

#ifdef TEST
void Find_Nuclei2()				// using iteration rather than recursion, some problems
{
   long		i, ii, j, jj, k;
   long		icell, jcell;
   long		id, idold, idnew;

   for (i=0; i<NPARTS; i++) {
      part[i].nuclid2	=	-1;
   }
   for (i=0; i<NPARTS+1; i++) {		
      sizeofnucl2[i]	=	0;		// initialize nuclei sizes
      sizedist2[i]	=	0;		// initialize nuclei size distribution
   }
   id	=	1;				// nucleus index starts from 1, rather than 0

   for (i=0; i<NPARTS; i++) {			//search over the particles
      
      if (part[i].nconnect >= critconnect) {		//this one is solid-like

         if (part[i].nuclid2 == -1) { 			//not identified yet
            part[i].nuclid2	=	id;		//assign new nucleus index #id
	    idnew		=	id;
            sizeofnucl2[id]	++;			//nucleus #id size increases by one
         }
	 else {
	    idnew	=	part[i].nuclid2;	//been identified before, pick up its nucleus index
	 }

#ifdef VERLET_LIST   
         for (jj=0; jj<part[i].nverlet; jj++) {		//check its neighbors
            j	=	part[i].vlist[jj];
            {
#elif CELL_LIST
         icell	=	part[i].icell;
	 for (ii=0; ii<Cell[icell].nneigh; ii++) {
	    jcell	=	Cell[icell].neigh[ii];

	    for (jj=0; jj<Cell[jcell].sites; jj++) {
	       j	=	Cell[jcell].list[jj];

#else
	 for (j=0; j<NPARTS; j++) {
	 {
#endif
	    if (j>i && part[j].nconnect >=critconnect && part[j].nuclid2!=idnew && DistSQ(part[i].p, part[j].p, part[k].box) < Rb2) {	
							//neighbor is solid-like, not same nucleus, and close enough

	       if (part[j].nuclid2 == -1) {		//if neighbor hasn't been identified
                  part[j].nuclid2	=	idnew;	//identify this neighbor
		  sizeofnucl2[id]	++;		//nucleus #id size increases by one
	       }
               else {
		  idold	=	part[j].nuclid2;		//if has been identified as #idold

		  for (k=0; k<NPARTS; k++) {		//change all particles with #idold to #id
		     if (part[k].nuclid2 == idold) {
			part[k].nuclid2	=	idnew;
			sizeofnucl2[idold]	--;
			sizeofnucl2[idnew]	++;
		     }
		  }
               }
	    }
	 }}

	 id	++;		//nucleus index increase by one
      }
   }

   MAXSIZE2	=	0;
   Nnucl2	=	0;		//number of Xtal nuclei
   Xtal2	=	0;		//number of Xtal-like particles
   for (id=1; id<NPARTS+1; id++) {		//determine the max nucleus size
      if (sizeofnucl2[id] != 0) {

         sizedist2[sizeofnucl2[id]]	++;
         Nnucl2	++;
         Xtal2	+=	sizeofnucl2[id];

         if (MAXSIZE2 < sizeofnucl2[id]){
	    MAXSIZE2	=	sizeofnucl2[id];
         }
      }
   }

   return;
}
#endif
/**********************************************************************************************************/
                                                                                                                                                                                            src/random.c                                                                                        0000600 0143352 0000144 00000002337 10714730275 012513  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    Name:       random.h
    Author:     Peng Yi at MIT
    Date:       October 23, 2006
    Purpose:    Using code from 'Numerical Recipes', chapter 7.
*/
#define __RANDOM_MODULE
#include "random.h"

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) 
    {
      if (-(*idum)< 1) *idum=1;
      else *idum = -(*idum);
      for (j=NTAB+7; j>=0; j--)
	{
	  k=(*idum)/IQ;
	  *idum=IA*(*idum-k*IQ) - IR*k;
	  if(*idum < 0) *idum += IM;
	  if (j < NTAB) iv[j] = *idum;
	}
      iy=iv[0];
    }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j]=*idum;
  if ((temp=AM*iy)> RNMX) return RNMX;
  else return temp;
}


double gauss(double std, double mean)		// Frenkel and Smit, Algorithm 44
{
   double	r, x, v1, v2;
   
   r	=	2.0;
   while (r >= 1.0) {
      v1	=	2.0 * ran1(seed) - 1.0;
      v2	=	2.0 * ran1(seed) - 1.0;
      r		=	v1*v1 + v2*v2;
   }
   x	=	v1 * sqrt(-2.0*log(r)/r);
   x	=	mean + std * x;
   return	x;
}




                                                                                                                                                                                                                                                                                                 src/rebridge.c                                                                                      0000644 0143352 0000144 00000032024 11306565150 013016  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	rebridge.c
    author:	Pieter J. in 't Veld for MIT Boston
    date:	May 4, 2001, January 28, 2002.
    purpose:	Trimer rebridging module, using Mavrantzas et al.,
    		Macromolecules 32, 5072 (1999).

    Notes:
      20010514	Check if spherical molecular notation corresponds with
		Mavrantas notation.
      20020128	Fixed angular incongruities in spherical coordinates.
      20090421  Mavrantas sphercal coordinates notation is different
		from the one used in Allen & Tiledesley, and the latter
		is adopted in our code.  Therefore, we need to adjust
		the implementation of Mavrantas's formulas in our coding.
		The difference is this, Mavrantas specified the bond angle
		theta_i as the angle at bead i, while A&T specified the bond
		angle theta_i as the angle at bead i-1. (PY)
		Therefore in function RebridgeSetup, A&T notation is used.
		Thus the Mavrantas's formula needs to be adjusted, e.g., the
		calculation of r_P in function Feasible().
      20091205  Added Jacobian calculation J(IV->I).

    Reference:  1. Macromolecules 1999, v32, 5072 
	        2. J.Chem.Phys. 2002, v117, 5465 
*/
#define __REBRIDGE_MODULE
#include "rebridge.h"


long Frame(
	vector *p0, vector *p1, vector *r_ref, 
	vector *x1, vector *x2, vector *x3)
{
  vector		t1, t2;
  double		f;
    
  t1			= V_Subtr(r_ref, p0);
  t2			= V_Subtr(p1, p0);
  if (!(f = sqrt(V_Dot(&t2, &t2)))) return 1;
  *x1			= V_Mult(1.0/f, &t2);
  *x2			= V_Mult(V_Dot(&t1, &t2), x1);
  *x3			= V_Cross(&t1, &t2);
  t2			= V_Mult(-f, &t1);
  *x2			= V_Add(x2, &t2);
  if (!(f = sqrt(V_Dot(x3, x3)))) return 1;
  *x2			= V_Mult(1.0/f, x2); 
  *x3			= V_Mult(1.0/f, x3);
  return 0;
}


double			R; // radius of C3
vector			u1, u2, u3,
			v1, v2, v3,
			w1, w2, w3;  // coordinate frame
vector			rN;

long Feasibility(
	vector *r, double l, double a, vector *u1, 
	double *phi_low, double *phi_high) // returns n_entries in phi interval
{
  static long		n_roots, n_real, i, j;
  static double		f, f2,
    			A, B, C, D, E, F,
  			c0, c1, c2, c3, c4, c5,
			c[5], phi[4];
  static double complex	z[4];
  static vector		t1;
    
  t1			= V_Subtr(r, &rN);
  c0			= V_Dot(&t1, &w2);
  c1			= a*V_Dot(&w2, u1);
  f			= 2/R;
  f			*= f;
  A			= f*(c0*c0+c1*c1);
  c2			= V_Dot(&t1, &w3);
  c3			= a*V_Dot(&w3, u1);
  B			= f*(c2*c2+c3*c3);
  C			= 2.0*f*(c1*c3+c0*c2);
  c4			= V_Dot(&t1, &t1);
  c5			= 2.0*a*V_Dot(&t1, u1);
  f			/= R;
  f2			= c4+R*R;
  F			= 0.25*f/R*(f2*(f2-2.0*(a*a+l*l))+
    				(a*a-l*l)*(a*a-l*l)+c5*c5);
  f2			-= a*a+l*l;
  D			= f*(f2*c0+c1*c5);
  E			= f*(f2*c2+c3*c5);
  
  c[0]			= -F-A+D;
  c[1]			= 2.0*(E-C);
  c[2]			= 2.0*(A-2.0*B-F);
  c[3]			= 2.0*(C+E);
  c[4]			= -(F+A+D);
  n_roots		= R_SolveFourthOrder(c, z);

  // Transcribe real roots
  
  n_real		= 0;
  for (i=0; i<n_roots; ++i)
  {
    f			= creal(z[i]);
    f2			= (((c[4]*f+c[3])*f+c[2])*f+c[1])*f+c[0];
    if (fabs(cimag(z[i]))<1e-9)
      phi[n_real++]	= 2.0*atan(creal(z[i]));
  }
  if (n_real&1)					// only even number of roots
    return -1;

  // Sort roots

  for (i=0; i<n_real-1; ++i)
    for (j=i+1; j<n_real; ++j)
      if (phi[j]<phi[i])
      {
	f		= phi[i];
	phi[i]		= phi[j];
	phi[j]		= f;
      }
	
  // Create feasibility intervals
    
  if (c[4]<0.0)
  {
    phi_low[0]		= phi[0];
    phi_high[0]		= phi[1];
    phi_low[1]		= phi[2];
    phi_high[1]		= phi[3];
    return (long) (n_real/2);
  }
  else
  {
    phi_low[0]		= -M_PI;
    phi_high[0]		= phi[0];
    phi_low[1]		= phi[1];
    phi_high[1]		= phi[2];
    phi_low[2]		= phi[3];
    phi_high[(long) (n_real/2)]		
    			= M_PI;
    return (long) (n_real/2)+1;
  }
}    


vector			rM,
 			rP, rQ,
  			r2, r3, r4;
double			a2, a4, phiL, phiR,
  			l3, cosa3, l4;
long			branch, R_FFAIL;
 
double Fx(double phi)
{ // revamped (vector manipulations written out)
  static double		f1, f2, fx, fy, fz, d;

  // r3
  
  r3.x			= (fx=R*cos(phi))*w2.x+(fy=R*sin(phi))*w3.x+rN.x;
  r3.y			= fx*w2.y+fy*w3.y+rN.y;
  r3.z			= fx*w2.z+fy*w3.z+rN.z;
  
  // r2
  
  f1			= (fx=rP.x-r3.x)*fx+(fy=rP.y-r3.y)*fy+(fz=rP.z-r3.z)*fz
  			    -l3*l3+a2*a2;
  f2			= 2.0*a2*(fx*u2.x+fy*u2.y+fz*u2.z);
  fy			= 4.0*a2*(fx*u3.x+fy*u3.y+fz*u3.z);
  fx			= f1-f2;
  fz			= f1+f2;
  d			= fy*fy-4.0*fx*fz;
  if ((d<0.0)||((fx==0.0)&&(fy==0.0)))
  {
    R_FFAIL		= TRUE;
    return 1e5;
  }
  fx			= a2*cos(phiL=2.0*atan(fx ? 0.5*(-fy+(branch&2 ?
  			  sqrt(d) : -sqrt(d)))/fx : -fz/fy));
  fy			= a2*sin(phiL);
  r2.x			= fx*u2.x+fy*u3.x+rP.x;
  r2.y			= fx*u2.y+fy*u3.y+rP.y;
  r2.z			= fx*u2.z+fy*u3.z+rP.z;

  // r4
  
  f1			= (fx=rQ.x-r3.x)*fx+(fy=rQ.y-r3.y)*fy+(fz=rQ.z-r3.z)*fz
  			    -l4*l4+a4*a4;
  f2			= 2.0*a4*(fx*v2.x+fy*v2.y+fz*v2.z);
  fy			= 4.0*a4*(fx*v3.x+fy*v3.y+fz*v3.z);
  fx			= f1-f2;
  fz			= f1+f2;
  d			= fy*fy-4.0*fx*fz;
  if ((d<0.0)||((fx==0.0)&&(fy==0.0)))
  {
    R_FFAIL		= TRUE;
    return 1e5;
  }
  fx			= a4*cos(phiR=2.0*atan(fx ? 0.5*(-fy+(branch&1 ?
  			  sqrt(d) : -sqrt(d)))/fx : -fz/fy));
  fy			= a4*sin(phiR);
  r4.x			= fx*v2.x+fy*v3.x+rQ.x;
  r4.y			= fx*v2.y+fy*v3.y+rQ.y;
  r4.z			= fx*v2.z+fy*v3.z+rQ.z;
  
  // Epilogue

  return ((fx=r4.x-r2.x)*fx+(fy=r4.y-r2.y)*fy+(fz=r4.z-r2.z)*fz
    -l3*l3-l4*(l4+2.0*l3*cosa3))/l3/l4;
}


double			phi_low[6], phi_high[6];

long Feasible(vector *p, sphere *s)
{
  long			i, j, k, nL, nR, ni;
  double		phiL_low[3], phiL_high[3],
  			phiR_low[3], phiR_high[3],
			f51, f31, f53;
  vector		rL, rR;
  register double	f;
  
  // Initialize frames and variables
  
  rM.x			= (p[0].x+p[6].x)/2.0;
  rM.y			= (p[0].y+p[6].y)/2.0;
  rM.z			= (p[0].z+p[6].z)/2.0;
  rL.x			= (p[0].x+p[1].x+p[2].x)/3.0;
  rL.y			= (p[0].y+p[1].y+p[2].y)/3.0;
  rL.z			= (p[0].z+p[1].z+p[2].z)/3.0;
  rR.x			= (p[4].x+p[5].x+p[6].x)/3.0;
  rR.y			= (p[4].y+p[5].y+p[6].y)/3.0;
  rR.z			= (p[4].z+p[5].z+p[6].z)/3.0;
 
  if (Frame(p, p+1, &rL, &u1, &u2, &u3)) return 0;
  if (Frame(p+6, p+5, &rR, &v1, &v2, &v3)) return 0;
  if (Frame(p+1, p+5, &rM, &w1, &w2, &w3)) return 0;
  w2.x			= -w2.x;
  w2.y			= -w2.y;
  w2.z			= -w2.z;
  w3.x			= -w3.x;
  w3.y			= -w3.y;
  w3.z			= -w3.z;
  
  rP.x			= (f=s[2].d*cos(s[2].alpha))*u1.x+p[1].x;
  rP.y			= f*u1.y+p[1].y;
  rP.z			= f*u1.z+p[1].z;
  rQ.x			= (f=s[5].d*cos(s[6].alpha))*v1.x+p[5].x;
  rQ.y			= f*v1.y+p[5].y;
  rQ.z			= f*v1.z+p[5].z;
  f31			= s[2].d*s[2].d+s[3].d*
  			    (s[3].d+2.0*s[2].d*cos(s[3].alpha));
  f53			= s[4].d*s[4].d+s[5].d*
  			    (s[5].d+2.0*s[4].d*cos(s[5].alpha));
  if (!(f51=(f=p[1].x-p[5].x)*f+(f=p[1].y-p[5].y)*f+(f=p[1].z-p[5].z)*f))
    return 0;
  f			= (f53-f31)/f51;
  rN.x			= 0.5*((1+f)*p[1].x+(1-f)*p[5].x);
  rN.y			= 0.5*((1+f)*p[1].y+(1-f)*p[5].y);
  rN.z			= 0.5*((1+f)*p[1].z+(1-f)*p[5].z);
  f			= 0.5*(f51+f31-f53);
  if ((f = floor((f31-f*f/f51)*1e15+0.5)/1e15)<=0.0)
    return 0;
  R			= sqrt(f);
 
  a2			= s[2].d*sin(s[2].alpha);
  a4			= s[5].d*sin(s[6].alpha);
  l3			= s[3].d;
  l4			= s[4].d;
  cosa3			= cos(s[4].alpha);
  
  // Feasibility intervals

  if ((nL = Feasibility(&rP, l3, a2, &u1, phiL_low, phiL_high))<0)
    return 0;
  if ((nR = Feasibility(&rQ, l4, a4, &v1, phiR_low, phiR_high))<0)
    return 0;
  
  // Construct phi interval from separate phiL and phiR intervals

  ni			= 0;
  for (i=0; i<nL; ++i)
    for (j=0; j<nR; ++j)
      if ((phiL_low[i]>=phiR_low[j])&&(phiL_low[i]<=phiR_high[j]))
      {
        phi_low[ni]	= phiL_low[i];
	if (phiL_high[i]<phiR_high[j])
	  phi_high[ni++]= phiL_high[i];
	else
	  phi_high[ni++]= phiR_high[j];
      }
  for (i=0; i<nR; ++i)
    for (j=0; j<nL; ++j)
      if ((phiR_low[i]>=phiL_low[j])&&(phiR_low[i]<=phiL_high[j]))
      {
        for (k=0; (k<ni)&&(phiR_low[i]!=phi_low[k]); ++k);
	if (k>=ni)
	{
	  phi_low[ni]	= phiR_low[i];
	  if (phiR_high[i]<phiL_high[j])
	    phi_high[ni++]= phiR_high[i];
	  else
	    phi_high[ni++]= phiL_high[j];
        }
      }

  return ni;
}


long Rebridge(vector *p, sphere *s)
{
  long			i, j, n, ni, nr, nr_old;
  double		phi[B_MAXNROOTS];
  vector		t1, r[B_MAXNROOTS][3];
  
  if (!(ni = Feasible(p, s)))
    return 0;

  // Regression:
  //   Search all feasible intervals over all branches
  
  nr			= 0;
  n			= 0;
  R_FFAIL		= FALSE;
  for (branch=0; branch<4; ++branch)
    for (i=0; (i<ni)&&(nr<B_MAXNROOTS); ++i)
    {
      nr_old		= nr;
      nr		+= R_SolveNumerical(Fx, phi_low[i], phi_high[i], 
      				20, 4, phi+nr, B_MAXNROOTS-nr);
      if (R_FFAIL)
        return 0;

      if (nr>=B_MAXNROOTS)
        return 0;
     
      // Chain segment reconstruction
      
      for (j=nr_old; (j<nr); ++j)
      {
	//fprintf(stderr, "F%d[%g] = %g\n", j, phi[j], Fx(phi[j]));
	Fx(phi[j]);
	t1		= V_Subtr(&r3, p+3);
        if (V_Dot(&t1, &t1)>1e-15)
	{
	  r[n][0]	= r2;
	  r[n][1]	= r3;
	  r[n][2]	= r4;
	  ++n;
	}
      }
    }
  
  // Wrap-up
 
  if (nr>1)
    if (nr%2) 					// Geometric fail
      return 0;
 
  i			= (long) ((double) n*ran1(seed));
  p[2]			= r[i][0];
  p[3]			= r[i][1];
  p[4]			= r[i][2];
  return n;
}


double Jacobian(vector *p, sphere *s)	// J(III->I), see ref.1 Appendix
{
  long			i;
  static double		l[7],
  			b1xb2b3, b1xb2b4, b1xb2b5,
			b2xb3b4, b2xb5b6,
			b3xb4b5, b3xb5b6,
			b4xb5b6,
  			J;
  static vector		b[7], v;

  for (i=1; i<7; ++i)
  {
    l[i]		= s[i].d;
    b[i]		= V_Subtr(p+i, p+i-1); 
    b[i]		= V_Mult(1.0/l[i], b+i);
  }
  v			= V_Cross(b+1, b+2);
  b1xb2b3		= V_Dot(&v, b+3);
  b1xb2b4		= V_Dot(&v, b+4);
  b1xb2b5		= V_Dot(&v, b+5);
  v			= V_Cross(b+2, b+3);
  b2xb3b4		= V_Dot(&v, b+4);
  v			= V_Cross(b+2, b+5);
  b2xb5b6		= V_Dot(&v, b+6);
  v			= V_Cross(b+3, b+4);
  b3xb4b5		= V_Dot(&v, b+5);
  v			= V_Cross(b+3, b+5);
  b3xb5b6		= V_Dot(&v, b+6);
  v			= V_Cross(b+4, b+5);
  b4xb5b6		= V_Dot(&v, b+6);

  J			= 1.0/(l[1]*l[6]*fabs(
  			    -l[3]*l[5]*b1xb2b3*b3xb4b5*(
			      l[2]*b2xb5b6+l[3]*b3xb5b6+l[4]*b4xb5b6)+
			    l[2]*l[4]*b2xb3b4*b4xb5b6*(
			      l[3]*b1xb2b3+l[4]*b1xb2b4+l[5]*b1xb2b5)));
  for (i=2; i<7; ++i)
    J			*= sqrt(l[i-1]*l[i-1]+l[i]*
			     (l[i]+2.0*l[i-1]*cos(s[i].alpha)));
  
  return J;
}

double new_Jacobian(vector *p, sphere *s)	// J(IV->I), see ref.2 Appendix
{
  long			i;
  static double		l[7],
  			b1xb2b3, b1xb2b4, b1xb2b5,
			b2xb3b4, b2xb5b6,
			b3xb4b5, b3xb5b6,
			b4xb5b6,
  			J;
  static vector		b[7], v;

  for (i=1; i<7; ++i)
  {
    l[i]		= s[i].d;
    b[i]		= V_Subtr(p+i, p+i-1); 
    b[i]		= V_Mult(1.0/l[i], b+i);
  }
  v			= V_Cross(b+1, b+2);
  b1xb2b3		= V_Dot(&v, b+3);
  b1xb2b4		= V_Dot(&v, b+4);
  b1xb2b5		= V_Dot(&v, b+5);
  v			= V_Cross(b+2, b+3);
  b2xb3b4		= V_Dot(&v, b+4);
  v			= V_Cross(b+2, b+5);
  b2xb5b6		= V_Dot(&v, b+6);
  v			= V_Cross(b+3, b+4);
  b3xb4b5		= V_Dot(&v, b+5);
  v			= V_Cross(b+3, b+5);
  b3xb5b6		= V_Dot(&v, b+6);
  v			= V_Cross(b+4, b+5);
  b4xb5b6		= V_Dot(&v, b+6);

  J			= 1.0/(l[1]*l[6]*fabs(
  			    -l[3]*l[5]*b1xb2b3*b3xb4b5*(
			      l[2]*b2xb5b6+l[3]*b3xb5b6+l[4]*b4xb5b6)+
			    l[2]*l[4]*b2xb3b4*b4xb5b6*(
			      l[3]*b1xb2b3+l[4]*b1xb2b4+l[5]*b1xb2b5)));
  for (i=1; i<=5; i++)
    J	*=	l[i];
  for (i=2; i<=6; i++)
    J	*=	l[i];
  for (i=2; i<=6; i++)
    J	*=	sin(s[i].alpha);
  
  return J;
}
// Setup input parameters for 7 vectors p[] and spherical coordinates s[]
// used by Feasible(), Rebridge(), and Jacobian() through the molecule
// parent list.  Used in conjunction with NeighborList(), NextEndBridge(),
// and NextRebridge().  
//
// Input:
//   molm:	Target molecule.
//   end:	Where to start in the parent list.
//   reverse:	Flag for transcription direction; 0: forwards, 1: reverse.
//
// Output:
//   p:		Transcribed positions.
//   s:		Calculated bond lengths and angles (no torsions!).
/*
long RebridgeSetup(
	molstruct *molm, long end, long reverse, vector *p, sphere *s)
{
   long			i, ip, flag = 0;
   sphere		*s1;
   vector		dr, drp;

   if (end >= molm->nsites) return 1;		// Exit on failure

   if (reverse)	{	p+=0;	s+=1;}		// end bead as p[0]
   else		{ 	p+=6;	s+=6;}		// end bead as p[6]

   ip	=	end;
   *p	=	molm->p[ip];

   for (i=0; (i<6)&&((ip = molm->parent[ip])>=0); ++i) {
      if (reverse)	p++;
      else		p--;

      dr	=	V_Subtr(p, molm->p+ip);		// bond vector
      *p	=	molm->p[ip];			// update p

      s->d	=	sqrt(V_Dot(&dr, &dr));		// bond length
      dr	=	V_Mult(1.0/s->d, &dr);		// unit bond vector

      if (i) {
         s1->alpha	=	acos(V_Dot(&dr, &drp));	// bond angle
      }
      drp	=	dr;			// previous dr

      if (reverse) {	s++;	s1=s;}		// or s1 = reverse ? ++s : s--;
      else	   {	s1=s;	s--;}
   }
   return (ip<0);	// return false if NO error occur
}
*/

long RebridgeSetup(
	molstruct *molm, long end, long reverse, vector *p, sphere *s)
{
  long			i, ip = end, flag = 0;
  sphere		*sl;
  vector		dr, drp, *pis, *pip;

  if (end>=molm->nsites) return 1;		// Exit on failure
  pis			= molm->p+ip;
  p			+= reverse ? 0 : 6;
  s			+= reverse ? 1 : 6;
  *p			= *pis;			// Transcribe positions
  for (i=0; (i<6)&&((ip = molm->parent[ip])>=0); ++i)
  {
    if (reverse) 
      ++p;
    else 
      --p;
    pip			= molm->p+ip;		// Parent position pointer
    *p			= *pip;			
    drp			= V_Subtr(pis, pip);	// Determine direction
    drp			= V_Mult(1.0/(s->d = sqrt(V_Dot(&drp, &drp))), &drp);
    if (i)
    {
      sl->alpha		= acos(V_Dot(&drp, &dr));
    }
    pis			= pip;
    dr			= drp;
    sl	 		= reverse ? ++s : s--;
  }
  return (ip<0);
}
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            src/rigidforce.c                                                                                    0000600 0143352 0000144 00000060542 10757337336 013362  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    forcefield.c
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Calculation of the potential of a system

    Modified:	Sept 20, 2007
		Long range correction no longer calculated for each particle,
		but calculated for the whole system at one time.
		
*/

#define __FORCEFIELD_MODULE
#include "forcefield.h"

/* Some basic operations to vstruct and wstruct */

void vstructNull(vstruct *V)
{
   V->lj	=	0.0;
   V->ljcorr	=	0.0;
   V->hs	=	0.0;
   V->stretch	=	0.0;
   V->bending	=	0.0;
   V->torsion	=	0.0;
   V->corr	=	0.0;
   V->bonded	=	0.0;
   V->nonbonded	=	0.0;
   V->tot	=	0.0;
}

void wstructNull(wstruct *VIR)
{
   VIR->lj	=	0.0;
   VIR->stretch	=	0.0;
   VIR->torsion	=	0.0;
   VIR->tot	=	0.0;
}

void vstructNegate(vstruct *V)
{
   V->lj	*=	-1.0;
   V->ljcorr	*=	-1.0;
   V->hs	*=	-1.0;
   V->stretch	*=	-1.0;
   V->bending	*=	-1.0;
   V->torsion	*=	-1.0;
   V->corr	*=	-1.0;
   V->bonded	*=	-1.0;
   V->nonbonded	*=	-1.0;
   V->tot	*=	-1.0;
}

void wstructNegate(wstruct *VIR)
{
   VIR->lj	*=	-1.0;
   VIR->stretch *=	-1.0;
   VIR->torsion	*=	-1.0;
   VIR->tot	*=	-1.0;
}

vstruct vstructSum(vstruct *V1, vstruct *V2)
{
   vstruct	V;
  
   V.lj		=	V1->lj + V2->lj;
   V.ljcorr	=	V1->ljcorr + V2->ljcorr;
   V.hs		=	V1->hs + V2->hs;
   V.stretch	=	V1->stretch + V2->stretch;
   V.bending	=	V1->bending + V2->bending;
   V.torsion	=	V1->torsion + V2->torsion;

   V.corr	=	V1->corr + V2->corr;
   V.bonded	=	V1->bonded + V2->bonded;
   V.nonbonded	=	V1->nonbonded + V2->nonbonded;
   V.tot	=	V1->tot + V2->tot;
   return	V;
}

wstruct wstructSum(wstruct *W1, wstruct *W2)
{
   wstruct	W;

   W.lj		=	W1->lj + W2->lj;
   W.stretch	=	W1->stretch + W2->stretch;	
   W.torsion	=	W1->torsion + W2->torsion;

   W.tot	=	W1->tot + W2->tot;
   return	W;
}

void Printvstruct(vstruct *V)
{
   printf("v.lj\t=\t%f\n", V->lj);
   printf("v.ljcorr\t=\t%f\n", V->ljcorr);
   printf("v.hs\t=\t%f\n", V->hs);
   printf("v.stretch\t=\t%f\n", V->stretch);
   printf("v.bending\t=\t%f\n", V->bending);
   printf("v.torsion\t=\t%f\n", V->torsion);
   printf("v.corr\t=\t%f\n", V->corr);
   printf("v.bonded\t=\t%f\n", V->bonded);
   printf("v.nonbonded\t=\t%f\n", V->nonbonded);
   printf("v.tot\t=\t%f\n", V->tot);
}


/* Lennard Jones combining rules */

/* SIGMA_eff 	= 0.5 * (SIGMA_1 + SIGMA_2) */
/* EPSILON_eff 	= sqrt( EPSILON_1 * EPSILON_2 ) */

void CalcMixSigma(long typei, long typej)
{
   if (typei == typej)
      type[typei].mix[typei].SIGMA	=	type[typei].SIGMA;
   else {
      type[typei].mix[typej].SIGMA	
		=	0.5 * (type[typei].SIGMA + type[typej].SIGMA);
      type[typej].mix[typei].SIGMA
		=	type[typei].mix[typej].SIGMA;
   }
}		

void CalcMixEpsilon(long typei, long typej)
{
   if (typei == typej)
      type[typei].mix[typej].EPSILON	=	type[typei].EPSILON;
   else {
      type[typei].mix[typej].EPSILON
		=	sqrt( type[typei].EPSILON * type[typej].EPSILON );
      type[typej].mix[typei].EPSILON
		=	type[typei].mix[typej].EPSILON;
   }
}

void InitForcefield()
{
   long		i, j;
   
   for (i=0; i<NTYPES; i++)
      for (j=i; j<NTYPES; j++) {	// symmetric combining rule
         CalcMixSigma(i, j);
         CalcMixEpsilon(i, j);
      }
}


/* Calculate different types of interactions */
/* V*Site() and V*Mol() do NOT initialize the virial */
/* in their parameter list, so we have to make sure  */
/* they do the right thing to the right virial.      */

double VHSSite(molstruct *molm, long site)
{
   static long		ibox, n;
   static molstruct	*moln;
   static double	r2;
   static vector	pm;
   static typestruct	*typema;
#ifdef CELL_LIST
   static long		i, j;
   static cellstruct	*cellm, *celli;
#endif

   if (!V_HS || 0==(molm->flags[site]) ) 
      return	0.0;
   
   ibox		=	molm->box;			// determine which box
   pm		=	molm->p[site];
   typema	=	type + molm->type[site];

#ifdef CELL_LIST
   cellm	=	molm->cell[site];
   for (j=0; j<cellm->nneigh; j++) {
      
      celli	=	cellm->neigh[j];

      for (i=0; i<celli->nsites; i++) {
         if (moln=celli->mol[i]) {

            n	=	celli->molsite[i];
#else
   for (moln=mol; moln < mol + NMOLS; moln++) {				// search through all mols

      if (moln->box	==	ibox) {					// check if two mols are in the same box

         for (n=0; n<moln->nsites; n++) {				// search through all sites
#endif /* CELL_LIST */
            
            if ( (moln->flags[n]>0) && (molm!=moln || n!=site) ) {	// interaction with itself is excluded

               r2	=	DistSQ(pm, moln->p[n], ibox);
               
               if (r2 < (typema->mix+moln->type[n])->SIGMA * (typema->mix+moln->type[n])->SIGMA)
		  return	1.0e4;
               else
		  return	0.0;
	    }
         }   
      }
   }
}


double VHSMol(molstruct *moli)
{
   long		i;
   double	vhs=0.0;

   if (moli->box <0 || !V_HS )
      return	0.0;

   for (i=0; i<moli->nsites; i++)
      if ( vhs=VHSSite(moli, i) > 1.0)
	 return	vhs;	
   
   return	vhs;
}


void CalcVHS()				// calculate total hard sphere interaction
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++)
      v[ibox].hs	=	0.0;

   for (moli=mol; moli<mol+NMOLS; moli++) {
      if ( (ibox=moli->box) >=0 && v[ibox].hs < 1.0) 
         v[ibox].hs	+=	VHSMol(moli);
   }
}


double VStretchSite(molstruct *molm, long site, double *w)	// V = 0.5 * k * (l-l0)^2
{
   long			i, k=0, system=molm->box;
   double		l0, cosa, v=0.0, f;
   vector		dr[2], *p0, *p1;
   typestruct		*t;

   if ( (!V_STRETCH) || (molm->flags[site]==0) || (site<0) || (site>=molm->nsites) )
      return	0.0;
   p0		=	molm->p+site;
   i		=	site;
   while ( (k<1) && ((i=molm->parent[i])>=0) && (molm->flags[i]>0) ) {
      p1	=	molm->p+i;
      dr[k]	=	V_Subtr(p0, p1);
      p0	=	p1;
      k	++;
   }
   if ( !k )
      return	0.0;
   t		=	type + molm->type[site];
   l0		=	sqrt( V_Dot(dr, dr) );
   f		=	l0 - t->LSTRETCH;
   v		=	0.5 * f * f * t->KSTRETCH;
   if (V_VIRIAL)
      *w	+=	l0 * f * t->KSTRETCH;
   return	v;
}


double VStretchMol(molstruct *molm, double *w)
{
   long		i;
   double	vstretch = 0.0;

   if (!V_STRETCH || molm->box<0) 
      return	0.0;
   for (i=0; i<molm->nsites; i++) 
      vstretch	+=	VStretchSite(molm, i, w);
   return 	vstretch;
}


void CalcVStretch()
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++) {
      v[ibox].stretch		=	0.0;
      vir[ibox].stretch		=	0.0;
   }
   if (!V_STRETCH) 	
      return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox=moli->box) >=0 ) 
         v[ibox].stretch	+=	VStretchMol(moli, &(vir[ibox].stretch));
}


double VLJSite(molstruct *molm, long site, double *w)	// calculate the LJ potential energy of one site
{
   static long		ibox, n, flags_m;
   static vector	dp;
   static double	vlj, wlj;			// lj potential, lj virial
   static molstruct	*moln;
   static double	Sigma, Epsilon, r2, r6i, r12i, rc2;
   static vector	pm;
   static double	ljcut;
   typestruct		*typem = type + molm->type[site];
   mixstruct		*typemix;
#ifdef CELL_LIST
   static long		i, j;
   static cellstruct	*cellm, *celli;
#endif

   if (!V_LJ || (flags_m=molm->flags[site]) == 0 || (typem->EPSILON==0.0) )
      return	0.0;

   vlj	=	0.0;					// MUST, because vlj is static variable
   wlj	=	0.0;
   
   ibox		=	molm->box;			// determine which box
   rc2		=	BOX[ibox].rc * BOX[ibox].rc;
   pm		=	molm->p[site];

#ifdef CELL_LIST
   cellm	=	molm->cell[site];

   for (j=0; j<cellm->nneigh; j++) {
      
      celli	=	cellm->neigh[j];

      for (i=0; i<celli->nsites; i++) {
         if (moln=celli->mol[i]) {

            n	=	celli->molsite[i];
#else
   for (moln=mol; moln < mol + NMOLS; moln++) {		// search through all mols

      if (moln->box	==	ibox) {			// check if two mols are in the same box

         for (n=0; n<moln->nsites; n++) {		// search through all sites
#endif /* CELL_LIST */

            if ( (moln->flags[n]>0) && (moln==molm ? abs(n-site)>=DLJ : 1) ) {	// LJ interaction condition
									// interaction with itself is excluded
									// automatically
               r2	=	DistSQ(pm, moln->p[n], ibox);	

               if (r2 < rc2) {
		  typemix	=	typem->mix + moln->type[n];
                  Sigma		=	typemix->SIGMA;
	          Epsilon	=	typemix->EPSILON;

         	  r6i 		= 	Sigma * Sigma/r2;
	          r6i 		= 	r6i * r6i * r6i;
		  r12i		=	r6i * r6i;
/*
		  if ( moln==molm && DLJ==abs(n-site) ) {	// LJ 1-4 pair scaling
         	     vlj 	+= 	2.0 * Epsilon * (r12i - r6i);
                     if (V_VIRIAL)
                        wlj	+=	-24.0 * Epsilon * (r12i - 0.5 * r6i);
  		  }
		  else { 
*/         	     vlj 	+= 	4.0 * Epsilon * (r12i - r6i);
                     if (V_VIRIAL)
                        wlj	+=	-48.0 * Epsilon * (r12i - 0.5 * r6i);
  //    		  }

		  if (V_LJSHIFT) {                
      		     ljcut	=	Sigma / BOX[ibox].rc;		// sigma/rc
		     ljcut	=	ljcut * ljcut * ljcut;		// (sigma/rc)^3
		     ljcut	=	ljcut * ljcut;			// (sigma/rc)^6
		     ljcut	=	4.0 * Epsilon * (ljcut * ljcut - ljcut);
//		     if ( moln==molm && DLJ==abs(n-site) ) 	// LJ 1-4 pair scaling
//		        vlj	-=	0.5 * ljcut;
//		     else
			vlj	-=	ljcut;
                  }
               }
	    }
         }   
      }
   }
   if (V_VIRIAL)
      *w	+=	wlj;

   return vlj;
}


double VLJMol(molstruct *molm, double *w)
{
   long		i;
   double	vlj = 0.0;

   if (!V_LJ || molm->box<0 )
      return	0.0;

   for (i=0; i<molm->nsites; i++)
      vlj	+=	VLJSite(molm, i, w);

   return	vlj; 
}


void CalcVLJ()
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++) {
      v[ibox].lj	=	0.0;
      vir[ibox].lj	=	0.0;
   }

   if (!V_LJ)
      return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox = moli->box) >=0)
         v[ibox].lj	+=	0.5 * VLJMol(moli, &(vir[ibox].lj));	// fix double-counting

   if (V_VIRIAL) {
      for (ibox=0; ibox<NBOX; ibox++)
         vir[ibox].lj	*=	0.5;			// fix double-counting of virial
   }
   return;
}


double VBendingSite(molstruct *molm, long site)		// Only calculate the bending
{							// energy on its parent side
   long		i, k=0;
   double	l0, l1, cosa, v=0.0, f;
   vector	dr[2], *p0, *p1;
   typestruct	*t;

   if ( (!V_BENDING) || (0==molm->flags[site]) || (site<0) || (site>=molm->nsites) )
      return	0.0;
   p0		=	molm->p+site;
   i		=	site;
   while ( (k<2) && ((i=molm->parent[i])>=0) && (molm->flags[i]>0)) {
      p1	=	molm->p+i;
      dr[k++]	=	V_Subtr(p0, p1);
      p0	=	p1;
   }
   if (k<2)
      return	0.0;
   t		=	type + molm->type[site];
   l0		=	sqrt(V_Dot(dr, dr));
   l1		=	sqrt(V_Dot(dr+1, dr+1));
   cosa		=	V_Dot(dr, dr+1)/(l0*l1);
//   f		=	acos(cosa) - t->THETA;		// OPLS model
   f		=	cosa - cos(t->THETA);		// Rigid model
   v		=	0.5 * t->KBENDING * f * f;
   return	v;
}


double VBendingMol(molstruct *molm)
{
   long		i;
   double	v_bending = 0.0;

   if (!V_BENDING || (molm->box<0))
      return	0.0;
   for (i=0; i<molm->nsites; i++)
      v_bending	+=	VBendingSite(molm, i);
   return	v_bending;
}


void CalcVBending()
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++) 
      v[ibox].bending	=	0.0;

   if (!V_BENDING)	return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox = moli->box) >=0 )
         v[ibox].bending	+=	VBendingMol(moli);	// no double-counting
}


double VTorsionSite(molstruct *molm, long site)		// no virial contribution
{							// only toward parents side !!
   long		i, k=0;
   double	cosb, b, l0, l1;
   vector	dr[3], n0, n1, *p0, *p1;
   typestruct	*t;

   if ( (!V_TORSION) || (0==molm->flags[site]) || (site<0) || (site>molm->nsites) ) 
      return	0.0;
   p0	=	molm->p+site;
   i	=	site;
   while ( (k<3) && ((i=molm->parent[i])>=0) && (molm->flags[i]>0)) {
      p1	=	molm->p+i;
      dr[k++]	=	V_Subtr(p0, p1);
      p0	=	p1;
   }
   if (k<3) {
      return	0.0;
   }
   n0	=	V_Cross(dr+1, dr);
   n1	=	V_Cross(dr+2, dr+1);
   l0	=	sqrt(V_Dot(&n0, &n0));
   l1	=	sqrt(V_Dot(&n1, &n1));
   cosb	=	-V_Dot(&n0, &n1)/(l0*l1);
   t	=	type + molm->type[site];
/*   
   if (fabs(cosb-1)<ZERO)
      b	=	0;
   else if (fabs(cosb+1) <ZERO)
      b	=	M_PI;
   else
      b	=	acos(cosb);

   return	0.5 * (t->TORSION[1]*(1-cosb) + t->TORSION[2]*(1-cos(2*b)) + t->TORSION[3]*(1-cos(3*b)) );
*/
   cosb	*=	-1;	// because now b=pi correspond to trans state

   return	t->TORSION[0] + cosb * t->TORSION[1] + cosb * cosb * t->TORSION[2] 
		+ cosb * cosb * cosb * t->TORSION[3] + cosb * cosb * cosb * cosb * t->TORSION[4] 
		+ cosb * cosb * cosb * cosb * cosb * t->TORSION[5]; 	// Rigid model
}


double VTorsionMol(molstruct *molm)
{
   long		i;
   double	v_tor = 0.0;

   if (!V_TORSION || (molm->box<0))
      return	0.0;
   for (i=0; i<molm->nsites; i++)
      v_tor	+=	VTorsionSite(molm, i);
   return	v_tor;
}


void CalcVTorsion()
{
   long		ibox;
   molstruct	*moli;

   for (ibox=0; ibox<NBOX; ibox++) 
      v[ibox].torsion	=	0.0;

   if (!V_TORSION)	return;

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (ibox = moli->box) >=0 )
         v[ibox].torsion	+=	VTorsionMol(moli);	// no double-counting
}


void CalcVLJCorr()		// recalculate LJ long-range energy correction
{
   long		n, m, flags_m, *nsites, system;
   double	Epsilon, Sigma, rc3, rc9, vol;
   molstruct	*molm, *moln;
   typestruct	*typem;
   mixstruct	*typemix;
   long		NActive[MAXNSYSTEMS];

   if (!V_LJLRC)	return;
   if (!(nsites = (long *) calloc(NSYSTEMS, sizeof(long))))
      Exit("forcefield", "CalcVLJCorr", "Out of memory");

   for (system=0; system<NSYSTEMS; system++) {
      v[system].ljcorr	=	0.0;
      nsites[system]	=	0;
      NActive[system]	=	0;
   }

   for (molm=mol; molm<mol+NMOLS; molm++)
      if ( (system=molm->box) >= 0)
         for (m=0; m<molm->nsites; m++)
            if (molm->flags[m]>0)
                NActive[system]	++;

   for (molm=mol; molm<mol+NMOLS; molm++)
      if ( (system=molm->box) >= 0)
         for (m=0; m<molm->nsites; m++) 
            if ( (molm->flags[m]>0) && ((typem=type+molm->type[m])->EPSILON) )
               for (moln=mol; moln<mol+NMOLS; moln++)
                  if (moln->box == system)
	             for (n=0; n<moln->nsites; n++)
			if ( (moln->flags[n]>0) && (moln==molm ? abs(n-m)>=DLJ : 1) 
				&& (Epsilon=(typemix=(typem->mix)+moln->type[n])->EPSILON) ) {

			   Sigma	=	typemix->SIGMA;
	                   rc3		=	BOX[system].rc / Sigma;
		           rc3		=	rc3 * rc3 * rc3;
		           rc9		=	rc3 * rc3 * rc3;
		           vol		=	BOX[system].vol;
			   if (V_LJSHIFT)
                              v[system].ljcorr	+=	Epsilon/vol * Sigma * Sigma * Sigma * (4.0/(3*rc9)-2.0/rc3);
			   else
                              v[system].ljcorr	+=	Epsilon/vol * Sigma * Sigma * Sigma * (1.0/(3*rc9)-1.0/rc3);
			   nsites[system]	++;
                        }
         
   for (system=0; system<NSYSTEMS; system++)
      v[system].ljcorr	*=	nsites[system] ? 8.0*pi/3.0*NSites[system]*NActive[system]/nsites[system] : 0.0;

   free(nsites);
}   


void CalcVCorr()
{
   long		ib;

   if (V_LJLRC)
      CalcVLJCorr();

   for (ib=0; ib<NBOX; ib++) 
      v[ib].corr	=	v[ib].ljcorr;
}


void CalcV()
{
   long		ib;

   for (ib=0; ib<NBOX; ib++) {
      vstructNull(v+ib);
      wstructNull(vir+ib);
   }

   if (V_LJ)	
      CalcVLJ();
   if (V_HS)
      CalcVHS();
   if (V_STRETCH)
      CalcVStretch();
   if (V_BENDING)
      CalcVBending();
   if (V_TORSION)
      CalcVTorsion();

   CalcVCorr();

   for (ib=0; ib<NBOX; ib++) {
      v[ib].bonded	=	v[ib].stretch + v[ib].bending + v[ib].torsion;
      v[ib].nonbonded	=	v[ib].lj + v[ib].hs + v[ib].corr;
      v[ib].tot		=	v[ib].bonded + v[ib].nonbonded;

      vir[ib].tot	=	vir[ib].lj + vir[ib].stretch + vir[ib].torsion;
   }
   return;
}


vstruct CalcVSite(molstruct *moli, long site, wstruct *VIRSite)	// energy and virial
{								// contribution of one
   long		i;						// single site to the system
   vstruct	V;

   vstructNull(&V);
   wstructNull(VIRSite);

   if (V_LJ) 
      V.lj		+=	VLJSite(moli, site, &(VIRSite->lj));
   if (V_HS)
      V.hs		+=	VHSSite(moli, site);
   if (V_STRETCH)
      for (i=0; i<2; i++) 
         V.stretch	+=	VStretchSite(moli, site+i, &(VIRSite->stretch));
   if (V_BENDING) 
      for (i=0; i<3; i++)
         V.bending	+=	VBendingSite(moli, site+i);
   if (V_TORSION)
      for (i=0; i<4; i++)
         V.torsion	+=	VTorsionSite(moli, site+i);

   // tail correction is calculated for the whole system, not for individual site or molecule
   V.bonded	=	V.stretch + V.torsion + V.bending;
   V.nonbonded	=	V.lj + V.hs;
   V.tot	=	V.bonded + V.nonbonded;

   if (V_VIRIAL)
      VIRSite->tot	=	VIRSite->lj + VIRSite->stretch + VIRSite->torsion;

   return	V;   
}

/* VDeleteSites, VAddSites, and grow update ALL energy and virial components */
/* bonded, nonbonded, tot, EXCEPT for long-range energy correction */

double VDeleteSites(molstruct *moli, long i_0, long i_n)
{
   long		i, j;
   vstruct	*v_new		=	v+moli->box;
   wstruct	*vir_new	=	vir+mol->box;
   double	v_old		=	v_new->tot;

   if (V_VIRIAL)
      wstructNegate(vir_new);		// need to be paired up in the end

   for (i=i_0; i<=i_n; i++) {
      
      if (V_STRETCH) 
         for (j=1; j<2; j++)
            v_new->stretch	-=	VStretchSite(moli, i+j, &(vir_new->stretch));
      if (V_BENDING)
         for (j=1; j<3; j++)
            v_new->bending	-=	VBendingSite(moli, i+j);
      if (V_TORSION)
         for (j=1; j<4; j++)
            v_new->torsion	-=	VTorsionSite(moli, i+j);

      if (moli->flags[i] > 0) {
         if (V_LJ) 
            v_new->lj		-=	VLJSite(moli, i, &(vir_new->lj));
         if (V_HS)
            v_new->hs		-=	VHSSite(moli, i);
         if (V_STRETCH) 
            v_new->stretch	-=	VStretchSite(moli, i, &(vir_new->stretch));
         if (V_BENDING)
	    v_new->bending	-=	VBendingSite(moli, i);
         if (V_TORSION)
            v_new->torsion	-=	VTorsionSite(moli, i);

#ifdef CELL_LIST
	 CL_Delete(moli, i);			// Unregister site.
#endif
         moli->flags[i]	=	0;		// Deactivate site. Important!! to avoid double-counting
      }
   }
   v_new->bonded	=	v_new->stretch + v_new->bending + v_new->torsion;
   v_new->nonbonded	=	v_new->lj + v_new->hs + v_new->corr;
   v_new->tot		=	v_new->bonded + v_new->nonbonded;
   if (V_VIRIAL) {				// update virial
      wstructNegate(vir_new);			// VIR = -VIR
      vir_new->tot	=	vir_new->lj + vir_new->stretch + vir_new->torsion;
   }
   return	v_new->tot - v_old;
}


double VAddSites(molstruct *moli, long i_0, long i_n)
{
   long		i, j;
   vstruct	*v_new		=	v+moli->box;
   wstruct	*vir_new	=	vir+moli->box;
   double	v_old		=	v_new->tot;

   for (i=i_0; i<=i_n; i++) {

      if (moli->flags[i] == 0) {
         moli->flags[i]		=	1;		// activate this site
#ifdef CELL_LIST
	 CL_Add(moli, i);
#endif
         if (V_LJ)
            v_new->lj		+=	VLJSite(moli, i, &(vir_new->lj));
         if (V_HS)
            v_new->hs		+=	VHSSite(moli, i);
         if (V_STRETCH)
            v_new->stretch	+=	VStretchSite(moli, i, &(vir_new->stretch));
         if (V_BENDING)
	    v_new->bending	+=	VBendingSite(moli, i);
         if (V_TORSION)
            v_new->torsion	+=	VTorsionSite(moli, i);
      }
      if (V_STRETCH)
         for (j=1; j<2; j++)
            v_new->stretch	+=	VStretchSite(moli, i+j, &(vir_new->stretch));
      if (V_BENDING)
         for (j=1; j<3; j++)
            v_new->bending	+=	VBendingSite(moli, i+j);
      if (V_TORSION)
         for (j=1; j<4; j++)
            v_new->torsion	+=	VTorsionSite(moli, i+j);
   }
   v_new->bonded	=	v_new->stretch + v_new->bending + v_new->torsion;
   v_new->nonbonded	=	v_new->lj + v_new->hs + v_new->corr;
   v_new->tot		=	v_new->bonded + v_new->nonbonded;
   if (V_VIRIAL)
      vir_new->tot	=	vir_new->lj + vir_new->stretch + vir_new->torsion;
   return	v_new->tot - v_old;
}


long Select(double *wt, double sumw)
{
   double	ws, cumw;
   long		n;

   ws		=	ran1(seed) * sumw;
   n		=	0;
   cumw		=	wt[0];

   while (cumw < ws) {
      n		++;
      cumw	+=	wt[n];
   }
   if (n>=NTRIALCONF)
      Exit("forcefield", "Select", "n>NTRIALCONF.");

   return	n;
}


double grow(char *s, molstruct *molm, long site)
{
   long		i, j, k, n, ib=molm->box;
   double	W, sumw, lbond;
   vector	dp;
   static long		init = 1;
   static double	*wt;
   static vector	*pt;
   static vstruct	*Vt;
   static wstruct	*VIRt;

   k = 		NTRIALCONF;

   if (init) {
      if (! (wt=(double*) calloc(k, sizeof(double))) )		// allocate for trial variables
         Exit("forcefield", "grow", "out of memory!");
      if (! (pt=(vector *) calloc(k, sizeof(vector))) )		// trial position
         Exit("forcefield", "grow", "out of memory!");
      if (! (Vt=(vstruct *) calloc(k, sizeof(vstruct))) )	// trial energy
         Exit("forcefield", "grow", "out of memory!");
      if (! (VIRt=(wstruct *) calloc(k, sizeof(wstruct))) )	// trial virial
         Exit("forcefield", "grow", "out of memory!");
      init	=	0;
   }

   if ( !strcmp(s, "old") ) {
      vstructNegate(v+ib);				// need to be paired up in the end
      wstructNegate(vir+ib);
   }

   for (i=site; i<molm->nsites; i++) {
      if (1 == molm->flags[i]) {			// if active site
         molm->flags[i]	=	0;			// deactivate sites
#ifdef CELL_LIST
         CL_Delete(molm, i);  				// remove these sites from cells
#endif
      }
   }

   W	=	0.0;		// Rosenbluth factor

   for (i=site; i<molm->nsites; i++) {

      if (0==i) {				// grow a whole chain, first atom
         if ( !strcmp(s, "new")) {		// new conf.
	    molm->p[0].x	=	(ran1(seed)-0.5) * BOX[ib].lx;
	    molm->p[0].y	=	(ran1(seed)-0.5) * BOX[ib].ly;
	    molm->p[0].z	=	(ran1(seed)-0.5) * BOX[ib].lz;
         }

#ifdef CELL_LIST
         CL_Add(molm, i);			// add into cell
#endif
         molm->flags[i]	=	1;		// activated this site 

         Vt[0]		=	CalcVSite(molm, i, VIRt);	// energy and virial contribution
         v[ib]		=	vstructSum(v+ib, Vt);		// of this site to the system

         if (V_VIRIAL)
            vir[ib]	=	wstructSum(vir+ib, VIRt);

         W	=	log(k) - Vt[0].nonbonded/BOX[ib].temp;
      }
      else {						// not the first atom

         for (j=0; j<k; j++) {				// search thru trial orient.

            molm->flags[i]	=	1;	// activate this site, before tors_bonda()

            if (!strcmp(s, "old") && (0==j))
               pt[0]	=	molm->p[i];
            else {
	       // generate trial conf. according to Vbonded
	       dp	=	tors_bonda(molm, i);
               lbond	=	bondl_g(type[0].LSTRETCH, BOX[ib].temp);
               dp	=	V_Mult(lbond, &dp);

	       pt[j]	=	V_Add(molm->p+i-1, &dp);	// trial position
               //MapInBox2(pt+j, PBC, BOX[ib].lbox);

               molm->p[i]	=	pt[j];
            }

#ifdef CELL_LIST
	    CL_Add(molm, i);				// add trial site into cell
#endif

            Vt[j]	=	CalcVSite(molm, i, VIRt+j);
            wt[j]	=	exp(-Vt[j].nonbonded/BOX[ib].temp);

            molm->flags[i]	=	0;		// deactivate this site
#ifdef CELL_LIST
	    CL_Delete(molm, i);				// delete trial site from cell list
#endif
            //printf("%s\ttrial %d\t of %d\t, vt[j].tot=%f\tvt[j].nonbonded=%f\n", s, j, i, Vt[j].tot, Vt[j].nonbonded);
         }

         sumw	=	0.0;				// sum of wt
         for (j=0; j<k; j++) 
	    sumw	+=	wt[j];
         
         W	+=	log(sumw);

	 if (!strcmp(s, "old"))				// grow old conf.
            n	=	0;				// pick old position
         else						// grow new conf.
            n	=	Select(wt, sumw);		// select one trial pos.

         molm->p[i]	=	pt[n];
         v[ib]		=	vstructSum(v+ib, Vt+n);
         if (V_VIRIAL)
            vir[ib]	=	wstructSum(vir+ib, VIRt+n);

#ifdef CELL_LIST
         CL_Add(molm, i);				// add trial site into cell
#endif
         molm->flags[i]	=	1;			

	 //printf("Adding %d done, trial #%d.\n", i, n);
      }		// not the atom #0 done
   }		// all atoms done
   if (!strcmp(s, "old")) {
      vstructNegate(v+ib);
      wstructNegate(vir+ib);
   }
   return	W;
}


void EnergyCheck()		// Check step-by-step energy update correctness
{
   double	*V, *VIR;
   long		i;

   V	=	(double *) calloc(NBOX, sizeof(double));
   VIR	=	(double *) calloc(NBOX, sizeof(double));

   for (i=0; i<NBOX; i++) {
      V[i]	=	v[i].tot;
      VIR[i]	=	vir[i].tot;
   }
   CalcV();

   for (i=0; i<NBOX; i++) {
      if ( fabs(V[i] - v[i].tot) > 1.0e-6)
          fprintf(foutput, "\n#############total energy problem in box %d.###############\n", i);
      if ( fabs(VIR[i] - vir[i].tot) > 1.0e-6)
          fprintf(foutput, "\n#############total virial problem in box %d.###############\n", i);

      fprintf(foutput, "\n\tBOX[%d]\n\n", i);
      fprintf(foutput, "\tTotal energy end of simulation:\t%f\n", V[i]);
      fprintf(foutput, "\t      energy recalculated:\t%f\n", v[i].tot);
      fprintf(foutput, "\tTotal virial end of simulation:\t%f\n", VIR[i]);
      fprintf(foutput, "\t      virial recalculated:\t%f\n", vir[i].tot);
      fprintf(foutput, "\n");
   }

   free(VIR);
   free(V);
   return;
}
                                                                                                                                                              src/roots.c                                                                                         0000644 0143352 0000144 00000022506 11200052666 012401  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	roots.c
    author:	Pieter J. in 't Veld
    date:	May 10, 2001
    purpose:	roots of upto fourth order polynomials
*/
#define __ROOTS_MODULE
#include "roots.h"
#include <stdio.h>
#define R_VERSION	"3.2"

#define R_ZERO		1e-14
#define R_CZERO		1e-8

long R_Round(long n, double complex *z)
{
  long			i;

  for (i=0; i<n; ++i)
  {
    if ((creal(z[i])<R_CZERO)&&(creal(z[i])>-R_CZERO)) 
       z[i]	= 0.0 + cimag(z[i]) * I;	// z[i].re = 0
    if ((cimag(z[i])<R_CZERO)&&(cimag(z[i])>-R_CZERO)) 
       z[i]	= creal(z[i]) + 0.0 * I;	// z[i].im = 0
  }
  return n;
}


long R_ComplexCount(long n, double complex *z)
{
  long			i, count = 0;

  for (i=0; i<n; ++i)
    if (cimag(z[i])) ++count;
  return count;
}


long R_Sort(long n, double complex *z)
{
  long			i, j;
  double complex		zt;
  
  for (i=0; i<n-1; ++i)
    for (j=i+1; j<n; ++j)
      if (creal(z[j]) < creal(z[i]))
      {
        zt		= z[i];
	z[i]		= z[j];
	z[j]		= zt;
      }
      else
        if ((creal(z[i])==creal(z[j]))&&(cimag(z[j])<cimag(z[i])))
	{
          zt		= z[i];
	  z[i]		= z[j];
	  z[j]		= zt;
	}
  return n;
}


long R_SolveFirstOrder(double *c, double complex *z)	// solve c[0] + c[1]*z = 0
{
  if (fabs(c[1])<R_ZERO)
    return 0;

  *z	=	-c[0]/c[1] + 0 * I;
  return 1;
}


long R_SolveSecondOrder(double *c, double complex *z)	// c[0] + c[1] * z + c[2] * z^2 = 0
{
  double complex	zR = c[1]*c[1]-4.0*c[0]*c[2] + 0.0 * I;

  if (fabs(c[2])<R_ZERO)
    return R_SolveFirstOrder(c, z);
  
  if (fabs(creal(zR))<R_ZERO)
    zR	=	0.0;	// zR.re = 0.0;

  zR	=	cpow(zR, 0.5);

  zR	/=	2.0 * c[2];
  z[0]	=	-c[1]/c[2]/2.0 + 0.0 * I;
  z[1]	=	z[0];

  z[0]	-=	zR;
  z[1]	+=	zR;

  return R_Round(2, z);
}


long R_SolveThirdOrder(double *c, double complex *z)	// c[0]+c[1]z+c[2]z^2+c[3]z^3=0
{
  double		f1, f2, f3;
  double complex	A, A1, B, z1, z2;

  if (fabs(c[3])<R_ZERO)
    return R_SolveSecondOrder(c, z);
  
  // real constants
  
  f1			= 3.0*c[1]*c[3]-c[2]*c[2];
  f2			= (c[2]*(9.0*c[1]*c[3]-2.0*c[2]*c[2])
  			  -27.0*c[0]*c[3]*c[3])/2.0;
  f3			= -c[2]/(3.0*c[3]);
 
  // complex constants
  A1		=	0.0 + 0.0 * I; 
  if (fabs(f1)<1e-3)					// Limit f1->0
  {
    A		=	(f2<0.0 ? pow(-2.0*f2, 1.0/3.0)/(3.0*c[3]) : 0.0) + 0.0 * I;
  } 
  else if (fabs(f2)<R_ZERO)				// Limit f2->0
  {
    A		=	sqrt(fabs(f1))/(3.0*c[3]) + 0.0 * I;
    if (f1 < 0.0) 
       A	=	creal(A)*sqrt(3.0)/2.0 - creal(A)/2.0 * I;
    else
       A	=	creal(A) + 0.0 * I;
  }
  else
  {
    A1	=	f1*f1*f1 + f2*f2 + 0.0 * I;
    A1	=	cpow(A1, 0.5);
    A1	+=	f2;
    A1	=	cpow(A1, 1.0/3.0);
    *z	=	f1/(3.0*c[3]) + 0.0 * I;
    A	=	(*z)/A1;
  }
  B	=	A/(3.0*c[3]);
  z1	=	1.0/2.0 + sqrt(3.0)/2.0 * I;
  z2	=	conj(z1);
  
  // roots
  z[0]	=	f3+creal(B)-creal(A) + (cimag(B)-cimag(A))*I;
  z[1]	=	z1 * A;
  A1	=	z2 * B;
  z[1]	+=	f3 - A1;
  z[2]	=	z2 * A;
  A1	=	z1 * B;
  z[2]	+=	f3 - A1;
  
  return R_Round(3, z);
}


long R_SolveFourthOrder(double *c, double complex *z)
{
  double		f1, f2, f3, f4 ;
  double complex	A, B, 
  			A1, A11, A111, B1, B2;

  if (fabs(c[4])<R_ZERO)
    return R_SolveThirdOrder(c, z);
  
  // real constant
  
  f1			= c[2]*c[2]-3.0*c[1]*c[3]+12.0*c[0]*c[4];
  f2			= (c[2]*(2.0*c[2]*c[2]-9.0*c[1]*c[3]-72.0*c[0]*c[4])+
  			  27.0*(c[0]*c[3]*c[3]+c[1]*c[1]*c[4]))/2.0;
  f3			= c[3]*c[3]/(4.0*c[4]*c[4])-2.0*c[2]/(3.0*c[4]);
  f4			= -c[3]*c[3]*c[3]/(c[4]*c[4]*c[4])
  		          +4.0*c[2]*c[3]/(c[4]*c[4])-8.0*c[1]/c[4];

  // complex constants
 
  A11	=	0.0 + 0.0 * I; 
  if (fabs(f1)<R_ZERO)				// limit f1->0
     A11	=	(f2>0.0 ? 2.0*f2 : 0.0) + 0.0 * I;
  else if (fabs(f2)<R_ZERO)			// limit f2->0
     A11	=	(f1>0.0 ? sqrt(f1/3.0) : 0.0) + 0.0 * I;
  else {
     A111	=	f2*f2 - f1*f1*f1 + 0.0 * I;
     A111	=	cpow(A111, 0.5);
     A111	+=	f2;
     A111	=	cpow(A111, 1.0/3.0);
     A		=	A111 * A111;
     A		+=	f1;
     A111	*=	3.0*c[4];
     A11	=	A/A111;
  }
  A1	=	A11+f3;
  A1	=	cpow(A1, 0.5);
  B1	=	2.0*f3-A11;

  B2	=	0.0 * I;
  if ((fabs(creal(A1))<5e-6)&&(fabs(cimag(A1))<5e-6)) {		// Limit A1->0
    A	=	c[0]/c[4] + cimag(A)*I;
    B2	=	creal(A)>=0.0 ? 0.0 : -4.0*sqrt(-creal(A));
  }
  else {
    A	=	f4/4.0 + 0.0 * I;
    B2	=	A/A1;
  }
  
  // roots

  A	=	-c[3]/(4.0*c[4])+A1/2.0;
  B	=	B1 + B2;
  B	=	cpow(B, 0.5);
  B	/=	2.0;
  z[0]	=	A+B;
  z[1]	=	A-B;
  
  A	=	-c[3]/(4.0*c[4])-A1/2.0;
  B	=	B1 - B2;
  B	=	cpow(B, 0.5);
  B	/=	2.0;
  z[2]	=	A+B;
  z[3]	=	A-B;
  return R_Round(4, z);
}


long R_SolveAnalytical(long n, double *c, double complex *z)
{
  switch (n)
  {
    case 1: return R_SolveFirstOrder(c, z);
    case 2: return R_SolveSecondOrder(c, z);
    case 3: return R_SolveThirdOrder(c, z);
    case 4: return R_SolveFourthOrder(c, z);
    default: return 0;
  }
}


#define R_MAXNROOTS	20
#define R_D2		1e-5
#define R_MAXITER	60
#define R_EPS		2.5e-13
#define sign(a)		((a<0) ? -1 : 1)

long R_NewtonRaphson(
  double y(double), double x_low, double x_high, double *root)
{
  static long		iter;
  static double		x0, x1, x2,
  			y0, y1, y2;

  x0			= x_low;
  y0			= y(x_low);
  x1			= x_high;
  y1			= y(x_high);
  iter			= 0;
  do
  {
    x2			= (x1*y0-x0*y1)/(y0-y1);
    y2			= y(x2);
    if (sign(y2)==sign(y1))
    {
      x0		= x1;
      y0		= y1;
      x1		= x2;
      y1		= y2;
    }
    else
    {
      x1		= x0;
      y1		= y0;
      x0		= x2;
      y0		= y2;
    }
  } while((fabs(x0-x1)>R_EPS)&&(++iter<R_MAXITER));
  *root			= x2;
  return (iter<R_MAXITER);
}

// Adapted from 
//   W.H. Press, Numerical Recipes in C, Chapter 9.3, Cambridge University
//   Press, 1997 (p. 361)

//long			R_NYCALL = 0;

long R_BrentInt(
  double y(double), double x_low, double x_high, 
  double y_low, double y_high, double tol, double *root)
{
  static long		iter;
  static double		a, b, c, d, e, min1,
  			min2, fa, fb, fc, p, q, r, s, tol1,
			xm;

  a			= x_low;
  b			= x_high;
  c			= x_high;
  fa			= y_low; 
  fb			= y_high;
  //if (sign(fa)!=sign(fb))
  {
    fc			= fb;
    for (iter=0; iter<R_MAXITER; ++iter)
    {
      if (((fb>0.0)&&(fc>0.0))||((fb<0.0)&&(fc<0.0)))
      {
        c		= a;
	fc		= fa;
	e		= d = b-a;
      }
      if (fabs(fc)<fabs(fb))
      {
        a		= b;
	b		= c;
	c		= a;
	fa		= fb;
	fb		= fc;
	fc		= fa;
      }
      tol1		= 2.0*R_EPS*fabs(b)+0.5*tol;
      xm		= 0.5*(c-b);
      if ((fabs(xm)<=tol1)||(fb==0.0))
      {
	//if (fabs(fb)<1e-10)
	{
	  *root		= b;
          return 1;
	}
	//else 
	//  return 0;
      }
      if ((fabs(e)>=tol1)&&(fabs(fa)>fabs(fb)))
      {
        s		= fb/fa;
	if (a==c)
	{
	  p		= 2.0*xm*s;
	  q		= 1.0-s;
	}
	else
	{
	  q		= fa/fc;
	  r		= fb/fc;
	  p		= s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
	  q		= (q-1.0)*(r-1.0)*(s-1.0);
	}  
	if (p>0.0)
	  q		= -q;
	p		= fabs(p);
	min1		= 3.0*xm*q-fabs(tol1*q);
	min2		= fabs(e*q);
	if (2.0*p<(min1<min2 ? min1 : min2))
	{
	  e		= d;
	  d		= p/q;
	}
	else
	{
	  d		= xm;
	  e		= d;
	}
      }
      else
      {
        d		= xm;
	e		= d;
      }
      a			= b;
      fa		= fb;
      if (fabs(d)>tol1)
        b		+= d;
      else
        b		+= sign(xm)*tol1;
      fb		= y(b);
      //++R_NYCALL;
    }
  }
  return 0;
}


long R_Brent(
  double y(double), double x_low, double x_high, double tol, double *root)
{
  return R_BrentInt(y, x_low, x_high, y(x_low), y(x_high), tol, root);
}


#define R_DX		1e-6
#define R_TOL		1e-12

long R_SolveNumericalInt(
  long recursive,
  double y(double), 
  double x_low, double x_high, 
  double y_low, double y_high, double dy_low, double dy_high,
  long n_mesh, long depth,
  double *roots, long n_max)
{
  long			i, nr, k, abort = 0;
  double	 	x, x_old, x_step, 
			yx, yx_old, 
			dx, dy, dy_old;
  
  if (depth<=0) 
    return 0;
  
  x_step		= (x_high-x_low)/n_mesh;
  dx			= R_DX*x_step;
  x			= x_low;
  yx			= y_low;
  dy			= dy_low;
  if (fabs(yx)<R_TOL)
  {
    *roots		= x;
    nr			= 1;
    x			+= x_step;
    yx			= y(x);
    dy			= y(x+dx)-yx;
    //R_NYCALL		+= 2;
  }
  else
    nr			= 0;
  for (i=nr ? 2 : 1; (i<=n_mesh)&&(nr<R_MAXNROOTS)&&(!abort); ++i)
  {
    // Function scan
    
    x_old		= x;
    yx_old		= yx;
    dy_old		= dy;
    if (i<n_mesh)
    {
      x			= x_low+i*x_step;
      yx		= y(x);
      dy		= y(x+dx)-yx;
      //R_NYCALL		+= 2;
    }
    else
    {
      x			= x_high;
      yx		= y_high;
      dy		= dy_high;
    }
    if (fabs(yx)>R_TOL)
    { 
      // Combined bisection and Newton-Raphson scheme

      if (sign(yx_old) != sign(yx))
      {
        //fprintf(stderr, "d = %d, x E [%g, %g], y E [%g, %g]",
	//  depth, x_old, x, yx_old, yx);
        if (R_BrentInt(y, x_old, x, yx_old, yx, R_TOL, roots+nr))
	{
	  //fprintf(stderr, ", x[%d] = %g", nr, roots[nr]);
	  ++nr;
	}
	//fprintf(stderr, "\n");
      }
      else
        if (sign(dy_old) != sign(dy))
	{
	  nr		+= k = R_SolveNumericalInt(
	  			1, y, x_old, x, yx_old, yx, dy_old, dy, 
				4, depth-1, roots+nr, n_max);
	  abort		=  k==0 ? recursive : 0;
	}
    }
    else
    {
      if (nr<n_max)
        roots[nr++]	= x;
      else
	return nr;
    }
  }
  return nr;
}


long R_SolveNumerical(
  double y(double), double x_low, double x_high, long n_mesh, long depth,
  double *roots, long n_max)
{
  long			nr;
  double		dx, y_low, y_high;

  //R_NYCALL		= 0;
  dx			= R_DX*(x_high-x_low)/n_mesh;
  x_low			+= dx;
  x_high		-= dx;
  y_low			= y(x_low);
  y_high		= y(x_high);
  dx			= R_DX*(x_high-x_low)/n_mesh;
  nr 			= R_SolveNumericalInt(0, y, x_low, x_high, 
    y_low, y_high, y(x_low+dx)-y_low, y_high-y(x_high-dx), 
    n_mesh, depth, roots, n_max);
  //fprintf(stderr, "R_NYCALL = %d\n", R_NYCALL);
  return nr;
}

                                                                                                                                                                                          src/sample.c                                                                                        0000600 0143352 0000144 00000245351 11555623455 012526  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    sample.c
    author:     Peng Yi at MIT
    date:       October 22, 2006
		October 13, 2007, adding Pieter's functions
    purpose:    Sample system properties
    notes:	April 18,2011	Add structure factor calculation based on definition
*/

#define __SAMPLE_MODULE
#include "sample.h"


long mod (long numb, long divisor)
{
  while (numb<0) {
    numb+=divisor;
  }
  return numb%divisor;
}


long factorial(long number)			//calculate factorial n! 
{
   long	result;

   switch(number){
	case 12:		result	=	479001600;		break;
	case 11:		result	=	39916800;		break;
	case 10:		result	=	3628800;		break;
	case 9:		result	=	362880;		break;
	case 8:		result	=	40320;		break;
	case 7:		result	=	5040;		break;
	case 6:		result	=	720;		break;
	case 5:		result	=	120;		break;
	case 4:		result	=	24;		break;
	case 3:		result	=	6;		break;
	case 2:		result	=	2;		break;
	case 1: case 0:		result	=	1;		break;
   }
   return	result;
/*
  if (number<=1)
    return 1;
  else
    return (number*factorial(number-1));
*/
}

long intpow(long base, long power)		//calculate integer power, base^power
{
   long		i, result;
   result	=	1;
   for (i=0; i<power; i++) {
      result	*= base;
   }
   return 	result;
}

long intlog(long base, long value)		//calculate log_base^value
{
   long		power, result;
   power	=	0;
   result	=	1;
   while (result<value) {
      result	*=	base;
      power	++;
   }
   return	power;
}

double plgndr(int l, int m, double x)  		// associated Legendre polynomial
{                                      		// m=0, 1, ..., l, x = [-1,1]
  double fact, pll, pmm, pmmp1, somx2;
  int i, ll;

  if (m<-l || m>l || fabs(x)>1.0) {
    printf("Bad arguments in routine plgndr.\n");
    return 0;
  }

//  if (m<0)
//    return ((mod(-m,2)==0)?1.0:-1.0) * factorial(l+m) / factorial(l-m) * plgndr(l, -m, x);

  pmm = 1.0;                           		//compute P_m^m
  if (m>0) {
    somx2 = sqrt((1.0-x)*(1.0+x));
    fact = 1.0;
    for (i=1; i<=m; i++) {
      pmm *= -fact*somx2;
      fact += 2.0;
    }
  }
  if (l==m)
    return pmm;
  else {                              		//compute P_{m+1}^m
    pmmp1 = x*(2*m+1)*pmm;
    if (l==(m+1))
      return pmmp1;
    else {                            		//compute P_l^m, l>m+1
      for (ll=m+2; ll<=l; ll++) {
	pll = (x*(2*ll-1)*pmmp1-(ll+m-1)*pmm)/(ll-m);
	pmm = pmmp1;
	pmmp1 = pll;
      }
      return pll;
    }
  }
}
/*
complex sfharmonics2(int l, int m, double costheta, double cosphi, double sinphi)	//spherical harmonics functions
{								//m = -l, -l+1, ..., l-1, l
  double 	Nlm, Plm;
  complex	Ylm;
  complex	eiphi, eimphi;
  int		i;

  if (abs(m)>l || fabs(costheta)>1.0) {
    printf("Bad arguments in routine sfharmonics.\n");
    return 0;
  }

  eiphi		=	cosphi + I*sinphi;
  eimphi	=	1.0;
  for (i=0; i<abs(m); i++) {
     eimphi	*=	(cosphi + I*sinphi);
  }

  Plm 	= 	plgndr(l, abs(m), costheta);
  Nlm	=	sqrt( (double) (2*l+1)/4/pi * factorial(l-abs(m)) / factorial(l+abs(m)));
//  Ylm	=	Nlm * Plm * cexp(I*abs(m)*phi);
  Ylm	=	Nlm * Plm * eimphi;

  if (m >= 0)
    return Ylm;
  else
//    return pow(-1,abs(m)) * conj(Ylm);
    return (mod(abs(m),2) ? -1 : 1) * conj(Ylm);
}
*/

double complex sfharmonics2(int l, int m, double costheta, double cosphi, double sinphi)
{
   double		SIN, SIN2, COS, COS2;
   double complex	eiphi, eimphi;		//exp(I*phi), exp(I*m*phi)
   int			i;
   double complex	Ylm;

   COS	=	costheta;
   COS2	=	COS * COS;
   SIN2	=	1-COS2;
   SIN	=	sqrt(SIN2);
   
   eiphi	=	cosphi + I*sinphi;
   eimphi	=	1.0;
   for (i=0; i<abs(m); i++) {
      eimphi	*=	eiphi;
   }

   switch (l) {
      case 6: {
         switch (abs(m)) {
            case 6:	Ylm	=	1.0/64 * sqrt(3003.0/pi) * eimphi * SIN2 * SIN2 * SIN2;
			break;
            case 5:	Ylm	=	-3.0/32 * sqrt(1001.0/pi) * eimphi * SIN2 * SIN2 * SIN * COS;
			break;
	    case 4:	Ylm	=	3.0/32 * sqrt(91.0/2/pi) * eimphi * SIN2 * SIN2 * (11*COS2-1);
			break;
            case 3:	Ylm	=	-1.0/32 * sqrt(1365.0/pi) * eimphi * SIN2 * SIN * (11*COS2-3) * COS;
			break;
	    case 2:	Ylm	=	1.0/64 * sqrt(1365.0/pi) * eimphi * SIN2 * (33*COS2*COS2 -18*COS2 +1);
			break;
	    case 1:	Ylm	=	-1.0/16 * sqrt(273.0/2/pi) * eimphi * SIN * (33*COS2*COS2-30*COS2+5) *COS;
			break;
	    case 0:	Ylm	=	1.0/32 * sqrt(13.0/pi) * (231*COS2*COS2*COS2-315*COS2*COS2+105*COS2-5);
			break;
	    default:	printf("l, m mismatch!\n");
			break;
	 }
         break;
      }
      case 4: {
         switch (abs(m)) {
            case 4:	Ylm	=	3.0/16 * sqrt(35.0/2/pi) * eimphi * SIN2 * SIN2;
			break;
	    case 3:	Ylm	=	-3.0/8 * sqrt(35.0/pi) * eimphi * SIN2 * SIN * COS;
			break;
	    case 2:	Ylm	=	3.0/8 * sqrt(5.0/2/pi) * eimphi * SIN2 * (7*COS2-1);
			break;
	    case 1:	Ylm	=	-3.0/8 * sqrt(5.0/pi) * eimphi * SIN * (7*COS2-3) * COS;
			break;
	    case 0:	Ylm	=	3.0/16 * sqrt(1.0/pi) * (35*COS2*COS2 - 30*COS2 +3);
			break;	
	    default:	printf("l, m mismatch!\n");
			break;
	 }
	 break;
      }
      default:	printf("could not find the l!\n");
		break;
   }

   if (m>=0)
      return	Ylm;
   else
      return	pow(-1, abs(m)) * conj(Ylm);
}

/*
complex sfharmonics(int l, int m, double costheta, double phi)	//spherical harmonics functions
{								//m = -l, -l+1, ..., l-1, l
  double 	Nlm, Plm;
  complex	Ylm;

  if (abs(m)>l || fabs(costheta)>1.0) {
    printf("Bad arguments in routine sfharmonics.\n");
    return 0;
  }

  Plm 	= 	plgndr(l, abs(m), costheta);
  Nlm	=	sqrt( (double) (2*l+1)/4/pi * factorial(l-abs(m)) / factorial(l+abs(m)));
  Ylm	=	Nlm * Plm * cexp(I*abs(m)*phi);

  if (m >= 0)
    return Ylm;
  else
//    return pow(-1,abs(m)) * conj(Ylm);
    return (mod(abs(m),2) ? -1 : 1) * conj(Ylm);
}
*/

double complex sfharmonics(int l, int m, double costheta, double phi)	//table from Wikipedia spherical harmonics entry
{
   double		SIN, SIN2, COS, COS2;
   double complex	Ylm;
   COS	=	costheta;
   COS2	=	COS * COS;
   SIN2	=	1-COS2;
   SIN	=	sqrt(SIN2);
   
   switch (l) {
      case 6: {
         switch (abs(m)) {
            case 6:	Ylm	=	1.0/64 * sqrt(3003.0/pi) * cexp(6*I*phi) * SIN2 * SIN2 * SIN2;
			break;
            case 5:	Ylm	=	-3.0/32 * sqrt(1001.0/pi) * cexp(5*I*phi) * SIN2 * SIN2 * SIN * COS;
			break;
	    case 4:	Ylm	=	3.0/32 * sqrt(91.0/2/pi) * cexp(4*I*phi) * SIN2 * SIN2 * (11*COS2-1);
			break;
            case 3:	Ylm	=	-1.0/32 * sqrt(1365.0/pi) * cexp(3*I*phi) * SIN2 * SIN * (11*COS2-3) * COS;
			break;
	    case 2:	Ylm	=	1.0/64 * sqrt(1365.0/pi) * cexp(2*I*phi) * SIN2 * (33*COS2*COS2 -18*COS2 +1);
			break;
	    case 1:	Ylm	=	-1.0/16 * sqrt(273.0/2/pi) * cexp(I*phi) * SIN * (33*COS2*COS2-30*COS2+5) *COS;
			break;
	    case 0:	Ylm	=	1.0/32 * sqrt(13.0/pi) * (231*COS2*COS2*COS2-315*COS2*COS2+105*COS2-5);
			break;
	    default:	printf("l, m mismatch!\n");
			break;
	 }
         break;
      }
      case 4: {
         switch (abs(m)) {
            case 4:	Ylm	=	3.0/16 * sqrt(35.0/2/pi) * cexp(4*I*phi) * SIN2 * SIN2;
			break;
	    case 3:	Ylm	=	-3.0/8 * sqrt(35.0/pi) * cexp(3*I*phi) * SIN2 * SIN * COS;
			break;
	    case 2:	Ylm	=	3.0/8 * sqrt(5.0/2/pi) * cexp(2*I*phi) * SIN2 * (7*COS2-1);
			break;
	    case 1:	Ylm	=	-3.0/8 * sqrt(5.0/pi) * cexp(I*phi) * SIN * (7*COS2-3) * COS;
			break;
	    case 0:	Ylm	=	3.0/16 * sqrt(1.0/pi) * (35*COS2*COS2 - 30*COS2 +3);
			break;	
	    default:	printf("l, m mismatch!\n");
			break;
	 }
	 break;
      }
      default:	printf("could not find the l!\n");
		break;
   }

   if (m>=0)
      return	Ylm;
   else
      return	pow(-1, abs(m)) * conj(Ylm);
}

/*************************************************************************************/

#ifdef TEST
complex aveqlm (int l, int m, long i)	//calculate average q_lm value for particle i
{					//also update part[i].nbond & part[i].qlm !!
  long		jj, k, nbond;
  vector	dp;
  double	r2, alpha, alphasum;
  double	costheta, phi;    	//theta=[0,PI], phi=[0,2*PI]
  complex	aveqlm;

  nbond		=	0;  			//neighbor bond number
  aveqlm	=	0;
  alphasum	=	0;
 
  for (jj=0; jj<part[i].nverlet; jj++) {

    k		=	part[i].vlist[jj];
    dp.x	=	part[k].p.x	-	part[i].p.x;
    dp.y	=	part[k].p.y	-	part[i].p.y;
    dp.z	=	part[k].p.z	-	part[i].p.z;
    MapInBox(&dp);
    r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    if (r2 <= Rb*Rb) {
	nbond 	++;

	costheta =  dp.z/sqrt(r2);     	//theta=[0,PI]
//	if (dp.x >= 0)                 	//atan=[-PI/2, PI/2]
	  phi	=	atan2( dp.y, dp.x ); 
//	else
//	  phi	=	atan2( dp.y, dp.x ) + pi;

	alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / 
			(sigma*sigma + Rb*Rb - 2*sigma*Rb);
	aveqlm 		+=	sfharmonics(l, m, costheta, phi) * alpha;
	alphasum	+=	alpha;
	//alpha is a damping factor, as explained in JCP v104 (1996) 9932
    }
  }
  if ( fabs(alphasum) < ZERO) {		//always be awared of divided-by-zero
    aveqlm		=	0;
  }
  else {
    aveqlm		=	aveqlm/alphasum;
  }
  part[i].nbond		=	nbond;		//update part[i].nbond
  part[i].qlm[m+l]	=	aveqlm;		//update part[i].qlm
  return aveqlm;
}

/*
complex tildeqlm(int l, int m, long i)	//calculate normalized q_lm for particle i
{
  int n;
  double sum;
  
  for (n=-l; n<=l; n++) {
    sum +=	pow(cabs(aveqlm(l,n,i)), 2); 
  }
  sum = sqrt(sum);
  
  if ( fabs(sum) < ZERO)
    return 0;
  else
    return aveqlm(l, m, i)/sum;
}
*/

/*
double ql(int l, long i)		//calculate q_l of particle i
{
  int 		m;
  double 	ql=0;
  complex 	term; 

  for (m=-l; m<=l; m++) {
    term	=	aveqlm(l, m, i);
    ql		+=	term * conj(term);
  }
  ql		*=	(double)4*pi/(2*l+1);
  return sqrt(ql);
}
*/

/******************************************************************************/

double ql(int l, long i)		//calculate q_l of particle i
{
  int		m;
  double	ql=0;

  for (m=-l; m<=l; m++) {
    ql	+=	part[i].qlm[m+l] * conj(part[i].qlm[m+l]);
  }
  ql	*=	(double)4*pi/(2*l+1);
  return	sqrt(ql);
}

/*
complex qlproduct(int l, long i, long j)	//calculate q_l(i) dot q_l(j)
{
  int		m;
  complex	termi, termj, term=0;
  double	sumi=0, sumj=0;

  for (m=-l; m<=l; m++) {
    termi	=	aveqlm(l, m, i);
    termj	=	aveqlm(l, m, j);
    term	+=	termi * conj(termj);
    sumi	+=	termi * conj(termi);
    sumj	+=	termj * conj(termj);
  }
  return	term/sqrt(sumi*sumj);		//need to check divided-by-zero
}
*/

/*********************************************************************************/

/**********************************************************************************/
/*
complex aveQlm1(int l, int m) 		//calculate average Q_lm
{
  long		i;
  long		sumbond	=	0;
  complex 	aveQlm	=	0;
  complex	term;

  for (i=0; i<NPARTS; i++)  {
    term	=	aveqlm(l, m, i);	//notice because aveqlm(l,m,i) updates part[i].nbond,
    aveQlm 	+=	part[i].nbond * term;	//then cannot appear in one expression at same time,
    sumbond	+=	part[i].nbond;		//e.g., aveQlm = part[i].bond * aveqlm(l,m,i).
  }
  
  if (sumbond == 0)
    return	0;
  else
    return aveQlm/sumbond;
}
*/

/****************************************************************************************/
/*
double CalcQl1(int l)    			//calculate Q_l
{
  int		m;
  double	Ql=0;
  complex	term;
 
  for (m=-l; m<=l; m++) {
    term	=	aveQlm1(l, m);
    Ql		+=	term * conj(term);    
  }
  Ql *= 4*pi/(2*l+1);
  return sqrt(Ql);
}
*/

/*************************************************************************************/

complex aveQlm(int l, int m) 		//calculate average Q_lm
{
  long		i;
  long		sumbond	=	0;
  complex 	aveQlm	=	0;

  for (i=0; i<NPARTS; i++)  {
    aveQlm 	+=	part[i].nbond * part[i].qlm[m+l];
    sumbond	+=	part[i].nbond;		
  }
  
  if (sumbond == 0)
    return	0;
  else
    return aveQlm/sumbond;
}

/******************************************************************************/

double CalcQl(int l)
{
  long		i;
  int		m;
  long		sumnbond;
  complex	aveQlm;
  double	Ql	=	0;
  
  for (m=-l; m<=l; m++) {
 
    sumnbond	=	0;	//notice where to set Ql=0 and where 
    aveQlm	=	0;	//to set aveQlm and sumnbond =0 !!

    for (i=0; i<NPARTS; i++) {		//calculate aveQlm
      aveQlm	+=	part[i].nbond * part[i].qlm[m+l];
      sumnbond	+=	part[i].nbond;
    }
    if (sumnbond	==	0)
      aveQlm	=	0;
    else
      aveQlm	=	aveQlm/sumnbond;

    Ql	+=	aveQlm * conj(aveQlm);
  }
  Ql	*=	(double) 4*pi/(2*l+1);
  return sqrt(Ql);
}

/****************************************************************************************/
/* void New_Qlm(int l)									*/
/*  											*/
/* Calculate aveqlm and ql of every particle, and also calculate			*/
/* aveQlm and Ql for the system.							*/
/* And since aveqlm = ylmalphasum/alphasum, we save two variables 			*/
/* ylmalphasum and alphasum, as well as aveqlm for each particle. 			*/
/* Same for Q.										*/
/*											*/
/* Input: configuration, l = l_of_Ylm, Verlet list.					*/
/****************************************************************************************/

void New_Qlmtemp(int l)
{
  int		m;
  long		i, jj, j;
  vector	dp;
  double	r2, alpha;	//alpha is a weighted # of bonds, see JCP V104 (1996) 9932
  double	costheta, phi;

   for (i=0; i<NPARTS; i++) {
      parttemp[i]	=	part[i];
   }

  AlphaSumtemp	=	0;			//initialize AlphaSum and YlmAlphaSum for
  for (m=-l; m<=l; m++) {			//the calculation of Qlm and Ql
    YlmAlphaSumtemp[m+l]	=	0;	
  }
  for (i=0; i<NPARTS; i++) {			//initialize alphasum and ylmalphasum for
    parttemp[i].alphasum	=	0;		//the calculation of qlm and ql
    for (m=-l; m<=l; m++) {			//for each particle
      parttemp[i].ylmalphasum[m+l]	=	0;
    }
  }

  for (i=0; i<NPARTS; i++) {			//Qlm is averaged among All particles
    for (jj=0; jj<part[i].nverlet; jj++) {	//search the Verlet list
			
      j		=	part[i].vlist[jj];
      if (j > i) {				//reduce theta and phi evaluation

        dp.x	=	part[j].p.x	-	part[i].p.x;
        dp.y	=	part[j].p.y	-	part[i].p.y;
        dp.z	=	part[j].p.z	-	part[i].p.z;
        MapInBox(&dp);
        r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    
        if (r2 < Rb*Rb) {
	  alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);
	  parttemp[i].alphasum	+=	alpha;
	  parttemp[j].alphasum	+=	alpha;

          costheta	= dp.z/sqrt(r2);     	//theta=[0,PI]
	  phi		= atan2( dp.y, dp.x ); 
	  for (m=-l; m<=l; m++) {
	    parttemp[i].ylmalphasum[m+l]	+=	sfharmonics(l, m, costheta, phi) * alpha;
	    parttemp[j].ylmalphasum[m+l]	+=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
	  }
        }
      }
    } 
  }

  for (i=0; i<NPARTS; i++) {			//calculate qlm and ql for each particle
    parttemp[i].ql		=	0; 
    if ( parttemp[i].alphasum < 1e-12) {	//always be aware of divided-by-zero
      if (parttemp[i].alphasum < 0)		//alphasum should always be positive
         printf("New_Qlm error, alphasum_i = %f\n", parttemp[i].alphasum);
      for (m=-l; m<=l; m++) {
        if (cabs(parttemp[i].ylmalphasum[m+l]) > ZERO) 
	  printf("Error, ylmalphasum[%d] = %f+%f*i\n", m, creal(parttemp[i].ylmalphasum[m+l]), 
			cimag(parttemp[i].ylmalphasum[m+l]));
        parttemp[i].qlm[m+l] 	=	0;	//update qlm for each particle and each m
      }
      parttemp[i].ql		=	0;
    }
    else {
      for (m=-l; m<=l; m++) {
        parttemp[i].qlm[m+l] 	=	parttemp[i].ylmalphasum[m+l] / parttemp[i].alphasum;
	parttemp[i].ql		+=	parttemp[i].qlm[m+l] * conj(parttemp[i].qlm[m+l]);
      }
    }
    parttemp[i].ql	=	sqrt(parttemp[i].ql*4*pi/(2*l+1));
    
    AlphaSumtemp	+=	parttemp[i].alphasum;	//sum over all particles to calculate Qlm
    for (m=-l; m<=l; m++) {
      YlmAlphaSumtemp[m+l]	+=	parttemp[i].ylmalphasum[m+l];
    } 
  }
  
  Qltemp		=	0;			//global variable
  if (AlphaSumtemp < 1e-12) {			//calculate aveQlm and Ql
      if (AlphaSumtemp < 0)
         printf("New_Qlm error, AlphaSum = %f\n", AlphaSumtemp);
    for (m=-l; m<=l; m++) {
	Qlmtemp[m+l]	=	0;
    }
    Qltemp	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
	Qlmtemp[m+l]	=	YlmAlphaSumtemp[m+l]/AlphaSumtemp;
	Qltemp		+=	Qlmtemp[m+l] * conj(Qlmtemp[m+l]);
    }
    Qltemp	=	sqrt(Qltemp*4*pi/(2*l+1));
  }

  return;	//variables updated are global, so no need to have a return value
}

/****************************************************************************************/

void New_Qlm(int l)		//calculate q and Q using verlet list	
{
  int		m;
  long		i, jj, j, k, icell, jcell;
  vector	dp;
  double	r2, rxy, alpha;	//alpha is a weighted # of bonds, see JCP V104 (1996) 9932
  double	costheta, phi, cosphi, sinphi;

  AlphaSum	=	0;			//initialize AlphaSum and YlmAlphaSum for
  for (m=-l; m<=l; m++) {			//the calculation of Qlm and Ql
    YlmAlphaSum[m+l]	=	0;	
  }
  for (i=0; i<NPARTS; i++) {			//initialize alphasum and ylmalphasum for
    part[i].alphasum	=	0;		//the calculation of qlm and ql
    for (m=-l; m<=l; m++) {			//for each particle
      part[i].ylmalphasum[m+l]	=	0;
    }
  }

#ifdef VERLET_LIST
  for (i=0; i<NPARTS-1; i++) {			//Qlm is averaged among All particles
    for (jj=0; jj<part[i].nverlet; jj++) {	//search the Verlet list
			
      j		=	part[i].vlist[jj];

      if (j > i) {				//reduce theta and phi evaluation
      {
#elif CELL_LIST
  for (i=0; i<NPARTS-1; i++) {
     icell	=	part[i].icell;

     for (jj=0; jj<Cell[icell].nneigh; jj++) {
        jcell	=	Cell[icell].neigh[jj];

        for (k=0; k<Cell[jcell].sites; k++) {
           j	=	Cell[jcell].list[k]; 

           if (j > i) {
#else
  for (i=0; i<NPARTS-1; i++) {
     for (j=i+1; j<NPARTS; j++) {
        {{					//match the # of loops

#endif
        dp.x	=	part[j].p.x	-	part[i].p.x;
        dp.y	=	part[j].p.y	-	part[i].p.y;
        dp.z	=	part[j].p.z	-	part[i].p.z;
        MapInBox(&dp);
        r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    
        if (r2 < Rb*Rb) {
	  alpha	= 1.0;
	  //alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);
	  part[i].alphasum	+=	alpha;
	  part[j].alphasum	+=	alpha;

          costheta	= dp.z/sqrt(r2);     	//theta=[0,PI]
//	  phi		= atan2( dp.y, dp.x ); 
	  rxy		=	sqrt(dp.x * dp.x + dp.y * dp.y);
	  cosphi	=	dp.x / rxy;
	  sinphi	=	dp.y / rxy;
	  if (rxy<ZERO) {
	     cosphi	=	0;
	     sinphi	=	0;
          }

	  for (m=-l; m<=l; m++) {
	    part[i].ylmalphasum[m+l]	+=	sfharmonics2(l, m, costheta, cosphi, sinphi) * alpha;
	    part[j].ylmalphasum[m+l]	+=	sfharmonics2(l, m, -costheta, -cosphi, -sinphi) * alpha;
	    
//	    part[i].ylmalphasum[m+l]	+=	sfharmonics(l, m, costheta, phi) * alpha;
//	    part[j].ylmalphasum[m+l]	+=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
	  }
        }
      }
    } 
  }}

  for (i=0; i<NPARTS; i++) {			//calculate qlm and ql for each particle
    part[i].ql		=	0; 
    if ( part[i].alphasum < 1e-12) {	//always be aware of divided-by-zero
      if (part[i].alphasum < 0)		//alphasum should always be positive
         printf("New_Qlm error, alphasum_i = %f\n", part[i].alphasum);
      for (m=-l; m<=l; m++) {
        if (cabs(part[i].ylmalphasum[m+l]) > ZERO) 
	  printf("Error, ylmalphasum[%d] = %f+%f*i\n", m, creal(part[i].ylmalphasum[m+l]), 
			cimag(part[i].ylmalphasum[m+l]));
        part[i].qlm[m+l] 	=	0;	//update qlm for each particle and each m
      }
      part[i].ql		=	0;
    }
    else {
      for (m=-l; m<=l; m++) {
        part[i].qlm[m+l] 	=	part[i].ylmalphasum[m+l] / part[i].alphasum;
	part[i].ql		+=	part[i].qlm[m+l] * conj(part[i].qlm[m+l]);
      }
    }
    part[i].ql	=	sqrt(part[i].ql*4*pi/(2*l+1));
    
    AlphaSum	+=	part[i].alphasum;	//sum over all particles to calculate Qlm
    for (m=-l; m<=l; m++) {
      YlmAlphaSum[m+l]	+=	part[i].ylmalphasum[m+l];
    } 
  }
  
  Ql		=	0;			//global variable
  if (AlphaSum < 1e-12) {			//calculate aveQlm and Ql
      if (AlphaSum < 0)
         printf("New_Qlm error, AlphaSum = %f\n", AlphaSum);
    for (m=-l; m<=l; m++) {
	Qlm[m+l]	=	0;
    }
    Ql	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
	Qlm[m+l]	=	YlmAlphaSum[m+l]/AlphaSum;
	Ql		+=	Qlm[m+l] * conj(Qlm[m+l]);
    }
    Ql	=	sqrt(Ql*4*pi/(2*l+1));
  }

  return;	//variables updated are global, so no need to have a return value
}

/****************************************************************************************/

void New_Qlm_NoVlist(int l)	//calculate New Qlm without using Verlet list	
{
  int		m;
  long		i, jj, j;
  vector	dp;
  double	r2, alpha;	//alpha is a weighted # of bonds, see JCP V104 (1996) 9932
  double	costheta, phi;

  AlphaSum	=	0;			//initialize AlphaSum and YlmAlphaSum for
  for (m=-l; m<=l; m++) {			//the calculation of Qlm and Ql
    YlmAlphaSum[m+l]	=	0;	
  }
  Ql		=	0;
  for (i=0; i<NPARTS; i++) {			//initialize alphasum and ylmalphasum for
    part[i].alphasum	=	0;		//the calculation of qlm and ql
    for (m=-l; m<=l; m++) {			//for each particle
      part[i].ylmalphasum[m+l]	=	0;
    }
    part[i].ql		=	0; 
  }

  for (i=0; i<NPARTS-1; i++) {			//Qlm is averaged among All particles
     for (j=i+1; j<NPARTS; j++) {		

        dp.x	=	part[j].p.x	-	part[i].p.x;
        dp.y	=	part[j].p.y	-	part[i].p.y;
        dp.z	=	part[j].p.z	-	part[i].p.z;
        MapInBox(&dp);
        r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    
        if (r2 < Rb*Rb) {
	  alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);
	  part[i].alphasum	+=	alpha;
	  part[j].alphasum	+=	alpha;

          costheta	= dp.z/sqrt(r2);     	//theta=[0,PI]
	  phi		= atan2( dp.y, dp.x ); 
	  for (m=-l; m<=l; m++) {
	    part[i].ylmalphasum[m+l]	+=	sfharmonics(l, m, costheta, phi) * alpha;
	    part[j].ylmalphasum[m+l]	+=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
	  }
        }
     } 
  }

  for (i=0; i<NPARTS; i++) {			//calculate qlm and ql for each particle
    if ( part[i].alphasum < 1e-12) {	//always be aware of divided-by-zero
      if (part[i].alphasum < 0)		//alphasum should always be positive
         printf("New_Qlm error, alphasum_i = %f\n", part[i].alphasum);
      for (m=-l; m<=l; m++) {
        if (cabs(part[i].ylmalphasum[m+l]) > ZERO) 
	  printf("Error, ylmalphasum[%d] = %f+%f*i\n", m, creal(part[i].ylmalphasum[m+l]), 
			cimag(part[i].ylmalphasum[m+l]));
        part[i].qlm[m+l] 	=	0;	//update qlm for each particle and each m
      }
      part[i].ql		=	0;
    }
    else {
      for (m=-l; m<=l; m++) {
        part[i].qlm[m+l] 	=	part[i].ylmalphasum[m+l] / part[i].alphasum;
	part[i].ql		+=	part[i].qlm[m+l] * conj(part[i].qlm[m+l]);
      }
    }
    part[i].ql	=	sqrt(part[i].ql*4*pi/(2*l+1));
    
    AlphaSum	+=	part[i].alphasum;	//sum over all particles to calculate Qlm
    for (m=-l; m<=l; m++) {
      YlmAlphaSum[m+l]	+=	part[i].ylmalphasum[m+l];
    } 
  }
  
  if (AlphaSum < 1e-12) {			//calculate aveQlm and Ql
      if (AlphaSum < 0)
         printf("New_Qlm error, AlphaSum = %f\n", AlphaSum);
    for (m=-l; m<=l; m++) {
	Qlm[m+l]	=	0;
    }
    Ql	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
	Qlm[m+l]	=	YlmAlphaSum[m+l]/AlphaSum;
	Ql		+=	Qlm[m+l] * conj(Qlm[m+l]);
    }
    Ql	=	sqrt(Ql*4*pi/(2*l+1));
  }

  return;	//variables updated are global, so no need to have a return value
}
/****************************************************************************************/

/*
void local_q_update()	//q update of all particles affected by a move
{
  long	jj;
  int	m;
  
  for (jj=0; jj<nverletplus; jj++) {	
    for (m=-l_of_Ylm; m<=l_of_Ylm; m++) {
      part[vlistplus[jj]].qlm[m+l_of_Ylm]	=	aveqlm(l_of_Ylm, m, vlistplus[jj]);
    }
  }
  return;
}
*/

/****************************************************************************************/
/* void Update_Qlm(int l, long i, vector p_old, vector p_new)				*/
/*											*/
/* Update Qlm and Ql due to displacement of particle i.  Update qlm also.		*/
/* The update of qlm is through update of ylmalphasum and alphasum!			*/
/* And the same is for Qlm.								*/
/*											*/
/* Input: configuration, verlet listplus						*/
/*        correct Q and q from last step.						*/
/* Caution: when use this function to recover Q and q, we need to be careful of		*/
/*	    the role of oldmol and part[n].						*/
/****************************************************************************************/

void Update_Qlm(int l, long i, vector p_old, vector p_new)
{			
  int		m;
  long		jj, j, k, jcell;
  vector	dp;
  double	r2, rxy, alpha, costheta, phi, cosphi, sinphi;
  complex	termi, termj;

#ifdef VERLET_LIST
  for (jj=0; jj<nverletplus-1; jj++) {	//the last element in vlistplus is itself

    j	=	vlistplus[jj];
    {{
#elif CELL_LIST
  for (jj=0; jj<nneighcellplus; jj++){
     jcell	=	neighcellplus[jj];

     for (k=0; k<Cell[jcell].sites; k++) {
        j	=	Cell[jcell].list[k];

        if (j!=i) {
#else
    {{{
#endif

    dp.x	=	part[j].p.x	-	p_old.x;	//old i, j distance
    dp.y	=	part[j].p.y	-	p_old.y;
    dp.z	=	part[j].p.z	-	p_old.z;
    MapInBox(&dp);
    r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;

    if (r2 < Rb*Rb) {
      //alpha	=	1.0;
      alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);
      part[i].alphasum	-=	alpha;		//substract contribution to qlm(i)
      part[j].alphasum	-=	alpha;		//substract contribution to qlm(j)
      AlphaSum		-=	2*alpha;	//substract contribution to Qlm

      costheta	=	dp.z/sqrt(r2);
//      phi	=	atan2( dp.y, dp.x ); 
      rxy		=	sqrt(dp.x*dp.x + dp.y*dp.y);
      cosphi	=	dp.x / rxy;
      sinphi	=	dp.y / rxy;
      if (rxy<ZERO) {
         cosphi	=	0;
         sinphi	=	0;
      }

      for (m=-l; m<=l; m++) {
        termi	=	sfharmonics2(l, m, costheta, cosphi, sinphi) * alpha;
	termj	=	sfharmonics2(l, m, -costheta, -cosphi, -sinphi) * alpha;
 
//	termi	=	sfharmonics(l, m, costheta, phi) * alpha;
//	termj	=	sfharmonics(l, m, -costheta, pi+phi) * alpha;	//reverse the vector direction
        part[i].ylmalphasum[m+l]	-=	termi;		//contribution to qlm(i)
	part[j].ylmalphasum[m+l]	-=	termj;		//contribution to qlm(j)
	YlmAlphaSum[m+l]	-=	(termi + termj);	//contribution to Qlm
      }
    }
   
    dp.x	=	part[j].p.x	-	p_new.x;	//new i, j distance
    dp.y	=	part[j].p.y	-	p_new.y;
    dp.z	=	part[j].p.z	-	p_new.z;
    MapInBox(&dp);
    r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;

    if (r2 < Rb*Rb) {
      //alpha	=	1.0;
      alpha	=	(r2 + Rb*Rb - 2*sqrt(r2)*Rb) / (sigma*sigma + Rb*Rb - 2*sigma*Rb);

      part[i].alphasum	+=	alpha;		//substract contribution to qlm(i)
      part[j].alphasum	+=	alpha;		//substract contribution to qlm(j)
      AlphaSum		+=	2*alpha;	//substract contribution to Qlm

      costheta	=	dp.z/sqrt(r2);
//      phi	=	atan2( dp.y, dp.x ); 
      rxy	=	sqrt(dp.x*dp.x + dp.y*dp.y);
      cosphi	=	dp.x / rxy;
      sinphi	=	dp.y / rxy;
      if (rxy<ZERO) {
         cosphi	=	0;
         sinphi	=	0;
      }
      for (m=-l; m<=l; m++) {
	termi	=	sfharmonics2(l, m, costheta, cosphi, sinphi) * alpha;
	termj	=	sfharmonics2(l, m, -costheta, -cosphi, -sinphi) * alpha;
//	termi	=	sfharmonics(l, m, costheta, phi) * alpha;
//	termj	=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
        part[i].ylmalphasum[m+l]	+=	termi;		
	part[j].ylmalphasum[m+l]	+=	termj;	
	YlmAlphaSum[m+l]	+=	(termi + termj);
      }
    }

    if (fabs(part[j].alphasum) < 1e-14) {		//update qlm(j) 
      if (part[j].alphasum < 0)
         printf("alpha_j = %10.8f\n", part[j].alphasum);
      for (m=-l; m<=l; m++) {			
        part[j].qlm[m+l]	=	0;
      }
      part[j].ql	=	0;		//update ql(j)
    }
    else {
      part[j].ql	=	0;		//don't forget initialization
      for (m=-l; m<=l; m++) {
        part[j].qlm[m+l]	=	part[j].ylmalphasum[m+l] / part[j].alphasum;
	part[j].ql		+=	part[j].qlm[m+l] * conj(part[j].qlm[m+l]);
      }
      part[j].ql	=	sqrt(part[j].ql * 4 * pi/(2*l+1));
    }
  }    
  }}

  if (fabs(part[i].alphasum) < 1e-14) {		//now it is time to update qlm(i)
      if (part[i].alphasum < 0)
         printf("alpha_i = %10.8f\n", part[i].alphasum);
    for (m=-l; m<=l; m++) {
      part[i].qlm[m+l]	=	0;
    }
    part[i].ql	=	0;
  }
  else {
    part[i].ql	=	0;
    for (m=-l; m<=l; m++) {
      part[i].qlm[m+l]	=	part[i].ylmalphasum[m+l] / part[i].alphasum;
      part[i].ql	+=	part[i].qlm[m+l] * conj(part[i].qlm[m+l]);
    }
    part[i].ql	=	sqrt(part[i].ql * 4 * pi/(2*l+1));
  }	

  Ql	=	0;			//update Qlm and Ql
  if (fabs(AlphaSum) < 1e-14) {
      if (AlphaSum < 0)
         printf("Alpha = %10.8f\n", AlphaSum);
    for (m=-l; m<=l; m++) {
      Qlm[m+l]	=	0;
    }
    Ql	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
      Qlm[m+l]	=	YlmAlphaSum[m+l] / AlphaSum;
      Ql	+=	Qlm[m+l] * conj(Qlm[m+l]);
    }
  }
  Ql	=	sqrt(Ql * 4 * pi/ (2*l+1));

  return;

/*
  //this alternative doesn't rely on the vlistplus, but actually due to the fact that
  //oldmol.vlist is changed when vlist is updated (because we transfer only pointer
  //rather than the whole set of array from part[n] to oldmol.

  for (jj=0; jj<oldmol.nverlet; jj++) {		//search the old vlist of particle i, 
						//substract their contribution 
    j	=	oldmol.vlist[jj];

    dp.x	=	part[j].p.x	-	oldmol.p.x;	//old i, j distance
    dp.y	=	part[j].p.y	-	oldmol.p.y;
    dp.z	=	part[j].p.z	-	oldmol.p.z;
    MapInBox(&dp);
    r2		=	dp.x*dp.x + dp.y*dp.y + dp.z*dp.z;
    
    if (r2 < BONDCUTSQ) {
      alpha	=	(r2 + BONDCUTSQ - 2*sqrt(r2*BONDCUTSQ)) / 
			(sigma*sigma + BONDCUTSQ - 2*sigma*BONDCUT);
      part[i].alphasum	-=	alpha;		//substract contribution to qlm(i)
      part[j].alphasum	-=	alpha;		//substract contribution to qlm(j)
      AlphaSum		-=	2*alpha;	//substract contribution to Qlm

      costheta	=	dp.z/sqrt(r2);
      phi	=	atan2( dp.y, dp.x ); 
      for (m=-l; m<=l; m++) {
	termi	=	sfharmonics(l, m, costheta, phi) * alpha;
	termj	=	sfharmonics(l, m, -costheta, pi+phi) * alpha;	//note the vector direction change
        part[i].ylmalphasum[m+l]	-=	termi;			//contribution to qlm(i)
	part[j].ylmalphasum[m+l]	-=	termj;			//contribution to qlm(j)
	YlmAlphaSum[m+l]		-=	(termi + termj);	//contribution to Qlm
      }
    }
  }

  for (jj=0; jj<part[i].nverlet; jj++) {	//search the new vlist of particle i, 
						//plus their contribution
    j	=	part[i].vlist[jj];

    dp.x	=	part[j].p.x	-	part[i].p.x;
    dp.y	=	part[j].p.y	-	part[i].p.y;
    dp.z	=	part[j].p.z	-	part[i].p.z;
    MapInBox(&dp);
    r2		=	dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;

    if (r2 < BONDCUTSQ) {
      alpha	=	(r2 + BONDCUTSQ - 2*sqrt(r2*BONDCUTSQ)) /
			(sigma * sigma + BONDCUTSQ - 2*sigma*BONDCUT);
      part[i].alphasum	+=	alpha;
      part[j].alphasum	+=	alpha;
      AlphaSum		+=	2*alpha;

      costheta	=	dp.z/sqrt(r2);
      phi	=	atan2( dp.y, dp.x);
      for (m=-l; m<=l; m++) {
	termi	=	sfharmonics(l, m, costheta, phi) * alpha;
	termj	=	sfharmonics(l, m, -costheta, pi+phi) * alpha;
	part[i].ylmalphasum[m+l]	+=	termi;
	part[j].ylmalphasum[m+l]	+=	termj;
	YlmAlphaSum[m+l]		+=	(termi + termj);
      }
    }
  }	

  for (jj=0; jj<oldmol.nverlet; jj++) {		//update qlm(j) in the old vlist of i
  						//do this after the update of new vlist 
    j	=	oldmol.vlist[jj];		//because some j's in the old vlist are 
    if (part[j].alphasum < ZERO) {		//also in the new vlist
      for (m=-l; m<=l; m++)
        part[j].qlm[m+l]	=	0;
    }
    else {
      for (m=-l; m<=l; m++) 
        part[j].qlm[m+l]	=	part[j].ylmalphasum[m+l] / part[j].alphasum;
    }
  }
  for (jj=0; jj<part[i].nverlet; jj++) {	//update qlm(j) in the new vlist of i
   
    j	=	part[i].vlist[jj];
    if (part[j].alphasum < ZERO) {
      for (m=-l; m<=l; m++)
        part[j].qlm[m+l]	=	0;
    }
    else {
      for (m=-l; m<=l; m++) 
        part[j].qlm[m+l]	=	part[j].ylmalphasum[m+l] / part[j].alphasum;
    }
  }

  if (part[i].alphasum < ZERO) {		//now it is time to update qlm(i)
    for (m=-l; m<=l; m++)
      part[i].qlm[m+l]	=	0;
  }
  else {
    for (m=-l; m<=l; m++)
      part[i].qlm[m+l]	=	part[i].ylmalphasum[m+l] / part[i].alphasum;
  }	

  Ql	=	0;				//update Qlm and Ql
  if (AlphaSum < ZERO) {
    for (m=-l; m<=l; m++) {
      Qlm[m+l]	=	0;
    }
    Ql	=	0;
  }
  else {
    for (m=-l; m<=l; m++) {
      Qlm[m+l]	=	YlmAlphaSum[m+l] / AlphaSum;
      Ql	+=	Qlm[m+l] * conj(Qlm[m+l]);
    }
  }
  Ql	=	sqrt(Ql * 4 * pi/ (2*l+1));

  return;
*/
}


int Qlbinfinder(double Ql)	//determine which Ql bin the sampled Ql value belongs to
{
   return	(int) ((Ql - Qlmiddle + Qlbins * Qlbinsize/2.0) / Qlbinsize); 
}

int NMAXbinfinder(int MAXSIZE)
{
   return	(int) ((MAXSIZE - NMAXmiddle + NMAXbins * NMAXbinsize/2.0) / NMAXbinsize);	
}

/*
double etaQl(double Ql)		//given etaQ[bins], find right eta using extrapolation and intrapolation
{
   int		i, bin;
   double 	eta;
   double	Qlmidpt[Qlbins];	//middle point of each Ql bin

   for (i=0; i<Qlbins; i++) {
      Qlmidpt[i]	=	(i+0.5-Qlbins/2.0) * Qlbinsize + Qlmiddle;
   }

   bin	=	Qlbinfinder(Ql);

   if ( fabs(Ql-Qlmidpt[bin]) < ZERO) {
      eta	=	etaQ[bin];
   }
   else if (Ql > Qlmidpt[bin]) {
      if (bin == Qlbins-1) {		//the rightmost bin, use extrapolation
	 eta	=	(Ql - Qlmidpt[bin-1])/Qlbinsize * (etaQ[bin] - etaQ[bin-1]) + etaQ[bin-1];
      }
      else {				//not the rightmost bin, use intrapolation
	 eta	=	(Ql - Qlmidpt[bin])/Qlbinsize * (etaQ[bin+1]-etaQ[bin]) + etaQ[bin];
      }
   }
   else {
      if (bin == 0) {			//the leftmost bin, use extrapolation
         eta	=	(Ql - Qlmidpt[bin+1])/Qlbinsize * (etaQ[bin+1]-etaQ[bin]) + etaQ[bin+1];
      }
      else {
	 eta	=	(Ql - Qlmidpt[bin-1])/Qlbinsize * (etaQ[bin]-etaQ[bin-1]) + etaQ[bin-1];
      }
   }
   return	eta;
}
*/

double etaQl(double Ql) 	//given Ql, use a fixed quadratic eta-Ql relation to calculate eta
{
   return	-0.5 * kQ * (Ql-Qlmiddle) * (Ql-Qlmiddle);	
}

double etaNMAX(long MAXSIZE)
{
   return	-0.5 * kN * ((double)MAXSIZE - NMAXmiddle) * ((double)MAXSIZE - NMAXmiddle);
}

double etaP2(double p2)
{
   return	-0.5 * kP2 * (p2 - P2middle) * (p2 - P2middle);
}

/*
double etaNMAX(int MAXSIZE)
{
   int		i, bin;
   double	etaN;
   double	NMAXmidpt[NMAXbins];

   for (i=0; i<NMAXbins; i++) {
      NMAXmidpt[i]	=	(i+0.5-NMAXbins/2.0) * NMAXbinsize + NMAXmiddle;
   }
   bin	=	NMAXbinfinder(MAXSIZE);   
   if ( fabs(MAXSIZE-NMAXmidpt[bin]) < ZERO) {
      etaN	=	eta[bin];
   }
   else if ( MAXSIZE > NMAXmidpt[bin]) {
      if (bin==NMAXbins-1) {
	 etaN	=	(MAXSIZE - NMAXmidpt[bin-1])/NMAXbinsize * (eta[bin] - eta[bin-1]) + eta[bin-1];	
      }
      else {
	 etaN	=	(MAXSIZE - NMAXmidpt[bin])/NMAXbinsize * (eta[bin+1] - eta[bin]) + eta[bin];	
      }
   }
   else {
      if (bin==0) {
	 etaN	=	(MAXSIZE - NMAXmidpt[bin+1])/NMAXbinsize * (eta[bin+1] - eta[bin]) + eta[bin+1];
      }
      else {
	 etaN	=	(MAXSIZE - NMAXmidpt[bin-1])/NMAXbinsize * (eta[bin] - eta[bin-1]) + eta[bin-1];
      }
   }
   return	etaN; 
}
*/
#endif	/* TEST */

void Calc_Qlm(long L)
{
   long		m;	// m=[-L, L];
   long		i, j, k, n, system;
   molstruct	*moli, *molj;
   vector	dp;
   double	r2, rxy, alpha, costheta, cosphi, sinphi,
  		AlphaSum[MAXNSYSTEMS];			// has to be MAXNSYSTEMS because NSYSTEMS unspecified
   double complex	YlmAlphaSum[MAXNSYSTEMS][2L+1], Qlm[MAXNSYSTEMS][2L+1];

#ifdef CELL_LIST
   cellstruct	*celli, *cellj, *celln;
#endif

   /* Initialization */

   for (system=0; system<NSYSTEMS; system++) {
      for (m=-L; m<=L; m++) {
         YlmAlphaSum[system][m+L]	=	0.0;
	 Qlm[system][m+L]		=	0.0;
      }
      AlphaSum[system]	=	0.0;
      Q6[system]	=	0.0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         for (m=-L; m<=L; m++) {
	    moli->qlm[i][m+L]		=	0.0;
	    moli->ylmalphasum[i][m+L]	=	0.0;
         }
         moli->alphasum[i]	=	0.0;
	 moli->q6[i]		=	0.0;
      }
   }
   /* Calculation */

   for (moli=mol; moli<mol+NMOLS; moli++) {
      system	=	moli->box;

      for (i=0; i<moli->nsites; i++) {

#ifdef CELL_LIST
	 celli	=	moli->cell[i];			// cell id of moli->i
	 for (n=0; n<celli->nneigh; n++) {
	    celln	=	celli->neigh[n];	// neighboring cell id

	    for (k=0; k<celln->nsites; k++) {		// sites in neighboring cell
	       if (molj=celln->mol[k]) {
		  j	=	celln->molsite[k];	// pick up one

                  if (moli!=molj || i!=j) {
#else
	 for (molj=mol; molj<mol+NMOLS; molj++) {	// note: molj starts from moli
	    if (molj->box == system) {
               for (j=0; j<molj->nsites; j++) {
                  if (molj>moli || j>i) {
#endif
                     dp	=	V_Subtr(molj->p+j, moli->p+i);
                     dp	=	MapInBox2(&dp, PBC, system);
		     r2	=	dp.x * dp.x + dp.y * dp.y + dp.z * dp.z;

		     if (r2 < Rb2) {
		        //alpha		=	1.0;
		        alpha	=	(r2 + Rb2 - 2*sqrt(r2)*Rb)/(1 + Rb2-2*Rb);
		        moli->alphasum[i]	+=	alpha;
		        molj->alphasum[j]	+=	alpha;

		        costheta	=	dp.z/sqrt(r2);
		        rxy	=	sqrt(dp.x*dp.x + dp.y*dp.y);
		        cosphi	=	(rxy<ZERO ? 0 : dp.x/rxy);
		        sinphi	=	(rxy<ZERO ? 0 : dp.y/rxy);
			
//			printf("costheta = %f cosphi = %f sinphi = %f\n", costheta, cosphi, sinphi);
                        for (m=-L; m<=L; m++) {
			   moli->ylmalphasum[i][m+L]	+=	
					sfharmonics2(L, m, costheta, cosphi, sinphi)*alpha;
			   molj->ylmalphasum[j][m+L]	+=	
					sfharmonics2(L, m, -costheta, -cosphi, -sinphi)*alpha;
		        }
                     }
		  }
         }  }  }
      }
   }

   /* NO need to correct double counting here because both numerator and denominator are double-counted */
   // only cell list case will result in double counting, no cell list will avoid double counting, see code

   /* Normalization */

   for (moli=mol; moli<mol+NMOLS; moli++) {
      system	=	moli->box;

      for (i=0; i<moli->nsites; i++) {
         if (moli->alphasum[i] < 0)			// alphasum should always > 0
	    Exit("sample.c", "Calc_Qlm", "alphasum < 0");
         else if (moli->alphasum[i] > ZERO) {		// always exclude divided-by-zero
	    for (m=-L; m<=L; m++) {
	       moli->qlm[i][m+L]	=	moli->ylmalphasum[i][m+L] / moli->alphasum[i];
       
               if (6==L)
                  moli->q6[i]		+=	moli->qlm[i][m+L] * conj(moli->qlm[i][m+L]);
               else if (4==L)
                  moli->q4[i]		+=	moli->qlm[i][m+L] * conj(moli->qlm[i][m+L]);
            }
            if (6==L)
               moli->q6[i]		=	sqrt(moli->q6[i] * 4 * pi/(2*L+1));
            else if (4==L)
               moli->q4[i]		=	sqrt(moli->q4[i] * 4 * pi/(2*L+1));

            AlphaSum[system]	+=	moli->alphasum[i];
            for (m=-L; m<=L; m++)
               YlmAlphaSum[system][m+L]	+=	moli->ylmalphasum[i][m+L];	// ylmalphasum=0 if alphasum=0
         }
      }
   }

   for (system=0; system<NSYSTEMS; system++) {
      if (AlphaSum[system] < 0)
	 Exit("sample.c", "Calc_Qlm", "AlphaSum < 0");
      else if (AlphaSum[system] > ZERO) {
         for (m=-L; m<=L; m++) {
            Qlm[system][m+L]	=	YlmAlphaSum[system][m+L] / AlphaSum[system];
	    
	    if (6==L)
	       Q6[system]	+=	Qlm[system][m+L]  * conj(Qlm[system][m+L]);
	    else if (4==L)
	       Q4[system]	+=	Qlm[system][m+L]  * conj(Qlm[system][m+L]);
         }
         if (6==L)
            Q6[system]		=	sqrt(Q6[system] * 4 * pi/(2*L+1));
         else if (4==L)
            Q4[system]		=	sqrt(Q4[system] * 4 * pi/(2*L+1));
      }
   }
}


double qlproductSQ(int l, molstruct *moli, long i, molstruct *molj, long j)	//calculate q_l(i) dot q_l(j) SQ
{						
  int 		m;
  float complex	qlproduct=0;
  double	sumi=0, sumj=0, sum;

  for (m=-l; m<=l; m++) {
    qlproduct	+=	moli->qlm[i][m+l] * conj(molj->qlm[j][m+l]);
    sumi	+=	moli->qlm[i][m+l] * conj(moli->qlm[i][m+l]);
    sumj	+=	molj->qlm[j][m+l] * conj(molj->qlm[j][m+l]);
  }
  sum		=	sumi * sumj;

  if (fabs(sum)<ZERO) 
    return	0;
  else 
    return 	qlproduct * qlproduct / sum;
}


vector CenterofMass(molstruct *moli)		// center of mass of a chain
{
   long		i;
   double	mx, my, mz, m, M;
   vector	p;

   m	=	0.0;			// mass of one site
   mx	=	0.0;			// mass weighted x coordinate
   my	=	0.0;
   mz	=	0.0;
   M	=	0.0;			// total mass of one molecule

   for (i=0; i<moli->nsites; i++) {
      m		=	type[moli->type[i]].M;
      mx	+=	m * (moli->p[i].x);
      my	+=	m * (moli->p[i].y);
      mz	+=	m * (moli->p[i].z);
      M		+=	m;
	//printf("%f\t%f\t%f\t%f\t%f\n", m, mx, my, mz, M);
   }
   p.x	=	mx/M;
   p.y	=	my/M;
   p.z	=	mz/M;

   return	p;
}

vector groupCoM(beadstruct *group, long nsites)	// June 6,2009
{
   long		i, site; 
   double	mx, my, mz, m, M;
   vector	p;
   molstruct	*moli;

   m		=	0.0;
   mx = my = mz =	0.0;
   M		=	0.0;

   for (i=0; i<nsites; i++) {
      moli	=	group[i].moli;
      site	=	group[i].site;
      m		=	type[moli->type[site]].M;
      mx	+=	m * (moli->p[i].x);
      my	+=	m * (moli->p[i].y);
      mz	+=	m * (moli->p[i].z);
      M		+=	m;
   }
   p.x	=	mx/M;
   p.y	=	my/M;
   p.z	=	mz/M;
   return	p;
}


double R2_gyration(molstruct *molm)	// squared radius of gyration of one molecule
{
   long		i;
   vector	rmean;
   double	Rg2;

   rmean.x	=	0;
   rmean.y	=	0;
   rmean.z	=	0;
   
   for (i=0; i<molm->nsites; i++) {
      rmean.x	+=	molm->p[i].x;
      rmean.y	+=	molm->p[i].y;
      rmean.z	+=	molm->p[i].z;
   }
   rmean.x	/=	molm->nsites;
   rmean.y	/=	molm->nsites;
   rmean.z	/=	molm->nsites;

   Rg2	=	0;

   for (i=0; i<molm->nsites; i++) {
      Rg2	+=	(molm->p[i].x - rmean.x) * (molm->p[i].x - rmean.x);
      Rg2	+=	(molm->p[i].y - rmean.y) * (molm->p[i].y - rmean.y);
      Rg2	+=	(molm->p[i].z - rmean.z) * (molm->p[i].z - rmean.z);
   }
   Rg2	/=	molm->nsites;

   return	Rg2;
}


matrix GyrationTensor(molstruct *molm)		// Gyration tensor of one molecule
{
   long         i;
   matrix       tensor;
   double       x, y, z;
   vector       rmean;

   M_Null(&tensor);

   rmean.x	=	0;
   rmean.y	=	0;
   rmean.z	=	0;
   
   for (i=0; i<molm->nsites; i++) {
      rmean.x	+=	molm->p[i].x;
      rmean.y	+=	molm->p[i].y;
      rmean.z	+=	molm->p[i].z;
   }
   rmean.x	/=	molm->nsites;
   rmean.y	/=	molm->nsites;
   rmean.z	/=	molm->nsites;

   for (i=0; i<molm->nsites; i++) {
      x         =       molm->p[i].x;
      y         =       molm->p[i].y;
      z         =       molm->p[i].z;

      tensor.x.x        +=      x * x;
      tensor.x.y        +=      x * y;
      tensor.x.z        +=      x * z;
      tensor.y.y        +=      y * y;
      tensor.y.z        +=      y * z;
      tensor.z.z        +=      z * z;
   }
   tensor	=	M_Mult(1.0/molm->nsites, &tensor);

   tensor.x.x   -=      rmean.x * rmean.x; 	// substract the center of mass contribution
   tensor.x.y   -=      rmean.x * rmean.y;
   tensor.x.z   -=      rmean.x * rmean.z;
   tensor.y.y   -=      rmean.y * rmean.y;
   tensor.y.z   -=      rmean.y * rmean.z;
   tensor.z.z   -=      rmean.z * rmean.z;
   tensor.y.x   =       tensor.x.y;             // a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}


matrix groupGyraTensor(beadstruct *group, long nsites)		// Gyration tensor of 
								// a group of particles
{
   long         i, site;
   matrix       tensor;
   double       x, y, z;
   vector       rmean;
   molstruct	*molm;

   M_Null(&tensor);

   rmean.x	=	0;
   rmean.y	=	0;
   rmean.z	=	0;

   for (i=0; i<nsites; i++) {		// center of mass, NO mass weighted
      molm	=	group[i].moli;   
      site	=	group[i].site;

      rmean.x	+=	molm->p[site].x;
      rmean.y	+=	molm->p[site].y;
      rmean.z	+=	molm->p[site].z;
   }
   rmean.x	/=	nsites;
   rmean.y	/=	nsites;
   rmean.z	/=	nsites;

   for (i=0; i<nsites; i++) {
      molm	=	group[i].moli;
      site	=	group[i].site;

      x         =       molm->p[site].x;
      y         =       molm->p[site].y;
      z         =       molm->p[site].z;

      tensor.x.x        +=      x * x;
      tensor.x.y        +=      x * y;
      tensor.x.z        +=      x * z;
      tensor.y.y        +=      y * y;
      tensor.y.z        +=      y * z;
      tensor.z.z        +=      z * z;
   }
   tensor	=	M_Mult(1.0/nsites, &tensor);

   tensor.x.x   -=      rmean.x * rmean.x; 	// substract the center of mass contribution
   tensor.x.y   -=      rmean.x * rmean.y;
   tensor.x.z   -=      rmean.x * rmean.z;
   tensor.y.y   -=      rmean.y * rmean.y;
   tensor.y.z   -=      rmean.y * rmean.z;
   tensor.z.z   -=      rmean.z * rmean.z;
   tensor.y.x   =       tensor.x.y;             // a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}

vector fshape(vector *eigvalue)		// shape descriptors, given eigvalues of gyration tensor
{
   vector	fshape;
   double	a, b, c, temp;

   a	=	eigvalue->x;
   b	=	eigvalue->y;
   c	=	eigvalue->z;

   if (a>b) {				// arrange the ordering so that a < b < c
      temp=b;	b=a;	a=temp;
   }
   if (b>c) {
      temp=c; 	c=b; 	b=temp;
   }
   if (a>b) {
      temp=b; 	b=a; 	a=temp;
   }
   fshape.x	=	c - 0.5* (a+b);		// asphericity
   fshape.y	=	b - a;			// acylindricity
   fshape.z	=	((fshape.x * fshape.x) + 0.75*(fshape.y * fshape.y)) / ((a+b+c)*(a+b+c));
						// relative shape anisotropy between 0 and 1
   return	fshape;
}

matrix InertiaTensor(molstruct *molm)	// added 1/30/08, calculate moment of inertia tensor
{					// inertia tensor is a mass weighted Gyration tensor
   long         i;
   matrix       tensor;
   double       x, y, z, m, mass;
   vector       com;

   M_Null(&tensor);
   com  =       CenterofMass(molm);

   mass	=	0;
   for (i=0; i<molm->nsites; i++) {
      mass	+=	type[molm->type[i]].M;
   }

   for (i=0; i<molm->nsites; i++) {
      x         =       molm->p[i].x - com.x;
      y         =       molm->p[i].y - com.y;
      z         =       molm->p[i].z - com.z;
      m         =       type[molm->type[i]].M;

      tensor.x.x        +=       m * (y*y + z*z);
      tensor.x.y        +=      -m * (x*y);
      tensor.x.z        +=      -m * (x*z);
      tensor.y.y        +=       m * (x*x + z*z);
      tensor.y.z        +=      -m * (y*z);
      tensor.z.z        +=       m * (x*x + y*y);
   }
   tensor.y.x   =       tensor.x.y;             //it is a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}


matrix groupInerTensor(beadstruct *group, long nsites)	// added June 6,2009, 
							// calculate moment of inertia tensor
							// of a group of beads		
{
   long         i, site;
   matrix       tensor;
   double       x, y, z, m, mass;
   vector       com;
   molstruct	*molm;

   M_Null(&tensor);
   com  =       groupCoM(group, nsites);

   mass	=	0;
   for (i=0; i<nsites; i++) {
      molm	=	group[i].moli;
      site	=	group[i].site;
      mass	+=	type[molm->type[site]].M;
   }

   for (i=0; i<nsites; i++) {
      molm	=	group[i].moli;
      site	=	group[i].site;

      x         =       molm->p[site].x;
      y         =       molm->p[site].y;
      z         =       molm->p[site].z;
      m         =       type[molm->type[site]].M;

      tensor.x.x        +=      m * ((y-com.y)*(y-com.y) + (z-com.z)*(z-com.z));
      tensor.x.y        +=      -m * (x-com.x) * (y-com.y);
      tensor.x.z        +=      -m * (x-com.x) * (z-com.z);
      tensor.y.y        +=      m * ((x-com.x)*(x-com.x) + (z-com.z)*(z-com.z));
      tensor.y.z        +=      -m * (y-com.y) * (z-com.z);
      tensor.z.z        +=      m * ((x-com.x)*(x-com.x) + (y-com.y)*(y-com.y));
   }
   tensor.y.x   =       tensor.x.y;             //it is a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}

double	R2_n2n(molstruct *moli)		// end-to-end distance square of one chain
{
   long		i, ib=moli->box;
   vector	p, dp;

   V_Null(&p);
   dp	=	V_Subtr(moli->p, moli->p+(moli->nsites)-1);

   return	V_Dot(&dp, &dp);
}


vector	CenterofNucleus(long nuclid, molstruct *molm) 	// center of specified nucleus, which contains molm
{
   vector	com, 
		rA,			// center of one molecule in nucleus
		rO,			// center of nucleus 
		rBA, 			// rB-rA, B and A belong to same nucleus
		rOA;			// rO-rA
   molstruct	*moli;
   long		system=0, n;

   // Step 1: find one molecule A belonging to the nucleus as a reference

   rA	=	CenterofMass(molm);
   system	=	molm->box;

   // Step 2: calc. the shift of center of nucleus to this molecule

   n	=	0;					// # of molecules in this nucleus
   V_Null(&rBA);
   V_Null(&rOA);

   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (moli->box == system && moli->nuclid[0] == nuclid) {
	 n	++;
	 com	=	CenterofMass(moli);
	 rBA	=	V_Subtr(&com, &rA);
	 rBA	=	MapInBox2(&rBA, PBC, system);
	 rOA	=	V_Add(&rOA, &rBA);		// assuming all molecules of same weight
      }
   }
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);			// center of nucleus
   return	rO;					// rO might be out of central box, but it is fine
}


void Sample_Pressure()
{
   long		i, j, k, n, nsites[2];
   double	vol, temp, virial, rc3, rc9;
   double	pres, dpres;
   double	Sigma, Epsilon;

   for (i=0; i<NBOX; i++) {
      n		=	NSites[i];			// # of atoms (not chains)
      vol	=	BOX[i].vol;
      temp	=	BOX[i].temp;
      virial	=	vir[i].tot;

      pres	=	( n*temp - virial/3) / vol;

      if (V_LJLRC) {
         if (2==NTYPES) {
            nsites[0]	=	2 * NMols[i];
            nsites[1]	=	NSites[i] - nsites[0];
            
            for (j=0; j<NTYPES; j++) {
               for (k=0; k<NTYPES; k++) {
                  Sigma		=	type[j].mix[k].SIGMA;
                  Epsilon	=	type[j].mix[k].EPSILON;

                  rc3	=	BOX[i].rc/Sigma;
                  rc3	=	rc3 * rc3 * rc3;
                  rc9	=	rc3 * rc3 * rc3;

                  dpres	=	2.0/(9*rc9) - 1.0/(3*rc3);
		  dpres	*=	nsites[j] * nsites[k] / (vol*vol);
                  dpres	*=	16.0 * pi * Epsilon * Sigma * Sigma * Sigma;

                  pres		+=	dpres;
               }
            }
         }
         if (1==NTYPES) {
            rc3	=	BOX[i].rc/type[0].SIGMA;
            rc3	=	rc3 * rc3 * rc3;
            rc9	=	rc3 * rc3 * rc3;

            dpres	=	16.0 * pi * n * n / (vol*vol) * (2.0/(9 * rc9) - 1.0/(3 * rc3));
            dpres	*=	type[0].EPSILON * type[0].SIGMA * type[0].SIGMA * type[0].SIGMA;
            pres	+=	dpres;
         }
      }
      BOX[i].pres	=	pres;
   }
}


void SampleDrift()
{
   molstruct	*moli;
   long		system;
   vector	p1, dp;

   for (system=0; system<NSYSTEMS; system++)
      drift2[system]	=	0.0;

   for (moli=mol; moli<mol+NMOLS; moli++) 
      if ( (system=moli->box)>= 0) {
         p1		=	CenterofMass(moli);
         dp		=	V_Subtr(&p1, &(moli->origin));
         drift2[system]	+=	V_Dot( &dp, &dp );
      }

   for (system=0; system<NSYSTEMS; system++)
      drift2[system]	/=	NMols[system];

/*
   if (D_DRIFT)				// temporary
      for (system=0; system<NSYSTEMS; system++) 
         for (moli=mol; moli<mol+NMOLS; moli++)
            if (moli->box == system) {
               drift2	=	V_Dot( &(moli->drift), &(moli->drift) );
               PutInDistribution(D_Drift+system, drift2, 1.0, 1.0);
            }
*/
}


void SampleSpherical()		// sample spherical coordinates distribution, 11/12/07
{				// given spherical coordinates updated
   molstruct	*moli;
   long		system, i, n[MAXNSYSTEMS];

   for (system=0; system<NSYSTEMS; system++) {
      n[system]		=	0;
      transfrac[system]	=	0.0;
   }

   AllSpherical();		// calculate spherical coordinates

   for (system=0; system<NSYSTEMS; system++)
      for (moli=mol; moli<mol+NMOLS; moli++)
         if (moli->box == system) {
            if (D_TORSION)
               for (i=3; i<moli->nsites; i++) {
                  //PutInDistribution(D_Torsion+system, AdjustAngle(moli->s[i].beta), 1.0, 1.0);
		  n[system]	++;
		  if (AdjustAngle(moli->s[i].beta) > 2.0*M_PI/3 || AdjustAngle(moli->s[i].beta) < -2.0*M_PI/3)
                     transfrac[system]	+=	1.0;
	       }
/*
            if (D_BONDA)
               for (i=2; i<moli->nsites; i++)
                  PutInDistribution(D_Bonda +system, AdjustAngle(moli->s[i].alpha), 1.0, 1.0);
            if (D_BONDL)
               for (i=1; i<moli->nsites; i++)
                  PutInDistribution(D_Bondl +system, moli->s[i].d, 1.0, 1.0);
*/
         }

   for (system=0; system<NSYSTEMS; system++)
      transfrac[system]	/=	n[system];
}


void Dist_Spherical()
{
   long		system, i;
   molstruct	*moli;

   for (system=0; system <NSYSTEMS; system++)
      for (moli=mol; moli<mol+NMOLS; moli++)
         if (moli->box == system) {
            if (D_TORSION)
               for (i=3; i<moli->nsites; i++) 
                  PutInDistribution(D_Torsion+system, AdjustAngle(moli->s[i].beta), 1.0, 1.0);
            if (D_BONDA)
               for (i=2; i<moli->nsites; i++)
                  PutInDistribution(D_Bonda +system, AdjustAngle(moli->s[i].alpha), 1.0, 1.0);
            if (D_BONDL)
               for (i=1; i<moli->nsites; i++)
                  PutInDistribution(D_Bondl +system, moli->s[i].d, 1.0, 1.0);
         }
}


void SampleRadial()
{
   long		i, j, system;
   molstruct	*moli, *molj;

   if (D_RADIAL)
      for (moli=mol; moli<mol+NMOLS; moli++)
         if ( (system=moli->box) >=0)
            for (molj=moli; molj<mol+NMOLS; molj++)
               if ( molj->box == system )
                  for (i=0; i<moli->nsites; i++)
                     for (j=0; j<molj->nsites; j++)
                        if ( moli==molj ? (j>i+2) : 1) {
                           PutInDistribution(D_Radial +system, sqrt(DistSQ(moli->p[i], molj->p[j], system)), 1.0, 1.0);
			}
}


void SampleP2()				// sample global orientational order
{
   molstruct	*moli, *molj;
   long		i, j, system, n[MAXNSYSTEMS], nz[MAXNSYSTEMS];
   vector	vi, vj;
   double	idotj;

   for (i=0; i<NSYSTEMS; i++) {		// global orientational order
      P2[i]	=	0.0;
      P2z[i]	=	0.0;
      n[i]	=	0;
      nz[i]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++)		// calculate P2
      if ( (system=moli->box)>=0 )
         for (i=1; i<moli->nsites-1; i++) {
            vi	=	V_Subtr(moli->p+i+1, moli->p+i-1);

            for (molj=moli; molj<mol+NMOLS; molj++)
               if ( molj->box == system )
                  for (j=1; j<molj->nsites-1; j++)
                     if ( (moli!=molj) ? 1 : j>i ) {
			vj		=	V_Subtr(molj->p+j+1, molj->p+j-1);			
                        idotj		=	V_Dot(&vi, &vj);
			idotj		*=	idotj;
			idotj		/=	V_Dot(&vi, &vi) * V_Dot(&vj, &vj);
                        P2[system]	+=	idotj;
			n[system]	++;
                     }
         }

   for (moli=mol; moli<mol+NMOLS; moli++)		// calculate P2z
      if ( (system=moli->box)>= 0)
         for (i=1; i<moli->nsites-1; i++) {
            vi		=	V_Subtr(moli->p+i+1, moli->p+i-1);
            P2z[system]	+=	(vi.z * vi.z) / V_Dot(&vi, &vi); 
            nz[system]	++;
         }
            
   for (i=0; i<NSYSTEMS; i++) {				// normalize
      P2[i]	/=	n[i];
      P2[i]	*=	1.5;
      P2[i]	-=	0.5;
      P2z[i]	/=	nz[i];
      P2z[i]	*=	1.5;
      P2z[i]	-=	0.5;
   }
}

//////////////////////////////////////////////////////////////////
/* Sample orientation order parameter, global P2 and local p2	*/
//////////////////////////////////////////////////////////////////
void SampleP2All()
{
   molstruct	*moli, *molj;
   long		i, j, system, n[MAXNSYSTEMS], nz[MAXNSYSTEMS];
   vector	vi, vj;
   double	idotj, r2;
#ifdef CELL_LIST
   long		k, l;
   static cellstruct	*celli, *cellk;
#endif

   // Initialization

   for (i=0; i<NSYSTEMS; i++) {		// global orientational order
      P2[i]	=	0.0;
      P2z[i]	=	0.0;
      P2M[i]	=	0.0;
      n[i]	=	0;
      nz[i]	=	0;
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {	// local orientational order
      for (i=0; i<moli->nsites; i++) {
         moli->p2[i]	=	0.0;
         moli->np2[i]	=	0;
      }
   }

   // Calculation

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (system=moli->box)>=0 )
         for (i=1; i<moli->nsites-1; i++) {			// exclude 2 end sites
            vi	=	V_Subtr(moli->p+i+1, moli->p+i-1);

#ifdef CELL_LIST
            celli	=	moli->cell[i];
            for (l=0; l<celli->nneigh; l++) {
               cellk	=	celli->neigh[l];
               for (k=0; k<cellk->nsites; k++) {
                  if (molj=cellk->mol[k]) {
		     j	=	cellk->molsite[k];
                     if (j!=0 && j!=molj->nsites-1 && (moli!=molj || i!=j)) {
#else
            for (molj=moli; molj<mol+NMOLS; molj++) {
               if ( molj->box == system ) {
                  for (j=1; j<molj->nsites-1; j++) {
                     if ( (moli!=molj) ? 1 : j>i ) {
#endif	/* CELL_LIST */

			vj		=	V_Subtr(molj->p+j+1, molj->p+j-1);			
                        idotj		=	V_Dot(&vi, &vj);
			idotj		*=	idotj;
			idotj		/=	V_Dot(&vi, &vi) * V_Dot(&vj, &vj);
                        P2[system]	+=	idotj;
			n[system]	++;
                        
                        if ((r2=DistSQ(moli->p[i], molj->p[j], system)) < Rp2*type[0].SIGMA*type[0].SIGMA) {	// calc. local p2
			   moli->p2[i]	+=	idotj;
			   molj->p2[j]	+=	idotj;
			   moli->np2[i]	++;
			   molj->np2[j]	++;
			}
            }  }  }  }
         }

   // Calc. the orientation along z-axis P2z
   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (system=moli->box)>= 0)
         for (i=1; i<moli->nsites-1; i++) {
            vi		=	V_Subtr(moli->p+i+1, moli->p+i-1);
            P2z[system]	+=	(vi.z * vi.z) / V_Dot(&vi, &vi); 
            nz[system]	++;
         }
            
   // Normalize P2 and P2z
   for (i=0; i<NSYSTEMS; i++) {
      if (n[i]) {
         P2[i]	/=	n[i];
         P2[i]	*=	1.5;
         P2[i]	-=	0.5;
      }
      else	P2[i]	=	0;
      if (nz[i]) {
         P2z[i]	/=	nz[i];
         P2z[i]	*=	1.5;
         P2z[i]	-=	0.5;
      }
      else	P2z[i]	=	0;
   }

   // calc. modified P2M using local p2(i)
   for (i=0; i<NSYSTEMS; i++) {
      n[i]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=1; i<moli->nsites-1; i++) {
         if (moli->np2[i]) {
            moli->p2[i]	/=	moli->np2[i];
            moli->p2[i]	*=	1.5;
	    moli->p2[i]	-=	0.5;
         }
	 else	moli->p2[i]	=	0;

         P2M[moli->box]	+=	moli->p2[i];
   	 n[moli->box]	++;
      }
      moli->p2[0]		=	moli->p2[1];
      moli->p2[moli->nsites-1]	=	moli->p2[moli->nsites-2];
   }
   for (i=0; i<NSYSTEMS; i++)
      P2M[i]	/=	n[i];
}

//////////////////////////////////////////////////////////
/* Sample number of connections based on chord vectors	*/
//////////////////////////////////////////////////////////
void SampleConnection()	
{
   molstruct	*moli, *molj;
   long		i, j, system, n[MAXNSYSTEMS], nz[MAXNSYSTEMS];
   vector	vi, vj;
   double	idotj, r2;
   static long		nconnect[1000];		// maximum 999 connections for one site
   static long		init=1;
#ifdef CELL_LIST
   static cellstruct	*cellm, *celli;
#endif
   if (init) {
      for (i=0; i<1000; i++)
          nconnect[i]	=	0;
      init	=	0;
   }

   for (i=0; i<NSYSTEMS; i++) {		// global orientational order
      P2[i]	=	0.0;
      P2z[i]	=	0.0;
      n[i]	=	0;
      nz[i]	=	0;
   }

   for (moli=mol; moli<mol+NMOLS; moli++)	// local orientational order
      for (i=0; i<moli->nsites; i++) {		// include end beads
         moli->p2[i]	=	0.0;
         moli->np2[i]	=	0;
      }

   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (system=moli->box)>=0 )
         for (i=0; i<moli->nsites; i++) {
            if (i==0)
               vi	=	V_Subtr(moli->p+2, moli->p);
            else if (i==moli->nsites-1)
               vi	=	V_Subtr(moli->p+moli->nsites-1, moli->p+moli->nsites-3);
            else
	       vi	=	V_Subtr(moli->p+i+1, moli->p+i-1);

            for (molj=moli; molj<mol+NMOLS; molj++)
               if ( molj->box == system )
                  for (j=0; j<molj->nsites; j++)
                     if ( (moli!=molj) ? 0 : j>i ) {
                        if (j==0)
                           vj	=	V_Subtr(molj->p+2, molj->p);
                        else if (j==molj->nsites-1)
                           vj	=	V_Subtr(molj->p+molj->nsites-1, molj->p+molj->nsites-3);
                        else
			   vj	=	V_Subtr(molj->p+j+1, molj->p+j-1);

                        idotj		=	V_Dot(&vi, &vj);
			idotj		*=	idotj;
			idotj		/=	V_Dot(&vi, &vi) * V_Dot(&vj, &vj);	//cos^2
			
                        if (idotj >= 0.933) {	// 0.933=cos^2 15degrees
//                        if (idotj >= 0.75) {	// 0.933=cos^2 30degrees
                        
                        if ((r2=DistSQ(moli->p[i], molj->p[j], system)) < Rconn2*type[0].SIGMA*type[0].SIGMA) {	// calc. distance between i and j 
			   moli->np2[i]	++;
			   molj->np2[j]	++;
			}

			}
                     }
         }
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         nconnect[moli->np2[i]]	++;

   for (i=0; i<200; i++)
      printf("%d  %d\n", i, nconnect[i]);
/*
   // calc. the orientation along z-axis
   for (moli=mol; moli<mol+NMOLS; moli++)
      if ( (system=moli->box)>= 0)
         for (i=1; i<moli->nsites-1; i++) {
            vi		=	V_Subtr(moli->p+i+1, moli->p+i-1);
            P2z[system]	+=	(vi.z * vi.z) / V_Dot(&vi, &vi); 
            nz[system]	++;
         }
            
   // normalize
   for (i=0; i<NSYSTEMS; i++) {
      if (n[i]) {
         P2[i]	/=	n[i];
         P2[i]	*=	1.5;
         P2[i]	-=	0.5;
      }
      else	P2[i]	=	0;
      if (nz[i]) {
         P2z[i]	/=	nz[i];
         P2z[i]	*=	1.5;
         P2z[i]	-=	0.5;
      }
      else	P2z[i]	=	0;
   }

   // calc. modified P2M using local p2(i)
   for (i=0; i<NSYSTEMS; i++) {
      P2M[i]	=	0.0;
      n[i]	=	0;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=1; i<moli->nsites-1; i++) {
         if (moli->np2[i]) {
            moli->p2[i]	/=	moli->np2[i];
            moli->p2[i]	*=	1.5;
	    moli->p2[i]	-=	0.5;
         }
	 else	moli->p2[i]	=	0;

         P2M[moli->box]	+=	moli->p2[i];
   	 n[moli->box]	++;
      }
      moli->p2[0]		=	moli->p2[1];
      moli->p2[moli->nsites-1]	=	moli->p2[moli->nsites-2];
   }
   for (i=0; i<NSYSTEMS; i++)
      P2M[i]	/=	n[i];
*/
}

//////////////////////////////////////////
/* Collect local p2 distribution	*/
//////////////////////////////////////////
void Dist_p2()		
{
   molstruct	*moli;
   long		i;

   if (D_LOCALP2)
      for (moli=mol; moli<mol+NMOLS; moli++)
         //for (i=1; i<moli->nsites-1; i++)
         for (i=0; i<moli->nsites; i++)
            PutInDistribution(D_Localp2+moli->box, moli->p2[i], 1.0, 1.0); 
}


//////////////////////////////////////////////////////////
/* Sample LC orientation matrix for all sites		*/
/* Principles of condensed matter physics, page 168	*/
//////////////////////////////////////////////////////////
void SampleM_Q()
{	
   molstruct		*moli, *molj, *molk;
   long			i, j, k, n, system, sitek;
   double		r2, ri;
   vector		evalue, evector, vi, vj;
   matrix		Q;
   static long		nconnect[1000];		// maximum 999 connections for one site
   static long		init=1;
   static double	cosangle;		// cosine of critical angle
#ifdef CELL_LIST
   cellstruct		*cellj, *cellk;
#endif

   if (init) {
      critangle	=	15.0/180*M_PI;
      cosangle	=	cos(critangle);
      for (i=0; i<200; i++)
         nconnect[i]	=	0;
      init	=	0;
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {
      //system=moli->box;
      system	=	0;				// as of Aug11,09, only one box
      for (i=1; i<moli->nsites-1; i++) {		// not for end sites

         moli->nconn[i]	=	0;

         for (molj=mol; molj<mol+NMOLS; molj++)		// list all chord vectors near i
            if ( molj->box == system )
               for (j=1; j<molj->nsites-1; j++)
                  if ((moli!=molj || i!=j) && (r2=DistSQ(moli->p[i], molj->p[j], system)) < Rp2*type[0].SIGMA*type[0].SIGMA) {

                     n		=	moli->nconn[i];
		     moli->connsite[i][n]	=	j;
		     moli->connmol[i][n]	=	molj;
		     moli->nconn[i]		++;
                  }

	 M_Null(&Q);

         vi	=	V_Subtr(moli->p+i+1, moli->p+i-1);	// self
         r2	=	V_Dot(&vi, &vi);
         vi	=	V_Mult(1.0/sqrt(r2), &vi);	// normalize
         Q.x.x	+=	vi.x * vi.x;
         Q.y.y	+=	vi.y * vi.y;
         Q.z.z	+=	vi.z * vi.z;
         Q.x.y	+=	vi.x * vi.y;
         Q.x.z	+=	vi.x * vi.z;
         Q.y.z	+=	vi.y * vi.z;

         for (k=0; k<moli->nconn[i]; k++) {		// neighbors
            molk	=	moli->connmol[i][k];
            sitek	=	moli->connsite[i][k];

            vi	=	V_Subtr(molk->p+sitek+1, molk->p+sitek-1);
            r2	=	V_Dot(&vi, &vi);
            vi	=	V_Mult(1.0/sqrt(r2), &vi);	// normalize
            Q.x.x	+=	vi.x * vi.x;
            Q.y.y	+=	vi.y * vi.y;
            Q.z.z	+=	vi.z * vi.z;
            Q.x.y	+=	vi.x * vi.y;
            Q.x.z	+=	vi.x * vi.z;
            Q.y.z	+=	vi.y * vi.z;
         }
         Q	=	M_Mult(1.0/(moli->nconn[i]+1), &Q);	// normalize
         Q.x.x	-=	0.333333;
         Q.y.y	-=	0.333333;
         Q.z.z	-=	0.333333;
         Q.y.x	=	Q.x.y;
         Q.z.x	=	Q.x.z;
	 Q.z.y	=	Q.y.z;

         evalue	=	M_eig(Q);				// find eigenvalue
         evector=	V_eig(Q, MAX(MAX(evalue.x, evalue.y),evalue.z));// find eigenvector
         r2	=	V_Dot(&evector, &evector);
         evector=	V_Mult(1.0/sqrt(r2), &evector);		// normalize eigenvector

         moli->vp2[i]	=	evector;		// assign value
//r2	=	V_Dot(&evector, &evector); 
//printf("%f %f %f %f ", evalue.x, evalue.y, evalue.z, evalue.x+evalue.y+evalue.z);
//printf("%f %f %f %f\n", evector.x, evector.y, evector.z, r2);
      }
      moli->vp2[0]	=	V_Add(moli->vp2+1, moli->vp2+2);	// end beads
      moli->vp2[0]	=	V_Mult(0.5, moli->vp2+0);

      moli->vp2[moli->nsites-1]	=	V_Add(moli->vp2+moli->nsites-2, moli->vp2+moli->nsites-3);
      moli->vp2[moli->nsites-1]	=	V_Mult(0.5, moli->vp2+moli->nsites-1);
   }

   // Now determine the number of connection for each site
   // Note, connectivity cutoff should be smaller than p2 cutoff !!!
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->np2[i]	=	0;			// NOTE: used to store connection #

         for (molj=mol; molj<mol+NMOLS; molj++) {
            for (j=0; j<molj->nsites; j++) {
               if (molj==moli && j==i)	break;
               r2	=	 DistSQ(moli->p[i], molj->p[j], moli->box);
               if (r2<Rconn2*type[0].SIGMA*type[0].SIGMA && fabs(V_Dot(moli->vp2+i, molj->vp2+j)) > cosangle)
                  moli->np2[i]	++;
/*
         for (k=0; k<moli->nconn[i]; k++) {		// search the neighbors
            molk	=	moli->connmol[i][k];
            sitek	=	moli->connsite[i][k];
	    r2		=	DistSQ(molk->p[sitek], moli->p[i], moli->box);
            if (r2 < Rconn2 && fabs(V_Dot(moli->vp2+i, molk->vp2+sitek)) > cosangle)
               moli->np2[i]	++;
         }            
*/	}}
         nconnect[moli->np2[i]]	++;
      }	    
   }
   printf("population distribution of number of connection\n");
   for (i=0; i<200; i++)	printf("%d\n", nconnect[i]);
}


//////////////////////////////////////////////////////////
/* sample the mean squared of end-to-end distance	*/
//////////////////////////////////////////////////////////
void SampleN2NSQ()	
{
   molstruct	*moli;
   double	R2 = 0.0;
   long		system;

   if (D_DENSITY)			// temporary
      for (system=0; system<NSYSTEMS; system++) 
         for (moli=mol; moli<mol+NMOLS; moli++)
            if (moli->box == system)
               PutInDistribution(D_Density+system, R2_n2n(moli), 1.0, 1.0);
}


//////////////////////////////////////////////////////////
/* sample energy					*/
//////////////////////////////////////////////////////////
void SampleEnergy()
{
   long		system;
   if (D_ENERGY)
      for (system=0; system<NSYSTEMS; system++)
         PutInDistribution(D_Energy+system, v[system].tot, 1.0, 1.0);
}

//////////////////////////////////////////////////////////
/* test function for 2D distribution			*/
//////////////////////////////////////////////////////////
void Sample2D()			
{
   long		system;
   double	dyn[2], y=1.0, weight=1.0;

   if (D_ENERGY)
      for (system=0; system<NSYSTEMS; system++) {
         dyn[0]	=	v[system].tot;
         dyn[1]	=	BOX[system].pres;
         D_Submit(D_Energy+system, dyn, &y, &weight);
      }
}

//////////////////////////////////////////////////////////
/* Print out distributions				*/
//////////////////////////////////////////////////////////
void S_PrintDistribution(diststruct **d)
{
   double		L = 1.0;				// change data units
   long			system, n = d-D_Distributions;

   if (*d) {
      for (system=0; system<NSYSTEMS; system++) {
         printf("%s\n", (*d+system)->header);
         D_Print(*d+system, 0, &L);
      }
   }
}

void S_PrintAll()
{
   long		i;
   for (i=0; i<D_NDIST; i++)
      S_PrintDistribution(D_Distributions+i);
}

void InitAllDistributions(diststruct **d, double binsize)	// add new distribution only to the end
{								// increase D_NDIST when adding
   long		i;
   double	dn[1], datom[1], denergy[1], ddensity[1], dpressure[1], dx[1];
   double	dangle[1];

   *dn		=	1.0;
   *datom	=	binsize;		// D_BINSIZE
/*   denergy[0]	=	2.000;
   denergy[1]	=	2.000;
*/
   *denergy	=	0.1;
   *ddensity	=	4.00;
   *dpressure	=	0.002;
   *dx		=	0.1;
   *dangle	=	M_PI / 180;

   for (i=0; i<D_NDIST; i++) {
      switch(i) {
         case 0: InitDist(D_DENSITY, d+i, "Density", 1, ddensity);	break;	// 1D distribution
	 case 1: InitDist(D_ENERGY, d+i, "Energy", 1, denergy);		break;	// 2D distribution
         case 2: InitDist(D_PRESSURE, d+i, "Pressure", 1, dpressure);	break;
         case 3: InitDist(D_DRIFT, d+i, "Drift", 1, dx);		break;
         case 4: InitDist(D_TORSION, d+i, "torsion angle", 1, dangle);	break;
         case 5: InitDist(D_BONDA, d+i, "bond angle", 1, dangle);	break;
         case 6: InitDist(D_BONDL, d+i, "bond length", 1, dpressure);	break;
         case 7: InitDist(D_RADIAL, d+i, "radial dist", 1, dx);		break;
         case 8: InitDist(D_LOCALP2, d+i, "local p2", 1, dpressure);	break;
         case 9: InitDist(D_XTALSIZE, d+i, "Xtal size", 1, dn);		break;
	 default:	break;
      }
   }
}

void LinkDistributions(diststruct **d)
{
   D_Density	=	*d++;			// i.e., D_Density = *d; d++;
   D_Energy	=	*d++;
   D_Pressure	=	*d++;
   D_Drift	=	*d++;
   D_Torsion	=	*d++;
   D_Bonda	=	*d++;
   D_Bondl	=	*d++;
   D_Radial	=	*d++;
   D_Localp2	=	*d++;
   D_Xtalsize	=	*d++;
}

void InitSample(char *argv[])
{
   static long		init = 1;
   char			name[256];
   long			n;
   double		binsize, d[1];

   if (init) {
      for (n=0; n<D_NDIST; n++) 
          D_Distributions[n]	=	NULL;
      init	=	0;
   }
   n		=	0;
   binsize	=	-1.0;

   if (binsize<=0.0) 
      binsize	=	D_BINSIZE;
   InitAllDistributions(D_Distributions, binsize);
   LinkDistributions(D_Distributions);
}

void ReinitSample(diststruct **d)	// only use when setting binsize afterwards
{
   InitAllDistributions(d, D_BINSIZE);
   if (d==D_Distributions)
      LinkDistributions(d);
}

void Sample_Histogram()
{
/*
   if (dynvar==1) {		//update the histogram
      p[MAXSIZE]	++;
      if (MAXSIZE == sizeright && mod(samplecycle, 2)==0)		samplecycle ++;
      else if (MAXSIZE == sizeleft && mod(samplecycle, 2)==1)	samplecycle ++;
   }
   else if (dynvar==2) {
      pQ[Qstatefinder(Ql)]	++;
      if (Qstatefinder(Ql) == Qrightid && mod(samplecycle, 2)==0)		samplecycle ++;
      else if (Qstatefinder(Ql) == Qleftid && mod(samplecycle, 2)==1)	samplecycle ++;
   }

   if (mod(counter*5, NCYCLE) == 0) {
      Print_Histogram();
   }
*/
}

vector Center_of_Mass(long *list, long number)		// center of mass of the molecules in a list
{						// need to be fixed to include different mass
   long         i, j;
   vector       p;

   V_Null(&p);

   for (i=0; i<number; i++) {
      for (j=0; j<Nsites; j++) {
         p.x       +=      mol[list[i]].p[j].x;
         p.y       +=      mol[list[i]].p[j].y;
         p.z       +=      mol[list[i]].p[j].z;
      }
   }
   V_Mult(1.0/(number*Nsites), &p);
   return       p;
}


matrix Gyration_Tensor(long *list, long number)
{
   long         i, j;
   matrix       tensor;
   double	x, y, z;
   vector	CoM;

   M_Null(&tensor);

   for (i=0; i<number; i++) {
      for (j=0; j<Nsites; j++) {

	 x	=	mol[list[i]].p[j].x;
	 y	=	mol[list[i]].p[j].y;
	 z	=	mol[list[i]].p[j].z;

         tensor.x.x        +=      x * x;
         tensor.x.y        +=      x * y;
         tensor.x.z        +=      x * z;
         tensor.y.y        +=      y * y;
         tensor.y.z        +=      y * z;
         tensor.z.z        +=      z * z;
      }
   }

   M_Mult(1.0/(number*Nsites), &tensor);

   CoM	=	Center_of_Mass(list, number);

   tensor.x.x	-=	CoM.x * CoM.x;		// substract the center of mass contribution
   tensor.x.y	-=	CoM.x * CoM.y;
   tensor.x.z	-=	CoM.x * CoM.z;
   tensor.y.y	-=	CoM.y * CoM.y;
   tensor.y.z	-=	CoM.y * CoM.z;
   tensor.z.z	-=	CoM.z * CoM.z;

   tensor.y.x   =       tensor.x.y;             //it is a symmetric matrix by definition
   tensor.z.x   =       tensor.x.z;
   tensor.z.y   =       tensor.y.z;

   return       tensor;
}


void Init_Sample()				// initialize sampling
{
   long		i;

   for (i=0; i<NBOX; i++) {
      chp[i]	=	0.0;
      cchp[i]	=	0.0;
      ichp[i]	=	0;
   }
   return;
}


void Sample_All()				// instant sampling during simulation
{
   if (V_VIRIAL)
      Sample_Pressure();
   SampleDrift();
   //SampleP2();
}


void Sample_Done()
{
   long		i;
/*
   for (i=0; i<NBOX; i++) {
      fprintf(foutput, "\n\tBOX[%d]:\n\n", i);
      fprintf(foutput, "\tChemical potential:\t%f\n", -log(cchp[i]/ichp[i]) * kT);
      fprintf(foutput, "\tTotal sampling:\t%d\n", ichp[i]);
   }
*/
   fflush(foutput);
   return;
}

/////////////////////////////////////////////////////////////////
/* Calculate radial distribution function and structure factor */
/////////////////////////////////////////////////////////////////

void radial(char *grswitch)		// calculate and print radial distribution function
{
   double	r, vb, nid, q;
   vector	rij, comi, comj, rij0;
   long		i, j, system, igr, i_q;
   molstruct	*moli, *molj;
   static long		init = 1, ngr=0;	// ngr: number of sampling
   static double	*binsize, *Lupper, 
			**gr, **grtot,			// g(r) of all beads
			**grcom, **grcomtot,		// g(r) of center of mass of chains
			**grcomxy, **grcomxytot,	// g(r) of com on z direction
			**grcomz, **grcomztot;		// g(r) of com on z direction
   static FILE		*frdf;

   // variables below for structure factor calculation, 06/25/09
   static double	**sq, **sqtot, dq, qmax;
   static long		nsq=0;			// nsq: # of sampling

   qmax	=	20.0;
   dq	=	qmax/NGRBINS;

   if (!strcmp(grswitch,"sample")) {		// do sampling

      // initialization
      if (init) {
         if (!(frdf=fopen("radial.dat", "w")))
            Exit("sample", "radial", "rdf output file failed to open");

         binsize	=	(double *) calloc(NSYSTEMS, sizeof(double));
         Lupper		=	(double *) calloc(NSYSTEMS, sizeof(double));

         gr		=	(double **) calloc(NSYSTEMS, sizeof(double *));
         grtot		=	(double **) calloc(NSYSTEMS, sizeof(double *));
         grcom		=	(double **) calloc(NSYSTEMS, sizeof(double *));
         grcomtot	=	(double **) calloc(NSYSTEMS, sizeof(double *));
         grcomz		=	(double **) calloc(NSYSTEMS, sizeof(double *));
         grcomztot	=	(double **) calloc(NSYSTEMS, sizeof(double *));
         grcomxy	=	(double **) calloc(NSYSTEMS, sizeof(double *));
         grcomxytot	=	(double **) calloc(NSYSTEMS, sizeof(double *));
 
         sq	=	(double **) calloc(NSYSTEMS, sizeof(double *));	// structure factor
         sqtot	=	(double **) calloc(NSYSTEMS, sizeof(double *));

         for (i=0; i<NSYSTEMS; i++) {
            gr[i]		=	(double *) calloc(NGRBINS, sizeof(double));
            grtot[i]		=	(double *) calloc(NGRBINS, sizeof(double));
            grcom[i]		=	(double *) calloc(NGRBINS, sizeof(double));
            grcomtot[i]		=	(double *) calloc(NGRBINS, sizeof(double));
            grcomz[i]		=	(double *) calloc(NGRBINS, sizeof(double));
            grcomztot[i]	=	(double *) calloc(NGRBINS, sizeof(double));
            grcomxy[i]		=	(double *) calloc(NGRBINS, sizeof(double));
            grcomxytot[i]	=	(double *) calloc(NGRBINS, sizeof(double));

            Lupper[i]	=	0.8*MIN(BOX[i].lx, MIN(BOX[i].ly, BOX[i].lz));	// once for all
	    binsize[i]	=	Lupper[i]/NGRBINS;				// once for all

	    sq[i]	=	(double *) calloc(NGRBINS, sizeof(double));
	    sqtot[i]	=	(double *) calloc(NGRBINS, sizeof(double));
         }
         init	=	0;
      }
      // Note: halfL and binsize are calculated once for all, thus we need to 
      //       make sure that the system is in equilibrium from the beginning!

      for (i=0; i<NSYSTEMS; i++) {
         for (j=0; j<NGRBINS; j++) {
            gr[i][j]		=	0.0;
	    grcom[i][j]		=	0.0;
            grcomz[i][j]	=	0.0;
            grcomxy[i][j]	=	0.0;
	 }
      }

      // calculate g(r) for all beads

      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box)>=0 ) {
            for (i=0; i<moli->nsites; i++) {

               for (molj=moli; molj<mol+NMOLS; molj++) {
                  if ( molj->box == system ) {
                     for (j=0; j<molj->nsites; j++) {
                        if ( (moli!=molj) ? 1 : j>=i+DLJ ) {
                        //if ( (moli!=molj) ? 1 : j>=i+1 ) {
	
                           rij0	=	V_Subtr(molj->p+j, moli->p+i);
			   rij	=	MapInBox2(&rij0, PBC, system);
			   r	=	sqrt(V_Dot(&rij, &rij));	// pair distance

			   if (r < Lupper[system]) {
			      igr	=	(int) (r/binsize[system]);
			      gr[system][igr]	+=	2;	// contribution from i and j
			   }
                        }
               }  }  }
      }  }  } 

      // calculate structure factor sq by doing FFT to rdf

      for (system=0; system<NSYSTEMS; system++) {
         for (i_q=0; i_q<NGRBINS; i_q++) {
            sq[system][i_q]	=	0.0;
  
            for (i=0; i<NGRBINS; i++) {
               if ((i+0.5)*binsize[system] >= 0.625*Lupper[system]) {
                  break;
               }			// upper limit of integral set to half box size
               sq[system][i_q]	+=	(gr[system][i]-1) * (i+0.5)*binsize[system]
					* sin((i_q+0.5)*dq*(i+0.5)*binsize[system])
					* binsize[system];
            }
            sq[system][i_q]	*=	4*M_PI*NSites[system]/BOX[system].vol/((i_q+0.5)*dq);
            sq[system][i_q]	+=	1.0;
            sqtot[system][i_q]	+=	sq[system][i_q];
         }
      }
      nsq	++;

      // calculate g(r) for center of mass of chains

      for (moli=mol; moli<mol+NMOLS-1; moli++) {
         if ( (system=moli->box)>=0 ) {
            comi	=	CenterofMass(moli); 

	    for (molj=moli+1; molj<mol+NMOLS; molj++) {
	       if (molj->box == system) {
		  comj	=	CenterofMass(molj);

		  rij0	=	V_Subtr(&comi, &comj);
		  rij	=	MapInBox2(&rij0, PBC, system);
		  r	=	sqrt(V_Dot(&rij, &rij));

		  if (r<Lupper[system]) {
	             igr	=	(int) (r/binsize[system]);
		     grcom[system][igr]	+=	2;

		     if (fabs(rij.z)*unit.LENGTH > 8.0)		// 8 angstrom apart in z direction
			grcomz[system][igr]	+=	2; 
		     else
			grcomxy[system][igr]	+=	2;
		  }
            }  }
      }  }

      // normalize g(r), do it every time because BOX.vol changes
      // In fact, it is not necessary as we assume the system in equilibrium

      for (i=0; i<NSYSTEMS; i++) {
         fprintf(frdf, "RDF in System %d (%d mols %d sites and vol = %f)\n", 
		system, NMols[system], NSites[system], BOX[system].vol);
         fprintf(frdf, "r \t gr \t grcom \t grcomz \t grcomxy\n");
         for (j=0; j<NGRBINS; j++) {
            vb	=	4.0*M_PI/3 
			* ((j+1)*(j+1)*(j+1)-j*j*j) * binsize[i] * binsize[i] * binsize[i];
            nid	=	vb * NSites[i]/BOX[i].vol; 
            gr[i][j]	/=	(NSites[i] * nid);
            grtot[i][j]	+=	gr[i][j];

	    nid =	vb * NMols[i]/BOX[i].vol;
            grcom[i][j]		/=	(NMols[i] * nid);
	    grcomz[i][j]	/=	(NMols[i] * nid);
	    grcomxy[i][j]	/=	(NMols[i] * nid);

	    grcomtot[i][j]	+=	grcom[i][j];
	    grcomztot[i][j]	+=	grcomz[i][j];
	    grcomxytot[i][j]	+=	grcomxy[i][j];

            fprintf(frdf, "%f %f %f %f %f\n", binsize[i]*(j+0.5)*unit.LENGTH, 
		gr[i][j], grcom[i][j], grcomz[i][j], grcomxy[i][j]);
         }
      }
      ngr	++;				// counter increase by 1
   }
   else if (!strcmp(grswitch, "print")) {
      for (i=0; i<NSYSTEMS; i++) {		// output rdf
         fprintf(frdf, "RDF in System %d, average of %d samples.\n", system, ngr);
         fprintf(frdf, "r \t grtot \t grcomtot \t grcomztot \t grcomxytot\n");
         for (j=0; j<NGRBINS; j++) {
            grtot[i][j]		/=	ngr;	// calculate ensemble average
            grcomtot[i][j]	/=	ngr;
	    grcomztot[i][j]	/=	ngr;
	    grcomxytot[i][j]	/=	ngr;

            fprintf(frdf, "%f %f %f %f %f\n", binsize[i]*(j+0.5)*unit.LENGTH, 
		grtot[i][j], grcomtot[i][j], grcomztot[i][j], grcomxytot[i][j]);
         }
      }
      for (i=0; i<NSYSTEMS; i++) {		// output structure factor
         fprintf(frdf, "Structure factor in System %d, average of %d samples.\n", system, nsq);
         fprintf(frdf, "q \t sqave\n");
         for (i_q=0; i_q<NGRBINS; i_q++) {
            fprintf(frdf, "%f \t %f\n", (i_q+0.5)*dq, sqtot[i][i_q]/nsq);
         }
      }
      fclose(frdf);
   }
   else {
      Exit("sample.c", "radial", "parameter not found");
   }
   return;
}


////////////////////////////////////////////////
/* Calculate structure factor from definition */
////////////////////////////////////////////////

void sq(FILE *fPtr, char *sqswitch)
{
   molstruct 		*moli;
   long			i, system, i_q, nq;		// nq: # of q points
   double		q, dq, qmax, argument;
   vector		p;

   static long		init=1, nsq=0; 			// nsq: # of samplings
   static double	*realx, *realy, *realz;
   static double	*imagx, *imagy, *imagz;
   static double	**sqx, **sqy, **sqz, **sq;	// sq: average of sqx,sqy,sqz

   nq	=	NGRBINS;
   qmax	=	20.0;		// C-C bond ~ 0.38sigma, q~2.6 sigma^-1
   dq	=	qmax/nq;

   if (!strcmp(sqswitch,"sample")) {		// do sampling
      if (init) {				// initialization
         sqx	=	(double **) calloc(NSYSTEMS, sizeof (double *));
         sqy	=	(double **) calloc(NSYSTEMS, sizeof (double *));
         sqz	=	(double **) calloc(NSYSTEMS, sizeof (double *));
         sq	=	(double **) calloc(NSYSTEMS, sizeof (double *));

         for (i=0; i<NSYSTEMS; i++) {
            sqx[i]	=	(double *) calloc(nq, sizeof(double));
            sqy[i]	=	(double *) calloc(nq, sizeof(double));
            sqz[i]	=	(double *) calloc(nq, sizeof(double));
            sq[i]	=	(double *) calloc(nq, sizeof(double));
         }

         realx	=	(double *) calloc (NSYSTEMS, sizeof(double));
         realy	=	(double *) calloc (NSYSTEMS, sizeof(double));
         realz	=	(double *) calloc (NSYSTEMS, sizeof(double));
         imagx	=	(double *) calloc (NSYSTEMS, sizeof(double));
         imagy	=	(double *) calloc (NSYSTEMS, sizeof(double));
         imagz	=	(double *) calloc (NSYSTEMS, sizeof(double));

         init	=	0;
      }

      // calculate structure factor sq

      for (i_q=0; i_q<nq; i_q++) {
         q	=	(i_q+0.5)*dq;		// center of the bin

         for (i=0; i<NSYSTEMS; i++) {
            realx[i]	=	0.0;	imagx[i]	=	0.0;
	    realy[i]	=	0.0;	imagy[i]	=	0.0;
	    realz[i]	=	0.0;	imagz[i]	=	0.0;	
         }

         for (moli=mol; moli<mol+NMOLS; moli++) {
            system	=	moli->box;
            for (i=0; i<moli->nsites; i++) {
               p	=	MapInBox2(moli->p+i, PBC, system);
               argument		=	q*p.x;
               realx[system]	+=	cos(argument);
               imagx[system]	+=	sin(argument);

               argument		=	q*p.y;
               realy[system]	+=	cos(argument);
               imagy[system]	+=	sin(argument);

               argument		=	q*p.z;
               realz[system]	+=	cos(argument);
               imagz[system]	+=	sin(argument);
            }
         }
         for (i=0; i<NSYSTEMS; i++) {
            sqx[i][i_q]		+=	(realx[i]*realx[i] + imagx[i]*imagx[i]);
            sqy[i][i_q]		+=	(realy[i]*realy[i] + imagy[i]*imagy[i]);
            sqz[i][i_q]		+=	(realz[i]*realz[i] + imagz[i]*imagz[i]);

/*            sqx[i][i_q]	=	(realx[i] * realx[i] + imagx[i] * imagx[i])/NSites[i];
            sqy[i][i_q]	=	(realy[i] * realy[i] + imagy[i] * imagy[i])/NSites[i];
            sqz[i][i_q]	=	(realz[i] * realz[i] + imagz[i] * imagz[i])/NSites[i];
            sq[i][i_q]	=	0.333333 * (sqx[i][i_q] + sqy[i][i_q] + sqz[i][i_q]);
*/
         }
      }
      nsq	++;
   }
   else if (!strcmp(sqswitch, "print")) {
      for (i_q=0; i_q<nq; i_q++) {
         for (i=0; i<NSYSTEMS; i++) {
            sqx[i][i_q]	/=	(NSites[i]*nsq);
            sqy[i][i_q]	/=	(NSites[i]*nsq);
            sqz[i][i_q]	/=	(NSites[i]*nsq);
            sq[i][i_q]	=	0.333333 * (sqx[i][i_q] + sqy[i][i_q] + sqz[i][i_q]);
         }
      }
      for (i=0; i<NSYSTEMS; i++) {
         fprintf(fPtr, "structure factor in system %d, average of %d samplings.\n", system, nsq);
         fprintf(fPtr, "q \t sqx \t sqy \t sqz \t sq\n");
         for (i_q=0; i_q<nq; i_q++) {
            fprintf(fPtr, "%f\t%f\t%f\t%f\t%f\n",
			 (i_q+0.5)*dq, sqx[i][i_q], sqy[i][i_q], sqz[i][i_q], sq[i][i_q]);
         }
      }
      fclose(fPtr);
   }
   return;
}
                                                                                                                                                                                                                                                                                       src/units.c                                                                                         0000600 0143352 0000144 00000011124 11046577643 012376  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	units.c
    author:	Peng Yi at MIT
    date:	Nov. 7, 2007
    purpose:	Units conversion from SI to program units
*/
#define __UNITS_MODULE
#include "units.h"

/*
			SI units used in input		Program units

	energy:		kJ/mol				epsilon[0]
	length:		Angstrom			sigma[0]
	temperature:	Kelvin				epsilon[0]
	pressure:	atm				epsilon[0]/sigma[0]^3
	angle:		degree				radian
	mass:		g/mol				1
	density:	g/cm^3				1/sigma[0]^3
*/

double CalcMass()		// calculate average mass of each site
{
   molstruct	*moli;
   long		i;
   double	avemass = 0.0;
/*
   for (moli=mol; moli<mol+NMOLS; moli++)
      for (i=0; i<moli->nsites; i++)
         avemass	+=	type[moli->type[i]].M;

   avemass	/=	NSITES;
   return	avemass;
*/

   // We don't use above calc. because we first need units to read in coord.
   // This is for now, because I don't want to assume that I know how many types are there

   if (1==NTYPES)
      return	type[0].M;

   if (NTYPES>=2) {
      avemass	=	2 * NMOLS * type[0].M + (NSITES - 2* NMOLS) * type[1].M;
      avemass	/=	NSITES;
      return	avemass;
   }
   return	0;
}


void CalcUnits(long flag)	// calc. the link b/w system units and SI units
{
   switch(flag){
      case 1: 
      	unit.MASS	=	CalcMass();	// ave mass of each site in g/mol
      	unit.LENGTH	=	type[0].SIGMA;
      	unit.ENERGY	=	type[0].EPSILON;
      	unit.TEMP	=	type[0].EPSILON / (N_KB * N_NAV * 1e-3);
      	unit.PRESSURE	=	type[0].EPSILON / pow(type[0].SIGMA, 3) / (N_ATM * N_NAV * 1e-33);
      	unit.ANGLE	=	180.0 / M_PI;
      	unit.DENSITY	=	unit.MASS / pow(unit.LENGTH, 3) / (N_NAV * 1e-24);
        break;

      case 0:
        unit.MASS	=	1.0;
  	unit.LENGTH	=	1.0;
      	unit.ENERGY	=	1.0;
      	unit.TEMP	=	1.0;
      	unit.PRESSURE	=	1.0;
      	unit.ANGLE	=	1.0;
      	unit.DENSITY	=	1.0;
        break;

      default:
        break;
   }
}


void PrintUnits()
{
   fprintf(foutput, "\nUnits:\n\n");
   fprintf(foutput, "\tMASS (g/mol):		\t%f\n", unit.MASS);
   fprintf(foutput, "\tLENGTH (Angstrom):	\t%f\n", unit.LENGTH);
   fprintf(foutput, "\tENERGY (kJ/mol): 	\t%f\n", unit.ENERGY);
   fprintf(foutput, "\tTEMP (K):		\t%f\n", unit.TEMP);
   fprintf(foutput, "\tPRESSURE (atm):		\t%f\n", unit.PRESSURE);
   fprintf(foutput, "\tANGLE (degree):		\t%f\n", unit.ANGLE);
   fprintf(foutput, "\tDENSITY (g/cm^3):	\t%f\n", unit.DENSITY);
}


void ConvertSIToSystem()		// convert read-in data to system units
{
   long		i, j;
  
   kT		/=	unit.TEMP;
   P		/=	unit.PRESSURE;
   Rho		/=	unit.DENSITY;

   for (i=0; i<NTYPES; i++) {
      type[i].M		/=	unit.MASS;
      type[i].SIGMA	/=	unit.LENGTH;
      type[i].EPSILON	/=	unit.ENERGY;
      type[i].LSTRETCH	/=	unit.LENGTH;
      type[i].KSTRETCH	/=	unit.ENERGY;	// need to be fixed/improved ...
      type[i].THETA	/=	unit.ANGLE;
      type[i].KBENDING	/=	unit.ENERGY;
      type[i].HS	/=	unit.LENGTH;

      for (j=0; j<6; j++)
         type[i].TORSION[j]	/=	unit.ENERGY;
   }
}


void ConvertSystemToSI()		// convert output data to SI units
{
   molstruct	*moli;
   long		i, system;

   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->p[i].x		*=	unit.LENGTH;
         moli->p[i].y		*=	unit.LENGTH;
         moli->p[i].z		*=	unit.LENGTH;
         moli->s[i].d		*=	unit.LENGTH;
         moli->s[i].alpha	*=	unit.ANGLE;
         moli->s[i].beta	*=	unit.ANGLE;
      }
   }
}


void CoorSI2System()		// box size and particle coordinates unit conversion
{
   molstruct	*moli;
   long		i;

   for (i=0; i<NSYSTEMS; i++) {
      BOX[i].lx		/=	unit.LENGTH;
      BOX[i].ly		/=	unit.LENGTH;
      BOX[i].lz		/=	unit.LENGTH;
      BOX[i].lbox	/=	unit.LENGTH;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->p[i].x	/=	unit.LENGTH;
         moli->p[i].y	/=	unit.LENGTH;
         moli->p[i].z	/=	unit.LENGTH;
         moli->s[i].d	/=	unit.LENGTH;
         moli->s[i].alpha	/=	unit.ANGLE;
         moli->s[i].beta	/=	unit.ANGLE;
      }
   }
}


void CoorSystem2SI()		// box size and particle coordinates unit conversion
{
   molstruct	*moli;
   long		i;

   for (i=0; i<NSYSTEMS; i++) {
      BOX[i].lx		*=	unit.LENGTH;
      BOX[i].ly		*=	unit.LENGTH;
      BOX[i].lz		*=	unit.LENGTH;
      BOX[i].lbox	*=	unit.LENGTH;
   }
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         moli->p[i].x	*=	unit.LENGTH;
         moli->p[i].y	*=	unit.LENGTH;
         moli->p[i].z	*=	unit.LENGTH;
         moli->s[i].d	*=	unit.LENGTH;
         moli->s[i].alpha	*=	unit.ANGLE;
         moli->s[i].beta	*=	unit.ANGLE;
      }
   }
}


void InitUnits()
{
   CalcUnits(ConvertUnits);
   ConvertSIToSystem();
}
                                                                                                                                                                                                                                                                                                                                                                                                                                            src/varbridge.c                                                                                     0000600 0143352 0000144 00000023676 10772306222 013204  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	varbridge.c
    author:	Peng Yi borrowed from Pieter J. in 't Veld
    date:	January 10, 2008
    purpose:	This module functions as a bridge between old style
    		system variables and new style system variables.
*/
#define __VARBRIDGE_MODULE
#include "varbridge.h"

void OrganizeMols()			// sort molecules according to their system id
{
  long			i, j, k = 0, nmols, nsites;
  molstruct		molt;
  
  for (i=0; i<NSYSTEMS; ++i)
  {
    nmols		= 0;
    nsites		= 0;
    for (j=0; j<NMOLS; ++j)
      if (mol[j].box==i)
      {
        molt		= mol[k];
	mol[k++]	= mol[j];
	mol[j]		= molt;
	nsites		+= mol[j].nsites;
        ++nmols;
      }
    NMols[i]		= nmols;
    NSites[i]		= nsites;
  }
}

bridgestruct		VB_BridgeMap;	// a collection of values and pointers pointing to system variables 
					// and molecules.  For binary history file output

bridgestruct *Bridge()
{
  return &VB_BridgeMap;
}


void BridgeSwap(long system1, long system2)
{
  molstruct		*mol_t;
  bridgestruct		*c = &VB_BridgeMap;

  mol_t			= c->mol[system1];
  c->mol[system1]	= c->mol[system2];
  c->mol[system2]	= mol_t;
}


bridgestruct *BridgeMap(long flag_mpi)
{
  long			i;
  bridgestruct		*c = &VB_BridgeMap;
  
  if (!flag_mpi)
    OrganizeMols();
    
  for (i=0; i<MAXNSYSTEMS; ++i)			// link to system and mol variables
  {
    c->system[i].n	= i;			// system.n is value, not pointer

    //GetBoxVectors(&(c->system[i].a), &(c->system[i].b), &(c->system[i].c));
    /*
    c->system[i].a.x	= BOX[i].lx;
    c->system[i].a.y	= 0;
    c->system[i].a.z	= 0;
    c->system[i].b.x	= 0;
    c->system[i].b.y	= BOX[i].ly;
    c->system[i].b.z	= 0;
    c->system[i].c.x	= 0;
    c->system[i].c.y	= 0;
    c->system[i].c.z	= BOX[i].lz;
    */
    c->system[i].lx	= &(BOX[i].lx);
    c->system[i].ly	= &(BOX[i].ly);
    c->system[i].lz	= &(BOX[i].lz);

    c->system[i].nmols	= NMols+i;		// pass pointer (address) rather than value
    c->system[i].nsites	= NSites+i;
    c->system[i].pres	= &(BOX[i].pres);
    c->system[i].vol	= &(BOX[i].vol);
    c->system[i].temp	= &(BOX[i].temp);
    c->system[i].drmax	= &(BOX[i].drmax);
    c->system[i].dlmax	= &(BOX[i].dlmax);
    c->system[i].damax	= &(BOX[i].damax);
    
    c->av[i]		= av+i;
    c->v[i]		= v+i;
    c->vir[i]		= vir+i;
    c->mol[i]		= flag_mpi ? mol : mol+i*NMols[i];
  } 
  for (i=0; i<D_NDIST; ++i)			// link to distributions
    c->dist[i]		= D_Distributions[i];

  //c->d_type		= &D_TYPE;
  c->d_type		= &NTYPES;
  c->type		= type;

  c->command.hs		= &V_HS;		// link to flags
  c->command.lj		= &V_LJ;
  c->command.ljshift	= &V_LJSHIFT;
  c->command.ljlrc	= &V_LJLRC;
  c->command.stretch	= &V_STRETCH;
  c->command.bending	= &V_BENDING;
  c->command.torsion	= &V_TORSION;
  c->command.virial	= &V_VIRIAL;
  c->command.scalecut	= &V_SCALECUTOFF;
  c->command.nvt	= &E_NVT;
  c->command.npt	= &E_NPT;
  c->command.gibbs	= &E_GIBBS;
  c->command.mpi	= &E_MPI;
  c->command.density	= &D_DENSITY;
  c->command.energy	= &D_ENERGY;
  c->command.pressure	= &D_PRESSURE;
  c->command.drift	= &D_DRIFT;
  c->command.torsional	= &D_TORSION;
  c->command.bonda	= &D_BONDA;
  c->command.bondl	= &D_BONDL;
  c->command.radial	= &D_RADIAL;
  c->command.localp2	= &D_LOCALP2;
  c->command.xtalsize	= &D_XTALSIZE;
  
/* 
  c->command.hs		= &P_HS;
  c->command.lj		= &P_LJ;
  c->command.stretch	= &P_STRETCH;
  c->command.ptorsion	= &P_TORSION;
  c->command.coulomb	= &P_COUL;
  c->command.virial	= &P_VIRIAL;
  c->command.polymer	= &E_POLYMER;
  c->command.mpi	= &E_MPI;
  c->command.async	= &S_ASYNC;
  c->command.monodisperse= &E_MONO;
  c->command.bias	= &E_BIAS;
  c->command.jacob	= &E_JACOB;
  c->command.nvt	= &E_NVT;
  c->command.npt	= &E_NPT;
  c->command.gibbs	= &E_GIBBS;
  c->command.insert	= &D_INSERT;
  c->command.widom	= &E_WIDOM;
  c->command.canonical	= &E_CANON;
  c->command.cavity	= &D_CAVITY;
  c->command.density	= &D_DENSITY;
  c->command.density3d	= &D_DENSITY3D;
  c->command.densfree	= &D_DENSITY_FREE;
  c->command.denstalobr	= &D_DENSITY_TALOBR;
  c->command.orient	= &D_ORIENTATION;
  c->command.orientfree = &D_ORIENT_FREE;
  c->command.orienttalobr= &D_ORIENT_TALOBR;
  c->command.orientcorr01 = &D_ORIENT_CORR_01;
  c->command.orientcorr02 = &D_ORIENT_CORR_02;
  c->command.orientcorr03 = &D_ORIENT_CORR_03;
  c->command.orientcorr04 = &D_ORIENT_CORR_04;
  c->command.orientcorr05 = &D_ORIENT_CORR_05;
  c->command.orientcorrCR = &D_ORIENT_CORR_CR;
  c->command.tails_etc	= &D_TAILS_ETC;
  c->command.e_n_function= &D_E_N_FUNCTION;
  c->command.radial	= &D_RADIAL;
  c->command.energy	= &D_ENERGY;
  c->command.torsion	= &D_TORSION;
  c->command.re_torsion	= &D_RE_TORSION;
  c->command.loopreentry= &D_LOOPREENTRY;
  c->command.b_length	= &D_B_LENGTH;
  c->command.b_angle	= &D_B_ANGLE;
  c->command.d_bridge	= &D_D_BRIDGE;
  c->command.e_profile	= &D_E_PROFILE;
  c->command.w_profile	= &D_W_PROFILE;
  c->command.n_profile	= &D_N_PROFILE;
  c->command.hs_dens	= &HS_DENS;
  c->command.temper	= &E_TEMPER;
  
  c->n.cycle		= &CYCLE;
  c->n.systems		= &NSYSTEMS;
  c->n.mols		= &NMOLS;
  c->n.sites		= &NSITES;
  c->n.types		= &NTYPES;
  c->n.volumes		= &NVOLUMES;
  c->n.swaps		= &NSWAPS;
  c->n.inserts		= &NINSERTS;
  c->n.cavinserts	= &NCAVINSERTS;
  c->n.blocks		= &NBLOCKS;
  c->n.cycles		= &NCYCLES;
  c->n.box		= &NBOX;
  c->n.rotation		= &NROTATION;
  c->n.reptation	= &NREPTATION;
  c->n.endbridge	= &NENDBRIDGE;
  c->n.rebridge		= &NREBRIDGE;
  c->n.bridges		= &NBRIDGES;
  c->n.fixed		= &NFIXED;
  c->n.semifixed	= &NSEMIFIXED;
  c->n.seed		= SEED;
  c->n.temper		= &NTEMPER;
  c->n.system		= &SYSTEM;
  c->n.sample		= &NSAMPLE;
  c->n.stretch		= &NSTRETCH;
  c->n.free		= &NFREE;
  c->n.ends		= &NENDS;
*/
  c->n.cycle		= &NCYCLE;		// link to numbers
  c->n.tape		= &ITAPE;
  c->n.systems		= &NSYSTEMS;
  c->n.mols		= &NMOLS;
  c->n.sites		= &NSITES;
  c->n.types		= &NTYPES;
  c->n.volchange	= &NVOLCHANGE;
  c->n.swap		= &NSWAP;
  c->n.displace		= &NDISPLACE;
  c->n.cbmc		= &NCBMC;

  return &VB_BridgeMap;
}

/*
#define VB_NLONG	2
#define VB_NDOUBLE	6
#define VB_NDIST	D_NDIST

long			**long_copy = NULL;
double			**double_copy = NULL;
diststruct		**dist_copy = NULL;
avstruct		*av_copy = NULL;
molstruct		*mol_copy = NULL;
vstruct			*v_copy = NULL;
wstruct			*w_copy = NULL;

void BridgeInitMemory(bridgestruct *bridge, long nsystems, long NMols)
{
  const char		module[32] = "varbridge",
			procedure[32] = "BridgeCreateCopy";
  long			i, nmols = 0;
  bridgestruct		*c = bridge;

  nmols			= nsystems*NMols;
  *c			= *Bridge();		// copy current bridge
  
  if (!long_copy)				// allocate copy space
  {
    if(!(long_copy = (long **) calloc(VB_NLONG, sizeof(long *))))
      Exit(module, procedure, "first long calloc error");
    for (i=0; i<VB_NLONG; ++i)
      if (!(long_copy[i] = (long *) calloc(nsystems, sizeof(long))))
        Exit(module, procedure, "second long calloc error");
  }
  if (!double_copy)
  {
    if(!(double_copy = (double **) calloc(VB_NDOUBLE, sizeof(double *))))
      Exit(module, procedure, "first double calloc error");
    for (i=0; i<VB_NDOUBLE; ++i)
      if (!(double_copy[i] = (double *) calloc(nsystems, sizeof(double))))
        Exit(module, procedure, "second double calloc error");
  }
  if (!dist_copy)
    if(!(dist_copy = (diststruct **) calloc(D_NDIST, sizeof(diststruct *))))
      Exit(module, procedure, "diststruct calloc error");
  if (!av_copy)
    if (!(av_copy = (avstruct *) calloc(nsystems, sizeof(avstruct))))
      Exit(module, procedure, "avstruct calloc error");
  if (!v_copy)
    if (!(v_copy = (vstruct *) calloc(nsystems, sizeof(vstruct))))
      Exit(module, procedure, "vstruct calloc error");
  if (!w_copy)
    if (!(w_copy = (wstruct *) calloc(nsystems, sizeof(wstruct))))
      Exit(module, procedure, "wstruct calloc error");
  if (!mol_copy)
    if (!(mol_copy = (molstruct *) calloc(nmols, sizeof(molstruct))))
      Exit(module, procedure, "molstruct calloc error");
  
  for (i=0; i<nsystems; ++i)			// set links to copy space
  {
    c->system[i].nmols	= long_copy[0]+i;
    c->system[i].nsites	= long_copy[1]+i;
    c->system[i].pres	= double_copy[0]+i;
    c->system[i].vol	= double_copy[1]+i;
    c->system[i].temp	= double_copy[2]+i;
    c->system[i].drmax	= double_copy[3]+i;
    c->system[i].dlmax	= double_copy[4]+i;
    c->system[i].damax	= double_copy[5]+i;

    c->av[i]		= av_copy+i;
    c->v[i]		= v_copy+i;
    c->vir[i]		= w_copy+i;
    c->mol[i]		= mol_copy+i*NMols;
  }
  ReinitSample(c->dist);
  for (i=0; i<D_NDIST; ++i)
    c->dist[i]		= dist_copy[i];
}


void BridgeResetVariable(bridgestruct *bridge)
{
  long			system;
  
  for (system=0; system<NSYSTEMS; ++system)
  {
    ResetAverages(bridge->av[system]);
    ResetAcceptance(bridge->av[system]);
  }
  ReinitSample(bridge->dist);
}


void BridgeCopyInstant(bridgestruct *to, bridgestruct *from, long system)
{
  long			i = system;
  
  *(to->system[i].nmols)	= *(from->system[i].nmols);
  *(to->system[i].nsites)	= *(from->system[i].nsites);
  *(to->system[i].pres)		= *(from->system[i].pres);
  *(to->system[i].vol)		= *(from->system[i].vol);
  *(to->system[i].temp)		= *(from->system[i].temp);
  *(to->system[i].drmax)	= *(from->system[i].drmax);
  *(to->system[i].dlmax)	= *(from->system[i].dlmax);
  *(to->system[i].damax)	= *(from->system[i].damax);
  to->system[i].a		= from->system[i].a;
  to->system[i].b		= from->system[i].b;
  to->system[i].c		= from->system[i].c;

  *(to->v[i])			= *(from->v[i]);
  *(to->vir[i])			= *(from->vir[i]);
  //to->av[i]->q_offset		= from->av[i]->q_offset;
  
  for (i=0; i<*(from->system[system].nmols); ++i)
    *(to->mol[system]+i)	= *(from->mol[system]+i);

  for (i=0; i<D_NDIST; ++i)
    if (from->dist[i]&&to->dist[i])
      for (system=0; system<*(from->n.systems); ++system)
        (to->dist[i]+system)->binsize = (from->dist[i]+system)->binsize;
}


void BridgeAddVariable(bridgestruct *to, bridgestruct *from)
{
  long			i, system;
  
  for (system=0; system<*(from->n.systems); ++system)
  {
    AddAvToAv(to->av[system], from->av[system]);
    for (i=0; i<D_NDIST; ++i)
      if (from->dist[i]&&to->dist[i])
        D_Add(to->dist[i]+system, from->dist[i]+system);
  }
}
*/
                                                                  src/vector.c                                                                                        0000600 0143352 0000144 00000016447 11213211077 012531  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	vector.c
    author:	Pieter J. in 't Veld for UT at Austin
    date:	April 1, 1999
    purpose:	3D vector algebra
*/
#define __VECTOR_MODULE
#include "vector.h"


void V_Null(vector *a)
{
  a->x			= 0.0;
  a->y			= 0.0;
  a->z			= 0.0;
}


void M_Null(matrix *a)
{
  V_Null(&(a->x));
  V_Null(&(a->y));
  V_Null(&(a->z));
}


void V_Unit(vector *a)
{
  a->x			= 1.0;
  a->y			= 1.0;
  a->z			= 1.0;
}


void V_Negate(vector *a)
{
  a->x			= -a->x;
  a->y			= -a->y;
  a->z			= -a->z;
}


void M_Unit(matrix *a)
{
  M_Null(a);
  a->x.x		= 1.0;
  a->y.y		= 1.0;
  a->z.z		= 1.0;
}


void M_Negate(matrix *a)
{
  a->x.x		= -a->x.x;
  a->x.y		= -a->x.y;
  a->x.z		= -a->x.z;
  a->y.x		= -a->y.x;
  a->y.y		= -a->y.y;
  a->y.z		= -a->y.z;
  a->z.x		= -a->z.x;
  a->z.y		= -a->z.y;
  a->z.z		= -a->z.z;
}


void V_PV(vector a)
{
  printf("{%g, %g, %g}", a.x, a.y, a.z);
//  fprintf(stderr, "{%g, %g, %g}", a.x, a.y, a.z);
}


void V_Print(vector a)
{
  V_PV(a);
  printf("\n");
//  fprintf(stderr, "\n");
}


void M_Print(matrix a)
{
  printf("{"); V_PV(a.x);
  printf(", "); V_PV(a.y);
  printf(", "); V_PV(a.z);
  printf("}\n");
/*
  fprintf(stderr, "{"); V_PV(a.x);
  fprintf(stderr, ", "); V_PV(a.y);
  fprintf(stderr, ", "); V_PV(a.z);
  fprintf(stderr, "}\n");
*/
}


double V_Dot(vector *a, vector *b)
{
  return a->x*b->x+a->y*b->y+a->z*b->z;
}

vector V_Add(vector *a, vector *b)
{
  vector		c;

  c.x			= a->x+b->x;
  c.y			= a->y+b->y;
  c.z			= a->z+b->z;
  return c;
}

vector V_Subtr(vector *a, vector *b)
{
  vector		c;

  c.x			= a->x-b->x;
  c.y			= a->y-b->y;
  c.z			= a->z-b->z;
  return c;
}

vector V_Cross(vector *a, vector *b)
{
  vector		c;

  c.x			= a->y*b->z - a->z*b->y;
  c.y			= a->z*b->x - a->x*b->z;
  c.z			= a->x*b->y - a->y*b->x;
  return c;
}


vector V_Mult(double f, vector *a)
{
  vector		b;

  b.x			= f*a->x;
  b.y			= f*a->y;
  b.z			= f*a->z;
  return b;
}


void V_Swap(vector *a, vector *b)
{
  vector		c = *a;

  *a			= *b;
  *b			= c;
}


vector MV_Dot(matrix *a, vector *b)
{
  vector		c;

  c.x			= a->x.x*b->x + a->y.x*b->y + a->z.x*b->z;
  c.y			= a->x.y*b->x + a->y.y*b->y + a->z.y*b->z;
  c.z			= a->x.z*b->x + a->y.z*b->y + a->z.z*b->z;
  return c;
}


matrix M_Rotation(double alpha, double beta)
{
  double		cosa = cos(alpha), cosb = cos(beta),
  			sina = sin(alpha), sinb = sin(beta);
  matrix		R;

  R.x.x			= cosa;
  R.x.y			= cosb*sina;
  R.x.z			= sinb*sina;
  R.y.x			= -sina;
  R.y.y			= cosb*cosa;
  R.y.z			= sinb*cosa;
  R.z.x			= 0.0;
  R.z.y			= -sinb;
  R.z.z			= cosb;
  return R;
}


vector V_Rotate(vector *a, double alpha, double beta)
{
  matrix		R = M_Rotation(alpha, beta);

  return MV_Dot(&R, a);
}


matrix M_Transpose(matrix *x)
{
  matrix		m;

  m			= *x;
  m.x.y			= x->y.x;
  m.x.z			= x->z.x;
  m.y.x			= x->x.y;
  m.y.z			= x->z.y;
  m.z.x			= x->x.z;
  m.z.y			= x->y.z;
  return m;
}


matrix M_Dot(matrix *x, matrix *y)
{
  matrix		m;

  m.x			= MV_Dot(x, &(y->x));
  m.y			= MV_Dot(x, &(y->y));
  m.z			= MV_Dot(x, &(y->z));
  return m;
}


double M_Det(matrix *x)
{
  vector		y = V_Cross(&(x->x), &(x->y));
  
  return V_Dot(&y, &(x->z));
}


int M_Inverse(matrix *a, matrix *b)
{
  double		f = M_Det(a);

  if (fabs(f)<1e-14) return V_E_SINGULAR;
  b->x			= V_Cross(&(a->y), &(a->z));
  b->y			= V_Cross(&(a->z), &(a->x));
  b->z			= V_Cross(&(a->x), &(a->y));
  b->x			= V_Mult(f = 1.0/f, &(b->x));
  b->y			= V_Mult(f, &(b->y));
  b->z			= V_Mult(f, &(b->z));
  *b			= M_Transpose(b);
  return 0;
}


matrix M_Mult(double f, matrix *x)
{
  matrix		y;

  y.x			= V_Mult(f, &(x->x));
  y.y			= V_Mult(f, &(x->y));
  y.z			= V_Mult(f, &(x->z));
  return y;
}


matrix M_Add(matrix *x, matrix *y)
{
  matrix		z;

  z.x			= V_Add(&(x->x), &(y->x));
  z.y			= V_Add(&(x->y), &(y->y));
  z.z			= V_Add(&(x->z), &(y->z));
  return z;
}


matrix M_Subtr(matrix *x, matrix *y)
{
  matrix		z;
  
  z.x			= V_Subtr(&(x->x), &(y->x));
  z.y			= V_Subtr(&(x->y), &(y->y));
  z.z			= V_Subtr(&(x->z), &(y->z));
  return z;
}


matrix M_XNormalRotate(vector *x)
{
  matrix		m;
  vector		r;
  double		f = V_Dot(x, x), 
  			sina, cosa, sinb, cosb;

  M_Unit(&m);
  if (f)
  {
    r			= V_Mult(1.0/sqrt(f), x);
    f			= (r.y<0) ? -sqrt(r.y*r.y+r.z*r.z) :
    				sqrt(r.y*r.y+r.z*r.z);
    if (f)
    {
      sina		= f;
      cosa		= r.x;
      sinb		= r.z/f;
      cosb		= -r.y/f;
      m.x.x		= cosa;
      m.x.y		= -sina*cosb;
      m.x.z		= sina*sinb;
      m.y.x		= sina;
      m.y.y		= cosa*cosb;
      m.y.z		= -cosa*sinb;
      m.z.x		= 0;
      m.z.y		= sinb;
      m.z.z		= cosb;
    }
  }
  return m;
}


// Determines the rotation matrix for rotating x to y

matrix M_Rotate(vector *x, vector *y)
{
  matrix		m1, m2;

  m1			= M_XNormalRotate(x);
  m2			= M_XNormalRotate(y);
  m1			= M_Transpose(&m1);
  m1			= M_Dot(&m2, &m1);
  return m1;
}


// Determines the orientation of a plane made up out of x and y

matrix M_Orientation(vector *x, vector *y)
{
  matrix		m;

  m.x			= V_Mult(1.0/sqrt(V_Dot(y, y)), y);
  m.z			= V_Cross(x, y);
  m.z			= V_Mult(1.0/sqrt(V_Dot(&(m.z),&(m.z))), &(m.z));
  m.y			= V_Cross(&(m.z), &(m.x));
  return m;
}


vector M_eig(matrix M)	// only for 3x3 REAL matrix, from Numerical Recipes
{			// note: real matrix can still have complex eigenvalues
   double	a, b, c, temp;
   double	Q, R, theta;
   double complex	A, B;
   vector	eigenvalue;

   a	=	-(M.x.x + M.y.y + M.z.z);
   b	=	M.x.x * M.y.y + M.x.x * M.z.z + M.y.y * M.z.z - M.z.y * M.y.z - M.x.y * M.y.x - M.x.z * M.z.x;
   c	=	-1.0* M.x.x * (M.y.y * M.z.z - M.y.z * M.z.y) 
		+ M.x.y * (M.y.x * M.z.z - M.y.z * M.z.x) 
		- M.x.z * (M.y.x * M.z.y - M.z.x * M.y.y);

   Q	=	(a*a - 3*b) /9;
   R	=	(2*a*a*a - 9*a*b + 27*c) / 54;

   if (R*R < Q*Q*Q) {
      theta	=	acos(R/sqrt(Q*Q*Q));
      eigenvalue.x	=	-2.0 * sqrt(Q) * cos(theta/3) - a/3;
      eigenvalue.y	=	-2.0 * sqrt(Q) * cos((theta+2*pi)/3) - a/3;	
      eigenvalue.z	=	-2.0 * sqrt(Q) * cos((theta-2*pi)/3) - a/3;

      //printf("%f\t%f\t%f\n", eigenvalue.x, eigenvalue.y, eigenvalue.z);
   }
   else {

      if (R<0)
         A	=	-cpow(R-sqrt(R*R-Q*Q*Q), 1.0/3);
      else
	 A	=	-cpow(R+sqrt(R*R-Q*Q*Q), 1.0/3);

      if (fabs(A)<ZERO)
         B	=	0;
      else
	 B	=	Q/A;

      //printf("%f\t%f\n", creal(A+B-a/3), cimag(A+B-a/3));
      //printf("%f\t%f\n",creal( -0.5*(A+B) - a/3 + 1.732050807569 * 0.5 * (A-B) * I ), cimag( -0.5*(A+B) - a/3 + 1.732050807569 * 0.5 * (A-B) * I ));
      //printf("%f\t%f\n",creal( -0.5*(A+B) - a/3 - 1.732050807569 * 0.5 * (A-B) * I ), cimag( -0.5*(A+B) - a/3 - 1.732050807569 * 0.5 * (A-B) * I ));

      eigenvalue.x	=	creal( A+B-a/3 );
      eigenvalue.y	=	creal( -0.5*(A+B) - a/3 + 1.732050807569 * 0.5 * (A-B) * I );
      eigenvalue.z	=	creal( -0.5*(A+B) - a/3 - 1.732050807569 * 0.5 * (A-B) * I );
   }
   return	eigenvalue;
}


vector V_eig(matrix M, double eigvalue) // the eigenvector with eigvalue
{
   vector       eigvector;
   double       xx, xy, xz, yx, yy, yz, zx, zy, zz;

   xx   =       M.x.x;          xy      =       M.x.y;          xz      =       M.x.z;
   yx   =       M.y.x;          yy      =       M.y.y;          yz      =       M.y.z;
   zx   =       M.z.x;          zy      =       M.z.y;          zz      =       M.z.z;

   eigvector.x  =       1.0;
   eigvector.y  =       -(zx * xz - (zz-eigvalue)*(xx-eigvalue))/(zy*xz-(zz-eigvalue)*xy);
   eigvector.z  =       -((xx-eigvalue)/xz + xy/xz*eigvector.y);

   return       eigvector;
}
                                                                                                                                                                                                                         src/correlation.h                                                                                   0000600 0143352 0000144 00000012256 11454720423 013556  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    correlation.h
    author:     Peng Yi at MIT
    date:       March 29, 2008
    purpose:    Correlation function calculation, only used in history file analysis
		so no need to compile these functions separately, just directly include
		into history file analysis source file hst.c
*/

#ifndef __CORR_HEADER
#define __CORR_HEADER

#define Lcorr		512			// maximum corr. separation b/w two records 

typedef struct{
   long		t;				// the index in the record file
   vector	com[MAXNMOLS];			// center of mass of each molecule
   vector	end2end[MAXNMOLS];		// end-to-end vector of each molecule
   vector	pos[MAXNMOLS*MAXNMOLSITES];	// position of each bead
   double	test;
} corr_data;					// data whose time correlation need to be calculated


typedef struct{
   double	com;				// mean square displacement of com of each chain
   double	end2end;			// end-to-end vector correlation
   double	pos;				// mean square displacement of each bead
   double	test;				// test
   long		norm;				// normalization factor
} corr_results;					// correlation function


corr_data	store[Lcorr+1];
corr_results	corr[Lcorr+1];


void corr_store(long nstore, long nlabel)	// do some calc. on each record and put in store
{						// nstore is the index in store and nlabel is the index in record
   long		i, j;
   molstruct 	*moli;
   vector	r;
   double	d;

   store[nstore].t	=	nlabel;

   store[nstore].test	=	(double)MAXSIZE[0];
   for (i=0; i<NMOLS; i++) {
      moli	=	mol+i;

      //////////

      store[nstore].com[i]		=	CenterofMass(moli);

      /////////

      if (moli->flip) {
         r	=	V_Subtr(moli->p, moli->p+(moli->nsites-1));
      }	
      else {
         r	=	V_Subtr(moli->p+(moli->nsites-1), moli->p);
      }
      d 	=	sqrt(V_Dot(&r, &r));
      store[nstore].end2end[i]	=	V_Mult(1.0/d, &r);

      /////////

      for (j=0; j<moli->nsites; j++) {
         store[nstore].pos[i*(NSITES/NMOLS)+j]	=	moli->p[j];
      }
   }
}


void corr_calc(long n1, long n2, long dt)	// calculate time correlation b/w two records, n1<=n2
{						// n1 and n2 are positions in store, dt is the real separation
   static long	init = 1;
   long		i, j;
   vector	dr;

   if (init) {					// initialization
      for (i=0; i<=Lcorr; i++) {
         corr[i].com		=	0;
         corr[i].end2end	=	0;
         corr[i].pos		=	0;
         corr[i].test		=	0;
         corr[i].norm		=	0;
      }
      init	=	0;
   }

   for (i=0; i<NMOLS; i++) {
      dr		=	V_Subtr(&(store[n1].com[i]), &(store[n2].com[i]));
      corr[dt].com	+=	V_Dot(&(dr), &(dr));

      corr[dt].end2end	+=	V_Dot(&(store[n1].end2end[i]),&(store[n2].end2end[i]));

      for (j=0; j<NSITES/NMOLS; j++) {
         dr		=	V_Subtr(&(store[n1].pos[i*NSITES/NMOLS+j]), 
					&(store[n2].pos[i*NSITES/NMOLS+j]));
         corr[dt].pos	+=	V_Dot(&(dr), &(dr));
      }
   }
   corr[dt].test	+=	store[n1].test * store[n2].test;
   corr[dt].norm	++;
}


void corr_norm()				// normalization
{
   long		i;

   corr[0].com		/=	corr[0].norm * NMOLS;		// m.s.d. = 0 if corr. length = 0
   corr[0].end2end	/=	corr[0].norm * NMOLS;
   corr[0].pos		/=	corr[0].norm * NSITES;
   corr[0].test		/=	corr[0].norm;

   for (i=1; i<=Lcorr; i++) {
      corr[i].com	/=	corr[i].norm * NMOLS;
      corr[i].end2end	/=	corr[i].norm * NMOLS * corr[0].end2end;
      corr[i].pos	/=	corr[i].norm * NSITES;
      corr[i].test	/=	corr[i].norm * corr[0].test;
   }

   //corr[0].com		=	1.0;	// m.s.d. = 0 if corr. length = 0
   corr[0].end2end	=	1.0;		// end2end correlation normalized to 1 at corr. length=0
   corr[0].test		=	1.0;
}


void corr_print()				// print out time correlation function
{
   long		i;

   printf("*****Time correlation function output*****\n");
   printf("Monte Carlo cycles between consecutive sampling = %d\n\n", ITAPE);

   printf("Column 1. Correlation time step.\n");
   printf("Column 2. Center of mass of chains.\n");
   printf("Column 3. Center of mass of bead.\n");
   printf("Column 4. End-to-end vector of chains.\n");
   printf("Column 5. Test.\n");

   for (i=0; i<=Lcorr; i++) {
      if (corr[i].norm > 0) {
         printf("%d\t", corr[i].norm);
         printf("%f\t", corr[i].com);
         printf("%f\t", corr[i].pos);
         printf("%f\t", corr[i].end2end);
         printf("%f\n", corr[i].test);
      }
   }
}


void correlation()
{
   long		i;
   static long	nstore = 0, 			// index in store
		nlabel = 0,			// index in record
		Ninit = 0; 			// index of the leftMOST record in corr calc.

   /* put record in store, whose length is max. correlation separation */

   corr_store(nstore, nlabel);			// record #nlabel is #nstore in store

   //printf("store[%d].test = %f\n", nstore, store[nstore].test);
   //printf("nstore = %d   nlabel = %d   Ninit = %d\n", nstore, nlabel, Ninit);

   /* calculate correlation */

   for (i=Ninit; i<=nlabel; i++) {
      //printf("mod(i, Lcorr+1) = %d  nstore = %d  nlabel-i = %d\n", mod(i, Lcorr+1), nstore, nlabel-i);
      corr_calc(mod(i, Lcorr+1), nstore, nlabel-i);	// corr. b/w current record and previous ones
   }

   nstore	++;				// store index increases by 1
   nlabel	++;				// record index increases by 1

   if (nlabel > Lcorr) {
      Ninit	++;
      if (mod(nlabel, Lcorr+1)==0)
         nstore	=	0;			// flush old store place for the new data
   }
}
#endif
                                                                                                                                                                                                                                                                                                                                                  src/distributions.h                                                                                 0000600 0143352 0000144 00000002123 10714730275 014133  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	distributions.h
    author:	Pieter J. in 't Veld for UT at Austin
    date:	June 12, 1998
    purpose:	header file for distributions.c
*/
#ifndef __DISTRIBUTIONS_HEADER
#define __DISTRIBUTIONS_HEADER

#include "header.h"

#define D_NAVERAGE	6		// do average to x, x^2, x^3 up to x^6

#ifndef __DISTRIBUTIONS_MODULE

extern void D_Allocate(diststruct *d, long n);
extern void D_Reset(diststruct *d);
extern void D_Copy(diststruct *dest, diststruct *src);
extern void D_Add(diststruct *dest, diststruct *src);
extern void D_Subtr(diststruct *dest, diststruct *src);
extern void D_Submit(diststruct *d, double *x, double *y, double *weight);
extern void D_Print(diststruct *d, int dist_type, double *L);
extern void D_PrintMath(diststruct *d, int dist_type, double *L);
extern void D_Init(
  diststruct **dist, char *header, long nlevels, double *binsize);

extern void PutInDistribution(diststruct *d, double x, double y, double weight);
extern void PrintDistribution(diststruct *d);
extern void InitDist(long flag, diststruct **dist, char *s, long nlevels, double *binsize);

#endif

#endif

                                                                                                                                                                                                                                                                                                                                                                                                                                             src/ensembles.h                                                                                     0000600 0143352 0000144 00000006514 11312052463 013205  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    ensembles.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Header file for ensembles.c.
*/
#ifndef __ENSEMBLES_HEADER
#define __ENSEMBLES_HEADER

#include "header.h"
#define	E_NMINFREE	3	

#ifdef __ENSEMBLES_MODULE


long		E_NVT,			// use NVT ensemble
		E_NPT,			// use NPT ensemble
		E_GIBBS,		// use GIBBS ensemble
		E_MPI;			// use MPI
long		NDISPLACE,		// % of trial displacement in each MC cycle among NVT moves
		NREPTATION,		// % of trial reptation in MC cycle among NVT moves
		NENDROT,		// end-mer rotation
		NCBMC,			// % of trial CBMC moves in each MC cycle among NVT moves
		NTRIALCONF,		// # of trial conf. in cmbc
		NTRIALFIRSTBEAD,	// # of trial positions of first bead in cmbc
		NVOLCHANGE,		// # of trial volume change in each MC cycle
		NGIBBSVOL,		// % of trial Gibbs volume change in each MC cycle
		NSWAP,			// # of trial Gibbs swap moves in each MC cycle
		NENDBR, NREBR,		// # of trial end-bridging and rebridging
		NDB, NIDR,		// # of trial double bridging and intra-double-rebridging
		NFLIP;			// # of trial one mer flip move
long            acc_move, rjc_move;  	// # of accepted/rejected trial move
long		acc_vol, rjc_vol;	// # of accepted/rejected trial volume change
long		acc_seq, rjc_seq;	// # of accepted/rejected trial move sequence
long		acc_gibbsvol, rjc_gibbsvol;
long		acc_swap, rjc_swap;

long		cbmcsucc[MAXNMOLSITES];
avstruct	av[MAXNBOX], av_past[MAXNBOX];
double		SUCC_DISP, SUCC_VOL;	// target accepted ratio of different kinds of move

#else

extern long	E_NVT, E_NPT, E_GIBBS, E_MPI;
extern long	NDISPLACE, NVOLCHANGE, NGIBBSVOL, NSWAP, NCBMC, NREPTATION, NENDROT, 
		NENDBR, NREBR, NDB, NIDR, NFLIP, NTRIALCONF, NTRIALFIRSTBEAD;

extern long     acc_move, rjc_move;  
extern long	acc_vol, rjc_vol;
extern long	acc_seq, rjc_seq;
extern long	acc_gibbsvol, rjc_gibbsvol;
extern long	acc_swap, rjc_swap;
extern long	cbmcsucc[MAXNMOLSITES];
extern avstruct	av[MAXNBOX], av_past[MAXNBOX];
extern double	SUCC_DISP, SUCC_VOL;

/* Move preparations */

extern void	StoreSystem();
extern void	RestoreSystem();
extern void	StoreOneMol(long);
extern void	RestoreOneMol(long);
extern void	StoreMols();
extern void	RestoreMols();

extern void	ResetAcceptance();	// reset move step sizes
extern void	Adjust_Stepsize();	// adjust move step sizes

extern void	MolInBox(molstruct *molm);
extern void	MolInBox2(molstruct *molm);
extern long	SiteSelect(molstruct **);
extern double	bondl_g(double, double);	// generate a bond length according to stretching energy
extern void	bonda_tors(double, double, double *, double *);
extern double	bonda_g(double, double);
extern double	tors(double);
extern vector	tors_bonda(molstruct *molm, long site);		// generate a unit vector according to bending and torsional energy

/* Moves */

extern long	rotation(long n);	// rotation of n sites, bond length and angle unchanged
extern long	movemol();		// move a molecule with bond length fixed
extern long	mcmove();   		// single particle displacement
extern long	mccbmc();		// conf.-biased Monte Carlo move
extern long	reptation();
extern long	mcvol();    		//system volume change
extern void	Gibbsmove();
extern void	NextMove();
extern void	Cycle();		//a Monte Carlo cycle
extern void	InitEnsemble();

extern void	ParentCheck();

extern void	Update_Eta(int);	//update weighting factors
extern int	CheckConstraint(char *);

#endif

#endif
                                                                                                                                                                                    src/forcefield.h                                                                                    0000600 0143352 0000144 00000007175 11241056715 013343  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    forcefield.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Calculation of the potential of a system
*/
#ifndef __FORCEFIELD_HEADER
#define __FORCEFIELD_HEADER


#ifdef __FORCEFIELD_MODULE
#include "header.h"

long		V_LJ;			// LJ interaaction V = 4*epsilon*( (sigma/r)^12 - (sigma/r)^6 )
long		V_HS;			// hard-sphere
long		V_RPL;			// repulsive spheres V = epsilon * (sigma/r)^12
long		V_LJSHIFT;		// LJ forcefield shift at cutoff
long		V_LJLRC;		// LJ forcefield long range correction
long		V_STRETCH;		// bond stretching
long		V_BENDING;
long		V_TORSION;
long		FIXBONDL, FIXBONDA;	// flag of fixed bond length/angle
long		V_VIRIAL;		// calculate virial function
long		DLJ;			// intrachain min. distance to have LJ interaction (>=DLJ)
double		Rc, Rclow,		// LJ interaction cutoff, and smaller cutoff
		Rv,			// verlet list cutoff
		Rb, 			// neighbor cutoff
		Rp, 			// local p2 cutoff
		Rconn;			// connectivity cutoff
long		V_SCALECUTOFF;		// scale cutoff radii with box dimension
typestruct	type[MAXNTYPES];	// force coefficient
char		TORTYPE[80];		// torsion energy type

#else

extern long	V_LJ, V_HS, V_RPL;
extern long	V_LJSHIFT, V_LJLRC;
extern long	V_STRETCH, V_BENDING, V_TORSION;
extern long	FIXBONDL, FIXBONDA;
extern long	V_VIRIAL;
extern long	DLJ;

extern double	Rc, Rclow, Rb, Rv, Rp, Rconn;
extern long	V_SCALECUTOFF;
extern typestruct	type[MAXNTYPES];	// force coefficient
extern char	TORTYPE[80];

/* Basic operations */

extern void	vstructNull(vstruct *);
extern void	wstructNull(wstruct *);
extern void	vstructNegate(vstruct *);			// v =  -v
extern void	wstructNegate(wstruct *);			// virial =  -virial
extern vstruct	vstructSum(vstruct *, vstruct *);
extern wstruct	wstructSum(wstruct *, wstruct *);
extern void	Printvstruct(vstruct *);

/* LJ mixing */

extern void	CalcMixSigma(long typei, long typej);
extern void	CalcMixEpsilon(long typei, long typej);

/* Initialization */

extern void	InitForcefield();

/* Calculate interaction energy */

extern double	OPLS(double phi, double k1, double k2, double k3);
extern double	OPLS2(double cosphi, double k1, double k2, double k3);
extern double	VLJSite(molstruct *molm, long site, double *w);
extern double	VLJMol(molstruct *, double *w);
extern double	VHSSite(molstruct *, long);
extern double	VHSMol(molstruct *);
extern double	VStretchSite(molstruct *, long, double *w);
extern double	VStretchMol(molstruct *, double *w);
extern double	VBendingSite(molstruct *, long);
extern double	VBendingMol(molstruct *);
extern double	VTorsionSite(molstruct *, long);
extern double	VTorsionMol(molstruct *);

extern vstruct	CalcVSite(molstruct *, long, wstruct *w);	// calculate energy and virial for single site
extern vstruct  CalcVSiteinner(molstruct *, long, wstruct *w);
extern vstruct	CalcVSiteouter(molstruct *, long, wstruct *w);

extern void	CalcVLJ();
extern void	CalcVHS();
extern void	CalcVStretch();
extern void	CalcVBending();
extern void	CalcVTorsion();
extern void	CalcVCorr();
extern void	CalcVLJCorr();
extern void	CalcV();
extern void	CalcV_mcvol(double volscale);

extern void	EnergyCheck();

/* Energy differece by moving molecule */

extern double	VDeleteSites(molstruct *, long, long);
extern double	VAddSites(molstruct *, long, long);
//extern double	grow(char *, molstruct *, long, vstruct *, wstruct *);
extern double	grow(char *, molstruct *, long);
extern long	Select(double *, double);

/* Old stuff */

extern vstruct 	VMol(long); 
extern double	Vtailco(double, double);
extern void 	Vtotal();
extern double	VDeleteMol(long);
extern double	VAddMol(long);
extern vstruct	VTestMol(molstruct *, long);
extern double   VAddTestMol(molstruct *, long);

#endif

#endif
                                                                                                                                                                                                                                                                                                                                                                                                   src/globals.h                                                                                       0000600 0143352 0000144 00000007675 11232631363 012667  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    globals.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Centralization of all global variables used in multiple
                programs
*/
#ifndef __MY_GLOBALS
#define __MY_GLOBALS

//#include "input.h"
#include "types.h"

#ifdef __MAIN_PROGRAM   
                        /* globals.h is included in header.h which in turn 
			   will be included in main.c, so we want the following
			   variables to be global variables.  But if any 
			   subroutine include header.h (thus global.h), since 
			   in them __MAIN_PROGRAM is not defined, only main.c
			   will define __MAIN_PROGRAM, they will only see like
			   extern long NPARTS , etc.
			   in this way we gurantee that these global variables
			   will not be defined multiple times.
			*/

//env variables
char			moltype[80];			// molecule type: LJ, monochain, polychain
long			NPARTS, NBOX, NSYSTEMS, NTYPES,
			NMOLS, NSITES,			// total mols, total sites
			Nsites,				// sites in one mol
   			Nequil, Nprod, NCYCLE, TRIALRUN, PROD,
			ITAPE, ICONF, NBLOCKS, NGSAMPLE, IRADIAL;
long			SEQUENCE;			// how many sequences in each MC cycle
long			TIMESTEP;			// used in i/o
long			NMols[MAXNSYSTEMS], NSites[MAXNSYSTEMS];
long			Ntest;
double			Dpoly, *mupoly;			// polydispersity variables

double			kT, P,
	                Rho, LBOX0,
			LBOX, VOL,
	                DRMAX, DLMAX, GDLMAX, DAMAX;

molstruct		*mol, *oldmol;
molstruct        	*part;			// define NPARTS particles in the system
molstruct		*parttemp;

wstruct			VIRMol;			// Virial for single molecule
wstruct			VIRSite;
vstruct			v[MAXNBOX];		// potential energy of each box
wstruct			vir[MAXNBOX];		// virial of each box
double			chp[MAXNBOX];		// chemical potential
double			cchp[MAXNBOX];		// cumulative chemical potential
long			ichp[MAXNBOX];		// chemical potential sampling times
boxstruct		BOX[MAXNBOX];		// enviro. parameter of each box
long			N[MAXNBOX];		// particle # in each box

//Monte Carlo step variables
long			Stage;			//which Stage it is in (melting, quenching, etc)
long			counter;   		//MC cycle counter

//block average
//double			ave_G[Qlbins][NGSAMPLE*2];
//double			ave_sqG[Qlbins];

//expanded ensemble variables

/*
double			*p;			//probability distribution 
double			*eta;			//weighting factor
double			*pt;			//target probability distribution
double			*G;			//Gibbs free energy

double			*pQ;
double			*etaQ;
double			*ptQ;
double			*GQ;
*/

int			trial;			//trial run counter

//random number generator variables
int                     *tim;  	
long int                *seed; 

//time and date variables
time_t			curtime;
struct tm 		*loctime;		
char			buffer[32];

//others
int			S_STOP;

#else

extern char		moltype[80];
extern long		NPARTS, NBOX, NSYSTEMS, NTYPES, NMOLS, NSITES, Nsites;
extern long		NMols[MAXNSYSTEMS], NSites[MAXNSYSTEMS];
extern double		Dpoly, *mupoly;	
extern double		Rho, LBOX0;
extern double		kT, P;
extern long		Nequil, Nprod, NCYCLE, TRIALRUN, PROD,
			ITAPE, ICONF, NBLOCKS, NGSAMPLE, IRADIAL;
extern long		SEQUENCE;
extern long		TIMESTEP;
extern long		Ntest;

extern molstruct	*mol, *oldmol;
extern molstruct        *part;
extern molstruct		*parttemp;
extern double           LBOX, VOL;
extern double           DRMAX, DLMAX, GDLMAX, DAMAX;

extern wstruct		VIRMol, VIRSite;
extern vstruct		v[MAXNBOX];
extern wstruct		vir[MAXNBOX];
extern double		chp[MAXNBOX];
extern double		cchp[MAXNBOX];
extern long		ichp[MAXNBOX];
extern boxstruct	BOX[MAXNBOX];
extern long		N[MAXNBOX];

extern long		Stage;
extern long 		counter;	

//extern double		ave_G[Qlbins][NGSAMPLE*2];
//extern double		ave_sqG[Qlbins];

/*
extern double		*p;
extern double		*eta;
extern double		*pt;
extern double		*G;

extern double		*pQ;
extern double		*etaQ;
extern double		*ptQ;
extern double		*GQ;
*/

extern int		trial;

extern int              *tim;  
extern long int         *seed;        

extern time_t		curtime;		//date and time
extern struct tm	*loctime;		
extern char		buffer[32];

extern int		S_STOP;
#endif

#endif
                                                                   src/header.h                                                                                        0000600 0143352 0000144 00000001235 11173374326 012465  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    header.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Header file for all programs.
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include "input.h"
#include "types.h"
#include "globals.h"
#include "distributions.h"
#include "init.h"
#include "random.h"
#include "sample.h"
#include "position.h"
#include "forcefield.h"
#include "ensembles.h"
#include "history.h"
#include "varbridge.h"
#include "lists.h"
#include "vector.h"
#include "io.h"
#include "motion.h"
#include "units.h"
#include "rebridge.h"
#include "roots.h"
                                                                                                                                                                                                                                                                                                                                                                   src/history.h                                                                                       0000600 0143352 0000144 00000001464 11006365314 012732  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
  program:	history.h
  author:	Peng Yi borrowed from Pieter J. in 't Veld
  date:		January 10, 2008
  purpose:	i/o for binary history files
*/
#ifndef __HISTORY_HEADER
#define __HISTORY_HEADER

#define HIST_VERSION	25
#define HIST_IDENT	"HIST"

#include "header.h"

#ifdef __HISTORY_MODULE

#else

extern int StartNewHistoryFile(char *name, long flag_mpi);
extern int H_StoreCurrent(char *name);
extern int H_StoreBridge(char *name, bridgestruct *bridge);
extern int H_GetHeader(char *name, long *version);
extern FILE *H_GetFirst(char *name, long *version, long flag_mpi);
extern int H_GetNext(FILE *fp, long version);
extern void H_InitSystem(char *argv[]);
extern bridgestruct *H_Bridge();

extern FILE *fcreate_hist(char *);		// for hstcomb
extern FILE *fread_hist(char *, long *);	// for hstcomb
#endif

#endif

                                                                                                                                                                                                            src/init.h                                                                                          0000600 0143352 0000144 00000000675 11002166627 012201  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    init.h
    author:     Peng Yi at MIT
    date:       October 23, 2006
    purpose:    Suggested initialization sequence
*/

#ifndef __INIT_HEADER
#define __INIT_HEADER

#ifdef __INIT_MODULE

#include "header.h"

char		INITCONF[256];		//initial configuration

#else

extern char	INITCONF[256];

extern void	GetCoordinates(char * filename);
extern void	InitMols(long, long);
extern void 	InitAll(char *argv[]);

#endif

#endif
                                                                   src/inline.h                                                                                        0000600 0143352 0000144 00000002467 10714730275 012522  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    file:	inline.h
    author:	Peng Yi
    date:	June 12, 2007.
    purpose:	Collection of inline functions.
*/
#ifndef __INLINE_HEADER
#define __INLINE_HEADER

#include "header.h"


#ifndef INLINE
# define INLINE extern inline
#endif

INLINE double	DistSQ(vector p, vector q)
{

//double DistSQ(vector p, vector q)	//distance square between vector p and q
//{
   double	xij, yij, zij;
   double	r2, r2new;
   vector	pimg;

/*   if (PBC==1) {
      xij	=	MIN(fabs(p.x-q.x), LBOX-fabs(p.x-q.x));
      yij	=	MIN(fabs(p.y-q.y), LBOX-fabs(p.y-q.y));
      zij	=	MIN(fabs(p.z-q.z), LBOX-fabs(p.z-q.z));
      r2	=	xij * xij + yij * yij + zij * zij;
   }
*/
//   if (PBC==2) {			//truncated octahedron box PBC
      xij	=	MIN(fabs(p.x-q.x), LBOX-fabs(p.x-q.x));		//min distance among the cubic images
      yij	=	MIN(fabs(p.y-q.y), LBOX-fabs(p.y-q.y));
      zij	=	MIN(fabs(p.z-q.z), LBOX-fabs(p.z-q.z));
      r2	=	xij * xij + yij * yij + zij * zij;

      pimg.x	=	p.x + ((q.x-p.x) > 0 ? 0.5 : -0.5)*LBOX;	//look around the octahedron images
      pimg.y	=	p.y + ((q.y-p.y) > 0 ? 0.5 : -0.5)*LBOX;
      pimg.z	=	p.z + ((q.z-p.z) > 0 ? 0.5 : -0.5)*LBOX; 
      r2new	=	(pimg.x-q.x) * (pimg.x-q.x) + (pimg.y-q.y) * (pimg.y-q.y) + (pimg.z-q.z)*(pimg.z-q.z); 

      if (r2new < r2)
	 r2	=	r2new;     
//   }
   return	r2;
//}

}

#endif
                                                                                                                                                                                                         src/input.h                                                                                         0000600 0143352 0000144 00000003164 11566346104 012376  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    input.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    System setup
*/

/*
	The difference between input.h and globals.h is that input.h defines CONSTANTS 
	used in the program, while globals.h defined VARIABLES.
*/

#ifndef __INPUT_HEADER
#define __INPUT_HEADER

#define VERSION		"2009_04_30"

#define pi		M_PI		// 3.1415926535897932
#define ZERO 		1E-14		// floating point number of real zero
#define TRUE		1
#define FALSE		0	

#define epsl            1.0   		// Lennard Jones parameters
#define sigma           1.0		// basis for reduced unit

#define Rc2		Rc*Rc
#define Rb2		Rb*Rb
#define Rv2		Rv*Rv
#define Rp2		Rp*Rp
#define Rconn2		Rconn*Rconn
#define critqlproductSQ	(critqlproduct * critqlproduct)
#define Rc3		(Rc*Rc*Rc)

#define l_of_Ylm	6		// self-explained

#define LINKLIST	0		// use linked list for verlet list
//#define VERLET_LIST	1
//#define CELL_LIST	1
#define	MAXNCELLSITES	9100		// max atom # in each cell 
#define MAXVERLETNEIGH	190		// max # of verlet neighbors of each particle
#define MAXCONNEIGH	30		// max # of connected neighbors of each particle

#define MAXNBOX		10		// at most 10 boxes
#define MAXNSYSTEMS	10		// at most 10 systems

#define MAXNMOLS	1000		// max # of molecules in ALL boxes
#define MAXNMOLSITES	160		// max # of sites on one single chain molecule

#define MAXNTYPES	4		// max # of different forcefield coeff. sets

#define MAXNDIST	16		// max of distribution we can handle

#define UMBRELLA	0
//#define MPI		1

#define MAX(a,b) (((a)>(b))?a:b)	// greater of two
#define MIN(a,b) (((a)<(b))?a:b)	// less of two

#define DEBUG	1
//#define TEST	1

#endif

                                                                                                                                                                                                                                                                                                                                                                                                            src/io.h                                                                                            0000600 0143352 0000144 00000002447 11176204760 011647  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	io.h
    author:     Peng Yi at MIT
    date:       October 24, 2006
    purpose:    Integral module for input/output handling
*/
#ifndef __IO_HEADER
#define __IO_HEADER


#ifdef __IO_MODULE

#include "header.h"

FILE		*foutput, *fhst, *fdump;	// output file handler
char		file_hst[256];		// name of binary history file
long		frame;			// visualization output frame #

#else

extern FILE	*foutput, *fhst, *fdump;
char		file_hst[256];
extern long	frame;

extern int	samestr(char *, char *);	// compare two strings, case insensitive
extern void	GetLVar(FILE *, long, long *);
extern void	GetDVar(FILE *, long, double *);
extern void	GetSetup(char * argv[]);	// read in setup info.
extern void 	InitFile();
extern void 	CloseFile();

extern void	Exit(char *module, char *procedure, char *error); 

extern void	PrintSetup();		// print out setup for record
extern void	Print_Header(FILE *);
extern void	Print_Histogram();
extern void	Print_Nuclei();
extern void	Print_Verlet();
extern void	Print_Verlet2();
extern void	Print_Clist();
extern void	Print_q();
extern void	Print_Q();
extern void	Print_qproduct();
extern void	Printout();
extern int 	Visualize(int);
extern void	Print_gr();

extern int	Write_Conf(long timestep);
extern int	Read_Conf(char *);
extern int	Read_MultiConf(FILE *fPtr);
#endif

#endif
                                                                                                                                                                                                                         src/lists.h                                                                                         0000600 0143352 0000144 00000005065 11311315621 012363  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    lists.h  
    author:     Peng Yi at MIT
    date:       October 22, 2006
    purpose:    header file for lists.c
*/
#ifndef __LISTS_HEADER
#define __LISTS_HEADER

#include "header.h"

#define MAXNNEIGHBORS	500

typedef struct {
   long			n;			// # of neighbors
   long			reverse[MAXNNEIGHBORS];	// whether reverse or not
   long			site[MAXNNEIGHBORS];	// id of site
   vector		dr[MAXNNEIGHBORS];	// distance to that neighbor
   molstruct		*mol[MAXNNEIGHBORS];	// id of mol
} neighborlist;					// neighbor list of one end bead

#ifdef __LISTS_MODULE


#ifdef CELL_LIST
cellstruct	*Cell;
lvector		M[MAXNBOX];		// number of cells along box length
long		NCELLS;
long		neighcellplus[54];	// neighbor cells before a move plus 
					// those after a move (no double count)
long		nneighcellplus;
#endif	/* CELL_LIST */

#ifdef VERLET_LIST
long		*vlistplus;		// the Verlet neighbors before a move PLUS 
					// those after a move (no double count)
					// PLUS itself
long		nverletplus;		// number of the vlistplus elements 
long		maxnverlet;		// maximum # of verlet neighbors of any particle
#endif	/* VERLET_LIST */


#else

extern long	NeighborList(molstruct *molm, neighborlist *list);
extern long	DB_NeighborList(molstruct *molm, long sitem, neighborlist *list);
extern long	IDR_NeighborList(molstruct *molm, long sitem, neighborlist *list);

#ifdef CELL_LIST
extern cellstruct	*Cell;
extern lvector	M[MAXNBOX];
extern long	NCELLS;
extern long	neighcellplus[54];
extern long	nneighcellplus;
extern void	New_CL();
extern void	CL_Update(long, long, long);

extern long	CL_Findcell(molstruct *, long, long ibox, long PBC);
extern long	CL_Neighbor(long, long, long);
extern void	CL_Init();
extern void	CL_Build();
extern void	CL_Destroy();
extern void	CL_Add(molstruct *, long);
extern void	CL_Delete(molstruct *, long);
extern void	CL_Relink(molstruct *moli);

#endif	/* CELL_LIST */

#ifdef VERLET_LIST
extern long	*vlistplus;
extern long	nverletplus;
extern long	maxnverlet;
extern void 	New_Vlist();
extern void	New_Vlist_LL();
extern void	Update_Vlist(long, vector, vector); 
extern void	Update_Vlist2(long, vector, vector); 
#endif	/* VERLET_LIST */

extern long	elem_index(long *, long, long);
extern long	Listlength(liststruct *);
extern int	List_Insert(liststruct **, long);
extern int	List_Remove(liststruct **, long);
extern void	Printlist(liststruct *);
extern void	Free_List(liststruct **);	//free the whole list
extern int	List_is_Empty(liststruct *);

extern void	New_ConnectList();
extern void	New_Clist();
extern void	Update_Clist(long, vector);

#endif	//ifdef __LIST_MODULE

#endif
                                                                                                                                                                                                                                                                                                                                                                                                                                                                           src/motion.h                                                                                        0000600 0143352 0000144 00000002051 11175453634 012541  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    motion.h
    author:     Peng Yi at MIT
    date:       October 26, 2007
    purpose:    move molecule, sites, functions borrowed and modified from Pieter's
*/
#ifndef __MOTION_HEADER
#define __MOTION_HEADER

#include "types.h"

#ifdef __MOTION_MODULE

#include "header.h"

#else

extern sphere	SiteSpherical(molstruct *, long);
extern vector	SiteCartesian(molstruct *, long, sphere);
extern void	MolSpherical(molstruct *);
extern void	MolCartesian(molstruct *);
extern void	AllSpherical();
extern void	AllCartesian();

// Molecule manipulators

extern void	SiteCopy(molstruct *molm, long m, molstruct *moln, long n, long p);
extern void	SiteSwap(molstruct *molm, long m, molstruct *moln, long n);
extern void	MolReverse(molstruct *mol1);	// reverse a molecule withOUT cell list update
extern void	MolFlip(molstruct *moli);	// reverse a molecule WITH cell list update
extern molstruct	MolAdd(molstruct *mol1, molstruct *mol2);
extern void	ChangeAxis(long system, vector scale);
extern void	ChangeVolume(long system, double scale);
#endif

#endif
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       src/mymath.h                                                                                        0000600 0143352 0000144 00000000377 10714730275 012541  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    mymath.h  
    author:     Peng Yi at MIT
    date:       October 22, 2006
    purpose:    header file for mymath.c
*/
#ifndef __MYMATH_HEADER
#define __MYMATH_HEADER

#ifdef __MYMATH_MODULE
#include "header.h"


#else


#endif

#endif
                                                                                                                                                                                                                                                                 src/position.h                                                                                      0000600 0143352 0000144 00000003413 11624275667 013112  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
  program:      position.c
  author:       Peng Yi at MIT
  date:         October 20, 2006
  purpose:      Position calculation, Xtal nuclei identification, and speed-up
*/
#ifndef __POSITION_HEADER
#define __POSITION_HEADER

#include "types.h"

#ifdef __POSITION_MODULE

#include "header.h"

long		PBC;	//periodic boundary condition variable
			//0: no pbc; 1: cubic pbc; 2: truncated octahedron pbc
long		critconnect;	// # of critical connection in LJ system

long		nsegment[MAXNMOLS];			// # of segments on a chain
long		seg_stat[MAXNMOLS][MAXNMOLSITES];	// segments stat for chains

#else

// variables

extern long	PBC;
extern long	critconnect;

extern long	nsegment[MAXNMOLS];
extern long	seg_stat[MAXNMOLS][MAXNMOLSITES];


// functions

extern double	AdjustAngle(double x);
extern void	MapInBox(vector *p);
extern vector 	MapInBox2(vector *p, long PBC, long system);
extern double	DistSQ(vector p, vector q, long system);

extern void	InitLattice(long, long, double, long);
extern vector	ranor();
extern void	Amorph(long nmols, long nmolsites, double LX, double LY, double LZ, long PBC);

extern void	sc_lattice(long, double, long);
extern void	bcc_lattice(long, double, long);
extern void	fcc_lattice(long, double, long);
extern void	randomconf(long, double, long);
extern long	crystal(molstruct *moli, long site);
extern long     getnuclsize(long, long);
extern void 	Find_Nuclei(long clusdef);
extern void	Find_Nuclei_LJ();
extern void	Find_Nuclei2();
extern void	Find_Nuclei1();
extern void	Find_Nuclei_p2(long clusdef);
extern vector	CoM_MaxNucleus(long system);
extern vector   cylindershape(beadstruct *nucleus, long size, long nuclid);
extern vector   cylinder(beadstruct *nucleus, long size, long nuclid);
extern void	Find_segments();
extern vector	Seg_smooth(long);

#endif

#endif
                                                                                                                                                                                                                                                     src/random.h                                                                                        0000600 0143352 0000144 00000000630 10714730275 012512  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    Name:       random.h
    Author:     Peng Yi at MIT
    Date:       October 23, 2006
    Purpose:    Using code from 'Numerical Recipes', chapter 7.
*/
#ifndef __RANDOM_HEADER
#define __RANDOM_HEADER

#ifdef __RANDOM_MODULE
/*
#include <stdio.h>
#include <math.h>
#include <time.h>
*/

#include "header.h"

#else

extern float 	ran1(long *idum);
extern double 	gauss(double, double);

#endif

#endif

                                                                                                        src/rebridge.h                                                                                      0000644 0143352 0000144 00000001720 11175735547 013037  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	rebridge.c
    author:	Pieter J. in 't Veld for MIT Boston
    date:	May 4, 2001
    purpose:	Trimer rebridging module, using Mavrantzas et al.,
    		Macromolecules 32, 5072 (1999).
*/
#ifndef __REBRIDGE_HEADER
#define __REBRIDGE_HEADER

#include "header.h"
/*
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "complex.h"
#include "roots.h"
#include "random.h"
#include "types.h"
#include "globals.h"
#include "move.h"
#ifdef MPI
#include "mpicomm.h"
#endif
*/

#define B_MAXNROOTS	17

typedef struct {
  vector		r;
  sphere		s;
} position;

#endif

#ifndef __REBRIDGE_MODULE

// enter 7 positions (both cartesian and spheric)

extern long Frame(vector *, vector *, vector *, vector *, vector *, vector *);
extern long Feasible(vector *p, sphere *s);
extern long Rebridge(vector *p, sphere *s);
extern double Jacobian(vector *p, sphere *s);
extern long RebridgeSetup(
	molstruct *molm, long end, long reverse, vector *p, sphere *s);

#endif

                                                src/roots.h                                                                                         0000644 0143352 0000144 00000001670 11173374411 012411  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	roots.h
    author:	Pieter J. in 't Veld
    date:	May 10, 2001
    purpose:	roots of upto fourth order polynomials
*/
#ifndef __ROOTS_HEADER
#define __ROOTS_HEADER

#include <math.h>
#include <complex.h>


#ifdef __ROOTS_MODULE

#else

extern long R_ComplexCount(long n, double complex *z);
extern long R_SolveFirstOrder(double *c, double complex *z);
extern long R_SolveSecondOrder(double *c, double complex *z);
extern long R_SolveThirdOrder(double *c, double complex *z);
extern long R_SolveFourthOrder(double *c, double complex *z);
extern long R_SolveAnalytical(long n, double *c, double complex *z);
extern long R_NewtonRaphson(
  double y(double), double x_low, double x_high, double *root);
extern long R_Brent(
  double y(double), double x_low, double x_high, double *root);
extern long R_SolveNumerical(
  double y(double), double x_low, double x_high, long n_mesh, long depth,
  double *roots, long n_max);

#endif

#endif

                                                                        src/sample.h                                                                                        0000600 0143352 0000144 00000016720 11553053603 012515  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    sample.h  
    author:     Peng Yi at MIT
    date:       October 22, 2006
    purpose:    header file for sample.c
*/
#ifndef __SAMPLE_HEADER
#define __SAMPLE_HEADER

#include "header.h"

#define D_NDIST		MAXNDIST	// MAXNDIST defined in input.h
#define D_BINSIZE	0.02		// not sure what is this, should make sense if reduced unit is used

#ifdef __SAMPLE_MODULE


long			D_DENSITY, D_ENERGY, D_PRESSURE, D_DRIFT,
			D_TORSION, D_BONDA, D_BONDL, D_RADIAL, 
			D_LOCALP2, D_XTALSIZE;
diststruct		*D_Density, *D_Energy, *D_Pressure, *D_Drift,
			*D_Torsion, *D_Bonda, *D_Bondl, *D_Radial, 
			*D_Localp2, *D_Xtalsize;
diststruct		*D_Distributions[D_NDIST];		// hook all sorts of distribution to this
								// adding convenience when operate to all
								// distributions

long			dynvar;			//(1=NMAX; 2=Ql)
long			NMAXmiddle, NMAXbinsize, NMAXbins, Qlbins;	//dynamic variable range
double			Qlmiddle, Qlbinsize, kQ, kN;
double			P2middle, kP2;

double			drift2[MAXNSYSTEMS];
//structure property variables
double			P2[MAXNSYSTEMS];	// global orientation order
double			P2M[MAXNSYSTEMS];	// global orientation order
double			P2z[MAXNSYSTEMS];	// P2 with respect to z-axis
double			transfrac[MAXNSYSTEMS];	// trans state fraction
double 			Q6[MAXNSYSTEMS];	//bond orientation order parameter Q_l
double 			Q4[MAXNSYSTEMS];
/*
complex			Qlm[2*l_of_Ylm+1];	//Qlm = YlmAlphaSum/AlphaSum
complex			YlmAlphaSum[2*l_of_Ylm+1];
double			AlphaSum;
*/ // no need to be global

long			cnndist[25];		//connected neighbors number distribution

double				Qltemp;
/*
double				AlphaSumtemp;
complex				YlmAlphaSumtemp[2*l_of_Ylm+1];
complex				Qlmtemp[2*l_of_Ylm+1];
*/  // no need to be global

//crystal nuclei related variables
long			*sizeofnucl, *sizeofnuclp2;		//size of certain crystal nucleus
long			*sizedist, *sizedistp2;		//how many nuclei of each size
long			MAXSIZE[MAXNSYSTEMS];	// maximum nuclei size in each system
long			Nnucl[MAXNSYSTEMS];	// # of Xtal nuclei
long			Xtal[MAXNSYSTEMS],	// total # of Xtal-like segments
			realXtal[MAXNSYSTEMS],
			secondNmax[MAXNSYSTEMS];
long			nmax[MAXNSYSTEMS][10];	// 10 biggest size, nmax[i][0]=MAXSIZE[i]
beadstruct		nucleus[MAXNMOLS*MAXNMOLSITES/8];	// beads in the biggest nucleus
long				*sizeofnucl2;		//size of certain crystal nucleus
long				*sizedist2;		//how many nuclei of each size
long				MAXSIZE2;		//maximum nuclei size in the system
long				Nnucl2;			//# of Xtal nuclei
long				Xtal2;			//total # of Xtal-like particles

double			critp2,		// critical local p2 value to be considered Xtal phase
			critqlproduct,	// critical qlproduct to be considered connected
			critangle;	// critical angle between two vectors
long			critconnect;	// critical connection number to distinguish phase
long			NGRBINS;		// bins used to measure radial distribution function
double			Alpha;			//damping factor of eta update
double			CRIT;			//criterion of uniform sampling

#else

extern long		D_DENSITY, D_ENERGY, D_PRESSURE, D_DRIFT,
			D_TORSION, D_BONDA, D_BONDL, D_RADIAL,
			D_LOCALP2, D_XTALSIZE;
extern diststruct	*D_Density, *D_Energy, *D_Pressure, *D_Drift,
			*D_Torsion, *D_Bonda, *D_Bondl, *D_Radial,
			*D_Localp2, *D_Xtalsize;
extern diststruct	*D_Distributions[D_NDIST];

extern long		dynvar;
extern long		NMAXmiddle, NMAXbinsize, NMAXbins, Qlbins;	//dynamic variable range
extern double		Qlmiddle, Qlbinsize, kQ, kN;
extern double		P2middle, kP2;

extern double		drift2[MAXNSYSTEMS];
//structure property variables
extern double		P2[MAXNSYSTEMS], P2M[MAXNSYSTEMS], P2z[MAXNSYSTEMS];
extern double		transfrac[MAXNSYSTEMS];
extern double		Q6[MAXNSYSTEMS], Q4[MAXNSYSTEMS];
/*
extern complex		Qlm[2*l_of_Ylm+1];	//Qlm = YlmAlphaSum/AlphaSum
extern complex		YlmAlphaSum[2*l_of_Ylm+1];
extern double		AlphaSum;
*/
extern long		cnndist[25];		//connected neighbors number distribution

extern double			Qltemp;
/*
extern double			AlphaSumtemp;
extern complex			YlmAlphaSumtemp[2*l_of_Ylm+1];
extern complex			Qlmtemp[2*l_of_Ylm+1];
*/

//crystal nuclei related variables
extern long		*sizeofnucl, *sizeofnuclp2;		//size of certain crystal nucleus
extern long		*sizedist, *sizedistp2;		//how many nuclei of each size
extern long		MAXSIZE[MAXNSYSTEMS];
extern long		Nnucl[MAXNSYSTEMS];		
extern long		Xtal[MAXNSYSTEMS],	
			realXtal[MAXNSYSTEMS],
			secondNmax[MAXNSYSTEMS];
extern long		nmax[MAXNSYSTEMS][10];		// size of 10 biggest nuclei
extern beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES/8];
extern long			*sizeofnucl2;		//size of certain crystal nucleus
extern long			*sizedist2;		//how many nuclei of each size
extern long			MAXSIZE2;		//maximum nuclei size in the system
extern long			Nnucl2;			//# of Xtal nuclei
extern long			Xtal2;			//total # of Xtal-like particles

extern double		critp2, critqlproduct, critangle;
extern long		critconnect;	
extern long		NGRBINS;
extern double		Alpha, CRIT;


extern long 	mod(long, long);
extern long 	factorial(long);
extern long	intpow(long, long);
extern long	intlog(long, long);
extern double 	plgndr(int, int, double);
extern double complex 	sfharmonics(int, int, double, double);
extern double complex	sfharmonics2(int, int, double, double, double);

extern double complex 	aveqlm(int, int, long);
//extern complex 	tildeqlm(int, int, long);
extern double complex	qlproduct(int, long, long);
extern double	qlproductSQ(int, molstruct *, long, molstruct *, long);
extern double  	ql(int, long);
extern double complex 	aveQlm(int, int);
extern double	CalcQl(int);
extern double complex 	aveQlm1(int, int);
extern double	CalcQl1(int);
extern void	local_q_update();
extern void	New_Qlm(int);
extern void	New_Qlmtemp(int);
extern void	New_Qlm_NoVlist(int);
extern void	Update_Qlm(int, long, vector, vector); 

extern void	Calc_Qlm(long L);

extern int	Qlbinfinder(double);
extern int	NMAXbinfinder(int);
extern double	etaQl(double Ql);
extern double   etaNMAX(long MAXSIZE);
extern double	etaP2(double);

extern vector	CenterofMass(molstruct *moli);			// center of mass of one chain
extern vector	groupCoM(beadstruct *group, long nsites);	// center of mass of a list
extern matrix	GyrationTensor(molstruct *moli);
extern matrix	groupGyraTensor(beadstruct *group, long nsites);
extern vector	fshape(vector *eigvalue);			// shape descriptors
extern matrix	InertiaTensor(molstruct *moli);
extern matrix	groupInerTensor(beadstruct *group, long nsites);
extern vector		CenterofNucleus(long nuclid, molstruct *molm);
extern void		InitSample();
extern void		ReinitSample(diststruct **d);
extern void		SampleEnergy();			// it will be called by outsider, so claimed here
extern void		SampleDistributions();
extern void		PrintDistributions();
extern void		S_PrintAll();
extern void		Sample2D();
extern void		SampleN2NSQ();
extern void		SampleDrift();
extern void		SampleP2();
extern void		SampleP2All();
extern void		Dist_p2();
extern void	SampleM_Q();
extern void	SampleConnection();
extern void		SampleSpherical();
extern void		Dist_Spherical();

extern void	Init_Sample(char *argv[]);
extern void	Sample_Energy();
extern void	Sample_Done();
extern void	Sample_Pressure();
extern void	Sample_All();
extern void	Sample_Histogram();
extern void	Sample_G(int);		//sample Gibbs free energy

extern double	R2_gyration(molstruct *moli);		// radius of gyration square
extern double	R2_n2n(molstruct *);			// end-to-end distance square
extern vector	Center_of_Mass(long *, long);		// center of mass of listed atoms
extern matrix   Gyration_Tensor(long *, long);		// gyration tensor of listed atoms
extern void	radial(char *);				// radial distribution function
extern void	sq(FILE *, char *);			// structure factor calculation

#endif

#endif
                                                src/types.h                                                                                         0000600 0143352 0000144 00000022001 11576550177 012402  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:    types.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Definition of types used in the NPT project
*/
#ifndef __TYPES_HEADER
#define __TYPES_HEADER

#include	<complex.h>
#include	"input.h"

typedef struct{   
   double	x, y, z;
} vector;

typedef struct{
   double	alpha,		// bond angle
		beta,		// torsion angle
		d;		// bond length
} sphere;

typedef struct{
  vector x, y, z;
} matrix;
  
typedef struct{
   long		x, y, z;
} lvector;

typedef struct{   	
  double 	lj, ljcorr,	// LJ potential, tail correction
		ljin, 		// LJ energy b/w sites on same chain
		ljout,		// LJ energy b/w sites on diff chains
		hs,		// hard sphere interaction
		stretch,	// bond stretching interaction
		bending,	// bending energy
		torsion,	// torsional energy
		corr,		// total correction
		bonded,		// bond related energy
		nonbonded,
		tot;
} vstruct;			// molecule potential energy structure

typedef struct{
  double	lj, 
		ljin,
		ljout,
		stretch,
		torsion,
		tot;
} wstruct;			//virial function structure

struct listNode{			//linked list structure
   long			neighbor;
   struct listNode * 	nextPtr;	//self-referential structure, typedef later
};
typedef	struct listNode	liststruct;
typedef liststruct *	liststrPtr;	//pointer pointing to a liststruct variable


/*
typedef struct{ 				// data structure for a single particle
  vector 	p,
		pv;				// pv is for Verlet list use
  long		box,				// which box this particle belongs to
  		nverlet,			// number of Verlet neighbors
 		vlist[MAXVERLETNEIGH],		// indexes of verlet neighbors
 		icell;				// which cell it belongs to, CELL_LIST

  // members above are essential for basic simulation
  // members below are for sampling, or more sophisticated simulation

  long 		nbond;		//number of nearest neighbor bonds
  complex	qlm[2*l_of_Ylm+1];	//aveq_lm for a certain l
  complex	ylmalphasum[2*l_of_Ylm+1];
  double	alphasum;	//qlm = ylmalphasum / alphasum
  double	ql;  		//ql value, but not used
  long		nconnect;	//how many particles are connected to this particle
  long		clist[MAXCONNEIGH];	
  long		nuclid;		//if type=1, which nucleus it belongs to, default -1
  long			nuclid2;	//for Find_Nuclei debug 6/9/2007
} molstruct;                
*/


typedef 
 struct cell_list {
   vector		center, p_min, p_max;	// center, lower corner, higher corner of this cell
   long			box,			// which box this cell belongs to
   			nsites,			// how many sites in this cell
			nempty,			// used as a stack
			empty[MAXNCELLSITES];
   struct mol		*mol[MAXNCELLSITES];	// pointer of molecules that have sites in this cell
   long			molsite[MAXNCELLSITES];	// id of sites (relative to mol) in this cell
   long			nneigh;			// # of neighboring cells
   struct cell_list	*neigh[27];		// neighboring cell id
} cellstruct;


typedef 
  struct mol {
   long		box, 			// which box this chain molecule belongs to
		nsites,			// how many sites on this chain molecule
		fix,			// whether this molecule is movable (1: fixed; 0: movable)
                flip,			// whether is a flip compared to the beginning
  		flags[MAXNMOLSITES],	// activated flag
		parent[MAXNMOLSITES],	// parent site
                type[MAXNMOLSITES];	// group type
   vector	p[MAXNMOLSITES];	// cartesian coord. of each site
   vector	origin;			// original position for drift calc.
   sphere	s[MAXNMOLSITES];	// spherical coord. of each site   
   cellstruct	*cell[MAXNMOLSITES];	// cell id of each site in CELL_LIST implementation
   long		cellsite[MAXNMOLSITES];	// position of this site in its corresponding cell

   // members above are essential for basic simulation
   // members below are for sampling, or more sophisticated simulation

   double	p2[MAXNMOLSITES];	// local orientational parameter
   long		np2[MAXNMOLSITES];	// # of neighbors to calculate p2
   long		nconn[MAXNMOLSITES],	// how many connection for each site
		connsite[MAXNMOLSITES][MAXCONNEIGH];	// site id of connected neighbors for each site
   struct mol   *connmol[MAXNMOLSITES][MAXCONNEIGH];	// mol id of connected neighbors for each site
   long		nuclid[MAXNMOLSITES];

   double complex	qlm[MAXNMOLSITES][2*l_of_Ylm+1],
		ylmalphasum[MAXNMOLSITES][2*l_of_Ylm+1];
   double	alphasum[MAXNMOLSITES],
		q6[MAXNMOLSITES],
		q4[MAXNMOLSITES];
   vector	vp2[MAXNMOLSITES];		// p2 vector added Aug11.09
} molstruct;


typedef
  struct site {
   long			system;
   struct mol		*mol;
   long			mol_site;
   struct cell_list	*cell;
   long			cellsite;
   struct type		*type;
   long			flags, parent;
   vector		p;
   sphere		s;
} sitestruct;				// to describe a site, like molstruct to a molecule

typedef struct {
   molstruct 		*moli;
   long	     		site;
} beadstruct;				// (Jun6,2009) a more efficient version of sitestruct

typedef struct {
   double	SIGMA, EPSILON;
} mixstruct;


typedef
  struct type {					// force coefficent for different type of atoms
  						// e.g., CH3 and CH2 in PE are different in some
						// forcefield models
   double	M, Q,				// mass and charge
   		SIGMA, EPSILON,			// LJ interaction
		LSTRETCH, KSTRETCH,		// stretching
		KBENDING, THETA,		// bending
		TORSION[6],			// torsion	
		HS;				// hard-sphere
   mixstruct	mix[MAXNTYPES];
} typestruct;					


typedef struct {
   long		move, acc_move,			// normal site displacement
		erot, acc_erot,			// end mer rotation
		flip, acc_flip,			// one mer flip move
		re, acc_re,			// rebridging
		eb, acc_eb,			// endbridging
                db, acc_db,			// doublebridging
                idr, acc_idr, 			// internal double rebridging
		rot, acc_rot,			// rotation move
		vol, acc_vol,			// volume change
		swap, acc_swap,			// gibbs swap
		gibbsvol, acc_gibbsvol,		// gibbs volume change
		cbmc, acc_cbmc,			// cbmc
                rep, acc_rep,			// reptation
		movemol, acc_movemol,		// move whole molecule
		seq, acc_seq;			// sequence
} avstruct;


typedef
  struct structdist {
   char			*header;
   long			level, 			// 0 (1D dist), 1 (2D dist), etc.
			nbins, 			// # of bins of this level
			n, 			// total counts
			startbin, 		// first bin
			ncount;
   double		*binsize;		// could be multi-dimension
   struct structdist	*dist;			// sub dist, for multi-dimension
   long			*bin;			// counts in each bin
   double		*data, 
			*cweight, 
			*average;		// average of x, x^2, ..., x^(D_NAVERAGE-1)
} diststruct;


typedef struct{
   double	lbox, lx, ly, lz,	// dimension
		xy, xz, yz,		// for triclinic box, definition see lammps manual
		vol, pres, temp,	// volume, pressure, temperature
		rv, rb, rc;		// cutoff
   double	drmax, dlmax, damax;	// move size: displacement, volume, angle
   long		maxsize, xtal, nnucl;
   double	Ql;
} boxstruct;


typedef struct
{ 
    long                *systems, *cycle, *tape, *mols, *sites, *types, *volchange,
                        *swap, *displace, *inserts, *cavinserts, *box,
                        *rotation, *reptation, *endbridge, *rebridge, *cbmc,
                        *bridges, *fixed, *semifixed, *seed, *temper, *sample,
                        *stretch, *free, *ends, *system;
} nstruct;

typedef struct
{
   long			*hs, *lj, *ljshift, *ljlrc, *stretch, *bending, 
			*torsion, *virial, *scalecut, 
			*nvt, *npt, *gibbs, *mpi,
			*density, *energy, *pressure, *drift, *torsional,
			*bonda, *bondl, *radial, *localp2, *xtalsize;
} commandstruct;			// used in bridgestruct for binary history i/o

/* Pieter's
typedef struct
{
    long                *hs, *lj, *stretch, *ptorsion, *coulomb, *virial,
                        *polymer, *mpi, *async, *monodisperse, *bias, 
                        *nvt, *npt, *gibbs, *insert, *widom, *canonical, 
                        *cavity, *density, *tails_etc, *radial, *energy, 
                        *torsion, *re_torsion, *hs_dens, *temper, *jacob,
                        *d_bridge, *e_profile, *w_profile, *n_profile,
                        *b_length, *b_angle, *bonded, *loopreentry,
                        *e_n_function, *orient, *density3d, *densfree,
                        *denstalobr, *orientfree, *orienttalobr,
                        *orientcorr01, *orientcorr02, *orientcorr03, 
                        *orientcorr04, *orientcorr05, *orientcorrCR;
} commandstruct;
*/

typedef struct
{
    long                n, 						// system id, value
			*nmols, *nsites;				// pointers
    double              *pres, *vol, *temp, *drmax, *dlmax, *damax;
    //vector              a, b, c;					// system dimen. values
    double		*lx, *ly, *lz;					// system dimen. pointers
} systemstruct;				// system info., for binary history file, used in bridgestruct


typedef struct
{
    long                *d_type;
    nstruct             n;
    commandstruct       command;		// flags
    typestruct          *type;			// type of sites
    systemstruct        system[MAXNSYSTEMS];	// system info., some pointers, some values
    avstruct            *av[MAXNSYSTEMS];
    vstruct             *v[MAXNSYSTEMS];	// potential energy, all pointers
    wstruct             *vir[MAXNSYSTEMS];	// virial, all pointers
    diststruct          *dist[MAXNDIST];	// distribution, all pointers
    molstruct           *mol[MAXNSYSTEMS];	// molecule info., all pointers
} bridgestruct;			// for binary history file i/o

#endif
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               src/units.h                                                                                         0000600 0143352 0000144 00000001572 10755407334 012404  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	units.h
    author:	Peng Yi at MIT
    date:	Nov. 7, 2007
    purpose:	Convert general SI units to program units
*/

#ifndef __UNITS_HEADER
#define __UNITS_HEADER

#define	N_NAV		6.0221367e23	/* Avogadro constant 1/mol */
#define	N_KB		1.380658e-23	/* Boltzmann constant J/K */
#define N_ATM		101325.024	/* pascals/atm */

typedef struct {
   double	LENGTH,			// SI length = system length * LENGTH
		ENERGY,
		PRESSURE,
		TEMP,
		ANGLE,
		MASS,
		DENSITY;
} unitstruct;

#ifdef __UNITS_MODULE

#include "header.h"

unitstruct		unit;
long			ConvertUnits;

#else

extern unitstruct	unit;
extern long		ConvertUnits;

extern double		CalcMass();
extern void		CalcUnits(long flag);
extern void		ConvertSIToSystem();
extern void		ConvertSystemToSI();
extern void		CoorSystem2SI();
extern void		CoorSI2System();
extern void		InitUnits();
extern void		PrintUnits();

#endif

#endif

                                                                                                                                      src/varbridge.h                                                                                     0000600 0143352 0000144 00000001335 10743214326 013176  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	varbridge.h
    author:	Peng Yi borrowed from Pieter J. in 't Veld
    date:	January 10, 2008
    purpose:	This module functions as a bridge between old style
    		system variables and new style system variables.
*/
#ifndef __VARBRIDGE_HEADER
#define __VARBRIDGE_HEADER

#ifdef __VARBRIDGE_MODULE

#include "header.h"


#else

extern bridgestruct *Bridge();
extern bridgestruct *BridgeMap(long flag_mpi);
/*
extern void BridgeInitMemory(bridgestruct *bridge, long nsystems, long NMols);
extern void BridgeResetVariable(bridgestruct *bridge);
extern void BridgeAddVariable(bridgestruct *to, bridgestruct *from);
extern void BridgeCopyInstant(
		bridgestruct *to, bridgestruct *from, long system);
*/

#endif

#endif

                                                                                                                                                                                                                                                                                                   src/vector.h                                                                                        0000600 0143352 0000144 00000003055 10750122066 012531  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  /*
    program:	vector.h
    author:	Pieter J. in 't Veld for UT at Austin
    date:	April 1, 1999
    purpose:	3D vector algebra
*/

#ifndef __VECTOR_HEADER
#define __VECTOR_HEADER

#define V_E_SINGULAR	1
/*
typedef struct {
  double	x, y, z; } vector;
typedef struct {
  vector	x, y, z; } matrix;
*/
#ifdef __VECTOR_MODULE

#include "header.h"
/*
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include "types.h"
*/
#else

// Vector operations

extern void V_Null(vector *);
extern void V_Unit(vector *);
extern void V_Negate(vector *);
extern double V_Dot(vector *, vector *);
extern vector V_Add(vector *, vector *);
extern vector V_Subtr(vector *, vector *);
extern vector V_Cross(vector *, vector *);
extern vector V_Mult(double, vector *);
extern vector V_Rotate(vector *, double, double);
extern void V_Print(vector);

// Matrix operations

extern void M_Null(matrix *);
extern void M_Unit(matrix *);
extern void M_Negate(matrix *);
extern matrix M_Dot(matrix *, matrix *);
extern matrix M_Add(matrix *, matrix *);
extern matrix M_Subtr(matrix *, matrix *);
extern matrix M_Mult(double, matrix *);
extern matrix M_Rotation(double, double);
extern matrix M_XNormalRotate(vector *);
extern matrix M_Rotate(vector *, vector *);
extern matrix M_Orientation(vector *, vector *);
extern matrix M_Transpose(matrix *);
extern int M_Inverse(matrix *, matrix *);
extern double M_Det(matrix *);
extern void M_Print(matrix);
extern vector M_eig(matrix);
extern vector V_eig(matrix, double);

// Hybrid operarions

extern vector MV_Dot(matrix *, vector *);


#endif

#endif

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   makefile                                                                                            0000700 0143352 0000144 00000011473 11564764045 012010  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  SHELL	=	/bin/sh
#
# Local definitions
#
#C      =       cc4 -O3 -ffast-math -mfpmath=387 -mieee-fp -ftree-vectorize
#C      =       cc4 -O3 -ffast-math -mfpmath=387 -mieee-fp -funroll-all-loops -ftree-vectorize
#C	=	mpicc4 -O3 
C	=	icc -O3 
#
# TARGET RULES
#
#
.f.o:;   
#
FLAGS  = 
#path for source files
SRCDIR = ./src
#path for header files	
VPATH = ./src		
#
LIB= -lm 
HEADERS = distributions.h ensembles.h forcefield.h history.h init.h io.h lists.h motion.h position.h random.h \
rebridge.h roots.h sample.h units.h varbridge.h vector.h input.h header.h globals.h types.h

OBJECTS = distributions.o ensembles.o forcefield.o history.o init.o io.o lists.o motion.o position.o random.o \
rebridge.o roots.o sample.o units.o varbridge.o vector.o
	
#

#

main: $(OBJECTS) main.o
	$(C) -o main  $(FLAGS) $(LIB) $(OBJECTS) main.o;
	rm $(OBJECTS) main.o;

block:	$(OBJECTS) block.o
	$(C) -o block $(FLAGS) $(LIB) block.o $(OBJECTS);
	rm $(OBJECTS) block.o

builder: $(OBJECTS) builder.o
	$(C) -o builder $(FLAGS) $(LIB) builder.o $(OBJECTS);
	rm $(OBJECTS) builder.o

conf2car: $(OBJECTS) conf2car.o
	$(C) -o conf2car $(FLAGS) $(LIB) conf2car.o $(OBJECTS);
	rm $(OBJECTS) conf2car.o

conf2lammps: $(OBJECTS) conf2lammps.o
	$(C) -o conf2lammps $(FLAGS) $(LIB) conf2lammps.o $(OBJECTS);
	rm $(OBJECTS) conf2lammps.o

conf2DL_POLY: $(OBJECTS) conf2DL_POLY.o
	$(C) -o conf2DL_POLY $(FLAGS) $(LIB) conf2DL_POLY.o $(OBJECTS);
	rm $(OBJECTS) conf2DL_POLY.o

conf2gro: $(OBJECTS) conf2gro.o
	$(C) -o conf2gro $(FLAGS) $(LIB) conf2gro.o $(OBJECTS);
	rm $(OBJECTS) conf2gro.o

config: $(OBJECTS) config.o
	$(C) -o config $(FLAGS) $(LIB) config.o $(OBJECTS);
	rm $(OBJECTS) config.o;

dumphst: $(OBJECTS) dumphst.o
	$(C) -o dumphst $(FLAGS) $(LIB) dumphst.o $(OBJECTS);
	rm $(OBJECTS) dumphst.o;

embed:	$(OBJECTS) embed.o
	$(C) -o embed $(FLAGS) $(LIB) embed.o $(OBJECTS);
	rm $(OBJECTS) embed.o

hist:	$(OBJECTS) histogram.o
	$(C) -o hist $(FLAGS) $(LIB) histogram.o $(OBJECTS);
	rm $(OBJECTS) histogram.o;

hst:	$(OBJECTS) hst.o
	$(C) -o hst $(FLAGS) $(LIB) hst.o $(OBJECTS);
	rm $(OBJECTS) hst.o

hstcomb:	$(OBJECTS) hstcomb.o
	$(C) -o hstcomb $(FLAGS) $(LIB) hstcomb.o $(OBJECTS);
	rm $(OBJECTS) hstcomb.o

lammps2car: $(OBJECTS) lammps2car.o
	$(C) -o lammps2car $(FLAGS) $(LIB) lammps2car.o $(OBJECTS);
	rm $(OBJECTS) lammps2car.o

lammps2hst: $(OBJECTS) lammps2hst.o
	$(C) -o lammps2hst $(FLAGS) $(LIB) lammps2hst.o $(OBJECTS);
	rm $(OBJECTS) lammps2hst.o

lammpstest: $(OBJECTS) lammpstest.o
	$(C) -o lammpstest $(FLAGS) $(LIB) lammpstest.o $(OBJECTS);
	rm $(OBJECTS) lammpstest.o

lmp_conf: $(OBJECTS) lmp_conf.o
	$(C) -o lmp_conf $(FLAGS) $(LIB) lmp_conf.o $(OBJECTS);
	rm $(OBJECTS) lmp_conf.o

extractlammps: $(OBJECTS) extractlammps.o
	$(C) -o extractlammps $(FLAGS) $(LIB) extractlammps.o $(OBJECTS);
	rm $(OBJECTS) extractlammps.o

#
#
main.o: main.c $(HEADERS)
	$(C) -c $(SRCDIR)/main.c
builder.o: builder.c $(HEADERS)
	$(C) -c $(SRCDIR)/builder.c
block.o:  block.c $(HEADERS)
	$(C) -c $(SRCDIR)/block.c
conf2car.o: conf2car.c $(HEADERS)
	$(C) -c $(SRCDIR)/conf2car.c
conf2lammps.o: conf2lammps.c $(HEADERS)
	$(C) -c $(SRCDIR)/conf2lammps.c
conf2DL_POLY.o: conf2DL_POLY.c $(HEADERS)
	$(C) -c $(SRCDIR)/conf2DL_POLY.c
conf2gro.o: conf2gro.c $(HEADERS)
	$(C) -c $(SRCDIR)/conf2gro.c
config.o: config.c $(HEADERS)
	$(C) -c $(SRCDIR)/config.c
dumphst.o: dumphst.c  $(HEADERS)
	$(C) -c $(SRCDIR)/dumphst.c
embed.o: embed.c $(HEADERS)
	$(C) -c $(SRCDIR)/embed.c
histogram.o:	histogram.c $(HEADERS)
	$(C) -c $(SRCDIR)/histogram.c
hst.o:	hst.c $(HEADERS)
	$(C) -c $(SRCDIR)/hst.c
hstcomb.o:	hstcomb.c $(HEADERS)
	$(C) -c $(SRCDIR)/hstcomb.c
lammps2car.o: lammps2car.c $(HEADERS)
	$(C) -c $(SRCDIR)/lammps2car.c
lammps2hst.o: lammps2hst.c $(HEADERS)
	$(C) -c $(SRCDIR)/lammps2hst.c
lmp_conf.o: lmp_conf.c $(HEADERS)
	$(C) -c $(SRCDIR)/lmp_conf.c
extractlammps.o: extractlammps.c $(HEADERS)
	$(C) -c $(SRCDIR)/extractlammps.c


distributions.o: distributions.c $(HEADERS)
	$(C) -c $(SRCDIR)/distributions.c
ensembles.o: ensembles.c $(HEADERS)
	$(C) -c $(SRCDIR)/ensembles.c
forcefield.o: forcefield.c $(HEADERS)
	$(C) -c $(SRCDIR)/forcefield.c
history.o: history.c $(HEADERS)
	$(C) -c $(SRCDIR)/history.c
init.o:	 init.c $(HEADERS)
	$(C) -c $(SRCDIR)/init.c
io.o:	io.c $(HEADERS)
	$(C) -c $(SRCDIR)/io.c
lists.o: lists.c $(HEADERS)
	$(C) -c $(SRCDIR)/lists.c
motion.o: motion.c $(HEADERS)
	$(C) -c $(SRCDIR)/motion.c
position.o: position.c $(HEADERS)
	$(C) -c $(SRCDIR)/position.c
random.o:  random.c $(HEADERS)
	$(C) -c $(SRCDIR)/random.c
rebridge.o: rebridge.c $(HEADERS)
	$(C) -c $(SRCDIR)/rebridge.c
roots.o: roots.c $(HEADERS)
	$(C) -c $(SRCDIR)/roots.c
sample.o: sample.c $(HEADERS)
	$(C) -c $(SRCDIR)/sample.c
units.o:  units.c $(HEADERS)
	$(C) -c $(SRCDIR)/units.c
varbridge.o: varbridge.c $(HEADERS)
	$(C) -c $(SRCDIR)/varbridge.c
vector.o: vector.c $(HEADERS)
	$(C) -c $(SRCDIR)/vector.c


clean:
	rm *~; rm output





                                                                                                                                                                                                     setup                                                                                               0000600 0143352 0000144 00000006322 11274610664 011362  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  Title:  Nucleation free energy barrier sampling,
        Monte Carlo, umbrella sampling, PYS forcefield
************************************************************************

/* system basic information */

moltype		=	chain
ConvertUnits	=	1		// convert SI to program units
Rho		=	0.635		// monomer number density
kT		=	450  
P		=	1.0
NBOX		=	1
NMOLS		=	40 		// for all boxes
NSITES		=	3120		// for all molecules in all boxes
NPARTS		=	3120		// total particle number

/* simulation specifics */

INITCONF	=	initconf 
PBC		=	1		// periodic boundary condition
Stage		=	3
Nequil		=	0     		// MC sequences for equilibration
Nprod		=	128		// MC sequences for production run
SEQUENCE	=	1		// number of MC cycles in each sequence
TRIALRUN	=	0
PROD		=	1
ITAPE		=	4
ICONF		=	2048
NBLOCKS		=	2048
NGSAMPLE	=	128
IRADIAL		=	256
NGRBINS		=	200
Alpha		=	0.4
CRIT		=	0.4

/* interaction model */

NTYPES		=	4
DLJ		=	4
TYPEMASS	=	{ 15  		14 		15		14 }     
SIGMA		=	{ 4.01		4.01		4.01		4.01 }       
EPSILON		=	{ 0.468608	0.468608	0.468608	0.468608 }
KSTRETCH	=	{ 47095.39688 	47095.39688 	47095.39688	47095.39688 }
LSTRETCH	=	{ 1.53 		1.53		1.53		1.53  }
KBENDING	=	{ 502.08    	502.08 		502.08		502.08  }
THETA		=	{ 70.5  	70.5 		70.5		70.5 }
HS		=	{ 0.0		0.0	 	0.0		0.0 }
TORTYPE		=	OPLS
TORSION0	=	{ 0.0		0.0		0.0		0.0 }
TORSION1	=	{ 6.6944	6.6944		6.6944		6.6944 }
TORSION2	=	{ -3.627528  	-3.627528 	-3.627528	-3.627528 }     
TORSION3	=	{ 13.55616   	13.55616  	13.55616	13.55616 }  
TORSION4	=	{ 0         	0		0.0		0.0 }           
TORSION5	=	{ 0         	0 		0.0		0.0 }                
V_LJ		=	1		// Lennard Jones force field
V_LJSHIFT	=	0		// LJ forcefield shift at cutoff
V_LJLRC		=	1		// LJ forcefield long range correction	
V_HS		=	0		// hard sphere force field
V_STRETCH	=	1
V_BENDING	=	1
V_TORSION	=	1
FIXBONDL	=	0
FIXBONDA	=	0
V_VIRIAL	=	1		// calculate virial function
Rc    	        =	2.5	   	// cutoff distance for LJ interaction
Rb       	=	1.5		// cutoff distance for nearest neighbor
Rv		=	2.9		// Verlet cutoff distance
Rp		=	2.0		// local p2 cutoff distance
Rclow		=	2.5
Rconn		=	1.5
SCALECUTOFF	=	0

/* ensemble setup */

E_NVT		=	0
E_NPT		=	1
E_GIBBS		=	0
NDISPLACE	=	0		// %, percentage of one site displacement
NREPTATION	=	6		// %
NENDROT		=	6		// %
NCBMC		=	0 		// %
NTRIALCONF	=	4		// cbmc parameter
NTRIALFIRSTBEAD	=	50		// cbmc parameter
NREBR		=	32		// %
NENDBR		=	50		// %
NFLIP		=	6		// %
NGIBBSVOL	=	10		// %
NSWAP		=	250		// %
NVOLCHANGE	=	15 		// #, number of volume change move per MC cycle
SUCC_DISP	=	0.25		// targeted acceptance ratio of one site displacement
SUCC_VOL	=	0.5		// targeted acceptance ratio of volume change
DRMAX		=	0.1		// max displacement distance
DLMAX		=	0.01		// max volume change scale
DAMAX		=	0.04

/* sampling */

D_DENSITY	=	0
D_ENERGY	=	0
D_PRESSURE	=	0
D_DRIFT		=	0
D_TORSION	=	1
D_BONDA		=	0
D_BONDL		=	0
D_RADIAL	=	0
D_LOCALP2	=	0
D_XTALSIZE	=	0

DYNVAR		=	NMAX		// dynamic variable
kP2		=	0.00
P2middle	=	40
NMAXmiddle	=	10		// center of sampling window
NMAXbinsize	=	1		// binsize of histogram
NMAXbins	=	21		// bins
kN		=	0.15		// bias potential depth
Qlmiddle	=	0.020
Qlbinsize	=	0.002
Qlbins		=	6
kQ		=	1e5
critqlproduct	=	0.5
critconnect	=	7
critp2		=	0.5
etaN		=	{ 1.0	2.0	3.0 }
etaQ		=	{ 4.0 	5.0	6.0 }

END_SETUP
                                                                                                                                                                                                                                                                                                              c24/in                                                                                              0000600 0143352 0000144 00000003024 11200067576 011211  0                                                                                                    ustar   pengyi                          users                                                                                                                                                                                                                  Title:  Nucleation free energy barrier sampling,
        Monte Carlo, umbrella sampling
************************************************************************

/* system basic information */

NBOX		=	1
NMOLS		=	32		// for all boxes
NSITES		=	768		// for all molecules in all boxes
Rho		=	0.86		// monomer number density
PBC		=	1		// periodic boundary condition
NX		=	4
NY		=	4
NZ		=	1
a		=	4.74		// 4.8 for octane
b		=	7.8		// 8.3 for octane
c		=	32.0		// 11.6 for octane, 30 for eicosane
angle		=	40.0

lattice constant at 300K by preliminary simulation

/* interaction model */

NTYPES		=	2 
DLJ		=	4
TYPEMASS	=	{ 15  14 }     
SIGMA		=	{ 4.01   4.01 }       
EPSILON		=	{ 0.468608   0.468608 }
KSTRETCH	=	{ 47095.39688 47095.39688 }
LSTRETCH	=	{ 1.53 	1.53  }
KBENDING	=	{ 502.08    502.08   }
THETA		=	{ 70.5  70.5 }
HS		=	{ 0.0	0.0 }
TORSION0	=	{ 0.0     0.0     }
TORSION1	=	{ 6.6944    6.6944    }
TORSION2	=	{ -3.627528  -3.627528 }     
TORSION3	=	{ 13.55616   13.55616  }  
TORSION4	=	{ 0         0 }           
TORSION5	=	{ 0         0 }                
V_LJ		=	1		// Lennard Jones force field
V_LJSHIFT	=	0		// LJ forcefield shift at cutoff
V_LJLRC		=	1		// LJ forcefield long range correction	
V_HS		=	0		// hard sphere force field
V_STRETCH	=	1
V_BENDING	=	1
V_TORSION	=	1
FIXBONDL	=	0
FIXBONDA	=	0
V_VIRIAL	=	1		// calculate virial function
Rc    	        =	2.5     	// cutoff distance for LJ interaction
Rb       	=	1.5		// cutoff distance for nearest neighbor
Rv		=	2.9		// Verlet cutoff distance
Rconn		=	2.0
SCALECUTOFF	=	0

END_SETUP
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            