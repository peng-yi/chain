/*
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


