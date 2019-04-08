/*
    program:    correlation.h
    author:     Peng Yi at MIT
    date:       March 29, 2008
    purpose:    Correlation function calculation, only used in history file analysis
		so no need to compile these functions separately, just directly include
		into history file analysis source file hst.c
*/

#ifndef __CORR_HEADER
#define __CORR_HEADER

#define Lcorr		2000		// maximum corr. separation b/w two records 
#define	NN		8

//**************************************//
//	data structure and variables	//
//**************************************//
typedef struct{
   int		t;				// the index in the record file
   vector	com[MAXNMOLS];			// center of mass of each molecule
   vector	end2end[MAXNMOLS];		// end-to-end vector of each molecule
   vector	pos[MAXNMOLS*MAXNMOLSITES];	// position of each bead
   vector	rhead[MAXNMOLS], 		// average of head NN beads
		rtail[MAXNMOLS],		// average of tail NN beads
		rctr[MAXNMOLS];			// average of center NN beads
   float	aveR2, PE;			// average n2n distance square, potential energy
   float	nonbond;
   float	test;
} corr_data;					// data whose time correlation need to be calculated

typedef struct{
   float	com;				// mean square displacement of com of each chain
   float	end2end;			// end-to-end vector correlation
   float	pos;				// mean square displacement of each bead
   float	rhead, rtail, rctr;
   float	test;				// test
   float	aveR2, PE, nonbond;

   int		norm;				// normalization factor
} corr_results;					// correlation function

corr_data	store[Lcorr+1];
corr_results	corr[Lcorr+1];

//******************************//
//	store variable values	//
//******************************//
void corr_store(int nstore, int nlabel)		// do some calc. on each record and store variable values
{						// nstore is the index in store and nlabel is the index in record
   molstruct 	*moli;
   int		i, j;
   float	d;
   vector	r, com_sys, rtemp;

   int		system=0;		// one system only
   float	n2n;			// n2n distance square

   store[nstore].t	=	nlabel;
   store[nstore].test	=	(float)nmax[system][0];
   store[nstore].PE	=	v[0].tot;
   store[nstore].nonbond	=	v[0].nonbonded;

   V_Null(&com_sys);			// set vector to zero
   n2n	=	0.0;

   for (i=0; i<NMOLS; i++) {
      moli	=	mol + i;

      // Position of center of mass of chains

      store[nstore].com[i]		=	CenterofMass(moli);
      com_sys				=	V_Add(&com_sys, store[nstore].com+i);

      // Average position of first NN beads
/*
      V_Null(&rtemp);
      for (j=0; j<NN; j++) {
          rtemp		=	V_Add(&rtemp, moli->p+j);
      }
      store[nstore].rhead[i]		=	V_Mult(1.0/NN, &rtemp);

      // Average position of last NN beads
      V_Null(&rtemp);
      for (j=moli->nsites-NN; j<moli->nsites; j++) {
          rtemp		=	V_Add(&rtemp, moli->p+j);
      }
      store[nstore].rtail[i]		=	V_Mult(1.0/NN, &rtemp);

      // Average position of center NN beads
      V_Null(&rtemp);
      for (j=NSITES/NMOLS/2-NN/2; j<NSITES/NMOLS/2+NN/2; j++) {
          rtemp		=	V_Add(&rtemp, moli->p+j);
      }
      store[nstore].rctr[i]		=	V_Mult(1.0/NN, &rtemp);
*/

      // End-to-end vector

      if (moli->flip) {
         r	=	V_Subtr(moli->p, moli->p+(moli->nsites-1));
         rtemp			=	store[nstore].rtail[i];
         store[nstore].rtail[i]	=	store[nstore].rhead[i];
         store[nstore].rhead[i]	=	rtemp;
      }	
      else {
         r	=	V_Subtr(moli->p+(moli->nsites-1), moli->p);
      }
      d 	=	sqrt(V_Dot(&r, &r));
      store[nstore].end2end[i]	=	V_Mult(1.0/d, &r);

      // End-to-end distance
      
      n2n	+=	d*d;		// square of end-to-end distance

      // Position of every bead

      if (moli->flip) {
         for (j=0; j<moli->nsites; j++)
            store[nstore].pos[i*(NSITES/NMOLS)+j]	=	moli->p[moli->nsites-1-j];
      }
      else {
         for (j=0; j<moli->nsites; j++)
            store[nstore].pos[i*(NSITES/NMOLS)+j]	=	moli->p[j];
      }
   }
   com_sys		=	V_Mult(1.0/NMOLS, &com_sys);
   store[nstore].aveR2	=	n2n/NMOLS;

   // Fix system drift (2/28/2012)
/*
   for (i=0; i<NMOLS; i++) {
      moli	=	mol+i;
      for (j=0; j<moli->nsites; j++) {
         store[nstore].pos[i*(NSITES/NMOLS)+j]	=	
		V_Subtr(store[nstore].pos+(i*(NSITES/NMOLS)+j), &com_sys);
      }
      store[nstore].com[i]	=	V_Subtr(store[nstore].com+i, &com_sys);

      //store[nstore].rhead[i]	=	V_Subtr(store[nstore].rhead+i, &com_sys);
      //store[nstore].rtail[i]	=	V_Subtr(store[nstore].rtail+i, &com_sys);
      //store[nstore].rctr[i]	=	V_Subtr(store[nstore].rctr+i, &com_sys);
   }
*/
   return;
}

//**************************************************************//
//	calculate correlation between time n1 and n2>=n1	//
//**************************************************************//
void corr_calc(int n1, int n2, int dt)	// n1 and n2 are positions in store, dt is the real separation
{
   static int	init = 1;
   int		i, j;
   vector	dr;

   if (init) {					// initialization
      for (i=0; i<=Lcorr; i++) {
         corr[i].com		=	0;
         corr[i].rhead		=	0;
         corr[i].rtail		=	0;
         corr[i].rctr		=	0;
         corr[i].end2end	=	0;
         corr[i].pos		=	0;
         corr[i].test		=	0;
	 corr[i].PE		=	0;
	 corr[i].aveR2		=	0;
	 corr[i].nonbond	=	0;
         corr[i].norm		=	0;
      }
      init	=	0;
   }

   for (i=0; i<NMOLS; i++) {
      dr		=	V_Subtr(&(store[n1].com[i]), &(store[n2].com[i]));
      corr[dt].com	+=	V_Dot(&(dr), &(dr));
/*
      dr		=	V_Subtr(&(store[n1].rhead[i]), &(store[n2].rhead[i]));
      corr[dt].rhead	+=	V_Dot(&(dr), &(dr));

      dr		=	V_Subtr(&(store[n1].rtail[i]), &(store[n2].rtail[i]));
      corr[dt].rtail	+=	V_Dot(&(dr), &(dr));

      dr		=	V_Subtr(&(store[n1].rctr[i]), &(store[n2].rctr[i]));
      corr[dt].rctr	+=	V_Dot(&(dr), &(dr));
*/
      corr[dt].end2end	+=	V_Dot(&(store[n1].end2end[i]),&(store[n2].end2end[i]));

      for (j=0; j<NSITES/NMOLS; j++) {
         dr		=	V_Subtr(&(store[n1].pos[i*NSITES/NMOLS+j]), 
					&(store[n2].pos[i*NSITES/NMOLS+j]));
         corr[dt].pos	+=	V_Dot(&(dr), &(dr));
      }

      for (j=0; j<NN; j++) {
         dr		=	V_Subtr(&(store[n1].pos[i*NSITES/NMOLS+j]), 
					&(store[n2].pos[i*NSITES/NMOLS+j]));
         corr[dt].rhead	+=	V_Dot(&(dr), &(dr));
      }
      for (j=NSITES/NMOLS-NN; j<NSITES/NMOLS; j++) {
         dr		=	V_Subtr(&(store[n1].pos[i*NSITES/NMOLS+j]), 
					&(store[n2].pos[i*NSITES/NMOLS+j]));
         corr[dt].rtail	+=	V_Dot(&(dr), &(dr));
      }
      for (j=NSITES/NMOLS/2-NN/2; j<NSITES/NMOLS/2+NN/2; j++) {
         dr		=	V_Subtr(&(store[n1].pos[i*NSITES/NMOLS+j]), 
					&(store[n2].pos[i*NSITES/NMOLS+j]));
         corr[dt].rctr	+=	V_Dot(&(dr), &(dr));
      }
   }
   corr[dt].test	+=	store[n1].test * store[n2].test;
   corr[dt].PE		+=	store[n1].PE * store[n2].PE;
   corr[dt].nonbond	+=	store[n1].nonbond * store[n2].nonbond;
   corr[dt].aveR2	+=	store[n1].aveR2 * store[n2].aveR2;

   corr[dt].norm	++;
   return;
}

//**********************************************//
//	correlation normalized by sample number	//
//**********************************************//
void corr_norm()
{
   int		i;

   corr[0].com		/=	corr[0].norm * NMOLS;		// m.s.d. = 0 if corr. length = 0
   corr[0].end2end	/=	corr[0].norm * NMOLS;
   corr[0].pos		/=	corr[0].norm * NSITES;
   corr[0].rhead	/=	corr[0].norm * NMOLS * NN;
   corr[0].rtail	/=	corr[0].norm * NMOLS * NN;
   corr[0].rctr		/=	corr[0].norm * NMOLS * NN;
   corr[0].test		/=	corr[0].norm;
   corr[0].PE		/=	corr[0].norm;
   corr[0].nonbond	/=	corr[0].norm;
   corr[0].aveR2	/=	corr[0].norm;

   //corr[0].com	=	1.0;	// m.s.d. = 0 if corr. length = 0
   corr[0].end2end	=	1.0;	// n2n correlation normalized to 1 at corr. length=0
   corr[0].test		=	1.0;

   for (i=1; i<=Lcorr; i++) {
      corr[i].com	/=	corr[i].norm * NMOLS;
      corr[i].end2end	/=	corr[i].norm * NMOLS * corr[0].end2end;
      corr[i].pos	/=	corr[i].norm * NSITES;
      corr[i].rhead	/=	corr[i].norm * NMOLS * NN;
      corr[i].rtail	/=	corr[i].norm * NMOLS * NN;
      corr[i].rctr	/=	corr[i].norm * NMOLS * NN;
      corr[i].test	/=	corr[i].norm * corr[0].test;
      corr[i].PE	/=	corr[i].norm;
      corr[i].nonbond	/=	corr[i].norm;
      corr[i].aveR2	/=	corr[i].norm;
   }
   //corr[0].PE		=	1.0;
   //corr[0].aveR2	=	1.0;
   return;
}

//**************************************//
//	print correlation results	//
//**************************************//
void corr_print()				// print out time correlation function
{
   int		i;

   printf("Monte Carlo cycles between consecutive sampling = %ld\n\n", ITAPE);

   printf("Column 1. Correlation time step.\n");
   printf("Column 2. Center of mass of chains.\n");
   printf("Column 3. Center of mass of bead.\n");
   printf("Column 4. Average position of head NN beads.\n");
   printf("Column 5. Average position of tail NN beads.\n");
   printf("Column 6. Average position of ctr NN beads.\n");
   printf("Column 7. End-to-end vector of chains.\n");
   printf("Column 8. potential energy.\n");
   //printf("Column 8. nonbonded energy.\n");
   printf("Column 9. End-to-end distance square.\n");
   printf("Column 10. Number of Data points for average.\n");
   //printf("Column 5. Test.\n");
   printf("\n");

   for (i=0; i<=Lcorr; i++) {
      if (corr[i].norm > 0) {
         printf("%-5d ", corr[0].norm-corr[i].norm);
         printf("%8.3f ", corr[i].com);
         printf("%8.3f ", corr[i].pos);
	 printf("%8.3f ", corr[i].rhead);
	 printf("%8.3f ", corr[i].rtail);
	 printf("%8.3f ", corr[i].rctr);
         printf("%6.3f ", corr[i].end2end);
         printf("%8.3f ", corr[i].PE);
         //printf("%f\t", corr[i].nonbond);
         printf("%8.3f ", corr[i].aveR2);
         printf("%5d\n", corr[i].norm);
         //printf("%f\n", corr[i].test);
      }
   }
   return;
}

//**************************************//
//	main correlation function	//
//**************************************//
void correlation()
{
   int		i;
   static int	nstore = 0, 			// index in store
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
   return;
}
#endif
