/*
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
