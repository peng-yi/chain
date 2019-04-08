/*
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

// OPLS() and OPLS2() are two equivalent forms, with same k1, k2, and k3. Difference is that OPLS2 
// require fewer cos operations, thus faster.

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
   double		result;

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
   dr[0]	=	MapInBox2(dr, PBC, system);
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

//==============================================================================//
//	CalcVLJ_MPI(): use multiple processes to calc. LJ energy and virial	//
//		       added on 5/4/2013					//
//==============================================================================//
void CalcVLJ_MPI()
{
   long		i, ibox;
   molstruct	*moli;

#ifdef myMPI
   int		proc;
   MPI_Status	status;
   double	vlj[2*MAXNSYSTEMS];			// LJ energy and virial for each box
   
   static int	molid[MAXNMOLS];			// id of first mol for each process
   static int	init = 1;
   int		ave, extra;

   if (init) {
      init	=	0;
      
      //------Assign molecules to processes, try to leave proc #0 has small load------//
      ave	=	NMOLS/MPI_numprocs;
      extra	=	NMOLS%MPI_numprocs;

      molid[0]	=	0;
      for (i=1; i<MPI_numprocs-extra+1; i++) {
	 molid[i]	=	molid[i-1] + ave;
      }
      for (i=MPI_numprocs-extra+1; i<MPI_numprocs+1; i++) {
	 molid[i]	=	molid[i-1] + ave + 1;
      }
      /*
      if (MPI_myid==0) {
	 printf("# Molecules assignment\n");
	 for (i=0; i<MPI_numprocs; i++) {
	    printf("# Process # %d: molecules # %d - %d\n", i, molid[i], molid[i+1]-1);
	 }
      }
      */
   }

   for (i=0; i<2*MAXNSYSTEMS; i++) {
      vlj[i]	=	0.0;
   }
#endif

   for (ibox=0; ibox<NBOX; ibox++) {
      v[ibox].lj	=	0.0;
      vir[ibox].lj	=	0.0;
   }
   if (!V_LJ)	return;

#ifdef myMPI 
   //------calculate v[box].lj and vir[box].lj 
   for (moli=mol+molid[MPI_myid]; moli<mol+molid[MPI_myid+1]; moli++) {
      if ( (ibox=moli->box) >=0) {
	 vlj[ibox]	+=	VLJMol(moli, &(vlj[MAXNSYSTEMS+ibox]));
      }
   }
  
   //------transfer data to process #0 and process data
   if (MPI_myid > 0) {					// send data to process #0
      MPI_Send(vlj, 2*MAXNSYSTEMS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
   }
   else {						// process #0
      for (ibox=0; ibox<NSYSTEMS; ibox++) {		// count results from own process
	 v[ibox].lj	+=	vlj[ibox];
	 vir[ibox].lj	+=	vlj[MAXNSYSTEMS+ibox];
      }
      for (proc=1; proc<MPI_numprocs; proc++) {		// receive results from other processes
	 MPI_Recv(vlj, 2*MAXNSYSTEMS, MPI_DOUBLE, proc, 0, MPI_COMM_WORLD, &status);

	 for (ibox=0; ibox<NSYSTEMS; ibox++) {		// count results from other processes
	    v[ibox].lj		+=	vlj[ibox];
	    vir[ibox].lj	+=	vlj[MAXNSYSTEMS+ibox];
	 }
      }
      for (ibox=0; ibox<NSYSTEMS; ibox++) {		// fix double-counting
	 v[ibox].lj	*=	0.5;
	 if (V_VIRIAL) {
	    vir[ibox].lj	*=	0.5;
	 }
      }
   }
#endif
  
   return;
}


//============================================================================//

void CalcVLJout()		// LJ energy b/w sites on different chains ONLY
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
   long		i, k=0, system=molm->box;
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
   dr[0]	=	MapInBox2(dr, PBC, system);
   dr[1]	=	MapInBox2(dr+1, PBC, system);

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
   long		i, k=0, system=molm->box;
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
   dr[0]=	MapInBox2(dr, PBC, system);
   dr[1]=	MapInBox2(dr+1, PBC, system);
   dr[2]=	MapInBox2(dr+2, PBC, system);

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
/*
int	i,j,k;
int		pair[MAXNTYPES][MAXNTYPES];		// # of pairs for each type combo
for (i=0; i<MAXNTYPES; i++) {
   for (j=0; j<MAXNTYPES; j++) {
      pair[i][j]	=	0;
   }
}
*/
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

//printf("NActive[0]=%d\n", NActive[0]);

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
			   if (V_LJSHIFT) {
                              v[system].ljcorr	+=	Epsilon/vol * Sigma * Sigma * Sigma * (4.0/(3*rc9)-2.0/rc3);
			   }
			   else {
                              v[system].ljcorr	+=	Epsilon/vol * Sigma * Sigma * Sigma * (1.0/(3*rc9)-1.0/rc3);
			   }
			   nsites[system]	++;

//pair[molm->type[m]][moln->type[n]]++;
                        }
/*
printf("nsites[system]=%d\n", nsites[0]);
for (j=0;j<NTYPES; j++){
   for (k=0; k<NTYPES; k++) {
      printf("pair[%d][%d]=%d\n", j, k, pair[j][k]);
   }
}
*/
   for (system=0; system<NSYSTEMS; system++) {
      v[system].ljcorr	*=	nsites[system] ? 8.0*pi/3.0*NSites[system]*NActive[system]/nsites[system] : 0.0;
   }
//v[0].ljcorr	=	8.0*pi/3.0*81000000*Epsilon/vol*Sigma*Sigma*Sigma*(1.0/(3*rc9)-1.0/rc3);

   free(nsites);
}   

//==============================================================================//
//	CalcVLJLRC(): LJ energy long range correction, a simple version		//
//		      added on 5/5/2013						//
//==============================================================================//
void CalcVLJLRC()
{
   molstruct	*moli;
   int		i, j, k, ibox;
   int		nsites[MAXNTYPES];
   double	vol, rc3, rc9;
   double	lrc, totlrc;
   double	Sigma, Epsilon;

   /*
   molstruct	*molj;
   int		pair[MAXNTYPES][MAXNTYPES];		// # of pairs for each type combo
   */

   for (ibox=0; ibox<NBOX; ibox++) {			// do for each box
      
      v[ibox].ljcorr	=	0.0;

      /*
      for (i=0; i<MAXNTYPES; i++) {
	 for (j=0; j<MAXNTYPES; j++) {
	    pair[i][j]	=	0;
	 }
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
	 for (i=0; i<moli->nsites; i++) {

	    for (molj=mol; molj<mol+NMOLS; molj++) {
	       if (moli!=molj) {
		  for (j=0; j<molj->nsites; j++) {
		     pair[moli->type[i]][molj->type[j]]	++;
		  }
	       }
	    }
	 }
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
	 for (i=0; i<moli->nsites; i++) {
	    for (j=0; j<moli->nsites; j++) {
	       if (labs(i-j)>=4) {
		  pair[moli->type[i]][moli->type[j]]	++;
	       }
	    }
	 }
      }
      
      vol		=	BOX[ibox].vol;
      totlrc		=	0.0;

      for (j=0; j<NTYPES; j++) {
	 for (k=0; k<NTYPES; k++) {
	    if ( pair[j][k]>0 ) {

	       
	       Sigma	=	type[j].mix[k].SIGMA;
	       Epsilon	=	type[j].mix[k].EPSILON;
	       
               //printf("pair[%d][%d] = %d\n", j, k, pair[j][k]);
               //printf("Sigma = %f  Epsilon = %f\n", Sigma, Epsilon);

	       rc3		=	BOX[ibox].rc/Sigma;
	       rc3		=	rc3 * rc3 * rc3;
	       rc9		=	rc3 * rc3 * rc3;

	       lrc		=	Epsilon/vol * Sigma * Sigma * Sigma * (1.0/(3*rc9) - 1.0/rc3);
	       lrc		*=	pair[j][k];
	       totlrc		+=	lrc;
	    }
	 }
      }
      v[ibox].ljcorr	=	totlrc * 8.0 *pi/3.0 * 81000000.0/80937108;
      */
	    
      for (i=0; i<MAXNTYPES; i++) {
	 nsites[i]	=	0;
      }

      //------Counting how many atoms of each type
      for (moli=mol; moli<mol+NMOLS; moli++) {
	 if (moli->box == ibox) {
	    for (i=0; i<moli->nsites; i++) {
	       nsites[moli->type[i]]	++;
	    }
	 }
      }

      vol		=	BOX[ibox].vol;
      
      for (j=0; j<NTYPES; j++) {
	 for (k=0; k<NTYPES; k++) {
	    if (nsites[j]>0 && nsites[k]>0) {
	       
	       //printf("nsites[%d]*nsites[%d] = %d\n", j, k, nsites[j]*nsites[k]);

	       Sigma	=	type[j].mix[k].SIGMA;
	       Epsilon	=	type[j].mix[k].EPSILON;

	       rc3		=	BOX[ibox].rc/Sigma;
	       rc3		=	rc3 * rc3 * rc3;
	       rc9		=	rc3 * rc3 * rc3;

	       lrc		=	8.0/(9*rc9) - 8.0/(3*rc3);
	       lrc		*=	nsites[j] * nsites[k] / vol * Sigma * Sigma * Sigma;
	       lrc		*=	pi * Epsilon;
	       
	       v[ibox].ljcorr	+=	lrc;	
	    }
	 }
      }
   }

   return;
}

//===========================================================================//

void CalcVCorr()
{
   long		ib;

   if (V_LJLRC) {
      //CalcVLJCorr();
      CalcVLJLRC();
   }

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

   CalcVCorr();
/*
   for (ib=0; ib<NBOX; ib++) {				// energy tail correction simple scaling
      v[ib].ljcorr	/=	volscale;		// if cutoff NOT scale with box size
      v[ib].corr	=	v[ib].ljcorr;
   }
*/
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
   int		ib;
   //time_t	start, end;
   
   for (ib=0; ib<NBOX; ib++) {
      vstructNull(v+ib);
      wstructNull(vir+ib);
   }
   //start	=	time(NULL);
#ifdef	myMPI
   if (V_LJ)		CalcVLJ_MPI();
#else
   if (V_LJ)		CalcVLJ();
#endif
   //end	=	time(NULL);
   //printf("\nCalcVLJ() running time was %lf seconds.\n", difftime(end, start));
   
   //start	=	time(NULL);
   if (V_HS)		CalcVHS();
   if (V_STRETCH)	CalcVStretch();
   if (V_BENDING)	CalcVBending();
   if (V_TORSION)	CalcVTorsion();
   //end	=	time(NULL);
   //printf("\nCalcVbond() running time was %lf seconds.\n", difftime(end, start));

   //start	=	time(NULL);
   CalcVCorr();
   //end	=	time(NULL);
   //printf("\nCalcVCorr() running time was %lf seconds.\n", difftime(end, start));

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
            //printf("%s\ttrial %ld\t of %ld\t, vt[j].tot=%f\tvt[j].nonbonded=%f\n", s, j, i, Vt[j].tot, Vt[j].nonbonded);
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

	 //printf("Adding %ld done, trial #%ld.\n", i, n);
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
          fprintf(foutput, "\n#############total energy problem in box %ld.###############\n", i);
      if ( fabs(VIR[i] - vir[i].tot) > 1.0e-6)
          fprintf(foutput, "\n#############total virial problem in box %ld.###############\n", i);

      fprintf(foutput, "\n\tBOX[%ld]\n\n", i);
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
