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
