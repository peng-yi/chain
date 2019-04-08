/* 
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
long		oldXtal[MAXNSYSTEMS];
long		oldrealXtal[MAXNSYSTEMS];
long		oldNnucl[MAXNSYSTEMS];
long		oldNmax[MAXNSYSTEMS][5];

void StoreSystem()				//store system information
{
   long		i;

   for (i=0; i<NBOX; i++) {
      oldv[i]		=	v[i];		// energy
      oldvir[i]		=	vir[i];		// virial
      oldBOX[i]		=	BOX[i];		// box dimensions
      oldNSites[i]	=	NSites[i];
      oldNMols[i]	=	NMols[i];

      oldP2[i]		=	P2[i];
      oldP2M[i]		=	P2M[i];
      oldP2z[i]		=	P2z[i];
      oldQ6[i]		=	Q6[i];
      oldtransfrac[i]	=	transfrac[i];      
      oldXtal[i]	=	Xtal[i];
      oldrealXtal[i]	=	realXtal[i];
      oldNnucl[i]	=	Nnucl[i];
      oldNmax[i][0]	=	nmax[i][0];
      oldNmax[i][1]	=	nmax[i][1];
      oldNmax[i][2]	=	nmax[i][2];
      oldNmax[i][3]	=	nmax[i][3];
   }
   return;
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
      Xtal[i]		=	oldXtal[i];
      realXtal[i]	=	oldrealXtal[i];
      Nnucl[i]		=	oldNnucl[i];
      nmax[i][0]	=	oldNmax[i][0];
      nmax[i][1]	=	oldNmax[i][1];
      nmax[i][2]	=	oldNmax[i][2];
      nmax[i][3]	=	oldNmax[i][3];
   }
   return;
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
   double	drbound;

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
	 drbound		=	2.0 * (BOX[ib].rv - BOX[ib].rc);
#else
	 drbound		=	0.25 * BOX[ib].lbox;
#endif	/* VERLET_LIST */

         if (BOX[ib].drmax > drbound)	BOX[ib].drmax = drbound;
         
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
	 drbound		=	2.0 * (BOX[ib].rv - BOX[ib].rc);
#else
	 drbound		=	0.25 * BOX[ib].lbox;
#endif	/* VERLET_LIST */

         if (BOX[ib].drmax > drbound)	BOX[ib].drmax = drbound;
         
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
        printf("molcheck flag error, molid %ld siteid %ld\n", moli-mol, i);
     if (moli->parent[i] != i-1)
        printf("molcheck parent error, molid %ld siteid %ld\n", moli-mol, i);
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

   // Need not do MapInBox to change particle position, just let the particle move.  
   // The PBC effect is taken care of in energy calc. and other places.

//printf("mcmove del dV = %f\n", dV);

   dV		+=	VAddSites(molm, site, site);
   arg		=	dV/(BOX[ibox].temp);

/*printf("mcmove add dV = %f\n", dV);
if (fabs(dV-(v[ibox].tot-v_old.tot)) > 0.000001)  {printf("mcmove_err\n"); exit(1);}
*/

   if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {	// accept
      av[ibox].acc_move	++; 

/*printf("%ld %ld\n", molm-mol, site);
printf("%f %f %f %f %f %f\n", mol_old.p[site].x, mol_old.p[site].y, mol_old.p[site].z, molm->p[site].x, molm->p[site].y, molm->p[site].z);
printf("displacement accept  old energy =%f  new energy =%f\n", v_old.tot, v[ibox].tot);
*/
/*
      for (i=0; i<4; i++)		// update spherical coord.
         if (site+i < molm->nsites)
            molm->s[site+i]	=	SiteSpherical(molm, site+i);
*/
   }
   else {	// reject
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


void displacemol()		// displace one molecule without changing its internal conf.
{
   long		i, system;
   molstruct	*moli;
   //	No content for now
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

//******************************//
//	Volume change move	//
//******************************//

long mcvol()
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

   // each dimension changes independently
   scale.x	=	exp( BOX[system].dlmax * (ran1(seed)-0.5) );
   scale.y	=	exp( BOX[system].dlmax * (ran1(seed)-0.5) );
   scale.z	=	exp( BOX[system].dlmax * (ran1(seed)-0.5) );
   volscale	=	scale.x * scale.y * scale.z;

   ChangeAxis(system, scale);	// scale coordinates and box sizes

#ifdef CELL_LIST
   CL_Destroy();		// need to remake cell list due to vol change
   CL_Init();
   CL_Build();
#endif

   //scale	=	exp( BOX[system].dlmax * (ran1(seed)-0.5) );
   //ChangeVolume(system, scale);
   
   //CalcV_mcvol(volscale);
   CalcV();
 
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
#ifdef CELL_LIST
      CL_Destroy();
      CL_Init();
      CL_Build();
#endif
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
printf("counter = %ld chainid = %ld site = %ld oldenergy = %f newenergy = %f \n", counter, molm-mol, site, v_old.tot, v[ib].tot);
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
   av[ib].rep	++;
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
//      fprintf(stderr, "MOVE = %ld\n");
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
printf("another endbridging move, site = %ld\n", site);
printf("Old %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
printf("molm\n");
for (i=molm->nsites-4; i<molm->nsites; i++)
   printf("%ld %ld %ld %f %f %f\n", i, molm->parent[i], molm->flags[i], molm->p[i].x, molm->p[i].y, molm->p[i].z);      
printf("moln\n");
for (i=site-7; i<=site+2; i++)
   printf("%ld %ld %ld %f %f %f\n", i, moln->parent[i], moln->flags[i], moln->p[i].x, moln->p[i].y, moln->p[i].z); 
*/
      VDeleteSites(moln, site-3-flag, site-1);	// Delete changed sites
      if (flag)
        VDeleteSites(molm, molm->nsites-1, molm->nsites-1);
/*
printf("Del done %f %f %f %f %f \n", v[system].tot, v[system].stretch, v[system].bending, v[system].torsion, v[system].nonbonded);
printf("molm\n");
for (i=molm->nsites-4; i<molm->nsites; i++)
   printf("%ld %ld %ld %f %f %f\n", i, molm->parent[i], molm->flags[i], molm->p[i].x, molm->p[i].y, molm->p[i].z);      
printf("moln\n");
for (i=site-7; i<=site+2; i++)
   printf("%ld %ld %ld %f %f %f\n", i, moln->parent[i], moln->flags[i], moln->p[i].x, moln->p[i].y, moln->p[i].z); 
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
   printf("%ld %ld %ld %f %f %f\n", i, molm->parent[i], molm->flags[i], molm->p[i].x, molm->p[i].y, molm->p[i].z);      
printf("moln\n");
for (i=moln->nsites-4; i<moln->nsites; i++)
   printf("%ld %ld %ld %f %f %f\n", i, moln->parent[i], moln->flags[i], moln->p[i].x, moln->p[i].y, moln->p[i].z); 
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
   printf("%ld %ld %ld %f %f %f\n", i, molm->parent[i], molm->flags[i], molm->p[i].x, molm->p[i].y, molm->p[i].z);      
printf("moln\n");
for (i=0; i<moln->nsites; i++)
   printf("%ld %ld %ld %f %f %f\n", i, moln->parent[i], moln->flags[i], moln->p[i].x, moln->p[i].y, moln->p[i].z); 
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
//printf("accept %ld\n",counter);
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
printf("counter = %ld chainid = %ld i_start = %ld oldenergy = %f newenergy = %f \n", counter, molm-mol, i_start, v_old.tot, v[ib].tot);
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

//   printf("counter=%ld\tNMols[out]=%ld\tNSites[out]=%ld\n", counter, NMols[out], NSites[out]);
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
            fprintf(foutput, "mol[%ld] site[%ld] parent site = %ld\n", moli-mol, i, moli->parent[i]);
}

//******************************//
//	Perform one MC move	//
//******************************//

#define NVFIX		10 	// first NVFIX moves with NO volume change
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
   return;
}

//******************************//
//	Monte Carlo cycles	//
//******************************//

void Cycle()			// one cycle contains one move per particle on average
{
   long		i, j, k, nmoves;
   double	arg=0.0, choice;
   double	ReCoor[MAXNSYSTEMS], oldReCoor[MAXNSYSTEMS];	// reaction coordinates

   // Test Area Begin
   time_t	start, end;   
   
/*   
   start	=	time(NULL);
   CL_Init();
   end		=	time(NULL);
   printf("cell_init() %lf\n", difftime(end, start));

   start	=	time(NULL);
   CL_Build();
   end		=	time(NULL);
   printf("cell_build() %lf\n", difftime(end, start));

   start	=	time(NULL);
   CalcV();
   end  	=	time(NULL);
   printf("calcv() %lf\n", difftime(end, start));
   
   start	=	time(NULL);
   for (i=0; i<NSITES; i++) {
      mcmove();
   }
   end  	=	time(NULL);
   printf("NSITES mcmoves %lf\n", difftime(end, start));
   
   start	=	time(NULL);
   mcvol();
   end  	=	time(NULL);
   printf("mcvol %lf\n", difftime(end, start));
   
   start	=	time(NULL);
   SampleSpherical();
   end  	=	time(NULL);
   printf("samplespherical %lf\n", difftime(end, start));
   
   start	=	time(NULL);
   Dist_Spherical();
   end  	=	time(NULL);
   printf("dist_sperical %lf\n", difftime(end, start));
   
   start	=	time(NULL);
   Sample_All();
   end  	=	time(NULL);
   printf("sampleall %lf\n", difftime(end, start));
   
   start	=	time(NULL);
   SampleP2All();
   end  	=	time(NULL);
   printf("sampleP2all %lf\n", difftime(end, start));
   
   start	=	time(NULL);
   Dist_p2();
   end  	=	time(NULL);
   printf("dist_p2 %lf\n", difftime(end, start));
   
   start	=	time(NULL);
   Find_Nuclei_p2(1);
   end  	=	time(NULL);
   printf("find nuclei %lf\n", difftime(end, start));
    
   start	=	time(NULL);
   CL_Destroy();
   end		=	time(NULL);
   printf("cell_destroy() %lf\n", difftime(end, start));
*/
   // Test Area End

   
   start	=	time(NULL);
   StoreMols();			// backup system and molecule info. before MC moves

   if (fabs(kP2) > ZERO) {	// if bias is turned on, backup r.c. value

      for (i=0; i<NSYSTEMS; i++) {
	 switch (dynvar) {					// different reaction coordinates
	    case 1:	oldReCoor[i]	=	nmax[i][0];	// nmax as r.c.
			break;
	    case 2:	oldReCoor[i]	=	Q6[i];		// Q6 as r.c.
			break;						
	    case 3:	oldReCoor[i]	=	P2M[i];		// P2_modified
			break;
	    case 4:	oldReCoor[i]	=	P2[i];		// P2
			break;
            case 5:	oldReCoor[i]	=	oldNSites[i]/oldBOX[i].vol;	// # density
			break;
	    case 6:	oldReCoor[i]	=	nmax[i][0];	// nmax_p2 as r.c.
			break;
	    default:	break;
         }
      }
   }
   end		=	time(NULL);
   printf("store old parameters %lf\n", difftime(end, start));

   /* Perform one sequence of MC moves */

   start	=	time(NULL);
   for (k=0; k<SEQUENCE; k++) {		// one sequence contains SEQUENCE MC cycles
      for (i=0; i<NSITES; i++) {
         NextMove();
      }
   }
   end		=	time(NULL);
   printf("one seq of nextmove() %lf\n", difftime(end, start));

   /* Basic sampling */

   // The real basic sampling is energy calculation which has been taken care of 
   // in NextMove().  Other samplings are basically for instant histogram output,
   // e.g., Ptrans and Pressure.

   start	=	time(NULL);
   SampleSpherical();
   Sample_All();
   end		=	time(NULL);
   printf("SampleSpherical and Sample_All %lf\n", difftime(end, start));

   /* Bias sampling */

   if (fabs(kP2) > ZERO) {		// if bias is turned on

      /* Sample the reaction coordinate */
   
   start	=	time(NULL);
      switch (dynvar) {
	 case 1:	Find_Nuclei(dynvar);	break;
	 //case 2:	Calc_Qlm(6);		break;		// calc. Q6
	 case 2: break;
	 case 3:	SampleP2All();		break;		// calc. global P2, P2M and local p2
	 case 4:	SampleP2();		break;
	 case 5:	break;
	 case 6:	SampleP2All();
			Find_Nuclei_p2(1);	break;
	 default:	break;
      }

      for (i=0; i<NSYSTEMS; i++) {
	 switch (dynvar) {
	    case 1:	ReCoor[i]	=	nmax[i][0];	break;
	    case 2:	ReCoor[i]	=	Q6[i];		break;
	    case 3:	ReCoor[i]	=	P2M[i];		break;
	    case 4:	ReCoor[i]	=	P2[i];		break;
	    case 5:	ReCoor[i]	=	NSites[i]/BOX[i].vol;	break;
	    case 6:	ReCoor[i]	=	nmax[i][0];		break;
	    default:	break;
	 }
      }
   end		=	time(NULL);
   printf("Sample biasing parameter %lf\n", difftime(end, start));

      /* Calculate bias potential */

      arg	=	(oldReCoor[0] - P2middle) * (oldReCoor[0] - P2middle) 
			- (ReCoor[0] - P2middle) * (ReCoor[0] - P2middle);
      arg	*=	(-0.5 * kP2);

      /*
      for (i=0; i<NSYSTEMS; i++)	// 5/1/08, allow only one nucleus in the system at any time
         if (dynvar==1)
	    if (sizedist[nmax[i][0]]>1 || nmax[i][1]>0)
	       arg	=	10000;
      */
   }
   else {				// bias is turned off
      arg	=	0.0;		// always accept sequence if NO bias
   }

   /* Determine acceptance of the sequence */

   if ((arg>0 ? (arg>75 ? 0 : (exp(-arg)<ran1(seed) ? 0 : 1)) : 1)) {
      av[0].acc_seq	++;
   }
   else {						// reject this sequence
      for (i=0; i<NSYSTEMS; i++) {
	 switch (dynvar) {
	    case 1:	nmax[i][0]	=	oldReCoor[i];		break;
	    case 2:	Q6[i]		=	oldReCoor[i];		break;
	    case 3:	P2M[i]		=	oldReCoor[i];		break;
	    case 4:	P2[i]		=	oldReCoor[i];		break;
	    case 5:	break;
	    case 6:	nmax[i][0]	=	oldReCoor[i];		break;
	    default:	break;
	 }
      }
      RestoreMols();			// recover system and molecules
#ifdef CELL_LIST
      CL_Destroy();
      CL_Init();
      CL_Build();
#endif
   }
   av[0].seq	++;
   //ParentCheck();

   if (!mod(counter, 1000))
      Adjust_Stepsize();			// Adjust MC move step sizes

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
