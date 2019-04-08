/*
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
   long		i1, i2, i3;
   vector	n1, d1, d2, d3 = {1.0, 0.0, 0.0};
   sphere	s;
   matrix	A;
   long		system	=	molm->box;		// added by PY

   if ((i1 = molm->parent[i0])>=0) {
      d1	=	V_Subtr(molm->p+i0, molm->p+i1);
      d1	=	MapInBox2(&d1, PBC, system);	// added by PY

      s.d                 = sqrt(V_Dot(&d1, &d1));

      if ((i2 = molm->parent[i1])>=0) {
         d2	=	V_Subtr(molm->p+i1, molm->p+i2);
      	 d2	=	MapInBox2(&d2, PBC, system);	// added by PY

         s.alpha	=	acos((d1.x*d2.x+d1.y*d2.y+d1.z*d2.z)/
                                s.d/sqrt(d2.x*d2.x+d2.y*d2.y+d2.z*d2.z));
         n1	=	V_Cross(&d2, &d1);
         n1	=	V_Mult(1.0/sqrt(V_Dot(&n1, &n1)), &n1);

         if ((i3 = molm->parent[i2])>=0) {
            d3	=	V_Subtr(molm->p+i2, molm->p+i3);
            d3	=	MapInBox2(&d3, PBC, system);	// added by PY
         }
         A	=	M_Orientation(&d3, &d2);
         d1.y	=	A.y.x*n1.x+A.y.y*n1.y+A.y.z*n1.z;
         d1.z	=	A.z.x*n1.x+A.z.y*n1.y+A.z.z*n1.z;
         s.beta	=	atan2(-d1.y, d1.z);
         return	s;
      }
      s.alpha	=	acos(d1.x/s.d);
      s.beta	=	atan2(d1.z, d1.y);
      return	s;
   }
   s.d		=	sqrt(V_Dot(molm->p+i0, molm->p+i0));
   s.alpha	=	acos(molm->p[i0].x/s.d);
   s.beta	=	atan2(molm->p[i0].z, molm->p[i0].y);
   return	s;
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
   d.x			=	scale.x - 1.0;
   d.y			=	scale.y - 1.0;
   d.z			=	scale.z - 1.0;

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
   return;
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


