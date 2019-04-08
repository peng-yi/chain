/*
  program:      position.c
  author:       Peng Yi at MIT
  date:         October 22, 2006
  purpose:      Position calculation, Xtal nuclei identification
  references:	Stoddard, J. Comput. Phys. v27, 291 (1978)
*/

#define __POSITION_MODULE
#include "position.h"

#define SIZECAP	3

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

/************************************************************************/
/*	DistSQ(): calculate distance square between vector p and q	*/
/************************************************************************/
double DistSQ(vector p, vector q, long system)
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

//==============================================================//
//	MapInBox2(): If particle escapes the box, map it back, 	//
//	and pick up the image due to PBC so that the distance 	//
//	is the minimum.						//
//								//
// 	Notice: here we do not change *p, unlike in MapInBox()	//
//==============================================================//
vector MapInBox2(vector *p, long PBC, long system)	
{
   static vector	l;
   double		zlo, zhi, ylo, yhi, xlo, xhi, yoff, xoff;
   vector		A, B, C;

   if (PBC==1) {		//cubic periodic boundary condition
      l.x		=	p->x / BOX[system].lx;	
      l.y		=	p->y / BOX[system].ly;
      l.z		=	p->z / BOX[system].lz;
      //
      // Notice: here we do not change *p, unlike in MapInBox()
      
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

      /* bring z coordinate to be b/w zlo and zhi */
      while (l.z < zlo) { 
         l.z += C.z;	l.y += C.y;	l.x += C.x;
      }
      while (l.z >= zhi) {
         l.z -= C.z;	l.y -= C.y;	l.x -= C.x;
      }
      
      yoff	=	BOX[system].yz * (l.z-zlo) / (zhi-zlo);
      xoff	=	BOX[system].xz * (l.z-zlo) / (zhi-zlo);
      //yoff = 0.0;
      //xoff = 0.0;

      /* bring y coordinate to be b/w ylo+yoff and yhi+yoff */
      while (l.y < ylo+yoff) {
         l.y += B.y;	l.x += B.x;
      }
      while (l.y >= yhi+yoff) {
         l.y -= B.y;	l.x -= B.x;
      }
      xoff	+=	BOX[system].xy * (l.y - (ylo+yoff))/(yhi-ylo);
      //xoff = 0.0;

      /* bring x coordinate to the center box */
      while (l.x < xlo+xoff) {
         l.x += A.x;
      }
      while (l.x >= xhi+xoff) {
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

/********************************************************************************/
/*	unfold(): unfold all chains so that they are continuous in space	*/
/*	Added:	Aug.17, 2013							*/
/********************************************************************************/
void unfold()
{
   molstruct	*moli;
   long		i, system;
   vector	dr, comold, com;

   for (moli=mol; moli<mol+NMOLS; moli++) {
      system	=	moli->box;
      for (i=1; i<moli->nsites; i++) {
         dr	=	V_Subtr(moli->p+i, moli->p+i-1);
         dr	=	MapInBox2(&dr, PBC, system);
         moli->p[i]	=	V_Add(moli->p+i-1, &dr);
      }
      //MolInBox2(moli);		// com of chains in center box
   }
   return;
}

/********************************************************************************/
/*	unfoldchain(): unfold one chain so that it is continuous in space	*/
/*	anchor site is provided as input parameter				*/
/********************************************************************************/
void unfoldchain(molstruct *moli, long site)
{
   long		i, system=moli->box;
   vector	dr;

   for (i=site+1; i<moli->nsites; i++) {
      dr	=	V_Subtr(moli->p+i, moli->p+i-1);
      dr	=	MapInBox2(&dr, PBC, system);
      moli->p[i]	=	V_Add(moli->p+i-1, &dr);
   }
   for (i=site-1; i>=0; i--) {
      dr	=	V_Subtr(moli->p+i, moli->p+i+1);
      dr	=	MapInBox2(&dr, PBC, system);
      moli->p[i]	=	V_Add(moli->p+i+1, &dr);
   }
   return;
}

/****************************************************************/
/*	MapInNucleus(): Map a nucleus into the center box	*/
/* 	Added on 8/31/2012					*/
/****************************************************************/
void MapInNucleus(int system, beadstruct *nucleus, int nsites, vector *rbead, float rcutoff)
{
   // The idea is to search the connected neighbors (Rconn) within the nucleus
   // and to adjust the positions of these neighbors so that their distance
   // is the least among pbc images.  The procedure will continue through all
   // sites in the nucleus.

   // Previously (MapInNucleus1() below), I did not use connectivity criterion, 
   // rather I simply adjust the positions of the sites so that the distance b/w
   // any two sites in the nucleus is the least among pbc images.  However, problem
   // arises if the nucleus is large.  It is possible that the distance between 
   // a pair of sites should be greater than the least among pbc images if this
   // pair wants to stay "connected" through a path.  For example, if an extended 
   // segment along box boundary and of length 3/4 of the box dimension L is part of
   // a nucleus, then the distance b/w the head and end sites of this segment 
   // must be 3/4*L, rather than 1/4*L, if they want to be connected through the 
   // middle sites.
   //
   // (1/27/2017) Note: the coordinates are not necessarily back in the central box.

   int 		i, j, k, temp;

   //int	*L;
   //L = (int *)calloc(nsites, sizeof(int));
   static int	init=1;
   static int   *L;
   //int		L[MAXNMOLS*MAXNMOLSITES];
   vector	ri, rj, rk, dr;
   float	dr2;
   float	rc2 = rcutoff * rcutoff;

   if (nsites==0) {
      printf("# MapInNucleus() error: There is no atom in the nucleus.\n");
      fflush(stdout);
   }
   
   if (init) {
      init  =  0;
      L     = (int *) calloc(MAXNMOLS*MAXNMOLSITES, sizeof(int));
   }

   for (i=0; i<nsites; i++) {
      L[i]	=	i;
   }

   //***** Bring all sites connected (by rcutoff)	*****//
   //***** Same algorithm used in Find_nuclei_p2()	*****//
   //***** Added on 8/31/2012 				*****//	

   // Create connection list
    
   for (i=0; i<nsites-1; i++) {
      if (L[i] == i) {					// Not scanned 
         j	=	i;

         for (k=i+1; k<nsites; k++) {			// search through the whole list

	    if (L[k] == k) {				// also not scanned
               rj	=	nucleus[j].moli->p[nucleus[j].site];
               rk	=	nucleus[k].moli->p[nucleus[k].site];
	       dr2	=	DistSQ(rk, rj, system);
	       
               //if (dr2 < Rconn2) {			// if j and k are connected
               if (dr2 < rc2) {			// if j and k are connected
                  temp	=	L[k];
		  L[k]	=	L[j];			// save the last particle in L[j]
		  L[j]	=	temp;
	       }
	    }
	 }
 
         j	=	L[j];				// start from the last site on the link

         while (j!=i) {					// connection link not complete yet

            for (k=i+1; k<nsites; k++) {		// search through the whole list again

               if (L[k] == k) {
		  rj	=	nucleus[j].moli->p[nucleus[j].site];
		  rk	=	nucleus[k].moli->p[nucleus[k].site];
		  dr2	=	DistSQ(rk, rj, system);
		  
		  //if (dr2 < Rconn2) {		
		  if (dr2 < rc2) {		
	             temp	=	L[k];
		     L[k]	=	L[j];
		     L[j]	=	temp;
		  }
	       }
	    }

	    j	=	L[j];				// try another particle in the same cluster
	 } 
      }
   }
   
   // Adjust position and store in rbead[]

   i	=	0;
   rbead[0]	=	nucleus[0].moli->p[nucleus[0].site];
   
   while (L[i]!=0) {
      j		=	L[i];
      rbead[j]	=	nucleus[j].moli->p[nucleus[j].site];

      dr	=	V_Subtr(rbead+j, rbead+i);
      dr	=	MapInBox2(&dr, PBC, system);
      rbead[j]	=	V_Add(rbead+i, &dr);
      
      i		=	L[i];
   }

   // Shift the center of nucleus to the center box
   /*
   V_Null(&ri);
   for (i=0; i<nsites; i++) {
      ri	=	V_Add(&ri, rbead+i);
   }
   ri	=	V_Mult(1.0/nsites, &ri);
   rj	=	MapInBox2(&ri, PBC, system);
   dr	=	V_Subtr(&rj, &ri);
   for (i=0; i<nsites; i++) {
      rbead[i]	=	V_Add(rbead+i, &dr);
   }
   */
   return;
}
   
//////////////////////////////////////////
/* Map a nucleus into the center box 	*/
/* Added on 4/16/2012			*/
//////////////////////////////////////////

void MapInNucleus1(long system, beadstruct *nucleus, long nsites, vector *rbead)
{
   long		i;
   double	x, y, z;

   for (i=0; i<nsites; i++) {
      rbead[i]		=	nucleus[i].moli->p[nucleus[i].site];
      rbead[i].x	/=	BOX[system].lx;
      rbead[i].y	/=	BOX[system].ly;
      rbead[i].z	/=	BOX[system].lz;
   }
  
   // Move the first bead into the center box
 
   if (PBC==1) {
      if (rbead[0].x < -0.5 || rbead[0].x >=0.5) {
         rbead[0].x	-=	floor(rbead[0].x + 0.5);
      }
      if (rbead[0].y < -0.5 || rbead[0].y >=0.5) {
         rbead[0].y	-=	floor(rbead[0].y + 0.5);
      }
      if (rbead[0].z < -0.5 || rbead[0].z >=0.5) {
         rbead[0].z	-=	floor(rbead[0].z + 0.5);
      }
   }

   // Move all other beads close to the first bead

   for (i=1; i<nsites; i++) {
      x		=	rbead[i].x - rbead[0].x;
      y		=	rbead[i].y - rbead[0].y;
      z		=	rbead[i].z - rbead[0].z;

      if (PBC==1) {
         if (x<-0.5 || x>=0.5) {
            x	-=	floor(x+0.5);
            rbead[i].x	=	rbead[0].x + x;
         }
         if (y<-0.5 || y>=0.5) {
            y	-=	floor(y+0.5);
            rbead[i].y	=	rbead[0].y + y;
         }
         if (z<-0.5 || z>=0.5) {
            z	-=	floor(z+0.5);
            rbead[i].z	=	rbead[0].z + z;
         }
      }
   }

   for (i=0; i<nsites; i++) {
      rbead[i].x	*=	BOX[system].lx;
      rbead[i].y	*=	BOX[system].ly;
      rbead[i].z	*=	BOX[system].lz;
   }
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
    fprintf(foutput,"%ld\t%ld\n", part[i].nconnect, part[i].nuclid);
  }
  fprintf(foutput,"\n");
   
  for (i=0; i<NPARTS+1; i++) {
    if (sizedist[i]!=0) {
	fprintf(foutput, "Find_Nuclei, Number of nuclei of size %ld\t=\t%ld\n", i, sizedist[i]);
    }
  }
  fprintf(foutput, "Maxsize=\t%ld\n",MAXSIZE); 
  */
  return;
}
#endif

//////////////////////////////////////////////
/* Decide whether one site is in xtal phase */
//////////////////////////////////////////////
int crystal(molstruct *moli, long i)		// determine a xtal site
{
   int	result;
   
   if (samestr(moltype, "LJ")) {			// Lennard Jones system
      result	=	moli->nconn[i]	> critconnect;
   }
   else if (samestr(moltype, "monochain")) {	// monodisperse chain molecules
      result	=	moli->p2[i] > critp2;
   }
   return 	result;
}
/////____________________________________////


//**************************************//
//	Find nuclei in LJ system	//
//**************************************//

void Find_Nuclei_LJ()				// deal with multiple system
{
   molstruct		*moli;
   long			i, j, k, jj, n, NIT, LIT, Lk, nuclid, system,
			nsites[MAXNSYSTEMS];	// # of sites in each system
   double		r2;
   vector		pj;
#ifdef CELL_LIST
   cellstruct		*cellj, *cellk;
#endif
   static long		init=1;
   static long		**L;
   static sitestruct	**site;

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
/*
      sizedist	=	(long *) calloc(NSITES+1, sizeof(long));// not for all systems, later
      sizeofnucl=	(long *) calloc(NSITES+1, sizeof(long));
*/
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
      Nnucl[i]		=	0;
      Xtal[i]		=	0;
      nmax[i][0]	=	0;
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

            if (nmax[system][0] < NIT)
	       nmax[system][0]	=	NIT;
            
	    nuclid	++;
         }
      }

      /* calculate more nuclei variables */

      for (i=1; i<10; i++)				// calc. the size of 10 biggest nuclei
         nmax[system][i]	=	0;

      k	=	0;			
      for (i=nmax[system][0]; i>0; i--) {		// start from largest size
         for (j=1; j<=sizedist[i]; j++) {
            nmax[system][k]	=	i;
	    k	++;
            if (k>=10)	    break;
         }
         if (k>=10)         break;
      }
   }	// done for all systems
   return;
}
///______Find_Nuclei_LJ_______//////

///////////////////////////////////////////////////
/* Decide whether chains j and k are connected   */
/* based on Esselink's definition (JCPv101p9033) */
///////////////////////////////////////////////////

long connect2(long j, long k, long clusdef)
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
         printf("mol %ld and %ld are connected\n", j, k);
         printf("mol %ld's coordinates:\n", j);
         for (i=0; i<(mol+j)->nsites; i++)
            printf("\t%f\t%f\t%f\t%f\n", (mol+j)->p[i].x, (mol+j)->p[i].y, (mol+j)->p[i].z, type[(mol+j)->type[i]].M);
         printf("mol %ld's coordinates:\n", k);
         for (i=0; i<(mol+k)->nsites; i++)
            printf("\t%f\t%f\t%f\t%f\n", (mol+k)->p[i].x, (mol+k)->p[i].y, (mol+k)->p[i].z, type[(mol+k)->type[i]].M);

         printf("mol %ld's center of mass\n", j);
	 V_Print(pj); 
         printf("mol %ld's center of mass\n", k);
	 V_Print(pk); 

         printf("mol %ld's coordinates:\n", j);
         for (i=0; i<(mol+j)->nsites; i++)
            printf("\t%f\t%f\t%f\t%f\n", (mol+j)->p[i].x-pj.x, (mol+j)->p[i].y-pj.y, (mol+j)->p[i].z-pj.z, type[(mol+j)->type[i]].M);
         printf("mol %ld's coordinates:\n", k);
         for (i=0; i<(mol+k)->nsites; i++)
            printf("\t%f\t%f\t%f\t%f\n", (mol+k)->p[i].x-pk.x, (mol+k)->p[i].y-pk.y, (mol+k)->p[i].z-pk.z, type[(mol+k)->type[i]].M);

         printf("mol %ld's gyration matrix\n", j);
	 M_Print(Mj); 
         printf("mol %ld's gyration matrix\n", k);
	 M_Print(Mk); 
   
         printf("mol %ld's eig\n", j);
	 V_Print(eigj);
         printf("mol %ld's eig\n", k);
	 V_Print(eigk);

         printf("mol %ld's main axis\n", j);
         V_Print(vj);
         printf("mol %ld's main axis\n", k);
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
//////____connect2()____________________//////


//======================================================================//
//	Find_Nuclei(): Use Esselink JCP 1994 method to define crystal	//
//		       nuclei.  					//
//		       Added on 12/18/2007.				//
//		       For one system only.				//
//======================================================================//
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
/*
      sizeofnucl	=	(long *) calloc(nsites+1, sizeof(long));
      sizedist		=	(long *) calloc(nsites+1, sizeof(long));
*/
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

               if (connect2(j, k, clusdef)) {		// if j and k are connected
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
		  if (connect2(j, k, clusdef)) {
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
      nmax[i][0]	=	0;
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

	 if (nmax[system][0] < NIT) {
	    nmax[system][0]	=	NIT;
         }
	 nuclid		++;
      }
   }

   /* calculate more nuclei variables */

   for (i=1; i<10; i++)				// calc. the size of 10 biggest nuclei
      nmax[system][i]	=	0;

   k	=	0;			
   for (i=nmax[system][0]; i>0; i--) {		// start from largest size
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
   if (sizedist[nmax[system][0]] > 1)
      secondNmax[system]	=	nmax[system][0];
   else {
      for (i=nmax[system][0]-1; i>SIZECAP; i--) {
         if (sizedist[i] > 0) {
	          secondNmax[system]	=	i;
	          break;
         }
      }
   }
   return;
}
/////_____Find_Nuclei()______/////


//======================================================================//
//	connect_p2(): determine the connectivity between two beads. 	//
//		      Added on 2/23/2009.  Ref: Yi, JCP, 2011.		//
//		      For one system only.				//
//======================================================================//
int connect_p2(long j, long k, long clusdef, float rcutoff)
{
   int		result;
   long		sitej, sitek;
   long		sitespermol = NSITES/NMOLS;
   double	r2;
   vector	pj, pk;
   molstruct	*molj, *molk;
   float 	rc2 = rcutoff * rcutoff;

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
   pj		=	molj->p[sitej];
   pk		=	molk->p[sitek];
   //pj	=	mol[j/sitespermol].p[mod(j, sitespermol)];
   //pk	=	mol[k/sitespermol].p[mod(k, sitespermol)];

   r2		=	DistSQ(pj, pk, molj->box);
   //if (r2 < Rconn2) {
   if (r2 < rc2) {
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
      result	=	1;
   }
   else {
      result	=	0;
   }
   return	result;
}

//////////////////////////////////////////////////////////////////////
/* Add on 8/18/2011, find out all the segments, xtal and amorphous  */
//////////////////////////////////////////////////////////////////////

void Find_segments()
{
   molstruct	*moli;
   int 		i, j, k, nseg, previous;

   // Initialization

   for (i=0; i<MAXNMOLS; i++) {
      nsegment[i]	=	0;	// # of total segs on each chain
      for (j=0; j<3*MAXNMOLSITES; j++) {
         seg_stat[i][j]	=	-1;	// seg stat, each seg needs 3 elements:
					// -- head, id, and tail
      }
   }

   // Find all segments (xtal and amorphous) for all chains

   for (moli=mol; moli<mol+NMOLS; moli++) {
      j		=	moli-mol;
      previous	=	moli->nuclid[0];
      nseg	=	0;
      seg_stat[j][nseg*3]	=	0;

      for (i=0; i<moli->nsites; i++) {
         if (moli->nuclid[i]  != previous) {
            seg_stat[j][nseg*3+2]	=	i-1;		// id of last bead of segment
            seg_stat[j][nseg*3+1]	=	previous;	// segment id
            nseg	++;

            seg_stat[j][nseg*3]	=	i;	// head site of the segment

            if (i==moli->nsites-1) {		// head site is the last of the chain
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
   return;
}

//==============================================//
//	Xseg_tilt(): calculate xtal segment	//
//		     tilt vector		//
//		     added 5/15/13		//
//==============================================//
vector	Xseg_tilt(int nuclid)
{
   molstruct	*moli;
   int		j, k, n;
   int		length;
   int		flag;
   vector	r1A, r1B, r1V, rV, rtemp;
   
   //------Find a long xtal segment for reference------//
   flag	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      k	=	moli-mol;
      
      for (j=0; j<nsegment[k]; j++) {
	 length	=	seg_stat[k][j*3+2] - seg_stat[k][j*3] + 1;
	 
	 if (seg_stat[k][3*j+1] == nuclid && length >= 20) {	// long xseg in nuclid

	    r1A		=	moli->p[seg_stat[k][j*3]];	// head of xseg
	    r1B		=	moli->p[seg_stat[k][j*3+2]];	// tail of xseg
	    rtemp	=	V_Subtr(&r1A, &r1B);
	    flag	=	1;
	    break;
	 }
      }
      if (flag) {
	 break;
      }
   }

   //------Average xtal stem vectors------//
   V_Null(&rV);
   n	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      k	=	moli-mol;

      for (j=0; j<nsegment[k]; j++) {
	 if (seg_stat[k][3*j+1] == nuclid) {

	    r1A	=	moli->p[seg_stat[k][j*3]];		// head of xseg
	    r1B	=	moli->p[seg_stat[k][j*3+2]];		// tail of xseg
	    r1V	=	V_Subtr(&r1A, &r1B);

	    if (V_Dot(&rtemp, &r1V) <0) {
	       r1V	=	V_Mult(-1.0, &r1V);
	    }
	    rV	=	V_Add(&rV, &r1V);
	    n	++;
	 }
      }
      rV	=	V_Mult(1.0/n, &rV);
   }
   return	rV;
}

//======================================================//
//	Lamella_norm(): calculate vector normal to 	//
//			lamella plane			//
//======================================================//
vector Lamella_norm(int nuclid)
{
   molstruct	*moli;
   int		i, m, size;
   int		system = 0;
   float	arg;
   vector	rbead[MAXNMOLS*MAXNMOLSITES];
   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];		// group beads in the same nucleus
   
   matrix	tensor;
   matrix	evec;
   vector	eval;
   vector	norm;
   vector	ev[3];
   vector	rV;
   /* 
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
   size	=	m;
   
   MapInNucleus(system, nucleus, size, rbead);		// make nucleus beads continuous in space
   */
   /*
   tensor	=	grpgytensor(size, rbead); 
   eigensol(&tensor, &eval, &evec);			// solve eigenvalue problem
   
   ev[0].x	=	evec.x.x;			// evector of the smallest abs(evalue)
   ev[0].y	=	evec.x.y;
   ev[0].z	=	evec.x.z;

   ev[1].x	=	evec.y.x;
   ev[1].y	=	evec.y.y;
   ev[1].z	=	evec.y.z;

   ev[2].x	=	evec.z.x;			// evector of the largest abs(evalue)
   ev[2].y	=	evec.z.y;
   ev[2].z	=	evec.z.z;
   */
   ev[0].x	=	1.0;
   ev[0].y	=	0.0;
   ev[0].z	=	0.0;

   ev[1].x	=	0.0;
   ev[1].y	=	1.0;
   ev[1].z	=	0.0;

   ev[2].x	=	0.0;
   ev[2].y	=	0.0;
   ev[2].z	=	1.0;

   //------Xtal segment tilt vector------//
   rV	=	Xseg_tilt(nuclid);

   for (i=0; i<3; i++) {
      arg	=	V_Dot(ev+i, &rV)/sqrt(V_Dot(&rV, &rV));

      if (fabs(arg) > 0.707) {					// cos(45 degrees)
	 norm	=	ev[i];
	 break;
      }
   }
   return	norm;
}

//======================================================//
//	Rho_profile(): calculate density profile	//
//		       along lamella normal vector	//
//		       added 5/15/13			//
//======================================================//
void Rho_profile()
{
   molstruct	*moli;
   int		i, j, k;
   int		system = 0;
   int		nuclid;
   int		flag;
   vector	r, norm;
   float	max, min, value, temp, dl;
   int		bin;
   int		Nbin=10;
   float	rho[100];
   float	p2[100];

   for (i=0; i<Nbin; i++) {
      rho[i]	=	0.0;
      p2[i]	=	0.0;
   }
   
   //------find nuclid for the lamella------//
   flag	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 if (sizeofnucl[moli->nuclid[i]] == nmax[system][0]) {
	    nuclid	=	moli->nuclid[i];
	    flag	=	1;
	    break;
	 }
      }
      if (flag) {
	 break;
      }
   }

   //------unit vector normal to lamella plane------//
   norm	=	Lamella_norm(nuclid);
   temp	=	sqrt(V_Dot(&norm, &norm));
   norm	=	V_Mult(1.0/temp, &norm);

   //------decide range and binsize------//
   if (fabs(norm.x-1) < FLT_EPSILON) {
      max	=	0.5 * BOX[system].lx;
   }
   else if (fabs(norm.y-1) <FLT_EPSILON) {
      max	=	0.5 * BOX[system].ly;
   }
   else if (fabs(norm.z-1) <FLT_EPSILON) {
      max	=	0.5 * BOX[system].lz;
   }
   //max	=	MAX(MAX(BOX[system].lx, BOX[system].ly), BOX[system].lz);
   min	=      -max;

   dl	=	(max-min)/Nbin;
   
   //------collect density count------//
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 r	=	MapInBox2(moli->p+i, PBC, system);
	 value	=	V_Dot(&r, &norm);
	 bin	=	(int) ((value-min)/dl);
	 //bin	=	(int) (value-min)/dl;
	 rho[bin]	+=	1.0;
	 p2[bin]	+=	moli->p2[i];
      }
   }
   printf("# Density profile along lamella normal direction (%5.3f, %5.3f, %5.3f)\n", norm.x, norm.y, norm.z);
   printf("# max = %f min = %f Nbin = %d binsize = %f (system unit)\n", max, min, Nbin, dl);
   for (i=0; i<Nbin; i++) {
      printf("%f  %f  %f\n", (i+0.5)*dl+min, rho[i], p2[i]/rho[i]);
   }
   printf("\n");
   
   return; 
}
   
//======================================================//
//	Seg_type(): determine the type of segments	//
//		    added 5/9/13			//	
//======================================================//
void Seg_type(int nuclid)
{
   molstruct	*moli;
   char		segtype;
   int		i, j, k, n;
   int		m, size;
   int		head, tail;
   int		neigh;
   int		flag;
   int		system = 0;
   float	r2, arg;
   vector	r1A, r1B, r2A, r2B, r1O, r1V, r2O, r2V, rV, dr;
   vector	rtemp;
   
   vector	rbead[MAXNMOLS*MAXNMOLSITES];
   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];		// group beads in the same nucleus

   matrix	tensor;
   matrix	evec;
   vector	eval;
   vector	norm;
   vector	ev[3];
   float	tilt;

   int		length;
   int		ntail, nloop, nbrdg, npbrdg, nxseg;	// # of tails, etc
   float	ltail, lloop, lbrdg, lpbrdg, lxseg;	// average leng of tails, etc
   int		ndfct;
   float	ldfct;

   rV	=	Xseg_tilt(nuclid);
   norm	=	Lamella_norm(nuclid);
   tilt	=	V_Dot(&norm, &rV)/sqrt(V_Dot(&rV, &rV));
   
   //------Determine segment types------//
   
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 moli->segtype[i]	=	'0';
      }
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {		// loop over all chains
      k		=	moli-mol;
      system	=	moli->box;

      if (nsegment[k] >=2 ) {				// chain containing at least 1 xtal segment
	 for (j=0; j<nsegment[k]; j++) {		// loop all segments on one chain

	    if (seg_stat[k][j*3+1]<0) {			// amorphous segment

	       head	=	(j==0 ? -1 : seg_stat[k][(j-1)*3+1]);		// nuclei type on one end
	       tail	=	(j==nsegment[k]-1 ? -1 : seg_stat[k][(j+1)*3+1]);	// the other end

	       if ((head<0 && tail==nuclid) || (head==nuclid && tail<0)) {	// a tail segment
		  segtype		=	'2';
		  seg_stat[k][j*3+1]	=	-2;
	       }
	       else if (head==nuclid && tail==nuclid) {				// a loop or PBC bridge

		  r1A	=	moli->p[seg_stat[k][(j-1)*3]];			// head of xseg 1 on one end
		  r1B	=	moli->p[seg_stat[k][(j-1)*3+2]];		// tail of xseg 1 on one end
		  r2A	=	moli->p[seg_stat[k][(j+1)*3]];
		  r2B	=	moli->p[seg_stat[k][(j+1)*3+2]];

		  r1V	=	V_Subtr(&r1A, &r1B);				// vector of xseg 1
		  r2V	=	V_Subtr(&r2A, &r2B);				// vector of xseg 2

		  if (V_Dot(&r1V, &r2V) < 0) {
		     rtemp	=	r2B;
		     r2B	=	r2A;
		     r2A	=	rtemp;
		  }

		  r1O	=	V_Subtr(&r1A, &r2B);
		  r2O	=	V_Subtr(&r1B, &r2A);

		  arg	=	V_Dot(&r1O, &norm) * V_Dot(&r2O, &norm);

		  if ( arg<0) {
		     segtype		=	'3';				// loop
		     seg_stat[k][j*3+1]	=	-3;
		  }
		  else {
		     if (seg_stat[k][j*3+2]-seg_stat[k][j*3]+1 >= 8) {		// pbc bridge
			segtype			=	'5';
			seg_stat[k][j*3+1]	=	-5;
		     }
		     else {							// defect
			segtype			=	'7';
			seg_stat[k][j*3+1]	=	-7;
		     }
		  }
		  /*
		  r1A	=	moli->p[seg_stat[k][(j-1)*3]];			// head of xseg 1 on one end
		  r1B	=	moli->p[seg_stat[k][(j-1)*3+2]];		// tail of xseg 1 on one end
		  r2A	=	moli->p[seg_stat[k][(j+1)*3]];
		  r2B	=	moli->p[seg_stat[k][(j+1)*3+2]];
		  r1O	=	V_Add(&r1A, &r1B);				// center of xseg 1 on one end
		  r1O	=	V_Mult(0.5, &r1O); 
		  r2O	=	V_Add(&r2A, &r2B);				// center of xseg 2 on one end
		  r2O	=	V_Mult(0.5, &r2O); 
		  dr	=	V_Subtr(&r1O, &r2O); 
		  r1V	=	V_Subtr(&r1A, &r1B);				// vector of xseg 1
		  r2V	=	V_Subtr(&r2A, &r2B);				// vector of xseg 2
		  if (V_Dot(&r1V, &r2V) < 0) {
		     rV	=	V_Subtr(&r1V, &r2V);
		  }
		  else {
		     rV	=	V_Add(&r1V, &r2V);
		  }

                  arg	=	V_Dot(&dr, &rV)/sqrt( V_Dot(&rV, &rV)*V_Dot(&dr, &dr) );
		  //if ( fabs(arg)<0.173648 ) {					// 10 degrees
		  //if ( fabs(arg)<0.342020 ) {					// 20 degrees
		  //if ( fabs(arg)<0.64 ) {					// 40 degrees
		  if ( fabs(arg)<0.707) {					// 45 degrees
		     segtype		=	'3';				// loop
		     seg_stat[k][j*3+1]	=	-3;
		  }
		  else {
		     segtype		=	'5';				// pbc bridge
		     seg_stat[k][j*3+1]	=	-5;
		  }
		  */
		     
		  /*
		  for (m=0; m<size; m++) {
		     if (nucleus[m].moli == moli && nucleus[m].site == seg_stat[k][j*3]) {
			break;
		     }
		  }
		  rA	=	rbead[m];					// position of head anchor A
		  dr	=	V_Subtr(moli->p+seg_stat[k][j*3], &rA);
		  
		  rB	=	moli->p[seg_stat[k][j*3+2]];			// position of tail anchor B
		  rB	=	V_Subtr(&rB, &dr);

		  neigh	=	0;
		  for (m=0; m<size; m++) {
		     r2	=	(rB.x - rbead[m].x)*(rB.x - rbead[m].x)
		              + (rB.y - rbead[m].y)*(rB.y - rbead[m].y)
		     	      + (rB.z - rbead[m].z)*(rB.z - rbead[m].z);
		     if (r2 < Rconn2) {						// connected to one
			neigh	=	1;
			break;
		     }
		  }
		  if (neigh == 1) {
		     segtype		=	'3';				// loop
		     seg_stat[k][j*3+1]	=	-3;
		  }
		  else {
		     segtype		=	'5';				// pbc bridge
		     seg_stat[k][j*3+1]	=	-5;
		  }
		  */
	       }
	       else if ( (head==nuclid && tail!=nuclid) || (tail==nuclid && head!=nuclid)) {
		  segtype		=	'4';				// bridge
		  seg_stat[k][j*3+1]	=	-4;
	       }
	       else {							// not related to nuclid
		  //printf("Seg_Type() error: amorphous segment not identified.\n");
	       }
	    }
	    else if (seg_stat[k][j*3+1]==nuclid) {				// xtal segment
	       segtype			=	'6';
	    }
	    
	    // label beads
	    for (i=seg_stat[k][j*3]; i<=seg_stat[k][j*3+2]; i++) {
	       moli->segtype[i]	=	segtype;
	       //printf("%c %c\n", segtype, moli->segtype[i]);
	    }
	 }	// loop all segments
      }
      else {				// chain contains NO xtal segment
	 segtype	=	'1';
	 for (i=0; i<moli->nsites; i++) {
	    moli->segtype[i]	=	segtype;
	 }
      }
   }	// loop over all chains

   
   //collect segment statistics
   ntail	=	0;
   nloop	=	0;
   nbrdg	=	0;
   npbrdg	=	0;
   nxseg	=	0;
   ndfct	=	0;
   ltail	=	0.0;
   lloop	=	0.0;
   lbrdg	=	0.0;
   lpbrdg	=	0.0;
   lxseg	=	0.0;
   ldfct	=	0.0;

   for (i=0; i<NMOLS; i++) {
      for (j=0; j<nsegment[i]; j++) {
	 length	=	seg_stat[i][j*3+2] - seg_stat[i][j*3] + 1;

	 switch(seg_stat[i][j*3+1]) {
	    case	-2:	ntail	++;
	    			ltail	+=	length;
				break;
	    case	-3:	nloop	++;
	    			lloop	+=	length;
				break;
	    case	-4:	nbrdg	++;
	    			lbrdg	+=	length;
				break;
	    case	-5:	npbrdg	++;
	    			lpbrdg	+=	length;
				break;
	    case	-7:	ndfct	++;
	    			ldfct	+=	length;
				break;
	    default:		break;
	 }
	 if (seg_stat[i][j*3+1] > 0) {
	    nxseg	++;
	    lxseg	+=	length;
	 }
      }
   }
   printf("%4d %4d ", nmax[0][0], nmax[0][1]);
   printf("___ ");
   printf("%4d %4d %4d %4d %4d %4d ", ntail, nloop, nbrdg, npbrdg, nxseg, ndfct);
   printf("___ ");
   printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
       ltail/ntail, lloop/nloop, lbrdg/nbrdg, lpbrdg/npbrdg, lxseg/nxseg, ldfct/ndfct, tilt);	 
   
   return;
}

//======================================================//
//	Seg_length(): determine the length distribution //
//		      of segments related to nuclid	//
//		      added 5/15/13			//	
//======================================================//
void Seg_length(int nuclid, int lsegdist[][MAXNMOLSITES])
{
   int		i, j, k;
   int		length;
   int		ntail, nloop, nbrdg, npbrdg, nxseg;	// # of tails, etc
   float	ltail, lloop, lbrdg, lpbrdg, lxseg;	// average leng of tails, etc
   int		ndfct;
   float	ldfct;

   //collect segment statistics
   ntail	=	0;
   nloop	=	0;
   nbrdg	=	0;
   npbrdg	=	0;
   nxseg	=	0;
   ndfct	=	0;
   ltail	=	0.0;
   lloop	=	0.0;
   lbrdg	=	0.0;
   lpbrdg	=	0.0;
   lxseg	=	0.0;
   ldfct	=	0.0;

   for (i=0; i<NMOLS; i++) {
      for (j=0; j<nsegment[i]; j++) {
	 length	=	seg_stat[i][j*3+2] - seg_stat[i][j*3] + 1;

	 switch(seg_stat[i][j*3+1]) {
	    case	-2:	ntail	++;
	    			ltail	+=	length;
				lsegdist[2][length]	++;
				break;
	    case	-3:	nloop	++;
	    			lloop	+=	length;
				lsegdist[3][length]	++;
				break;
	    case	-4:	nbrdg	++;
	    			lbrdg	+=	length;
				lsegdist[4][length]	++;
				break;
	    case	-5:	npbrdg	++;
	    			lpbrdg	+=	length;
				lsegdist[5][length]	++;
				break;
	    case	-7:	ndfct	++;
	    			ldfct	+=	length;
				lsegdist[7][length]	++;
				break;
	    default:		break;
	 }
	 if (seg_stat[i][j*3+1] == nuclid) {
	    nxseg	++;
	    lxseg	+=	length;
	    lsegdist[6][length]	++;
	 }
      }
   }
   /*
   for (i=1; i<NSITES/NMOLS; i++) {
      printf("%d", i);
      for (j=2; j<8; j++) {
	 printf(" %d", lsegdist[j][i]);
      }
      printf("\n");
   }
   */
   printf("%4d %4d ", nmax[0][0], nmax[0][1]);
   printf("___ ");
   printf("%4d %4d %4d %4d %4d %4d ", ntail, nloop, nbrdg, npbrdg, nxseg, ndfct);
   printf("___ ");
   printf("%6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n",
       ltail/ntail, lloop/nloop, lbrdg/nbrdg, lpbrdg/npbrdg, lxseg/nxseg, ldfct/ndfct);	 
   
   return;
}

//*******************************//
// Deal with segment fluctuation //
//*******************************//

vector Seg_smooth(long nucleusid) {

   molstruct	*moli;
   long		i, j, k, nseg, nsegtot, nsegtotnew, previous;
   long		seg_start, seg_O, seg_A, seg_B; 
   long		length, len_A, len_B, typ_O, typ_A, typ_B;
   long		max, maxid, xtal, xtalnew;
   long		done, round;
   double	temp;
   static long	Lmin = 4;
   vector	result;

   max		=	0;			// look for chain with most segments
   maxid	=	0;
   nsegtot	=	0;
   for (j=0; j<NMOLS; j++) {
      if (nsegment[j] > max) {
         max	=	nsegment[j];
         maxid	=	j;
      }
      nsegtot	+=	nsegment[j];
   }

//   printf("max = %ld maxid = %ld\n", max, maxid);

   // Calculate xtal beads in chain with most segments.  Reason to do so it to see
   // -- if fluctuation reduction is even

   xtal	=	0;
   for (i=0; i<nsegment[maxid]; i++) {
      if (seg_stat[maxid][i*3+1] > 0)
         xtal	+=	seg_stat[maxid][i*3+2]-seg_stat[maxid][i*3]+1;
//      printf("%6d%6d%6d\n", seg_stat[maxid][i*3], seg_stat[maxid][i*3+1], seg_stat[maxid][i*3+2]);
   }
//   printf("total xtal segment = %ld\n", xtal);

   //***** Perform smoothing for all chains *****//
   for (j=0; j<NMOLS; j++) {
      done	=	0;			// Check smoothing result
      round	=	0;			// Rounds of smoothing
      seg_start	=	-1;

      while (done==0 && round<3) {		// 3 rounds or less
      //while (done==0 && nsegment[j]>1) {	
         //printf("B chain %ld nsegment %ld done = %ld\n", j, nsegment[j], done);

         for (i=0; i<nsegment[j]; i++) {	// Find a xtal segment at least Lmin long
            //printf("seg %ld %ld %ld %ld\n", i, seg_stat[j][i*3], seg_stat[j][i*3+1], seg_stat[j][i*3+2]);

            //if ( seg_stat[j][i*3+1]>0 && (seg_stat[j][i*3+2] - seg_stat[j][i*3]+1 >= Lmin) ) {
            if ( seg_stat[j][i*3+2] - seg_stat[j][i*3]+1 >= Lmin)  {
            //if ( seg_stat[j][i*3+1]>0 ) {
               seg_start	=	i;
               break;
            }
         }  	// There should be at least one such segment on any chain

         //***** Proceed toward the tail of the chain *****//
         seg_O	=	seg_start;
         while (seg_O < nsegment[j] - 1) {
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
                  len_B		=	seg_stat[j][seg_B*3+2]-seg_stat[j][seg_B*3]+1;
                  seg_stat[j][seg_A*3+2]	=	seg_stat[j][seg_A*3] + len_B-1;
                  seg_stat[j][seg_B*3]	=	seg_stat[j][seg_B*3+2] - len_A+1;

                  seg_stat[j][seg_A*3+1]	=	seg_stat[j][seg_B*3+1];
                  seg_stat[j][seg_B*3+1]	=	typ_A;
               }
            } 	// if O, A different and len_A < Lmin
            seg_O	++;
         }	// proceed toward the tail of the chain

         //***** Proceed toward the head of the chain *****//
         seg_O	=	seg_start;
         while (seg_O > 0) {
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
                  len_B		=	seg_stat[j][seg_B*3+2]-seg_stat[j][seg_B*3]+1;
                  seg_stat[j][seg_A*3]	=	seg_stat[j][seg_A*3+2] - len_B+1;                
                  seg_stat[j][seg_B*3+2]	=	seg_stat[j][seg_B*3] + len_A-1;

                  seg_stat[j][seg_A*3+1]	=	seg_stat[j][seg_B*3+1];
                  seg_stat[j][seg_B*3+1]	=	typ_A;
               }
            }
            seg_O	--;
         } // proceed toward the head of the chain

         //***** Resorting, basically combining same type of segments *****//
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

         //***** Check smoothing result *****//
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
         //printf("chain %ld nsegment %ld done = %ld\n", j, nsegment[j], done);
         round	++;
      }
      //printf("\n"); 
   }	// Smoothing for all chains ends

//   printf("After resorting\n");
//   printf("max = %ld maxid = %ld\n", max, maxid);

   //***** Report changes due to smoothing *****//

   // Difference in # of xtal beads on a chain of max # of segments due to smoothing

   xtalnew	=	0;
   for (i=0; i<nsegment[maxid]; i++) {
      if (seg_stat[maxid][i*3+1] > 0)
         xtalnew	+=	seg_stat[maxid][i*3+2]-seg_stat[maxid][i*3]+1;
         //printf("%6d%6d%6d\n", seg_stat[maxid][i*3], seg_stat[maxid][i*3+1], seg_stat[maxid][i*3+2]);
   }
   //printf("total xtal segment = %ld\n", xtalnew);
   result.x	=	xtalnew - xtal;	

   // Difference in size of the nucleus with nuclid as input due to smoothing

   xtalnew	=	0;
   for (j=0; j<NMOLS; j++) {
      for (i=0; i<nsegment[j]; i++) {
         if (seg_stat[j][i*3+1] == nucleusid ) {
            xtalnew	+=	seg_stat[j][i*3+2]-seg_stat[j][i*3]+1;
   }  }  }
   result.y	=	xtalnew - sizeofnucl[nucleusid];

   // Difference in total number of xtal beads due to smoothing

   xtalnew	=	0;
   for (j=0; j<NMOLS; j++) {
      for (i=0; i<nsegment[j]; i++) {
         if (seg_stat[j][i*3+1] > 0) {
            xtalnew	+=	seg_stat[j][i*3+2]-seg_stat[j][i*3]+1;
   }  }  }
   result.z	=	xtalnew - Xtal[0];

   // Difference in total number of segments due to smoothing

   nsegtotnew	=	0;
   for (j=0; j<NMOLS; j++) {
      nsegtotnew	+=	nsegment[j];
   }
   result.x	=	nsegtotnew - nsegtot;
 
   return result;
}

///////////////////////////////////////////////////////////////////
/* 8/20/2011. Calculate cylinder shape based on segment analysis */
///////////////////////////////////////////////////////////////////

vector cylinder(beadstruct *nucleus, long size, long nuclid)
{
   molstruct		*moli;
   long			i, j, k, m, index, exist, max;
   long			first, firstid, nseg, sizenew;
   vector		rAB, rpre, rtemp, shape, director;
   double		LAB, Ltot, thickness, radius, rho;

   static long		init=1;
   static molstruct	**molx;
   static vector	*rtot;

   // Initialization

   if (init) {
      init	=	0;
      molx	=	(molstruct **) calloc(NMOLS, sizeof(molstruct *));
      rtot	=	(vector *) calloc(NMOLS, sizeof(vector));
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
                  printf("cylinder error: two nearby segments pointing same direction. chainid=%ld segment=%ld\n", j, i);
                  for (m=0; m<nsegment[j]; m++) {
                      printf("%ld %ld %ld\n", seg_stat[j][m*3], seg_stat[j][m*3+1], seg_stat[j][m*3+2]);
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
   radius	=	((double)sizenew)/ (M_PI*thickness*rho);	//size??
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
   static float		*l;
   float		lave, lmax, radius, rho;
   static vector	*rAB;
   vector		rABave;
   static long		init=1;
   vector		shape;

   if (init) {
      molx	=	(molstruct **) calloc(NMOLS, sizeof(molstruct *));
      l		=	(float *) calloc(NMOLS, sizeof(float));
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

   rho		=	((float)NSITES) / BOX[0].vol;
   radius	=	((float)size)/ (M_PI*lave*rho);
   radius	=	sqrt(radius);
/*
  printf("nmax=%ld\n", nmax);
  printf("index=%ld\n", index);
  for (i=0; i<index; i++) { printf("l[%ld]=%f  A=%ld  B=%ld\n",i, l[i], A[i], B[i]); }
  printf("maxid=%ld\n", maxid);
  printf("rho=%f\n", rho);
  printf("lave=%f\n", lave);
  printf("radius=%f\n", radius); 
*/
   shape.x	=	(double)index;
   shape.y	=	lave;
   shape.z	=	radius;

   return	shape;
}

//**********************************************************************//
//	Crystal site smoothing using Multi-point average (10/1/12)	//
//**********************************************************************//

float xtal_smooth() {
   molstruct	*moli;
   int		i, j, oldxtal, newxtal;
   int		old[MAXNMOLSITES], new[MAXNMOLSITES];
   int		Lave = 7, Lhalf = (Lave-1)/2;			// 7-point average
   float	ave;
   
   // how many xtal beads before smoothing
   
   oldxtal	=	0;				
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 if (moli->p2[i] > critp2) {
	    oldxtal	++;
	 }
      }
   }

   // do multi-point smoothing

   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 old[i]	=	(moli->p2[i] >critp2 ? 1 : -1);		// 1 (crystal) or -1 (liquid)
      }
      
      // treating the middle portion of each chain
      for (i = Lhalf; i < moli->nsites - Lhalf; i++) {
	 ave	=	0;
	 
	 for (j = i - Lhalf; j <= i + Lhalf; j++) {
	    ave	+=	old[j];
	 }
	 new[i]	=	(ave > 0 ? 1 : -1);
      }
      
      // treating two ends of each chain
      for (i=0; i<Lhalf; i++) {
	 new[i]	=	new[Lhalf];
      }
      for (i=moli->nsites-1; i>=moli->nsites-Lhalf; i--) {
	 new[i]	=	new[moli->nsites-1-Lhalf];
      }

      // update moli->p2[i]
      for (i=0; i<moli->nsites; i++) {
	 moli->p2[i]	=	critp2 + new[i] * 0.1;		// not related to old value
      }
   }
   
   // how many xtal beads after smoothing
   
   newxtal	=	0;	
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 if (moli->p2[i] > critp2) {
	    newxtal	++;
	 }
      }
   }

   // return
   return	(1.0*newxtal)/oldxtal;
}

//======================================================================//
//	Find_Nuclei_p2(): Use local p2 to define crystal nuclei. 	//
//		          Added on 2/23/2009.  Ref: Yi, JCP, 2011.	//
//		          For one system only.				//
//======================================================================//
void Find_Nuclei_p2(long clusdef, float rcutoff)	// Stoddard, J. Comput. Phys. v27, 291 (1978)
{
   molstruct	*moli;
   long		i, j, k, m, n, jj, kk, NIT, LIT, Lk, nuclid;
   long		system = 0, nsites, sitespermol = NSITES/NMOLS, size, sum;
   float	r2, temp;
   vector	pj;
   vector	nreduce;
   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];		// group beads in the same nucleus

   static int	init=1;
   static int	*L;
   static int	*id;

   // Variables for shape analysis

   static FILE		*fPtr;
   vector		shape;
   vector		evalue;			// eigenvalue of gytensor
   matrix		gytensor;		// tensor of gyration of one nucleus
   static int		*number;		// count for average
   static float		*nchain;		// # of chains participate in one nucleus
   static float		*thickness, *radius;	// cylindrical model parameters
   static float		*asphericity;		// asphericity, calc. from evalue

   nsites	=	NSITES;			// definition based on beads

   if (init) {
      init	=	0;

      L		=	(int *) calloc(nsites, sizeof(int));	// only for one system now
      if (L==NULL)
         Exit("position", "Find_Nuclei", "out of memory");
/* 8/29/2012
      fPtr	=	fopen("shape.out", "w");
      fprintf(fPtr, "##### Nucleus shape analysis output from Find_Nuclei_p2 #####\n");
*/
      // Allocate cylindrical model variables
      number		=	(int *) calloc(NSITES, sizeof(int));
      nchain		=	(float *) calloc(NSITES, sizeof(float));
      thickness		=	(float *) calloc(NSITES, sizeof(float));
      radius		=	(float *) calloc(NSITES, sizeof(float));
      asphericity	=	(float *) calloc(NSITES, sizeof(float));

      id	=	(int *) calloc(NSITES, sizeof(int));
   }

   for (i=0; i<nsites; i++) {
      if (mol[i/sitespermol].p2[mod(i, sitespermol)] > critp2)
         L[i]		=	i;				// crystal-like particle
      else
	 L[i]		=	-1;				// liquid-like particle

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

               if (connect_p2(j, k, clusdef, rcutoff)) {		// if j and k are connected
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
		  if (connect_p2(j, k, clusdef, rcutoff)) {
		     L[k]	=	L[j];
		     L[j]	=	Lk;
		  }
	       }
	    }

	    j	=	L[j];				// try another particle in the same cluster
	 } 
      }
   }

   /***** Collect nuclei size distribution *****/

   for (i=0; i<NSYSTEMS; i++) {
      nmax[i][0]	=	0;
      Nnucl[i]		=	0;			// # of Xtal nuclei
      Xtal[i]		=	0;			// # of Xtal-like particles
   }
   for (i=0; i<nsites+1; i++) {		
      sizeofnucl[i]	=	0;			// initialize nuclei sizes
      sizedist[i]	=	0;			// initialize nuclei size distribution
   }

   nuclid	=	1;				// nuclei index starts from 1

   for (i=0; i<nsites; i++) {
      if (L[i] >= 0) {
	 NIT	=	1;
         LIT	=	L[i];
	 L[i]	=	-1 * (L[i]+2);	// clear this particle, but they can be easily recovered

//         for (n=0; n<mol[LIT].nsites; n++)
//            mol[LIT].nuclid[n]	=	nuclid;
         mol[LIT/sitespermol].nuclid[mod(LIT, sitespermol)]	=	nuclid;	// should be more general

         while (LIT != i) {
            NIT		+=	1;
            j		=	LIT;
            LIT		=	L[LIT];
            L[j]	=	-1 * (L[j]+2);
            
//            for (n=0; n<mol[LIT].nsites; n++)
//               mol[LIT].nuclid[n]	=	nuclid;
            mol[LIT/sitespermol].nuclid[mod(LIT, sitespermol)]	=	nuclid;
         }

         if (D_XTALSIZE)    
            PutInDistribution(D_Xtalsize+system, NIT, 1.0, 1.0);

	 Xtal[system]		+=	NIT;
//	 if (NIT>SIZECAP)
	    Nnucl[system]		++;

	 sizeofnucl[nuclid]	=	NIT;
	 sizedist[NIT]		++;

	 if (nmax[system][0] < NIT) {
	    nmax[system][0]	=	NIT;
         }
	 nuclid		++;
      }
   }

   /***** Calculate more nuclei variables *****/

   for (i=1; i<10; i++)				// nmax[system][1-9]
      nmax[system][i]	=	0;

   k	=	0;
   for (i=nmax[system][0]; i>0; i--) {		// start from largest size
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
   if (sizedist[nmax[system][0]] > 1)
      secondNmax[system]	=	nmax[system][0];
   else {
      for (i=nmax[system][0]-1; i>SIZECAP; i--) {
         if (sizedist[i] > 0) {
	          secondNmax[system]	=	i;
	          break;
         }
      }
   }

   /***** Sort nucleus id based on their size, from small to big *****/

   for (i=0; i<NSITES; i++) {
      id[i]	=	0;
   }
   for (i=1; i<=Nnucl[system]; i++) {		// search for every nucleus
      size	=	sizeofnucl[i];	// find its size

      sum	=	0;
      for (j=0; j<size; j++) {
         sum	+=	sizedist[j];	// how many nuclei are smaller than size
      }
      while (id[sum] !=0) {		// one same size nucleus already been found
         sum	++;
      }
      id[sum]	=	i;    
   }
/*
   printf("Nuclid size (in the order of increasing nucleus size)\n");
   for (i=0; i<Nnucl[system]; i++) {
      printf("%ld  %ld\n", id[i], sizeofnucl[id[i]]);
   }
*/

   /***** Store all beads that belong to the biggest nucleus in a list *****/
   /***** --If there are two nuclei both of nmax, then only store one  *****/

   //Find_segments();
   //nreduce = Seg_smooth(id[Nnucl[system]-1]);

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

      if (size != sizeofnucl[nuclid] ) {
         printf("Error, size %3ld != sizeofnucl[nuclid%ld] %3ld\n",size, nuclid, sizeofnucl[nuclid]);
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

/* 8/29/2012
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
*/
   }

/* 8/29/2012
   fprintf(fPtr, "\n********** Cylinder Model Analysis (%ld < nmax < %ld)**********\n", 
		(nmax[0][0]/100)*100, (nmax[0][0]/100+1)*100);
   fprintf(fPtr, "Rp = %f\tRconn = %f\tcritp2 = %f\n", Rp, Rconn, critp2);
   fprintf(fPtr, "nmax = %ld nreduce = %ld nmaxreduce = %ld\n", nmax[0][0], (long)nreduce.x, (long)nreduce.y);
   fprintf(fPtr, "size\t count\t P(n)\t nchains\t thickness\t radius\t asphericity\n");

   for (size=1; size<=nmax[0][0]; size++) {	// increase by one
      if (number[size]==0) {
         fprintf(fPtr, "%ld\t %ld\t %f\t nan\t nan\t nan\t nan\n", size, number[size], (double)number[size]/number[1]);
      }
      else {
         temp	=	1.0/number[size];
         fprintf(fPtr,"%ld\t %ld\t %f\t %f\t ", size, number[size], (double)number[size]/number[1], nchain[size]*temp);
         fprintf(fPtr,"%lf\t %lf\t %lf\n", thickness[size]*temp, radius[size]*temp, asphericity[size]*temp);
      }
   }
*/

/*
   fprintf(fPtr, "----------increase by three-----------\n");
   for (size=3; size<=nmax[0][0]; size+=3) {	// increase by three
      if (number[size]==0) {
         fprintf(fPtr, "%ld\t %ld\t nan\t nan\t nan\t nan\n", size, number[size]);
      }
      else {
         temp	=	1.0/number[size];
         fprintf(fPtr,"%ld\t %ld\t %f\t %f\t %f\t %f\n", size, number[size], nchain[size]*temp, thickness[size]*temp, radius[size]*temp, asphericity[size]*temp);
      }
   }
*/
/* 8/29/2012
   fflush(fPtr);
*/
/*
   nuclid	=	-1;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         if (sizeofnucl[moli->nuclid[i]] == nmax[0][0]) {	// find one bead in the biggest
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
}	// End of Find_Nuclei_p2()

//======================================================================//
//	connect_cnp(): determine the connectivity between two atoms	//
//======================================================================//
int connect_cnp(long j, long k, float rcutoff)
{
   molstruct	*molj, *molk;
   int		result;
   long		sitej, sitek;
   long		sitespermol = NSITES/NMOLS;
   double	r2;
   float	Rnn2;
   vector	pj, pk;
   //Rnn2 = 14.44;	// 3.8^2
   //Rnn2 = 33.64;	// 5.8 ^2
   //Rnn2 = 10.24;	// 3.2^2
   //Rnn2 = 12.25;	// 3.5^2
   //Rnn2 = 225.0;	// 15^2
   float        rc2 = rcutoff * rcutoff;

   if (mol[j/(NSITES/NMOLS)].box != mol[k/(NSITES/NMOLS)].box)
      Exit("position", "connect", "not in the same box");	// only one box for now

   molj		=	mol + j/sitespermol;
   molk		=	mol + k/sitespermol;
   sitej	=	j % sitespermol;
   sitek	=	k % sitespermol;
   pj		=	molj->p[sitej];
   pk		=	molk->p[sitek];
   //pj	=	mol[j/sitespermol].p[mod(j, sitespermol)];
   //pk	=	mol[k/sitespermol].p[mod(k, sitespermol)];

   r2		=	DistSQ(pj, pk, molj->box);
   if (r2 < rc2) {
      result	=	1;
   }
   else {
      result	=	0;
   }
   return	result;
}

/*---------------------------------------------------------------------------*/
/*	connect_1d(): if the two atoms have close z coordinates	 	     */
/*	added on: 9/22/18 						     */
/*---------------------------------------------------------------------------*/
int connect_1d(long j, long k, char *neighbormethod, float rcutoff)
{
   molstruct	*molj, *molk;
   int		result;
   long		sitej, sitek;
   long		sitespermol = NSITES/NMOLS;
   vector	pj, pk;
   double       r;

   if (mol[j/(NSITES/NMOLS)].box != mol[k/(NSITES/NMOLS)].box)
      Exit("position", "connect", "not in the same box");	// only one box for now

   molj		=	mol + j/sitespermol;
   molk		=	mol + k/sitespermol;

   sitej	=	j % sitespermol;
   sitek	=	k % sitespermol;

   pj		=	molj->p[sitej];
   pk		=	molk->p[sitek];

   if (strcmp(neighbormethod, "x")==0) {
      r = pj.x - pk.x;
   }
   else if (strcmp(neighbormethod, "y")==0) {
      r = pj.y - pk.y;
   }
   else if (strcmp(neighbormethod, "z")==0) {
      r = pj.z - pk.z;
   }
   else {
      printf("Error: No direction is defined in connect_1D()\n");
      exit(1);
   }

   if (fabs(r) < rcutoff) {
      result	=	1;
   }
   else {
      result	=	0;
   }
   return	result;
}


//======================================================================//
//	connect_general(): determine the connection between two atoms,	//
//			based on method specified in parameter list.	//
//			Added on 2/23/2009. 				//
//		      	For one system only.				//
//======================================================================//
int connect_general(long j, long k, char *method, char *neighbormethod, float rcutoff) 
{
   int		clusdef;

   if (!strcmp(neighbormethod, "r")) {
      if (strcmp(method, "p2")==0) {
	 clusdef	=	1;
      }
      else if (strcmp(method, "cnp")==0) {
	 clusdef	=	2;
      }
      else if (strcmp(method, "cna")==0) {
	 clusdef	=	3;
      }
      else if (strcmp(method, "cna2")==0) {
	 clusdef	=	32;
      }
      else if (strcmp(method, "nneigh")==0) {
	 clusdef	=	4;
      }
      else if (strcmp(method, "custom")==0) {
	 clusdef	=	5;
      }
   }
   else if (!strcmp(neighbormethod, "x") || 
            !strcmp(neighbormethod, "y") ||
	    !strcmp(neighbormethod, "z") ) {
      clusdef 	= 	9;
   }
   else {
      printf("connect_general() error: cluster definition missing.\n");
      exit(0);
   }

   switch(clusdef) {
      case	1:	return	connect_p2(j, k, clusdef, rcutoff);
      case	2:	return	connect_cnp(j, k, rcutoff);
      case	3:	return	connect_cnp(j, k, rcutoff);	// same as cnp
      case	4:	return	connect_cnp(j, k, rcutoff);	// same as cnp
      case	5:	return	connect_cnp(j, k, rcutoff);	// same as cnp
      case 	9:	return  connect_1d(j, k, neighbormethod, rcutoff);
      case  	32:	return	connect_cnp(j, k, rcutoff);	// same as cnp
   }
}

//======================================================================//
//	find_nuclei_general(): find nuclei based on method specified . 	//
//		          in the parameter list.  Added on 3/10/2014. 	//
//		          For one system only.				//
//======================================================================//
void find_nuclei_general(char *method, char *modifier, char *neighbormethod, float rcutoff)	// Stoddard, J. Comput. Phys. v27, 291 (1978)
{
   molstruct	*moli;
   long		i, j, k, m, n, jj, kk, NIT, LIT, Lk, nuclid;
   long		system = 0, nsites, sitespermol = NSITES/NMOLS, size, sum;
   float	r2, temp;
   vector	pj;
   vector	nreduce;
   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];		// group beads in the same nucleus

   int		clusdef;
   int		itmp, itmp1;
   double	crit_p2 = critp2;
   double	crit_cnp = 24.0;
   int		crit_nneigh = 10;

   static int	init=1;
   static int	*L;
   static int	*id;

   // Variables for shape analysis

   static FILE		*fPtr;
   vector		shape;
   vector		evalue;			// eigenvalue of gytensor
   matrix		gytensor;		// tensor of gyration of one nucleus
   static int		*number;		// count for average
   static float		*nchain;		// # of chains participate in one nucleus
   static float		*thickness, *radius;	// cylindrical model parameters
   static float		*asphericity;		// asphericity, calc. from evalue

   nsites	=	NSITES;			// definition based on beads

   if (init) {
      init	=	0;

      L		=	(int *) calloc(nsites, sizeof(int));	// only for one system now
      if (L==NULL)
         Exit("position", "find_nuclei_general", "out of memory");

      // Allocate cylindrical model variables
      number		=	(int *) calloc(NSITES, sizeof(int));
      nchain		=	(float *) calloc(NSITES, sizeof(float));
      thickness		=	(float *) calloc(NSITES, sizeof(float));
      radius		=	(float *) calloc(NSITES, sizeof(float));
      asphericity	=	(float *) calloc(NSITES, sizeof(float));

      id	=	(int *) calloc(nsites, sizeof(int));
   }

   if (strcmp(method, "p2")==0) {
      clusdef	=	1;
   }
   else if (strcmp(method, "cnp")==0) {
      clusdef	=	2;
   }
   else if (strcmp(method, "cna")==0) {
      clusdef	=	3;
   }
   else if (strcmp(method, "cna2")==0) {
      clusdef	=	32;
   }
   else if (strcmp(method, "nneigh")==0) {
      clusdef	=	4;
   }
   else if (strcmp(method, "custom")==0) {
      clusdef	=	5;
   }
   else {
      printf("find_nuclei_general() error: cluster definition missing.\n");
      fflush(stdout);
      exit(0);
   }

   for (i=0; i<nsites; i++) {

      // identify particles to count
      switch (clusdef) {
	 case	1: 
	    itmp	=	mol[i/sitespermol].p2[mod(i, sitespermol)] > crit_p2;
	    break;

	 case	2:
	    itmp	=	(mol[i/sitespermol].cnp[i%sitespermol] > crit_cnp) &&
	    			(mol[i/sitespermol].type[i%sitespermol] ==3);
	    break;

	 case	3:
	    if (!strcmp(modifier,"nonhcp")) {
	    	itmp	=	(mol[i/sitespermol].cna[i%sitespermol] !=2) &&
	    			(mol[i/sitespermol].type[i%sitespermol] !=1 &&
	    			mol[i/sitespermol].type[i%sitespermol] !=2);
	    } 
	    else if (!strcmp(modifier,"partial")) {	
	    	itmp	=	(mol[i/sitespermol].cna[i%sitespermol] ==5) &&
	    			(mol[i/sitespermol].type[i%sitespermol] !=1 &&
	    			mol[i/sitespermol].type[i%sitespermol] !=2);
	    }
	    else if (!strcmp(modifier,"fcc")) {	
	    	itmp	=	(mol[i/sitespermol].cna[i%sitespermol] ==1) &&
	    			(mol[i/sitespermol].type[i%sitespermol] !=1 &&
	    			mol[i/sitespermol].type[i%sitespermol] !=2);
	    }
	    else if (!strcmp(modifier, "hcp")) {
	    	itmp	=	(mol[i/sitespermol].cna[i%sitespermol] == 2);
	    }
	    else if (!strcmp(modifier, "bcc")) {
	    	itmp	=	(mol[i/sitespermol].cna[i%sitespermol] == 3);
	    }
	    else if (!strcmp(modifier, "fcc/hcp")) {
	    	itmp	=	(mol[i/sitespermol].cna[i%sitespermol] == 1) ||
                                (mol[i/sitespermol].cna[i%sitespermol] == 2) ;
	    }
	    else if (!strcmp(modifier, "unknown")) {
	    	itmp	=	(mol[i/sitespermol].cna[i%sitespermol] == 5);
	    }
	    break;

         case  32:
	    if (!strcmp(modifier, "fcc/hcp")) {	// fcc or hcp
	    	itmp	=	(mol[i/sitespermol].cna2[i%sitespermol] == 1) ||
                                (mol[i/sitespermol].cna2[i%sitespermol] == 2) ;
	    }
	    else if (!strcmp(modifier, "bcc")) {
	    	itmp	=	(mol[i/sitespermol].cna2[i%sitespermol] ==3);
	    }
	    else if (!strcmp(modifier, "fcc")) {
	    	itmp	=	(mol[i/sitespermol].cna2[i%sitespermol] == 1);
	    }
	    else if (!strcmp(modifier, "unknown")) {
	    	itmp	=	(mol[i/sitespermol].cna2[i%sitespermol] == 5);
	    }
	    break;

	 case	4:
	    itmp1	=	atoi(modifier);
	    itmp	=	(mol[i/sitespermol].nneigh[i%sitespermol] <= itmp1) &&
	    			(mol[i/sitespermol].tmp[i%sitespermol]==0) &&
	    			(mol[i/sitespermol].type[i%sitespermol] !=1 &&
	    			mol[i/sitespermol].type[i%sitespermol] !=2);
	    break;
	 case	5: 
	    itmp	=	mol[i/sitespermol].tmp[mod(i, sitespermol)];
	    break;
	 default: break;
      }

      // set counting array
      if (itmp) {
	 L[i]	=	i;
      }
       else {
         L[i]	=	-1;
      }

      // initialize nucleus id associated with each atom
      mol[i/sitespermol].nuclid[mod(i,sitespermol)]	=	-1;
   }

   for (i=0; i<nsites-1; i++) {
      if (L[i] == i) {					// Not scanned crystal particle
         j	=	i;

         for (k=i+1; k<nsites; k++) {			// search through the whole list

            Lk	=	L[k];

	    if (Lk == k) {				// also not scanned crystal particle
               if (connect_general(j, k, method, neighbormethod, rcutoff)) {	// if j and k are connected
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
		  if (connect_general(j, k, method, neighbormethod, rcutoff)) {
		     L[k]	=	L[j];
		     L[j]	=	Lk;
		  }
	       }
	    }

	    j	=	L[j];				// try another particle in the same cluster
	 } 
      }
   }

   /***** Collect nuclei size distribution *****/

   for (i=0; i<NSYSTEMS; i++) {
      nmax[i][0]	=	0;
      Nnucl[i]		=	0;			// # of Xtal nuclei
      Xtal[i]		=	0;			// # of Xtal-like particles
   }
   for (i=0; i<nsites+1; i++) {		
      sizeofnucl[i]	=	0;			// initialize nuclei sizes
      sizedist[i]	=	0;			// initialize nuclei size distribution
   }

   nuclid	=	1;				// nuclei index starts from 1

   for (i=0; i<nsites; i++) {
      if (L[i] >= 0) {
	 NIT	=	1;
         LIT	=	L[i];
	 L[i]	=	-1 * (L[i]+2);	// clear this particle, but they can be easily recovered

//         for (n=0; n<mol[LIT].nsites; n++)
//            mol[LIT].nuclid[n]	=	nuclid;
         mol[LIT/sitespermol].nuclid[mod(LIT, sitespermol)]	=	nuclid;	// should be more general

         while (LIT != i) {
            NIT		+=	1;
            j		=	LIT;
            LIT		=	L[LIT];
            L[j]	=	-1 * (L[j]+2);
            
//            for (n=0; n<mol[LIT].nsites; n++)
//               mol[LIT].nuclid[n]	=	nuclid;
            mol[LIT/sitespermol].nuclid[mod(LIT, sitespermol)]	=	nuclid;
         }

         if (D_XTALSIZE)    
            PutInDistribution(D_Xtalsize+system, NIT, 1.0, 1.0);

	 Xtal[system]		+=	NIT;
//	 if (NIT>SIZECAP)
	    Nnucl[system]		++;

	 sizeofnucl[nuclid]	=	NIT;
	 sizedist[NIT]		++;

	 if (nmax[system][0] < NIT) {
	    nmax[system][0]	=	NIT;
         }
	 nuclid		++;
      }
   }

   /***** Calculate more nuclei variables *****/

   for (i=1; i<10; i++)				// nmax[system][1-9]
      nmax[system][i]	=	0;

   k	=	0;
   for (i=nmax[system][0]; i>0; i--) {		// start from largest size
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
   realNnucl[system]  	=	Nnucl[system];
   for (i=1; i<=SIZECAP; i++) {
      realXtal[system]	-=	sizedist[i] * i;
      realNnucl[system] -=	sizedist[i];
   }

   secondNmax[system]	=	0;
   if (sizedist[nmax[system][0]] > 1)
      secondNmax[system]	=	nmax[system][0];
   else {
      for (i=nmax[system][0]-1; i>SIZECAP; i--) {
         if (sizedist[i] > 0) {
	          secondNmax[system]	=	i;
	          break;
         }
      }
   }

   /***** Sort nucleus id based on their size, from small to big *****/

   for (i=0; i<NSITES; i++) {
      id[i]	=	0;
   }
   for (i=1; i<=Nnucl[system]; i++) {		// search for every nucleus
      size	=	sizeofnucl[i];	// find its size

      sum	=	0;
      for (j=0; j<size; j++) {
         sum	+=	sizedist[j];	// how many nuclei are smaller than size
      }
      while (id[sum] !=0) {		// one same size nucleus already been found
         sum	++;
      }
      id[sum]	=	i;    
   }
/*
   printf("Nuclid size (in the order of increasing nucleus size)\n");
   for (i=0; i<Nnucl[system]; i++) {
      printf("%ld  %ld\n", id[i], sizeofnucl[id[i]]);
   }
*/

   /***** Store all beads that belong to the biggest nucleus in a list *****/
   /***** --If there are two nuclei both of nmax, then only store one  *****/

   //Find_segments();
   //nreduce = Seg_smooth(id[Nnucl[system]-1]);

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

      if (size != sizeofnucl[nuclid] ) {
         printf("Error, size %3ld != sizeofnucl[nuclid%ld] %3ld\n",size, nuclid, sizeofnucl[nuclid]);
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

/* 8/29/2012
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
*/
   }

/* 8/29/2012
   fprintf(fPtr, "\n********** Cylinder Model Analysis (%ld < nmax < %ld)**********\n", 
		(nmax[0][0]/100)*100, (nmax[0][0]/100+1)*100);
   fprintf(fPtr, "Rp = %f\tRconn = %f\tcritp2 = %f\n", Rp, Rconn, critp2);
   fprintf(fPtr, "nmax = %ld nreduce = %ld nmaxreduce = %ld\n", nmax[0][0], (long)nreduce.x, (long)nreduce.y);
   fprintf(fPtr, "size\t count\t P(n)\t nchains\t thickness\t radius\t asphericity\n");

   for (size=1; size<=nmax[0][0]; size++) {	// increase by one
      if (number[size]==0) {
         fprintf(fPtr, "%ld\t %ld\t %f\t nan\t nan\t nan\t nan\n", size, number[size], (double)number[size]/number[1]);
      }
      else {
         temp	=	1.0/number[size];
         fprintf(fPtr,"%ld\t %ld\t %f\t %f\t ", size, number[size], (double)number[size]/number[1], nchain[size]*temp);
         fprintf(fPtr,"%lf\t %lf\t %lf\n", thickness[size]*temp, radius[size]*temp, asphericity[size]*temp);
      }
   }
*/

/*
   fprintf(fPtr, "----------increase by three-----------\n");
   for (size=3; size<=nmax[0][0]; size+=3) {	// increase by three
      if (number[size]==0) {
         fprintf(fPtr, "%ld\t %ld\t nan\t nan\t nan\t nan\n", size, number[size]);
      }
      else {
         temp	=	1.0/number[size];
         fprintf(fPtr,"%ld\t %ld\t %f\t %f\t %f\t %f\n", size, number[size], nchain[size]*temp, thickness[size]*temp, radius[size]*temp, asphericity[size]*temp);
      }
   }
*/
/* 8/29/2012
   fflush(fPtr);
*/
/*
   nuclid	=	-1;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
         if (sizeofnucl[moli->nuclid[i]] == nmax[0][0]) {	// find one bead in the biggest
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
}	// End of find_nuclei_general()


//==============================================================//
//	find_nuclei_lmp(): 					//
//		find nuclei based on LAMMPS clustering results	//
//		stored in moli->clus[i]				//
//	 	Added on 2/2/2018. 				//
//	        For one system only.				//
//==============================================================//
void find_nuclei_lmp()
{
   molstruct *moli;
   int 	i, j, k;
   int  nsites = 0;
   static int init = 1;
   static int *size, *id;

   // Initialization
   //
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 nsites ++;
      }
   }
   if (init) {
      init = 0;
      size = (int *) calloc(nsites+1, sizeof(int));
      id   = (int *) calloc(nsites+1, sizeof(int));
   }

   for (i=0; i<MAXNMOLS*MAXNMOLSITES; i++) {
      sizeofnucl[i] = 0;
      sizedist[i]   = 0;
   }
   for (i=0; i<NSYSTEMS; i++) {
      nmax[i][0] = 0;
      Nnucl[i]	 = 0;			// # of Xtal nuclei
      Xtal[i]	 = 0;			// # of Xtal-like particles
   }

   // Find the greatest index of clusters given by LAMMPS
   //
   int  N = 0; 		               	// largest cluster id assigned by lammps
   					// but doesn't mean there are N clusters
					
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 sizeofnucl[moli->clus[i]] ++;
	 if (moli->clus[i] > N) {
	    N = moli->clus[i];
	 }
      }
   }
   
   // Rearrange sizeofnucl
   //
   int realid = 0;

   for (i=1; i<=N; i++) {
      if (sizeofnucl[i] > 0) {		// a real nucleus
	 realid ++;			// realid starts from 1
	 id[realid] = i;
	 size[realid] = sizeofnucl[i];
      }
   }

   // Collect clusters statistics
   //
   int Nreal = realid;

   for (i=1; i<=Nreal; i++) {
      sizeofnucl[i] = size[i];
   }

   int total = 0;
   int max = 0, maxid = 0;

   for (i=1; i<=Nreal; i++) {
      if (sizeofnucl[i] > 0) {		// a real nucleus
	 total += sizeofnucl[i];
      }
      if (sizeofnucl[i] > max) {
	 max = sizeofnucl[i];
	 maxid = i;
      }
   }

   int secondmax = 0;

   for (i=1; i<=Nreal; i++) {
      if (sizeofnucl[i] > secondmax && i != maxid) {
	 secondmax = sizeofnucl[i];
      }
   }

   int system = 0;

   Nnucl[system] = Nreal;
   Xtal[system] = total;
   nmax[system][0] = max;
   nmax[system][1] = secondmax;

   // Update moli->nuclid[i]
   //
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 moli->nuclid[i] = -1;

	 for (j=1; j<=Nreal; j++) {
	    if (id[j] == moli->clus[i]) {
	       moli->nuclid[i] = j;
	       break;
	    }
	 }
      }
   }

   return;
}

//////////////////////////////////////////////////////////////////
/* Added on 10/7/09, find connections based on chord vectors	*/
//////////////////////////////////////////////////////////////////

long	connect_p2_new(long j, long k, long clusdef, float rcutoff)
{
   vector	pj, pk;
   double	r2, costheta;
   long		sitespermol = NSITES/NMOLS;
   float	rc2 = rcutoff * rcutoff;

   if (mol[j/(NSITES/NMOLS)].box != mol[k/(NSITES/NMOLS)].box)
      Exit("position", "connect", "not in the same box");	// only one box for now

   pj	=	mol[j/sitespermol].p[mod(j, sitespermol)];
   pk	=	mol[k/sitespermol].p[mod(k, sitespermol)];

   r2	=	DistSQ(pj, pk, mol[j/sitespermol].box);
    
   //if (r2 < Rconn2)
   if (r2 < rc2)
      return	1;
   else
      return	0;
}

/*
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

      if (L==NULL)
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
*/

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
      if (moli->box == system && sizeofnucl[moli->nuclid[0]] == nmax[system][0]) {
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
