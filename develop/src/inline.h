/*
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
