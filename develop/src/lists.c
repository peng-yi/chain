/*
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
         printf("%ld-->", currentPtr->neighbor);
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
   l	=	V_Mult(0.999999, &l);		// make sure atoms do not sit on the upper bound
						// this will not affect the real coordinate
						// it only helps cell assignment
   celli	=	Cell;
   while ( (celli->box != ib) && celli < Cell+NCELLS )
      celli	++;

   celli	+=	((int) ((l.x / BOX[ib].lx +0.5)* M[ib].x))
		+	((int) ((l.y / BOX[ib].ly +0.5)* M[ib].y)) * M[ib].x
		+	((int) ((l.z / BOX[ib].lz +0.5)* M[ib].z)) * M[ib].x * M[ib].y;

   if (celli-Cell >= NCELLS) {
      printf("site coordinates: %10.8f\t%10.8f\t%10.8f\n", l.x, l.y, l.z);
      printf("box dimension: %10.8f\t%10.8f\t%10.8f\n", BOX[ib].lx, BOX[ib].ly, BOX[ib].lz);
      printf("(Mx, My, Mz) = (%d, %d, %d)\n", M[ib].x, M[ib].y, M[ib].z);
      printf("NCELLS=%ld\t, celli-Cell=%ld\n", NCELLS, celli-Cell);
      if (l.x/BOX[ib].lx >= 0.5)     printf("lx/Bx=%f\n", l.x/BOX[ib].lx);
      if (l.y/BOX[ib].ly >= 0.5)     printf("ly/By=%f\n", l.y/BOX[ib].ly);
      if (l.z/BOX[ib].lz >= 0.5)     printf("lz/Bz=%f\n", l.z/BOX[ib].lz);
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
   return;
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
