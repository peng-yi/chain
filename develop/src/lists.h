/*
    program:    lists.h  
    author:     Peng Yi at MIT
    date:       October 22, 2006
    purpose:    header file for lists.c
*/
#ifndef __LISTS_HEADER
#define __LISTS_HEADER

#include "header.h"

#define MAXNNEIGHBORS	500

typedef struct {
   long			n;			// # of neighbors
   long			reverse[MAXNNEIGHBORS];	// whether reverse or not
   long			site[MAXNNEIGHBORS];	// id of site
   vector		dr[MAXNNEIGHBORS];	// distance to that neighbor
   molstruct		*mol[MAXNNEIGHBORS];	// id of mol
} neighborlist;					// neighbor list of one end bead

#ifdef __LISTS_MODULE


#ifdef CELL_LIST
cellstruct	*Cell;
ivector		M[MAXNBOX];		// number of cells along box length
long		NCELLS;
long		neighcellplus[54];	// neighbor cells before a move plus 
					// those after a move (no double count)
long		nneighcellplus;
#endif	/* CELL_LIST */

#ifdef VERLET_LIST
long		*vlistplus;		// the Verlet neighbors before a move PLUS 
					// those after a move (no double count)
					// PLUS itself
long		nverletplus;		// number of the vlistplus elements 
long		maxnverlet;		// maximum # of verlet neighbors of any particle
#endif	/* VERLET_LIST */


#else

extern long	NeighborList(molstruct *molm, neighborlist *list);
extern long	DB_NeighborList(molstruct *molm, long sitem, neighborlist *list);
extern long	IDR_NeighborList(molstruct *molm, long sitem, neighborlist *list);

#ifdef CELL_LIST
extern cellstruct	*Cell;
extern ivector	M[MAXNBOX];
extern long	NCELLS;
extern long	neighcellplus[54];
extern long	nneighcellplus;
extern void	New_CL();
extern void	CL_Update(long, long, long);

extern long	CL_Findcell(molstruct *, long, long ibox, long PBC);
extern long	CL_Neighbor(long, long, long);
extern void	CL_Init();
extern void	CL_Build();
extern void	CL_Destroy();
extern void	CL_Add(molstruct *, long);
extern void	CL_Delete(molstruct *, long);
extern void	CL_Relink(molstruct *moli);

#endif	/* CELL_LIST */

#ifdef VERLET_LIST
extern long	*vlistplus;
extern long	nverletplus;
extern long	maxnverlet;
extern void 	New_Vlist();
extern void	New_Vlist_LL();
extern void	Update_Vlist(long, vector, vector); 
extern void	Update_Vlist2(long, vector, vector); 
#endif	/* VERLET_LIST */

extern long	elem_index(long *, long, long);
extern long	Listlength(liststruct *);
extern int	List_Insert(liststruct **, long);
extern int	List_Remove(liststruct **, long);
extern void	Printlist(liststruct *);
extern void	Free_List(liststruct **);	//free the whole list
extern int	List_is_Empty(liststruct *);

extern void	New_ConnectList();
extern void	New_Clist();
extern void	Update_Clist(long, vector);

#endif	//ifdef __LIST_MODULE

#endif
