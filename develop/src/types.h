/*
    program:    types.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Definition of types used in the NPT project
*/
#ifndef __TYPES_HEADER
#define __TYPES_HEADER

#include	<complex.h>
#include	"input.h"

typedef struct{   
   double	x, y, z;
} vector;

typedef struct{
   double	alpha,		// bond angle
		beta,		// torsion angle
		d;		// bond length
} sphere;

typedef struct{
  vector x, y, z;
} matrix;
  
typedef struct{
   long		x, y, z;
} lvector;

typedef struct{
   int		x, y, z;
} ivector;

typedef struct{   	
  double 	lj, ljcorr,	// LJ potential, tail correction
		ljin, 		// LJ energy b/w sites on same chain
		ljout,		// LJ energy b/w sites on diff chains
		hs,		// hard sphere interaction
		stretch,	// bond stretching interaction
		bending,	// bending energy
		torsion,	// torsional energy
		corr,		// total correction
		bonded,		// bond related energy
		nonbonded,
		tot;
} vstruct;			// molecule potential energy structure

typedef struct{
  double	lj, 
		ljin,
		ljout,
		stretch,
		torsion,
		tot;
} wstruct;			//virial function structure

struct listNode{			//linked list structure
   long			neighbor;
   struct listNode * 	nextPtr;	//self-referential structure, typedef later
};
typedef	struct listNode	liststruct;
typedef liststruct *	liststrPtr;	//pointer pointing to a liststruct variable


/*
typedef struct{ 				// data structure for a single particle
  vector 	p,
		pv;				// pv is for Verlet list use
  long		box,				// which box this particle belongs to
  		nverlet,			// number of Verlet neighbors
 		vlist[MAXVERLETNEIGH],		// indexes of verlet neighbors
 		icell;				// which cell it belongs to, CELL_LIST

  // members above are essential for basic simulation
  // members below are for sampling, or more sophisticated simulation

  long 		nbond;		//number of nearest neighbor bonds
  complex	qlm[2*l_of_Ylm+1];	//aveq_lm for a certain l
  complex	ylmalphasum[2*l_of_Ylm+1];
  double	alphasum;	//qlm = ylmalphasum / alphasum
  double	ql;  		//ql value, but not used
  long		nconnect;	//how many particles are connected to this particle
  long		clist[MAXCONNEIGH];	
  long		nuclid;		//if type=1, which nucleus it belongs to, default -1
  long			nuclid2;	//for Find_Nuclei debug 6/9/2007
} molstruct;                
*/


typedef 
 struct cell_list {
   vector		center, p_min, p_max;	// center, lower corner, higher corner of this cell
   long			box,			// which box this cell belongs to
   			nsites,			// how many sites in this cell
			nempty,			// used as a stack
			empty[MAXNCELLSITES];
   struct mol		*mol[MAXNCELLSITES];	// pointer of molecules that have sites in this cell
   long			molsite[MAXNCELLSITES];	// id of sites (relative to mol) in this cell
   long			nneigh;			// # of neighboring cells
   struct cell_list	*neigh[27];		// neighboring cell id
} cellstruct;


typedef 
  struct mol {
   int		box, 			// which box this chain molecule belongs to
		nsites,			// how many sites on this chain molecule
		fix,			// whether this molecule is movable (1: fixed; 0: movable)
                flip,			// whether is a flip compared to the beginning
  		flags[MAXNMOLSITES],	// activated flag
		parent[MAXNMOLSITES],	// parent site
                type[MAXNMOLSITES];	// group type
   vector	p[MAXNMOLSITES];	// cartesian coord. of each site
   vector	origin;			// original position for drift calc.
   sphere	s[MAXNMOLSITES];	// spherical coord. of each site   
   cellstruct * cell[MAXNMOLSITES];	// cell id of each site in CELL_LIST implementation
   long		cellsite[MAXNMOLSITES];	// position of this site in its corresponding cell

   // members above are essential for basic simulation
   // members below are for sampling, or more sophisticated simulation
   float	velx[MAXNMOLSITES],
   		vely[MAXNMOLSITES],
   		velz[MAXNMOLSITES];
   float	p2[MAXNMOLSITES];	// local orientational parameter
   int		np2[MAXNMOLSITES];	// # of neighbors to calculate p2

   int		nuclid[MAXNMOLSITES];	// id of the crystal nucleus this site belongs to
   char		segtype[MAXNMOLSITES];	// segtype of this site (tail, loop, bridge, pbcbridge, xseg)

   // 150515, reduce size by not allocating array
   //int		nconn[MAXNMOLSITES],	// how many connection for each site
   //		connsite[MAXNMOLSITES][MAXCONNEIGH];	// site id of connected neighbors for each site
   //struct mol * connmol[MAXNMOLSITES][MAXCONNEIGH];	// mol id of connected neighbors for each site
   int *nconn;		
   int **connsite;
   struct mol ***connmol;
    
   // 150515, reduce size by not allocating array
   //float complex	qlm[MAXNMOLSITES][2*l_of_Ylm+1],
   //			ylmalphasum[MAXNMOLSITES][2*l_of_Ylm+1];
   //  float		alphasum[MAXNMOLSITES],
   //			q6[MAXNMOLSITES],
   //			q4[MAXNMOLSITES];
   float complex **qlm;		
   float complex **ylmalphasum;
   float *alphasum;
   float *q6;
   float *q4;

   vector	vp2[MAXNMOLSITES];	// p2 vector added Aug11.09

   int		nneigh[MAXNMOLSITES];		// number of nearest neighbors
   int		nneigh2[MAXNMOLSITES];

   int		neigh[MAXNMOLSITES][18];	// maximum 24 nearest neighbors

   short	cna[MAXNMOLSITES];		// value range 1-5
   short	cna2[MAXNMOLSITES];
   short	cna3[MAXNMOLSITES];
   int		clus[MAXNMOLSITES];  		// cluster id

   float	cnp[MAXNMOLSITES];		// common neighborhood parameter
   float	centro[MAXNMOLSITES];		// central symmetry parameter

   short	tmp[MAXNMOLSITES];

   float	pe[MAXNMOLSITES];		// potential energy per atom
   float	stress[MAXNMOLSITES][6];	// stress tensor per atom
   float        ql[MAXNMOLSITES][2]; 		// bond orientation order parameter ql

   float	rmin[MAXNMOLSITES];		// shortest nearest neighbor distance
   float	rmax[MAXNMOLSITES];		// longest nearest neighbor distance
   float	rave[MAXNMOLSITES];		// average nearest neighbor distance
} molstruct;

typedef struct {
   int		id;
   int		system;
   int		type;
   float	x, y, z;
   float	pe;
   float	stress[6];
} atomstruct;

typedef
  struct site {
   long			system;
   struct mol		*mol;
   long			mol_site;
   struct cell_list	*cell;
   long			cellsite;
   struct type		*type;
   long			flags, parent;
   vector		p;
   sphere		s;
} sitestruct;				// to describe a site, like molstruct to a molecule

typedef struct {
   molstruct 		*moli;
   int	     		site;
} beadstruct;				// (Jun6,2009) a more efficient version of sitestruct

typedef struct {
   double	SIGMA, EPSILON;
} mixstruct;


typedef
  struct type {					// force coefficent for different type of atoms
  						// e.g., CH3 and CH2 in PE are different in some
						// forcefield models
   double	M, Q,				// mass and charge
   		SIGMA, EPSILON,			// LJ interaction
		LSTRETCH, KSTRETCH,		// stretching
		KBENDING, THETA,		// bending
		TORSION[6],			// torsion	
		HS;				// hard-sphere
   mixstruct	mix[MAXNTYPES];
} typestruct;					


typedef struct {
   long		move, acc_move,			// normal site displacement
		erot, acc_erot,			// end mer rotation
		flip, acc_flip,			// one mer flip move
		re, acc_re,			// rebridging
		eb, acc_eb,			// endbridging
                db, acc_db,			// doublebridging
                idr, acc_idr, 			// internal double rebridging
		rot, acc_rot,			// rotation move
		vol, acc_vol,			// volume change
		swap, acc_swap,			// gibbs swap
		gibbsvol, acc_gibbsvol,		// gibbs volume change
		cbmc, acc_cbmc,			// cbmc
                rep, acc_rep,			// reptation
		movemol, acc_movemol,		// move whole molecule
		seq, acc_seq;			// sequence
} avstruct;


typedef
  struct structdist {
   char			*header;
   long			level, 			// 0 (1D dist), 1 (2D dist), etc.
			nbins, 			// # of bins of this level
			n, 			// total counts
			startbin, 		// first bin
			ncount;
   double		*binsize;		// could be multi-dimension
   struct structdist	*dist;			// sub dist, for multi-dimension
   long			*bin;			// counts in each bin
   double		*data, 
			*cweight, 
			*average;		// average of x, x^2, ..., x^(D_NAVERAGE-1)
} diststruct;


typedef struct{
   double	xlo,ylo,zlo,xhi,yhi,zhi;
   double	lbox, lx, ly, lz,	// dimension
		xy, xz, yz,		// for triclinic box, definition see lammps manual
		vol, pres, temp,	// volume, pressure, temperature
		rv, rb, rc;		// cutoff
   double	drmax, dlmax, damax;	// move size: displacement, volume, angle
   long		maxsize, xtal, nnucl;
   double	Ql;
} boxstruct;


typedef struct
{ 
    long                *systems, *cycle, *tape, *mols, *sites, *types, *volchange,
                        *swap, *displace, *inserts, *cavinserts, *box,
                        *rotation, *reptation, *endbridge, *rebridge, *cbmc,
                        *bridges, *fixed, *semifixed, *seed, *temper, *sample,
                        *stretch, *free, *ends, *system;
} nstruct;

typedef struct
{
   long			*hs, *lj, *ljshift, *ljlrc, *stretch, *bending, 
			*torsion, *virial, *scalecut, 
			*nvt, *npt, *gibbs, *mpi,
			*density, *energy, *pressure, *drift, *torsional,
			*bonda, *bondl, *radial, *localp2, *xtalsize;
} commandstruct;			// used in bridgestruct for binary history i/o

/* Pieter's
typedef struct
{
    long                *hs, *lj, *stretch, *ptorsion, *coulomb, *virial,
                        *polymer, *mpi, *async, *monodisperse, *bias, 
                        *nvt, *npt, *gibbs, *insert, *widom, *canonical, 
                        *cavity, *density, *tails_etc, *radial, *energy, 
                        *torsion, *re_torsion, *hs_dens, *temper, *jacob,
                        *d_bridge, *e_profile, *w_profile, *n_profile,
                        *b_length, *b_angle, *bonded, *loopreentry,
                        *e_n_function, *orient, *density3d, *densfree,
                        *denstalobr, *orientfree, *orienttalobr,
                        *orientcorr01, *orientcorr02, *orientcorr03, 
                        *orientcorr04, *orientcorr05, *orientcorrCR;
} commandstruct;
*/

typedef struct
{
    long                n, 						// system id, value
			*nmols, *nsites;				// pointers
    double              *pres, *vol, *temp, *drmax, *dlmax, *damax;
    //vector              a, b, c;					// system dimen. values
    double		*lx, *ly, *lz;					// system dimen. pointers
} systemstruct;				// system info., for binary history file, used in bridgestruct


typedef struct
{
    long                *d_type;
    nstruct             n;
    commandstruct       command;		// flags
    typestruct          *type;			// type of sites
    systemstruct        system[MAXNSYSTEMS];	// system info., some pointers, some values
    avstruct            *av[MAXNSYSTEMS];
    vstruct             *v[MAXNSYSTEMS];	// potential energy, all pointers
    wstruct             *vir[MAXNSYSTEMS];	// virial, all pointers
    diststruct          *dist[MAXNDIST];	// distribution, all pointers
    molstruct           *mol[MAXNSYSTEMS];	// molecule info., all pointers
} bridgestruct;			// for binary history file i/o

#endif
