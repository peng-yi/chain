/*
  program:      position.c
  author:       Peng Yi at MIT
  date:         October 20, 2006
  purpose:      Position calculation, Xtal nuclei identification, and speed-up
*/
#ifndef __POSITION_HEADER
#define __POSITION_HEADER

#include "types.h"

#ifdef __POSITION_MODULE

#include "header.h"

long		PBC;		//periodic boundary condition variable
				//0: no pbc; 1: cubic pbc; 2: truncated octahedron pbc
long		critconnect;	// # of critical connection in LJ system

int		nsegment[MAXNMOLS];			// # of segments on each chain
int		seg_stat[MAXNMOLS][3*MAXNMOLSITES];	// segment stat. for chains
vector	r;

#else

// variables

extern long	PBC;
extern long	critconnect;
extern int	nsegment[MAXNMOLS];
extern int	seg_stat[MAXNMOLS][3*MAXNMOLSITES];
extern vector r;

// functions

extern double	AdjustAngle(double x);
extern void	MapInBox(vector *p);
extern vector 	MapInBox2(vector *p, long PBC, long system);
extern void	unfold();
extern void	unfoldchain(molstruct *, long);
extern void 	MapInNucleus(int system, beadstruct *nucleus, int nsites, vector *rbead, float rcutoff);

extern double	DistSQ(vector p, vector q, long system);

extern void	InitLattice(long, long, double, long);
extern vector	ranor();
extern void	Amorph(long nmols, long nmolsites, double LX, double LY, double LZ, long PBC);

extern void	sc_lattice(long, double, long);
extern void	bcc_lattice(long, double, long);
extern void	fcc_lattice(long, double, long);
extern void	randomconf(long, double, long);
extern int	crystal(molstruct *moli, long site);
extern long     getnuclsize(long, long);
extern void 	Find_Nuclei(long clusdef);
extern void	Find_Nuclei_LJ();
extern void	Find_Nuclei2();
extern void	Find_Nuclei1();
extern void	Find_Nuclei_p2(long clusdef);
extern void	find_nuclei_general(char *method, char *modifier, char *neighbormethod, float rcutoff);
extern vector	CoM_MaxNucleus(long system);
extern vector   cylindershape(beadstruct *nucleus, long size, long nuclid);
extern vector   cylinder(beadstruct *nucleus, long size, long nuclid);
extern void	Find_segments();
extern vector	Xtal_tilt(int);
extern vector	Lamella_norm(int);
extern void	Rho_profile();
extern void	Seg_type(int);
extern void	Seg_length(int, int array[][MAXNMOLSITES]);
extern vector	Seg_smooth(long);
extern float	xtal_smooth();

#endif

#endif
