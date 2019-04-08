/*
    program:    forcefield.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Calculation of the potential of a system
*/
#ifndef __FORCEFIELD_HEADER
#define __FORCEFIELD_HEADER


#ifdef __FORCEFIELD_MODULE
#include "header.h"

long		V_LJ;			// LJ interaaction V = 4*epsilon*( (sigma/r)^12 - (sigma/r)^6 )
long		V_HS;			// hard-sphere
long		V_RPL;			// repulsive spheres V = epsilon * (sigma/r)^12
long		V_LJSHIFT;		// LJ forcefield shift at cutoff
long		V_LJLRC;		// LJ forcefield long range correction
long		V_STRETCH;		// bond stretching
long		V_BENDING;
long		V_TORSION;
long		FIXBONDL, FIXBONDA;	// flag of fixed bond length/angle
long		V_VIRIAL;		// calculate virial function
long		DLJ;			// intrachain min. distance to have LJ interaction (>=DLJ)
double		Rc, Rclow,		// LJ interaction cutoff, and smaller cutoff
		Rv,			// verlet list cutoff
		Rb, 			// neighbor cutoff
		Rp, 			// local p2 cutoff
		Rconn;			// connectivity cutoff
long		V_SCALECUTOFF;		// scale cutoff radii with box dimension
typestruct	type[MAXNTYPES];	// force coefficient
char		TORTYPE[80];		// torsion energy type

#else

extern long	V_LJ, V_HS, V_RPL;
extern long	V_LJSHIFT, V_LJLRC;
extern long	V_STRETCH, V_BENDING, V_TORSION;
extern long	FIXBONDL, FIXBONDA;
extern long	V_VIRIAL;
extern long	DLJ;

extern double	Rc, Rclow, Rb, Rv, Rp, Rconn;
extern long	V_SCALECUTOFF;
extern typestruct	type[MAXNTYPES];	// force coefficient
extern char	TORTYPE[80];

/* Basic operations */

extern void	vstructNull(vstruct *);
extern void	wstructNull(wstruct *);
extern void	vstructNegate(vstruct *);			// v =  -v
extern void	wstructNegate(wstruct *);			// virial =  -virial
extern vstruct	vstructSum(vstruct *, vstruct *);
extern wstruct	wstructSum(wstruct *, wstruct *);
extern void	Printvstruct(vstruct *);

/* LJ mixing */

extern void	CalcMixSigma(long typei, long typej);
extern void	CalcMixEpsilon(long typei, long typej);

/* Initialization */

extern void	InitForcefield();

/* Calculate interaction energy */

extern double	OPLS(double phi, double k1, double k2, double k3);
extern double	OPLS2(double cosphi, double k1, double k2, double k3);
extern double	VLJSite(molstruct *molm, long site, double *w);
extern double	VLJMol(molstruct *, double *w);
extern double	VHSSite(molstruct *, long);
extern double	VHSMol(molstruct *);
extern double	VStretchSite(molstruct *, long, double *w);
extern double	VStretchMol(molstruct *, double *w);
extern double	VBendingSite(molstruct *, long);
extern double	VBendingMol(molstruct *);
extern double	VTorsionSite(molstruct *, long);
extern double	VTorsionMol(molstruct *);

extern vstruct	CalcVSite(molstruct *, long, wstruct *w);	// calculate energy and virial for single site
extern vstruct  CalcVSiteinner(molstruct *, long, wstruct *w);
extern vstruct	CalcVSiteouter(molstruct *, long, wstruct *w);

extern void	CalcVLJ();
extern void 	CalcVLJ_MPI();
extern void	CalcVHS();
extern void	CalcVStretch();
extern void	CalcVBending();
extern void	CalcVTorsion();
extern void	CalcVCorr();
extern void	CalcVLJCorr();
extern void	CalcV();
extern void	CalcV_mcvol(double volscale);

extern void	EnergyCheck();

/* Energy differece by moving molecule */

extern double	VDeleteSites(molstruct *, long, long);
extern double	VAddSites(molstruct *, long, long);
//extern double	grow(char *, molstruct *, long, vstruct *, wstruct *);
extern double	grow(char *, molstruct *, long);
extern long	Select(double *, double);

/* Old stuff */

extern vstruct 	VMol(long); 
extern double	Vtailco(double, double);
extern void 	Vtotal();
extern double	VDeleteMol(long);
extern double	VAddMol(long);
extern vstruct	VTestMol(molstruct *, long);
extern double   VAddTestMol(molstruct *, long);

#endif

#endif
