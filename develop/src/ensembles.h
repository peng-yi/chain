/*
    program:    ensembles.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Header file for ensembles.c.
*/
#ifndef __ENSEMBLES_HEADER
#define __ENSEMBLES_HEADER

#include "header.h"
#define	E_NMINFREE	3	

#ifdef __ENSEMBLES_MODULE


long		E_NVT,			// use NVT ensemble
		E_NPT,			// use NPT ensemble
		E_GIBBS,		// use GIBBS ensemble
		E_MPI;			// use MPI
long		NDISPLACE,		// % of trial displacement in each MC cycle among NVT moves
		NREPTATION,		// % of trial reptation in MC cycle among NVT moves
		NENDROT,		// end-mer rotation
		NCBMC,			// % of trial CBMC moves in each MC cycle among NVT moves
		NTRIALCONF,		// # of trial conf. in cmbc
		NTRIALFIRSTBEAD,	// # of trial positions of first bead in cmbc
		NVOLCHANGE,		// # of trial volume change in each MC cycle
		NGIBBSVOL,		// % of trial Gibbs volume change in each MC cycle
		NSWAP,			// # of trial Gibbs swap moves in each MC cycle
		NENDBR, NREBR,		// # of trial end-bridging and rebridging
		NDB, NIDR,		// # of trial double bridging and intra-double-rebridging
		NFLIP;			// # of trial one mer flip move
long            acc_move, rjc_move;  	// # of accepted/rejected trial move
long		acc_vol, rjc_vol;	// # of accepted/rejected trial volume change
long		acc_seq, rjc_seq;	// # of accepted/rejected trial move sequence
long		acc_gibbsvol, rjc_gibbsvol;
long		acc_swap, rjc_swap;

long		cbmcsucc[MAXNMOLSITES];
avstruct	av[MAXNBOX], av_past[MAXNBOX];
double		SUCC_DISP, SUCC_VOL;	// target accepted ratio of different kinds of move

#else

extern long	E_NVT, E_NPT, E_GIBBS, E_MPI;
extern long	NDISPLACE, NVOLCHANGE, NGIBBSVOL, NSWAP, NCBMC, NREPTATION, NENDROT, 
		NENDBR, NREBR, NDB, NIDR, NFLIP, NTRIALCONF, NTRIALFIRSTBEAD;

extern long     acc_move, rjc_move;  
extern long	acc_vol, rjc_vol;
extern long	acc_seq, rjc_seq;
extern long	acc_gibbsvol, rjc_gibbsvol;
extern long	acc_swap, rjc_swap;
extern long	cbmcsucc[MAXNMOLSITES];
extern avstruct	av[MAXNBOX], av_past[MAXNBOX];
extern double	SUCC_DISP, SUCC_VOL;

/* Move preparations */

extern void	StoreSystem();
extern void	RestoreSystem();
extern void	StoreOneMol(long);
extern void	RestoreOneMol(long);
extern void	StoreMols();
extern void	RestoreMols();

extern void	ResetAcceptance();	// reset move step sizes
extern void	Adjust_Stepsize();	// adjust move step sizes

extern void	MolInBox(molstruct *molm);
extern void	MolInBox2(molstruct *molm);
extern long	SiteSelect(molstruct **);
extern double	bondl_g(double, double);	// generate a bond length according to stretching energy
extern void	bonda_tors(double, double, double *, double *);
extern double	bonda_g(double, double);
extern double	tors(double);
extern vector	tors_bonda(molstruct *molm, long site);		// generate a unit vector according to bending and torsional energy

/* Moves */

extern long	rotation(long n);	// rotation of n sites, bond length and angle unchanged
extern long	movemol();		// move a molecule with bond length fixed
extern long	mcmove();   		// single particle displacement
extern long	mccbmc();		// conf.-biased Monte Carlo move
extern long	reptation();
extern long	mcvol();    		//system volume change
extern void	Gibbsmove();
extern void	NextMove();
extern void	Cycle();		//a Monte Carlo cycle
extern void	InitEnsemble();

extern void	ParentCheck();

extern void	Update_Eta(int);	//update weighting factors
extern int	CheckConstraint(char *);

#endif

#endif
