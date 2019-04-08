/*
    program:    input.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    System setup
*/

/*
	The difference between input.h and globals.h is that input.h defines CONSTANTS 
	used in the program, while globals.h defined VARIABLES.
*/

#ifndef __INPUT_HEADER
#define __INPUT_HEADER

#define VERSION		"2014_01_20"

#define pi		M_PI		// 3.1415926535897932
#define ZERO 		1E-14		// floating point number of real zero
#define TRUE		1
#define FALSE		0	

#define epsl            1.0   		// Lennard Jones parameters
#define sigma           1.0		// basis for reduced unit

#define Rc2		Rc*Rc
#define Rb2		Rb*Rb
#define Rv2		Rv*Rv
#define Rp2		Rp*Rp
#define Rconn2		Rconn*Rconn
#define critqlproductSQ	(critqlproduct * critqlproduct)
#define Rc3		(Rc*Rc*Rc)

#define l_of_Ylm	6		// order of spherical harmonics

#define LINKLIST	0		// use linked list for verlet list
//#define VERLET_LIST	1
//#define CELL_LIST	1
#define	MAXNCELLSITES	3000		// max # of atoms in each cell 
#define MAXVERLETNEIGH	190		// max # of verlet neighbors of each particle
#define MAXCONNEIGH	14		// max # of connected neighbors of each particle

#define MAXNBOX		1		// max # of boxes
#define MAXNSYSTEMS	1		// max # of systems

#define MAXNMOLS	2000000     	// max # of molecules in ALL boxes
#define MAXNMOLSITES	1      		// max # of sites on one single chain molecule

#define MAXNTYPES	4		// max # of different forcefield coeff. sets

#define MAXNDIST	16		// max of distribution we can handle

#define UMBRELLA	0
//#define MPI		1
#define myMPI		1		// 5/4/2013, to speed up lammps2hst
					// not using MPI portion of code 
					// previously inherited from Pieter

#define MAX(a,b) (((a)>(b))?a:b)	// greater of two
#define MIN(a,b) (((a)<(b))?a:b)	// less of two

#define DEBUG	1
//#define TEST	1

#endif

