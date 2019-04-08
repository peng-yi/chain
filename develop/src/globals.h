/*
    program:    globals.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Centralization of all global variables used in multiple
                programs
*/
#ifndef __MY_GLOBALS
#define __MY_GLOBALS

//#include "input.h"
#include "types.h"

#ifdef __MAIN_PROGRAM   
                        /* globals.h is included in header.h which in turn 
			   will be included in main.c, so we want the following
			   variables to be global variables.  But if any 
			   subroutine include header.h (thus global.h), since 
			   in them __MAIN_PROGRAM is not defined, only main.c
			   will define __MAIN_PROGRAM, they will only see like
			   extern long NPARTS , etc.
			   in this way we gurantee that these global variables
			   will not be defined multiple times.
			*/

//env variables
char			moltype[80];			// molecule type: LJ, monochain, polychain
int			NSITES;				// int to be compatable with MPI_INT
long			NPARTS, NBOX, NSYSTEMS, NTYPES,
			NMOLS,				// total mols, total sites
			Nsites,				// sites in one mol
   			Nequil, Nprod, NCYCLE, TRIALRUN, PROD,
			ITAPE, ICONF, NBLOCKS, NGSAMPLE, IRADIAL;
long			SEQUENCE;			// how many sequences in each MC cycle
long			TIMESTEP;			// used in i/o

long			NMols[MAXNSYSTEMS], NSites[MAXNSYSTEMS];
long			Ntest;
double			Dpoly, *mupoly;			// polydispersity variables

double			kT, P,
	                Rho, LBOX0,
			LBOX, VOL,
	                DRMAX, DLMAX, GDLMAX, DAMAX;

molstruct		*mol, *oldmol;
molstruct        	*part;			// define NPARTS particles in the system
molstruct		*parttemp;

atomstruct		*atom;

wstruct			VIRMol;			// Virial for single molecule
wstruct			VIRSite;
vstruct			v[MAXNBOX];		// potential energy of each box
wstruct			vir[MAXNBOX];		// virial of each box
double			chp[MAXNBOX];		// chemical potential
double			cchp[MAXNBOX];		// cumulative chemical potential
long			ichp[MAXNBOX];		// chemical potential sampling times
boxstruct		BOX[MAXNBOX];		// enviro. parameter of each box
long			N[MAXNBOX];		// particle # in each box

//Monte Carlo step variables
long			Stage;			// which Stage it is in (melting, quenching, etc)
long			counter;   		// MC cycle counter
long			FLAG_restart;		// flag to restart an unfinished run
long			TIME_restart;		// restart timepoint

//block average
//double			ave_G[Qlbins][NGSAMPLE*2];
//double			ave_sqG[Qlbins];

//expanded ensemble variables

/*
double			*p;			//probability distribution 
double			*eta;			//weighting factor
double			*pt;			//target probability distribution
double			*G;			//Gibbs free energy

double			*pQ;
double			*etaQ;
double			*ptQ;
double			*GQ;
*/

int			trial;			//trial run counter

//random number generator variables
int                     *tim;  	
long int                *seed; 

//time and date variables
time_t			curtime;
struct tm 		*loctime;		
char			buffer[32];

//MPI variables
#ifdef myMPI
MPI_Comm		comm_world, comm_lammps;
MPI_Group		group_world, group_lammps;
MPI_Request		MPIrequest;
MPI_Status 		MPIstatus;
#endif
int			MPI_numprocs;			// number of processors
int			MPI_myid;			// process id
int			MPI_lammpsid;
int			nprocs_lammps;

//others
int			S_STOP;

#else

extern char		moltype[80];
extern int		NSITES;
extern long		NPARTS, NBOX, NSYSTEMS, NTYPES, NMOLS, Nsites;
extern long		NMols[MAXNSYSTEMS], NSites[MAXNSYSTEMS];
extern double		Dpoly, *mupoly;	
extern double		Rho, LBOX0;
extern double		kT, P;
extern long		Nequil, Nprod, NCYCLE, TRIALRUN, PROD,
			ITAPE, ICONF, NBLOCKS, NGSAMPLE, IRADIAL;
extern long		SEQUENCE;
extern long		TIMESTEP;
extern long		Ntest;

extern molstruct	*mol, *oldmol;
extern molstruct        *part;
extern molstruct		*parttemp;

extern atomstruct	*atom;

extern double           LBOX, VOL;
extern double           DRMAX, DLMAX, GDLMAX, DAMAX;

extern wstruct		VIRMol, VIRSite;
extern vstruct		v[MAXNBOX];
extern wstruct		vir[MAXNBOX];
extern double		chp[MAXNBOX];
extern double		cchp[MAXNBOX];
extern long		ichp[MAXNBOX];
extern boxstruct	BOX[MAXNBOX];
extern long		N[MAXNBOX];

extern long		Stage;
extern long 		counter;	
extern long		FLAG_restart;
extern long		TIME_restart;

//extern double		ave_G[Qlbins][NGSAMPLE*2];
//extern double		ave_sqG[Qlbins];

/*
extern double		*p;
extern double		*eta;
extern double		*pt;
extern double		*G;

extern double		*pQ;
extern double		*etaQ;
extern double		*ptQ;
extern double		*GQ;
*/

extern int		trial;

extern int              *tim;  
extern long int         *seed;        

extern time_t		curtime;		//date and time
extern struct tm	*loctime;		
extern char		buffer[32];

#ifdef myMPI
extern MPI_Comm		comm_world, comm_lammps;
extern MPI_Group	group_world, group_lammps;
extern MPI_Status 	MPIstatus;
extern MPI_Request	MPIrequest;
#endif
extern int		MPI_numprocs;		// number of processors
extern int		MPI_myid;		// proc id
extern int		MPI_lammpsid;		// proc id in lammps communicator
extern int		nprocs_lammps;

extern int		S_STOP;
#endif

#endif
