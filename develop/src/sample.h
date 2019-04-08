/*
    program:    sample.h  
    author:     Peng Yi at MIT
    date:       October 22, 2006
    purpose:    header file for sample.c
*/
#ifndef __SAMPLE_HEADER
#define __SAMPLE_HEADER

#include "header.h"

#define D_NDIST		MAXNDIST	// MAXNDIST defined in input.h
#define D_BINSIZE	0.02		// not sure what is this, should make sense if reduced unit is used

#ifdef __SAMPLE_MODULE

#define	SAMPLE_CODE
long			D_DENSITY, D_ENERGY, D_PRESSURE, D_DRIFT,
			D_TORSION, D_BONDA, D_BONDL, D_RADIAL, 
			D_LOCALP2, D_XTALSIZE;
diststruct		*D_Density, *D_Energy, *D_Pressure, *D_Drift,
			*D_Torsion, *D_Bonda, *D_Bondl, *D_Radial, 
			*D_Localp2, *D_Xtalsize;
diststruct		*D_Distributions[D_NDIST];	// hook all sorts of distribution to this
							// -- adding convenience when operate to 
							// -- all distributions

long			dynvar;				//(1=NMAX; 2=Ql)
long			NMAXmiddle, NMAXbinsize, NMAXbins, Qlbins;	//dynamic variable range
double			Qlmiddle, Qlbinsize, kQ, kN;
double			P2middle, kP2;

double			drift2[MAXNSYSTEMS];
//structure property variables
double			P2[MAXNSYSTEMS];	// global orientation order
double			P2M[MAXNSYSTEMS];	// global orientation order modified
double			P2z[MAXNSYSTEMS];	// P2 with respect to z-axis
double			transfrac[MAXNSYSTEMS];	// trans state fraction
double			Ptors2[MAXNSYSTEMS][3][3];	// Ptt, Ptg+, Ptg-, Pg+g+, Pg+g-, Pg-g-
double			Ptors3[MAXNSYSTEMS][3][3][3];	// Pttt, etc

double 			Q6[MAXNSYSTEMS];	//bond orientation order parameter Q_l
double 			Q4[MAXNSYSTEMS];
/*
complex			Qlm[2*l_of_Ylm+1];	//Qlm = YlmAlphaSum/AlphaSum
complex			YlmAlphaSum[2*l_of_Ylm+1];
double			AlphaSum;
*/ // no need to be global variables

long			cnndist[25];		//connected neighbors number distribution

double				Qltemp;
/*
double				AlphaSumtemp;
complex				YlmAlphaSumtemp[2*l_of_Ylm+1];
complex				Qlmtemp[2*l_of_Ylm+1];
*/  // no need to be global variables

//crystal nuclei related variables
int			sizeofnucl[MAXNMOLS*MAXNMOLSITES];	// size of nucleus
int			sizedist[MAXNMOLS*MAXNMOLSITES];	// nucl size distri.
//int			sizeofnuclp2[MAXNMOLS*MAXNMOLSITES];
//int			sizedistp2[MAXNMOLS*MAXNMOLSITES];
//int * sizeofnucl;
//int * sizedist;
int * sizeofnuclp2;
int * sizedistp2;

long			MAXSIZE[MAXNSYSTEMS];	// maximum nuclei size in each system
long			Nnucl[MAXNSYSTEMS];	// # of Xtal nuclei
long			Xtal[MAXNSYSTEMS],	// total # of Xtal-like segments
			realXtal[MAXNSYSTEMS],
			realNnucl[MAXNSYSTEMS],
			secondNmax[MAXNSYSTEMS];
long			nmax[MAXNSYSTEMS][10];	// 10 biggest size, nmax[i][0]=MAXSIZE[i]
//beadstruct		nucleus[MAXNMOLS*MAXNMOLSITES/8];	// beads in the biggest nucleus
long				*sizeofnucl2;		//size of certain crystal nucleus
long				*sizedist2;		//how many nuclei of each size
long				MAXSIZE2;		//maximum nuclei size in the system
long				Nnucl2;			//# of Xtal nuclei
long				Xtal2;			//total # of Xtal-like particles

double			critp2,		// critical local p2 value to be considered Xtal phase
			critqlproduct,	// critical qlproduct to be considered connected
			critangle;	// critical angle between two vectors
long			critconnect;	// critical connection number to distinguish phase
long			NGRBINS;		// bins used to measure radial distribution function
double			Alpha;			//damping factor of eta update
double			CRIT;			//criterion of uniform sampling

#else

extern long		D_DENSITY, D_ENERGY, D_PRESSURE, D_DRIFT,
			D_TORSION, D_BONDA, D_BONDL, D_RADIAL,
			D_LOCALP2, D_XTALSIZE;
extern diststruct	*D_Density, *D_Energy, *D_Pressure, *D_Drift,
			*D_Torsion, *D_Bonda, *D_Bondl, *D_Radial,
			*D_Localp2, *D_Xtalsize;
extern diststruct	*D_Distributions[D_NDIST];

extern long		dynvar;
extern long		NMAXmiddle, NMAXbinsize, NMAXbins, Qlbins;	//dynamic variable range
extern double		Qlmiddle, Qlbinsize, kQ, kN;
extern double		P2middle, kP2;

extern double		drift2[MAXNSYSTEMS];
//structure property variables
extern double		P2[MAXNSYSTEMS], P2M[MAXNSYSTEMS], P2z[MAXNSYSTEMS];
extern double		transfrac[MAXNSYSTEMS];
extern double		Ptors2[MAXNSYSTEMS][3][3];
extern double		Ptors3[MAXNSYSTEMS][3][3][3];
extern double		Q6[MAXNSYSTEMS], Q4[MAXNSYSTEMS];
/*
extern complex		Qlm[2*l_of_Ylm+1];	//Qlm = YlmAlphaSum/AlphaSum
extern complex		YlmAlphaSum[2*l_of_Ylm+1];
extern double		AlphaSum;
*/
extern long		cnndist[25];		//connected neighbors number distribution

extern double			Qltemp;
/*
extern double			AlphaSumtemp;
extern complex			YlmAlphaSumtemp[2*l_of_Ylm+1];
extern complex			Qlmtemp[2*l_of_Ylm+1];
*/

//crystal nuclei related variables
extern int		sizeofnucl[MAXNMOLS*MAXNMOLSITES];	// size of nucleus
extern int		sizedist[MAXNMOLS*MAXNMOLSITES];	// nucl size distri.
//extern int		sizeofnuclp2[MAXNMOLS*MAXNMOLSITES];
//extern int		sizedistp2[MAXNMOLS*MAXNMOLSITES];

//extern int  *sizeofnucl;
//extern int  *sizedist;
extern int  *sizeofnuclp2;
extern int  *sizedistp2;

extern long		MAXSIZE[MAXNSYSTEMS];
extern long		Nnucl[MAXNSYSTEMS];		
extern long		Xtal[MAXNSYSTEMS],	
			realXtal[MAXNSYSTEMS],
			realNnucl[MAXNSYSTEMS],
			secondNmax[MAXNSYSTEMS];
extern long		nmax[MAXNSYSTEMS][10];		// size of 10 biggest nuclei
//extern beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES/8];
extern long			*sizeofnucl2;		//size of certain crystal nucleus
extern long			*sizedist2;		//how many nuclei of each size
extern long			MAXSIZE2;		//maximum nuclei size in the system
extern long			Nnucl2;			//# of Xtal nuclei
extern long			Xtal2;			//total # of Xtal-like particles

extern double		critp2, critqlproduct, critangle;
extern long		critconnect;	
extern long		NGRBINS;
extern double		Alpha, CRIT;


extern long 	mod(long, long);
extern long 	factorial(long);
extern long	intpow(long, long);
extern long	intlog(long, long);
extern double 	plgndr(int, int, double);
#ifdef SAMPLE_CODE
extern double complex 	sfharmonics(int, int, double, double);
extern double complex	sfharmonics2(int, int, double, double, double);

extern double complex 	aveqlm(int, int, long);
//extern complex 	tildeqlm(int, int, long);
extern double complex	qlproduct(int, long, long);
extern double	qlproductSQ(int, molstruct *, long, molstruct *, long);
extern double  	ql(int, long);
extern double complex 	aveQlm(int, int);
extern double	CalcQl(int);
extern double complex 	aveQlm1(int, int);
extern double	CalcQl1(int);
extern void	local_q_update();
extern void	New_Qlm(int);
extern void	New_Qlmtemp(int);
extern void	New_Qlm_NoVlist(int);
extern void	Update_Qlm(int, long, vector, vector); 

extern void	Calc_Qlm(long L);

extern int	Qlbinfinder(double);
extern int	NMAXbinfinder(int);
extern double	etaQl(double Ql);
extern double   etaNMAX(long MAXSIZE);
extern double	etaP2(double);
#endif //SAMPLE_CODE

extern vector	CenterofMass(molstruct *moli);			// center of mass of one chain
extern vector	groupCoM(beadstruct *group, long nsites);	// center of mass of a list
extern matrix	GyrationTensor(molstruct *moli);
extern matrix	groupGyraTensor(beadstruct *group, long nsites);
extern vector	fshape(vector *eigvalue);			// shape descriptors
extern matrix	InertiaTensor(molstruct *moli);
extern matrix	groupInerTensor(beadstruct *group, long nsites);
extern matrix	grpgytensor(int size, vector *rbead);
extern vector		CenterofNucleus(long nuclid, molstruct *molm);
extern void		InitSample();
extern void		ReinitSample(diststruct **d);
extern void		SampleEnergy();			// it will be called by outsider, so claimed here
extern void		SampleDistributions();
extern void		PrintDistributions();
extern void		S_PrintAll();
extern void		Sample2D();
extern void		SampleN2NSQ();
extern void		SampleDrift();
extern void		SampleP2();
extern void		SampleP2All();
extern void		SampleP2All_MPI();
extern void		Dist_p2();
extern void	SampleM_Q();
extern void	SampleConnection();
extern void		SampleSpherical();
extern void		Dist_Spherical();

extern void	Init_Sample(char *argv[]);
extern void	Sample_Energy();
extern void	Sample_Done();
extern void	Sample_Pressure();
extern void	Sample_All();
extern void	Sample_Histogram();
extern void	Sample_G(int);		//sample Gibbs free energy

extern float	R2_gyration(molstruct *moli);		// radius of gyration square
extern float	rx2_gyration(molstruct *moli);
extern float	ry2_gyration(molstruct *moli);
extern float	rz2_gyration(molstruct *moli);
extern float	group_Rg2(vector *rbead, long nsites);	// R2 gyration of one group
extern vector	group_com(vector *rbead, long nsites);	// geometrical center of one group
extern float	R2_n2n(molstruct *);			// end-to-end distance square
extern vector	Center_of_Mass(long *, long);		// center of mass of listed atoms
extern matrix   Gyration_Tensor(long *, long);		// gyration tensor of listed atoms
extern void	radial(char *);				// radial distribution function
extern void	sq(FILE *, char *);			// structure factor calculation
extern void	sq1(char *);				// sq calc. by definition
extern void	sq2(char *);				// sq calc. by definition
extern void	sq3(char *);				// sq calc. by definition
extern void	sq4(char *);				// sq calc. by definition
extern void	atom_neighbor();			// build nearest neighbor list
extern void	comm_neigh_para();			// calculate common neighborhood parameter
extern float	hyper_distance(int,float *, float *);	// calculate hyper-distance
extern float	slip_distance(int,float *, float *);	// calculate slip-distance
extern float 	noncubic();				// noncubic factor of orthorgonal simulation box
extern void     composition_profile(long timestep, char, vector, double *, double *, double);
extern void     profile_2D(long timestep, char *, char *, double);
extern void     ave_property(char *property, float cutoff);    	// compute local average with a cutoff
#endif

#endif
