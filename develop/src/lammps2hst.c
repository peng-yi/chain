/*
	program:	lammps2hst.c
	author:		Peng Yi at MIT
	date:		Feb. 14, 2008
	purpose:	read lammps dump file, do some sampling 
			and output histogram file
	note:		require setup file
	update:
			Jul. 28, 2009	Add average center of mass calculation
 			Mar.  7, 2011	Add segments statistics
			Mar.  9, 2011	Allow skipping frames
			Apr. 28, 2011	Allow multiple input data files 
			June 15, 2011	Allow triclinic box, but energy and pressure 
					calculation not fixed for triclinic box yet.
			Nov. 18, 2011	Move car, conf and pdf output out of main program
			Nov. 19, 2011	Create output for pre-analysis results, e.g., 
					nuclid info.  Output file .pre
			Mar. 16, 2012	Add timestep range
			Aug. 27, 2012	Add unfold()
			May, 2013	MPI implementation
*/

#define __MAIN_PROGRAM
#include "header.h"

#define L2HVERSION	"1/20/2014"
#define CAP		5		// largest nuclei considered still the melt
#define L2HDEBUG	0
#define SIZECAP		5
#define DEN_BIN		12
#define MAXNINPUTFILE	32

#include "correlation.h"		// autocorrelation calculation module

//==============================//
//	Global variables	//
//==============================//
long		timestep;
int 		nmaxp2_1[10];		// nmax using p2 nucleus definition
int 		nmaxp2_2[10];
int 		nmaxp2_3[10];
int 		nmaxp2_4[10];
float		rshift2;		// shift of the biggest nucleus

float		R02;
float		Rg2[MAXNMOLS];		// Radius of gyration square
float		aveRg2;			// average Rg2 for the system
float		avexRg2, avemRg2;	// average Rg2 for xtal and melt phases
int		nRg2, nxRg2;
vector		smoothresult;
float		rgx2[MAXNMOLS], rgy2[MAXNMOLS], rgz2[MAXNMOLS];
float		avergx2, avergy2, avergz2;

vector		com_sys,		// com of system
		com_sysinit,		// com of initial system
		com_syslast,		// com of previous system
		com_init[MAXNMOLS],	// com of each chain in the beginning
		com_last[MAXNMOLS];	// com of each chain in the previous dump
vector		dr_sys;			// vector of system drift
float		d2initave,		// average displacement w.r.t. to the beginning
		d2lastave,		// average displacement w.r.t. to previous dump
                d2sys, d2syslast;
float		asphericity;		// asphericity of the box

double		lam1, lam2, amo1, amo2;	// thickness of lamellae and interface
double		v2, vf2, vcom2;		// velocity and kinetic energy analysis

//==============================================//
//	Function Declarations and Definitions	//
//==============================================//
void		output_Rg();				// Output Rg analysis
void 		output_pre(FILE *);			// Output pre-analysis results
void 		read_pre(FILE *, long timestep);	// Input pre-analysis results
void 		output_conf(FILE *, long timestep);	// Output configuration file
void 		output_carfile(FILE *);			// Output conf. to .car file
void 		output_pdbfile(FILE *, FILE *, long, double, double, double);	// to .pdb file
void 		output_xyzfile(FILE *);
void		output_Zinput(FILE *);			// Output to Z code input file
void		Seg_Type(int);

char *Element(long t, char *s)
{
   double	mass = type[t].M;

   s[0]		=	0;
/*
   if (14==mass || 15==mass)	strcpy(s, "C");
   else if (1.01==mass)		strcpy(s, "H");
   else if (28.086==mass)	strcpy(s, "Si");
   else if (26.982==mass)	strcpy(s, "Al");
   else if (16==mass)		strcpy(s, "O");
*/
   strcpy(s, "C");
   return	s;
}

void Print_hst(FILE *fPtr)
{
   long		i;
   static long	init_hst = 1;

   if (init_hst) {
      init_hst	=	0;
      // Print variable names
      fprintf(fPtr, "Timestep ");
      fprintf(fPtr, "Epot ");
      fprintf(fPtr, "Volume ");
      fprintf(fPtr, "Pressure ");

      fprintf(fPtr, "Lx Ly Lz ");	// box dimensions
      fprintf(fPtr, "xy xz yz ");	// triclinic box parameters

      fprintf(fPtr, "P2 P2m P2z ");	// orientation order parameters
      fprintf(fPtr, "Transfrac ");
      fprintf(fPtr, "Xtal ");
      fprintf(fPtr, "RealXtal ");
      fprintf(fPtr, "Nnucl ");

      fprintf(fPtr, "Nmaxp2 2ndNmaxp2 ");
/*
      fprintf(fPtr, "Nmaxp2_1 2ndNmaxp2_1 ");
      fprintf(fPtr, "Nmaxp2_2 2ndNmaxp2_2 ");
      fprintf(fPtr, "Nmaxp2_3 2ndNmaxp2_3 ");
      fprintf(fPtr, "Nmaxp2_4 2ndNmaxp2_4 ");
*/
      fprintf(fPtr, "Q6 ");
      fprintf(fPtr, "rshift2 ");

      fprintf(fPtr, "E_bond ");
      fprintf(fPtr, "E_angle ");
      fprintf(fPtr, "E_tors ");
      fprintf(fPtr, "E_lj ");
      fprintf(fPtr, "E_ljcorr");

      fprintf(fPtr, "\n");
   }

   // Write variable values

   fprintf(fPtr, "%-6ld ", timestep);
   fprintf(fPtr, "%8.4f ", v[0].tot);
   fprintf(fPtr, "%8.4f ", BOX[0].vol);
   fprintf(fPtr, "%8.4f ", BOX[0].pres);
   fprintf(fPtr, "%5.3f %5.3f %5.3f ", BOX[0].lx, BOX[0].ly, BOX[0].lz);
   fprintf(fPtr, "%5.3f %5.3f %5.3f ", BOX[0].xy, BOX[0].xz, BOX[0].yz);

   fprintf(fPtr, "%6.4f %6.4f %6.4f ", P2[0], P2M[0], P2z[0]);
   fprintf(fPtr, "%6.4f ", transfrac[0]);
   fprintf(fPtr, "%4ld ", Xtal[0]);
   fprintf(fPtr, "%4ld ", realXtal[0]);
   fprintf(fPtr, "%4ld ", Nnucl[0]);

   for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_3[i]);
/*
   for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_1[i]);
   for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_2[i]);
   for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_3[i]);
   for (i=0; i<2; i++)      fprintf(fPtr, "%4d ", nmaxp2_4[i]);
*/

   fprintf(fPtr, "%6.4f ", Q6[0]);
   fprintf(fPtr, "%6.4f ", rshift2);

   fprintf(fPtr, " %8.4f", v[0].stretch);
   fprintf(fPtr, " %8.4f", v[0].bending);
   fprintf(fPtr, " %8.4f", v[0].torsion);
   fprintf(fPtr, " %8.4f", v[0].lj);
   fprintf(fPtr, " %8.4f", v[0].ljcorr);

   fprintf(fPtr, "\n");
   fflush(fPtr);
   return;
}

void shiftbox(long system, beadstruct *nucleus, long nsites) //copied from conf2car.c (07/13/09)
{
   molstruct	*moli;
   long		i, n, site;
   vector	rA, rO, rBA, rOA;

   // step 1: find one segment that belongs to the biggest nucleus
   
   rA	=	nucleus[0].moli->p[nucleus[0].site];
  
   // step 2: calc. the shift of com of nucleus to this segment
   
   V_Null(&rBA);
   V_Null(&rOA);

   for (n=0; n<nsites; n++) { 
      moli	=	nucleus[n].moli;
      site	=	nucleus[n].site;
      rBA	=	V_Subtr(moli->p+site, &rA);
      rBA	=	MapInBox2(&rBA, PBC, system);
      rOA	=	V_Add(&rOA, &rBA);
   }
   rOA	=	V_Mult(1.0/n, &rOA);
   rO	=	V_Add(&rA, &rOA);

   // step 3: every segment shift and move to central box
  
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (site=0; site<moli->nsites; site++) {
         moli->p[site]	=	V_Subtr(moli->p+site, &rO);
         moli->p[site]	=	MapInBox2(moli->p+site, PBC, system);
      }
   }
   return;
} 

//======================//
//	Main Program	//
//======================//
int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL), prog_start, prog_end, start, end;
   molstruct	*moli;
   long		i, j, k, system, m, n, id, maxid, max, size;
   long		siteid, molid, type;
   int		nx, ny, nz, nx1, ny1, nz1;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo,			// box upper boundary
		xy=0, xz=0, yz=0;		// triclinic parameter
   long		LENGTH, accum;			// chain length variables

   // dummy variables
   
   long		ltemp1, ltemp2, ltemp3;
   double	temp, temp1, temp2, temp3;
   vector	rtemp, rtemp1, rtemp2;

   // file variables
   
   char		infile[MAXNINPUTFILE][80], filein[80], filename[80], procid[16];
   char		s[80], ff[80], par[80], dummy[255];
   char		atomname;			// for visualization file output

   FILE		*fin, *fhst, *fpre;		// input and analysis output
   FILE		*fxyz;
   FILE		*fcolor;
   FILE		*fcar, *fconf, *fpdb, *fdat;	// configuration/visualization output
   FILE		*fZfile;			// Z-code input configuration file

   int		preflag		= 0;		// flags for input command line
   int		polydisperse	= 0;
   int		confflag	= 0;
   int		carflag		= 0;
   int		pdbflag		= 0;
   int		xyzflag		= 0;
   int		Zflag		= 0;
   int		nfiles=1, ifile;		// number of input files 
   int		nframe, dnframe=1;		// analyze only every dnframe
   int		previous_timestep=-1;
   int		starttime=-1;
   int		endtime=-1;

   static int	clistinit=1, cominit=1;	

   vector		con;			// center of nucleus
   static vector	rO;			// original position of nucleus

   // a group of beads
   
   vector	rbead[MAXNMOLS*MAXNMOLSITES];

   // variables to characterize nucleus, chains, etc at ONE timestep

   long		nsites;
   int 		nuclid;
   int		nmaxid;					// nuclid of the largest nucleus
   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];		// group beads in the same nucleus
   vector	com[MAXNMOLS];				// average center of mass
   long		ncom[MAXNMOLS];
   int		nchainmax;				// # of chains in the largest nucleus

   // chain rotational angle distribution

   vector	chainface;
   double	orient;
   long		orientdist[180];

   // density profile

   long		direct;
   double	density[100];
   vector	rsite[MAXNMOLS * MAXNMOLSITES];	

   // velocity statistics
 
   double	dvel=0.001;
   double	ave_vx, ave_vy, ave_vz;
   double 	velmax=0.0, velmin=0.0;
   double	vel_distx[101], vel_disty[101], vel_distz[101];
 
   // variables to AVERAGE over a number of timesteps

   long		*Nn;				// size dist. of nucleus size
   long		ntot, maxnmax;			// total number of nuclei for normalization
   long		xbeadpos[MAXNMOLSITES];		// xtal bead position on chain
   float	avePtt, avePttt;		// prob. of 2 or 3 consecutive trans states
   float	stdPtt, stdPttt;
   float	*Rg2nucl;			// Rg2 of nucleus as a function of size

   // segment statistics variables

   int 		head, tail, seg_on, seg_id, nseg, nxtal, length;
   int 		segment[MAXNMOLS][MAXNMOLSITES];	// segments identification on a chain
							// -- segment[molid][siteid]
   int		ntail, nloop, nbridge, npbrdg, nxseg;	// # of tails, etc
   float	ltail, lloop, lbridge, lpbrdg, lxseg;	// average leng of tails, etc
   float	lfree;	
   int		lsegdist[10][MAXNMOLSITES];		// segment length distribution
   int		nltail[MAXNMOLSITES],			// tail length distribution
   		nlloop[MAXNMOLSITES], 			// loop
		nlbridge[MAXNMOLSITES],			// bridge
                nlpbrdg[MAXNMOLSITES],			// p.b.c bridge
		nlxseg[MAXNMOLSITES],			// xtal stem
		nlfree[MAXNMOLSITES];			// free chain
   int		oldxtal=0, oldnnucl=0,			// xtal variables before segment smoothing
                oldnmax=0, old2nmax=0;
   vector	rhead, rtail;

   // pbc image variables
   
   long		imagen=0;
   double	imagex, imagey, imagez;	

   //------MPI initialization------//
   MPI_numprocs	=	1;			// default value
   MPI_myid	=	0;			// default value

#ifdef myMPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &MPI_numprocs);
   MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myid);
#endif

   if (argc<2 && MPI_myid==0) {
      printf("lammps2hst (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tlammps2hst [-option] [x= 0.1 y= 0.1 z= 0.1 n= 1] [dn= 2 nfiles= 2] lammpsdumpfile(s)\n\n");
      printf("Notes:\n");
      printf("\t* -option = \n");
      printf("\t* -pre: generate pre-analysis results\n");
      printf("\t* -poly: polydisperse system, has molid in dump file\n");
      printf("\t* -conf: configuration file output\n");
      printf("\t* -car: car file output\n");
      printf("\t* -pdb: pdb file output\n");
      printf("\t* -xyz: xyz file output\n");
      printf("\t* x= y= z=: duplicate the system and shift unit vector\n");
      printf("\t* n=: multiple of shift vector\n");
      printf("\t* start=: the timestep that analysis stops\n");
      printf("\t* end=: the timestep that analysis stops\n");
      printf("\t* dn=: only analyze every dn frames\n");
      printf("\t* nfiles=: number of input files if more than 1 (must be <=MAXNINPUTFILE)\n");
      printf("\t* \"=\" must immediately follow x or y or z or n or dn\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-poly"))	polydisperse	=	1;
      else if(samestr(par, "-pre"))	preflag		=	1;
      else if (samestr(par, "-conf"))	confflag	=	1;
      else if (samestr(par, "-car"))	carflag		=	1;
      else if (samestr(par, "-pdb"))	pdbflag		=	1;
      else if (samestr(par, "-xyz"))	xyzflag		=	1;
      else if (samestr(par, "-Z"))	Zflag		=	1;
      else if (samestr(par, "x="))	imagex		=	atof(argv[i+1]);
      else if (samestr(par, "y=")) 	imagey		=	atof(argv[i+1]);
      else if (samestr(par, "z=")) 	imagez		=	atof(argv[i+1]);
      else if (samestr(par, "n=")) 	imagen		=	atol(argv[i+1]);
      else if (samestr(par, "start=")) 	starttime	=	atol(argv[i+1]);
      else if (samestr(par, "end=")) 	endtime		=	atol(argv[i+1]);
      else if (samestr(par, "dn=")) 	dnframe		=	atol(argv[i+1]);
      else if (samestr(par, "nfiles="))	nfiles		=	atoi(argv[i+1]);
   }
   if (nfiles > MAXNINPUTFILE) {
      printf("Error! number of input files exceeds MAXNINPUTFILE\n");
   }
   for (i=0; i<nfiles; i++) {
      strcpy(infile[nfiles-1-i], argv[argc-1-i]);	// get input filenames
   }

   //---------------------------//
   //	Open output files	//
   //---------------------------//
 
   if (nfiles==1)	strcpy(filein, infile[0]);	// determine output filename
   else			strcpy(filein, "multi");

   sprintf(procid, ".%d", MPI_myid);			// MPI: proc id for output files

#ifdef	myMPI
   //if (MPI_myid==0) {					// ONLY process 0 do the following
#endif

   if (preflag) {
      strcpy(filename, filein);
      strcat(filename, ".hst");
      strcat(filename, procid);
      if ((fhst=fopen(filename, "w"))==NULL )		// .hst output file
         Exit("lammps2hst", "main", "open .hst file failed.");

      strcpy(filename, filein);
      strcat(filename, ".pre");
      strcat(filename, procid);
      if ((fpre=fopen(filename, "w"))==NULL )		// .pre output file
         Exit("lammps2hst", "main", "open .pre file failed.");
   }

   //------Open optional output files------//

   if (Zflag) {
      fZfile	=	fopen("Z_data", "w");		// Z_code file
   }

   if (xyzflag) {
      strcpy(filename, filein);
      strcat(filename, ".xyz");
      strcat(filename, procid);
      if ((fxyz=fopen(filename, "w"))==NULL )
	 Exit("lammps2hst", "main", "open .xyz file failed.");

      fcolor	=	fopen("color", "w");
   }

   if (confflag) {
      strcpy(filename, filein);
      strcat(filename, ".conf");
      strcat(filename, procid);
      if ((fconf=fopen(filename, "w"))==NULL )
	 Exit("lammps2hst", "main", "open .conf file failed.");
   }

   if (carflag) {
      strcpy(filename, filein);
      strcat(filename, ".car");
      strcat(filename, procid);
      if ( (fcar=fopen(filename, "w"))==NULL )
	 Exit("lammps2hst", "main", "open .car file failed.");
   }

   if (pdbflag) {
      strcpy(filename, filein);
      strcat(filename, ".pdb");
      strcat(filename, procid);
      if ( ((fpdb=fopen(filename, "w"))==NULL || (fdat=fopen("vmd.dat","w"))==NULL))
	 Exit("lammps2hst", "main", "open .pdb/dat file failed.");
   }

#ifdef	myMPI
   //}
#endif

   if (!L2HDEBUG) {			
      strcpy(filename, filein);
      strcat(filename, ".out");			// every processor does output
      strcat(filename, procid);
      freopen(filename, "w", stdout);		// redirect standard output stream to a file
   }

   //------Write run info. to output file------//

   printf("# BEGINNING OF OUTPUT FILE\n");
   printf("%s\n", asctime(localtime(&t)));

   printf("# lammps2hst version: %s\n", L2HVERSION);
   printf("# command: ");
   for (i=0; i<argc; i++) {	
      printf("%s ", argv[i]);	
   }
   printf("\n");
   printf("# number of processor:  %d\n", MPI_numprocs);
   printf("# id of this processor: %d\n", MPI_myid);
   printf("\n");

   prog_start	=	time(NULL);			// to calculate running time

   printf("the memory needed for one mol is %d bytes\n", sizeof(molstruct));

   fflush(stdout);

   //-------------------//
   //	Initialization	//
   //-------------------//
   if (L2HDEBUG)	printf("DEBUG: Initialization ...\n");

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   InitUnits();				// initialize units for calculation
   InitForcefield();			// initialize Lennard-Jones potential mixing rule

   system	=	0;		// for now 2/14/08, only one system

   InitSample();			// initialize sampling

   nframe	=	-1;

   //------Initialize variable values------//

   for (i=0; i<180; i++) {
      orientdist[i]	=	0;	// chain orientation distribution
   }

   for (i=0; i<MAXNMOLSITES; i++) {
      xbeadpos[i]	=	0;	// crystal bead position
   }

   Nn		=	(long *) calloc(NSITES, sizeof(long));
   ntot		=	0;
   maxnmax	=	0;
   Rg2nucl	=	(float *) calloc(NSITES, sizeof(float));

   for (i=0; i<MAXNMOLSITES; i++) {	// segment length variables
      nlloop[i]		=	0;
      nltail[i]		=	0;
      nlbridge[i]	=	0;
      nlxseg[i]		=	0;
      nlfree[i]		=	0;
   }
   for (i=0; i<10; i++) {
      for (j=0; j<MAXNMOLSITES; j++) {
	 lsegdist[i][j]	=	0;
      }
   }

   //---------------------------//
   //	Start Data Processing	//
   //---------------------------//

   strcpy(filename, filein);
   strcat(filename, ".pre");
   if (!preflag) {
      if ((fpre=fopen(filename, "r"))==NULL )		// read .pre file
         Exit("lammps2hst", "main", "open .pre file failed.");
   }

   for (ifile=0; ifile<nfiles; ifile++) {		// multiple input dump files

      fin	=	fopen(infile[ifile], "r");

      if (L2HDEBUG) printf("# Current dump file: %s\n", infile[ifile]);

      while (!feof(fin)) {

	 if (L2HDEBUG)	printf("\nDEBUG: Reading one configuration ...\n");

	 //---------------------------------------------//
	 //	Read one configuration from dump file	//
	 //---------------------------------------------//
	 if (!fgets(dummy, sizeof(dummy), fin)) {	// line 1 of each conf. in dump file
	    break;					// end of file
	 }
	 fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);	// line 2
	 fgets(dummy, sizeof(dummy), fin);						// line 3
	 fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);	// line 4

	 fgets(dummy, sizeof(dummy), fin);						// line 5
	 if (strstr(dummy, "xy"))	PBC = 3;	// PBC=3 for triclinic box

	 fscanf(fin, "%lf%lf", &xlo, &xhi);	if (PBC==3) fscanf(fin, "%lf", &xy); 	// line 6
	 fgets(dummy, sizeof(dummy), fin);
	 fscanf(fin, "%lf%lf", &ylo, &yhi);	if (PBC==3) fscanf(fin, "%lf", &xz); 	// line 7
	 fgets(dummy, sizeof(dummy), fin);
	 fscanf(fin, "%lf%lf", &zlo, &zhi);	if (PBC==3) fscanf(fin, "%lf", &yz);	// line 8
	 fgets(dummy, sizeof(dummy), fin);
	 fgets(dummy, sizeof(dummy), fin);						// line 9

	 xlo	-=	MIN(0.0, MIN(xy, MIN(xz, xy+xz)));		// triclinic box in general
	 xhi	-=	MAX(0.0, MAX(xy, MAX(xz, xy+xz)));
	 ylo	-=	MIN(0.0, yz);
	 yhi	-=	MAX(0.0, yz);

	 system	=	0;
	 BOX[system].lx	=	xhi-xlo;
	 BOX[system].ly	=	yhi-ylo;
	 BOX[system].lz	=	zhi-zlo;
	 BOX[system].xy	=	xy;
	 BOX[system].xz	=	xz;
	 BOX[system].yz	=	yz;

	 temp1	=	MAX( MAX(BOX[system].lx, BOX[system].ly), BOX[system].lz);
	 temp3	=	MIN( MIN(BOX[system].lx, BOX[system].ly), BOX[system].lz);
	 if ( fabs(BOX[system].lx-temp1) > ZERO && fabs(BOX[system].lx-temp3) > ZERO) {
	    temp2	=	BOX[system].lx;
	 }
	 else if ( fabs(BOX[system].ly-temp1) > ZERO && fabs(BOX[system].ly-temp3) > ZERO) {
	    temp2	=	BOX[system].ly;
	 }
	 else if ( fabs(BOX[system].lz-temp1) > ZERO && fabs(BOX[system].lz-temp3) > ZERO) {
	    temp2	=	BOX[system].lz;
	 }
	 asphericity	=	temp1 * temp1 - 0.5*(temp2*temp2 + temp3*temp3);
	 asphericity	/=	(temp1 + temp2 + temp3) * (temp1 + temp2 + temp3) / 9;
    
	 LENGTH		=	NSITES/NMOLS;	// monodisperse system for now (4/26/2008)
	 accum		=	0;

	 for (i=0; i<nsites; i++) {
	    fscanf(fin, "%ld", &id);
	    fscanf(fin, "%ld", &molid);		// Need to check the lammps.dump file format
					   	// because some early lammps.dump file 
						// does not have molid output
						
	    // if (polydisperse)	fscanf(fin, "%ld", &molid);
	    fscanf(fin, "%ld", &type);
	    fscanf(fin, "%lf%lf%lf %lf%lf%lf %d%d%d", &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
	    fgets(dummy, sizeof(dummy), fin);

	    rsite[i].x	=	x;		// for density profile
	    rsite[i].y	=	y;
	    rsite[i].z	=	z;
    
	    if (polydisperse) {			// polydisperse
	       molid	--;
	       /*accum	=	0;
	       for (j=0; j<molid; j++)
		  accum	+=	(mol+j)->nsites;
	       siteid	=	id  - accum;
	       */
	    }
	    else {					// monodisperse
	       id	--; 				// -1 because lammps index starts from 1
	       molid	=	(long) (id/LENGTH);
	       siteid	=	id % LENGTH;
	    }
	    mol[molid].box	=	system;		// for now, only one system
	    mol[molid].nsites	=	LENGTH;		// for now, Jan/10/2010

	    // triclinic box in general
	    mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx) 
						  + ny*(BOX[system].xy) + nz*(BOX[system].xz);
	    mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly) + nz*(BOX[system].yz);
	    mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);

	    // velocity
	    mol[molid].velx[siteid]	=	vx;
	    mol[molid].vely[siteid]	=	vy;
	    mol[molid].velz[siteid]	=	vz;

	    // atom type 
	    mol[molid].type[siteid]=	type - 1;	// -1 because lammps index starts from 1
	 }	// One configuration read

	 //output_xyzfile(fxyz);		// output configurations using xyz format

	 unfold();				// unfold by pbc to make a chain continuous in space

	 for (i=0; i<NSYSTEMS; i++) {
	    NMols[i]	=	0;
	    NSites[i]	=	0;
	 }
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    if ( (i=moli->box) >= 0) {
	       NMols[i]	++;				// total # of mols in certain system
	       NSites[i]	+=	moli->nsites;	// total # of sites in certain system
	    }
	 }

	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    for (i=0; i<moli->nsites; i++)  {
	       moli->flags[i]	=	1;		// activate all the sites on this processor
	       moli->parent[i]	=	i-1;		// initialize parent site
	    }
	    moli->flip		=	0;		// flip to the original direction
	    moli->origin		=	CenterofMass(moli);
	 }

	 //-----------------------------------------------------//
	 //	Skip repeated frames from different input files	//
	 //-----------------------------------------------------//
	 if (ifile>=1 && timestep==previous_timestep) {
	    continue;
	 }
	 previous_timestep	=	timestep;	// check repeat frames

	 //-----------------------------//
	 //	Check timestep range	//
	 //-----------------------------//
	 if (timestep < starttime) {		// starttime=-1 by default
	    continue;
	 }
	 if (endtime > 0 && endtime < timestep) {	// endtime=-1 by default
	    break;
	 }

	 //-------------------------------------//
	 //	Skip (dnframe - 1) frames	//
	 //-------------------------------------//
	 nframe	++;
	 if (mod(nframe, dnframe))	continue;	// analyze every dnframe frames
	
         //-------------------------------------//
	 //	Convert to Z-code input file	//
         //-------------------------------------//
	 
	 if (Zflag && MPI_myid==0) {			// ONLY process 0 do the following
	    output_Zinput(fZfile);			// create init. conf. for Z code
	    continue;
	 }
	 

	 //-------------------------------------------------------------//
	 //	Convert coordinates and box size from SI to system unit	//
	 // 	-- For orthorombix box only as of 11/21/2011		//
	 //-------------------------------------------------------------//
	 CoorSI2System();		// Convert {x, y, z}'s, lx, ly, lz and lbox
	 
	 for (i=0; i<NSYSTEMS; i++) {	// Convert other lengths
	    BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
	    // repeatative, just in case .lbox was not converted before
	    BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
	    BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
	    BOX[i].rb		=	Rb;
	    BOX[i].rv		=	Rv;
	 } 

	 //-----------------------------//
	 //	Build cell list		//
	 //-----------------------------//
#ifdef CELL_LIST	

	 if (!clistinit) {
	    CL_Destroy();
	 }
	 clistinit=0;
	 CL_Init();		// need to init. cell list every time due to volume change
	 CL_Build();

#endif	/* CELL_LIST */

	 //-----------------------------------------------------//
	 //	Perform Analysis				//
	 //	-- Part 1: Coordinate and velocity analysis	//
	 // 	-- Part 2: Energy calculation			//
	 // 	-- Part 3: Crystal phase identification		//
	 // 	-- Part 4: Nucleus analysis			//
	 // 	-- Part 5: Correlation and s(k)			//
	 //-----------------------------------------------------//
	 if (L2HDEBUG)	printf("DEBUG: Analysis starts ...\n");

	 //-----------------------------------------------------------------------------//
	 //	1.1 - Calculate the average position of center of mass for each chain	//
	 //-----------------------------------------------------------------------------//
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	     rtemp	=	CenterofMass(moli);
	     if (fabs(rtemp.x) < 0.45 * BOX[moli->box].lx && 
		   fabs(rtemp.y) < 0.45 * BOX[moli->box].ly &&
		   fabs(rtemp.z) < 0.45 * BOX[moli->box].lz) {
		com[moli-mol]	=	V_Add(&rtemp, com+(moli-mol));
		ncom[moli-mol]	++;
	     }
	 }

	 //-----------------------------//
	 //	1.2 - Density profile	//	
	 //-----------------------------//
	 system	=	0;
	 temp	=	MAX(BOX[system].lx, MAX(BOX[system].ly, BOX[system].lz));
	 if (fabs(temp - BOX[system].lx) < ZERO)		direct	=	1;
	 else if (fabs(temp - BOX[system].ly) < ZERO)	direct	=	2;
	 else						direct	=	3;
	 /*
	 V_Null(&rtemp1);
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    for (i=0; i<moli->nsites; i++) {
	       rtemp	=	MapInBox2(moli->p+i, PBC, system);
	       rtemp1	=	V_Add(&rtemp1, &rtemp);
	    }
	 }
	 rtemp2	=	V_Mult(1.0/NSITES, &rtemp1);	// com of system
	 rtemp1	=	MapInBox2(&rtemp2, PBC, system);
	 */
	
	 /* 
	 for (i=0; i<NSITES; i++) {
	    rsite[i]	=	V_Mult(1/4.01, rsite+i);
	 }
	 V_Null(&rtemp1);
	 for (i=0; i<NSITES; i++) {
	    rtemp1	=	V_Add(&rtemp1, rsite+i);
	 }
	 rtemp1	=	V_Mult(1.0/NSITES, &rtemp1);
	 
	 printf("%lf\t%lf\t%lf\n", rtemp1,x, rtemp1,y, rtemp1.z);

	 for (i=0; i<DEN_BIN; i++)	density[i]	=	0;

	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    for (i=0; i<moli->nsites; i++) {
	       rtemp	=	MapInBox2(moli->p+i, PBC, system);
	       //rtemp	=	V_Subtr(&rtemp, &rtemp1);

	       switch(direct){
		  case 1:		temp1	=	rtemp.x;	break;
		  case 2:		temp1	=	rtemp.y;	break;
		  case 3:		temp1	=	rtemp.z;	break;
	       } 
	       density[(long) ((temp1/temp+0.5)*DEN_BIN)]	++;
	    }
	 }
	 */

	 /*
	 for (i=0; i<NSITES; i++) {
	    switch(direct){
	       case 1:		temp1	=	rsite[i].x-rtemp1.x;	break;
	       case 2:		temp1	=	rsite[i].y-rtemp1.y;	break;
	       case 3:		temp1	=	rsite[i].z-rtemp1.z;	break;
	    } 
	    density[(long) ((temp1/temp+0.5)*DEN_BIN)]	++;
	 }
	 */
	 
	 lam1	=	0;
	 amo1	=	0;
	 for (i=0; i<DEN_BIN; i++) {
	    if (density[i] > NSITES/DEN_BIN)	lam1	++;
	    else				amo1	++;
	 }
	 lam1	*=	temp/DEN_BIN;
	 amo1	*=	temp/DEN_BIN;

	 //-----------------------------------------------------------------------------//
	 //	1.3 - Compute chain orientation distribution, test the rotator phase	//
	 //-----------------------------------------------------------------------------//
	 if (L2HDEBUG)	printf("DEBUG: Computing chain orientation ...\n");

	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    V_Null(&chainface);

	    for (i=0; i<moli->nsites; i++) {
	       if (mod(i,2)) {
		  chainface	=	V_Add(&chainface, moli->p+i);
	       }
	       else {
		  chainface	=	V_Subtr(&chainface, moli->p+i);
	       }
	    }
	    orient	=	atan2(chainface.y, chainface.x)	+ M_PI;
	    //rtemp=CenterofMass(moli);
	    //if (rtemp.z < -0.25 * BOX[system].lz && rtemp.z <0) 
	    orientdist[(long)(orient/(2*M_PI)*36)]	++;
	 }

	 //---------------------------------------------//
	 //	1.4 - Spherical coordinate calculation	//
	 //---------------------------------------------//
	 if (L2HDEBUG)	printf("DEBUG: Sampling spherical coordinates ...\n");
	 SampleSpherical();			// sample spherical coordinates
	 Dist_Spherical();			// put spherical coord. into distribution
	 
	 // Average Ptt and Pttt

	 avePtt		+=	Ptors2[0][0][0];
	 avePttt	+=	Ptors3[0][0][0][0];
	 stdPtt		+=	Ptors2[0][0][0] * Ptors2[0][0][0];
	 stdPttt	+=	Ptors3[0][0][0][0] * Ptors3[0][0][0][0];

	 //-------------------------------------//
	 //	1.5 - Velocity analysis		//
	 //-------------------------------------//
	 ave_vx	=	0.0;
	 ave_vy	=	0.0;
	 ave_vz	=	0.0;
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    for (i=0; i<moli->nsites; i++) {
	       ave_vx	+=	moli->velx[i];  
	       ave_vy	+=	moli->vely[i];  
	       ave_vz	+=	moli->velz[i];  
	    }
	 }
	 ave_vx	/=	NSITES;
	 ave_vy	/=	NSITES;
	 ave_vz	/=	NSITES;
	 vcom2	=	ave_vx * ave_vx + ave_vy * ave_vy + ave_vz * ave_vz;

	 v2	=	0;
	 vf2	=	0;
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    for (i=0; i<moli->nsites; i++) {
	       vx	=	moli->velx[i];
	       vy	=	moli->vely[i];
	       vz	=	moli->velz[i];
	       v2	+=	vx*vx + vy*vy + vz*vz;

	       temp1	=	MAX(vx, MAX(vy, vz));
	       temp2	=	MIN(vx, MIN(vy, vz));

	       if (temp1 > velmax)	velmax	=	temp1;
	       if (temp2 < velmin)	velmin	=	temp2;
	 
	       ltemp1	=	(vx>=0 ? (long) (vx/dvel+0.5) : (long) (vx/dvel-0.5));
	       vel_distx[ltemp1+50]	++;
	       ltemp1	=	(vy>=0 ? (long) (vy/dvel+0.5) : (long) (vy/dvel-0.5));
	       vel_disty[ltemp1+50]	++;
	       ltemp1	=	(vz>=0 ? (long) (vz/dvel+0.5) : (long) (vz/dvel-0.5));
	       vel_distz[ltemp1+50]	++;

	       vx	=	moli->velx[i]-ave_vx;
	       vy	=	moli->vely[i]-ave_vy;
	       vz	=	moli->velz[i]-ave_vz;
	       vf2	+=	vx*vx + vy*vy + vz*vz;
	    }
	 }

	 //---------------------------------------------//
	 //	2.1 - Calculate energy and pressure	//
	 //---------------------------------------------//
	 // Only perform when creating .pre file because energy calc. is expensive

	 if (preflag) {					// perform pre-analysis
	    if (L2HDEBUG)	printf("DEBUG: Computing energy and virial ...\n");

	    CalcV();					// calc energy and virial
	    if (V_VIRIAL) {
	       Sample_Pressure();
	    }
	 }
	 
	 //-----------------------------//
	 //	3.1 - Calculate p2(i)	//
	 //-----------------------------//
	 // Only perform when creating .pre file because p2 calc. is expensive
	 
	 if (preflag) {					// generate pre-analysis results
							// Most time-consuming part is
							// energy and p2 calc. as they require
							// pair distance calc. for all pairs.

	    if (L2HDEBUG)	printf("DEBUG: Sampling p2 ...\n");
	    SampleP2All_MPI();				// sample P2 and P2m and local p2
	    Dist_p2();					// put local p2 into distribution
	    //SampleConnection();
	    //Calc_Qlm(6);				// calc. Qlm for LJ system, require CELL_LIST
	 }				
	 else {						// if pre-analysis exists
	    if (L2HDEBUG)	printf("DEBUG: Reading pre-analysis results ...\n");
	    read_pre(fpre, timestep);			// read pre-analysis results
	    
	    system	=	0;			// one system only for now
	    oldxtal	=	Xtal[system];		// store old nucleus variables
	    oldnnucl	=	Nnucl[system];		// ... because xtal_smooth() later
	    oldnmax	=	nmax[system][0]; 	// ... will change these quantities
	    old2nmax	=	nmax[system][1]; 

	    if (L2HDEBUG)	printf("DEBUG: Perform xtal smoothing ...\n");
	    for (i=0; i<4; i++) {			// perform xtal segment smoothing until convergence
	       xtal_smooth();				// ONLY perform when reading data from .pre file
	    }
	 }
	 
	 //-------------------------------------------------------------//
	 //	3.2 - Average xtal bead position on chain 3/19/2012	//
	 //-------------------------------------------------------------//
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    for (i=0; i<moli->nsites; i++) {
	       if (moli->p2[i] > critp2) {
		  xbeadpos[i]	++;
	 }  }  }

	 //-------------------------------------//
	 //	3.3 - End-to-end distance	//
	 //-------------------------------------//
	 R02	=	0.0;
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    R02	+=	R2_n2n(moli);
	 }
	 R02	/=	NMOLS;

	 //---------------------------------------------//
	 //	3.4 - Radius of gyration analysis	//
	 //---------------------------------------------//
	 aveRg2		=	0.0;	nRg2	=	0;
	 avexRg2	=	0.0;	nxRg2	=	0;

	 avergx2	=	0.0;
	 avergy2	=	0.0;
	 avergz2	=	0.0;

	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    j		=	moli-mol;
	    Rg2[j]	=	R2_gyration(moli);
	    aveRg2	+=	Rg2[j];
	    nRg2	++;

	    rgx2[j]	=	rx2_gyration(moli);
	    rgy2[j]	=	ry2_gyration(moli);
	    rgz2[j]	=	rz2_gyration(moli);

	    avergx2	+=	rgx2[j];
	    avergy2	+=	rgy2[j];
	    avergz2	+=	rgz2[j];

	    for (i=0; i<moli->nsites; i++) {
	       if (moli->p2[i] > critp2) {		// for chains have xtal segments
		  avexRg2	+=	Rg2[j];
		  nxRg2	++;
		  break;
	       }
	    }
	 }
	 avemRg2	=	(aveRg2 - avexRg2)/(nRg2 - nxRg2);
	 aveRg2		/=	nRg2;
	 avexRg2	/=	nxRg2;

	 avergx2	/=	nRg2;
	 avergy2	/=	nRg2;
	 avergz2	/=	nRg2;

	 printf("Rg %d %d %d %d %f %d %f  ", timestep, nmax[0][0], nmax[0][1], nRg2, aveRg2, nxRg2, avexRg2);
	 printf(" %f %f %f\n", avergx2, avergy2, avergz2);

	 //---------------------------------------------//
	 //	4.1 - Find nuclei based on p2(i)	//
	 //---------------------------------------------//
	 if (L2HDEBUG)	printf("DEBUG: Finding nuclei ...\n");
	 
	  
	 for (i=0; i<2; i++) {			// perform xtal segment smoothing until convergence
	    xtal_smooth();			// ONLY perform when reading data from .pre file
	    					// ... because it changes local p2 values
	 }
	 
	 Find_Nuclei_p2(1);				// always perform w or w/o .pre file

	 //********************************************************//
	 //	4.2 - Calculate the nucleus size distribution	//
	 //	      -- and the Rg of nuclei			//
	 //********************************************************//

	 system	=	0;
	 for (i=1; i<= nmax[system][0]; i++) {
	    Nn[i]	+=	sizedist[i];
	    ntot	+=	sizedist[i];
	 } 
	 maxnmax	=	MAX(maxnmax, nmax[system][0]);

	 for (nuclid=1; nuclid<=Nnucl[system]; nuclid++) {	// for all nuclei
	    size	=	sizeofnucl[nuclid];
	    m	=	0;
	    for (moli=mol; moli<mol+NMOLS; moli++) {
	       if (system == moli->box) {
		  for (i=0; i<moli->nsites; i++) {
		     if (moli->nuclid[i] == nuclid) {
			nucleus[m].moli	=	moli;
			nucleus[m].site	=	i;
			m	++;
	    }  }  }  }

	    if (size != m) {printf("Error! search nucleus size != m\n");}

	    MapInNucleus(system, nucleus, size, rbead);
	    Rg2nucl[size]	+=	group_Rg2(rbead, size);
	 }
    
	 //-------------------------------------------------------------//
	 //	4.3 - Find the nucleus id for the largest nucleus	//
	 //-------------------------------------------------------------//
	 if (L2HDEBUG)	printf("DEBUG: Finding nucleus id for the largest nucleus ...\n");

	 system	=	0;		// for now only one system, 4/16/2010
	 i	= 1;
	 while (sizeofnucl[i] != nmax[system][0]) {
	    i ++;
	 }
	 nmaxid	=	i;
      
	 //-------------------------------------------------------------//
	 //	4.4 - Group the beads belonging to the biggest nucleus	//
	 //-------------------------------------------------------------//
	 if (L2HDEBUG)	printf("DEBUG: Grouping the beads in the largest nucleus ...\n");
	 
	 m	=	0;
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    if (system == moli->box) {
	       for (i=0; i<moli->nsites; i++) {
		  if (moli->nuclid[i] == nmaxid) {
		     nucleus[m].moli	=	moli;
		     nucleus[m].site	=	i;
		     m	++;
		  }
	       }
	    }
	 }
	 if (m!= nmax[0][0]) {
	    printf("4.4 error!\n");
	 }

	 if (m>1) {
	    if (L2HDEBUG)	printf("DEBUG: Making the largest nucleus continuous in space ...\n");
	    MapInNucleus(system, nucleus, m, rbead);	// make the largest nucleus continuous in space
	 }

	 // cylindershape(nucleus, nmax[system][0], nmaxid); 

	 //------Calculate how many chains participating in the largest nucleus------//
	 if (L2HDEBUG)	printf("DEBUG: Calculating the number of chains in the largest nucleus ...\n");

	 nchainmax	=	0;
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    if (system == moli->box) {
	       for (i=0; i<moli->nsites; i++) {
		  if (moli->nuclid[i] == nmaxid) {
		     nchainmax	++;
		     break;
		  }
	       }
	    }
	 }

	 //********************************************************//
	 //	4.x - Calculate the drift of biggest nucleus	//
	 //********************************************************//

	 /*
	 if (timestep==0) {
	    rO		=	CoM_MaxNucleus(0);
	    rshift2	=	0;
	 }
	 else {
	    con	=	CoM_MaxNucleus(0);
	    con	=	V_Subtr(&con, &rO);
	    con	=	MapInBox2(&con, PBC, system);
	    rshift2	=	V_Dot(&con, &con);
	 }
	 */

	 //-------------------------------------//
	 //	4.5 - Segment analysis of nmax	//
	 //-------------------------------------//
	 if (L2HDEBUG)	printf("# Finding segments ...\n");
	 
	 Find_segments();				// Find all segments
	 //smoothresult	=	Seg_smooth(nmaxid);	// Reduce fluctuation for all segments 

	 if (L2HDEBUG)	printf("# Identifying type of each segment ...\n");
	 
	 printf("Seg %d ", timestep);
         Seg_type(nmaxid);
	 //Seg_length(nmaxid, lsegdist);
	 //Rho_profile();
         
	 /*
	 nloop		=	0;			// # of loop 
	 nbridge	=	0;			// # of bridge
	 ntail		=	0;			// # of tail
	 nxseg		=	0;			// # of xtal segment
	 npbrdg		=	0;			// # of pbc bridge

	 lloop		=	0.0;			// average length of loop
	 ltail		=	0.0;			// average length of tail
	 lbridge	=	0.0;			// average length of bridge
	 lxseg		=	0.0;			// average length of xtal segment
	 lpbrdg		=	0.0;			// average length of pbc bridge
	 
	 system	=	0;			// one system only for now

	 printf("Seg  %-6ld", timestep);
	 printf(" ___");
	 printf(" %4d %4d %4d %4d", oldxtal, oldnmax, old2nmax, oldnnucl);
	 printf(" ___");
	 printf(" %4ld %4ld %4ld %4ld %3d", Xtal[system], nmax[system][0], nmax[system][1], Nnucl[system], nchainmax);
	 printf(" ___");
	 printf(" %3d %3d %3d %3d %3d", ntail, nloop, nbridge, npbrdg, nxseg);
	 printf(" ___");
	 printf(" %5.1f %5.1f %5.1f %5.1f %5.1f\n", ltail/ntail, lloop/nloop, lbridge/nbridge, lpbrdg/npbrdg, lxseg/nxseg);
         */

	 //********************************************************//
	 //	x.x - Mean-squared displacement analysis	// 
	 //	      -- added on 2/14/2012			// 
	 //********************************************************//

	 if (cominit) {
	    V_Null(&com_sysinit);

	    for (moli=mol; moli<mol+NMOLS; moli++) {
	       j		=	moli-mol;
	       com_init[j]	=	CenterofMass(moli);
	       com_last[j]	=	com_init[j];
	       com_sysinit	=	V_Add(&com_sysinit, com_init+j);
	    }
	    com_sysinit	=	V_Mult(1.0/NMOLS, &com_sysinit);
	    com_syslast	=	com_sysinit;

	    cominit	=	0;
	 }

	 d2initave	=	0.0;			// squared displacement of com of chains
	 d2lastave	=	0.0;
	 d2sys	=	0.0; 			// squared displacement of com of system
	 d2syslast	=	0.0;

	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    j		=	moli-mol;
	    com[j]		=	CenterofMass(moli);
	    rtemp		=	V_Subtr(com+j, com_init+j);
	    d2initave	+=	V_Dot(&rtemp, &rtemp);
	    rtemp		=	V_Subtr(com+j, com_last+j);
	    d2lastave	+=	V_Dot(&rtemp, &rtemp);
	    com_last[j]	=	com[j];

	    com_sys	=	V_Add(&com_sys, com+j);
	 }
	 d2initave	/=	NMOLS;			// take average
	 d2lastave	/=	NMOLS;			// take average

	 com_sys	=	V_Mult(1.0/NMOLS, &com_sys);
	 dr_sys	=	V_Subtr(&com_sys, &com_sysinit);
	 d2sys	=	V_Dot(&dr_sys, &dr_sys);

	 rtemp	=	V_Subtr(&com_sys, &com_syslast);
	 d2syslast	=	V_Dot(&rtemp, &rtemp);
	 com_syslast	=	com_sys;

	 //************************************************************************//
	 //	x.x - Calc. nmax and 2ndnmax using other cutoff parameters	//
	 //************************************************************************//

	 if (L2HDEBUG)	printf("DEBUG: Use other cutoffs for nuclei definition ...\n");

	 nmaxp2_3[0]	=	nmax[0][0];	// Rp=2.5, Rconn=1.3, critp2=0.4
	 nmaxp2_3[1]	=	nmax[0][1];

	 /*
	 temp1	=	Rp;		// store original criterion
	 temp2	=	Rconn;
	 temp3	=	critp2;

	 Rp	=	1.5;		// new definition
	 Rconn	=	1.5;
	 critp2	=	0.6;
	 SampleP2All();			// sample P2 and P2m and local p2
	 Find_Nuclei_p2(1);
	 nmaxp2_1[0]	=	nmax[0][0];
	 nmaxp2_1[1]	=	nmax[0][1];

	 Rp	=	temp1;		// restore the original criterion
	 Rconn	=	temp2;
	 critp2	=	temp3;
	 */

	 //****************************************************************//
	 //	5.1 - Correlations and Radial distribution function	//
	 //****************************************************************//

	 if (L2HDEBUG)	printf("DEBUG: Calculate correlations ...\n");

	 correlation();			// calculate correlation, in system units

	 //if (!(timestep%IRADIAL)) {
	 //radial("sample");		// sample radial distribution function
	 //sq4("sample");
	 //         sq(stdout, "sample");
	 //}

	 //shiftbox(system, nucleus, n);

	 //********************************************************//
	 //	Output analysis results after each timestep	//
	 //********************************************************//
	 
	 if (L2HDEBUG)	printf("DEBUG: Output histogram / pre-analysis ...\n");

	 if (preflag) {
	    Print_hst(fhst);		// print out histgram
	    output_pre(fpre);		// must be after Print_hst
	 }
	 else {
	    output_Rg();			// output Rg analysis
	 }

	 //*********************************************************************//
	 //	Convert coordinates and box size from system to SI unit		//
	 //*********************************************************************//

	 if (L2HDEBUG)	printf("DEBUG: System2SI ...\n");

	 CoorSystem2SI();

	 //****************************************//
	 //	Output configuration files	//
	 //****************************************//

	 if (L2HDEBUG)	printf("DEBUG: Output configuration files ...\n");
	 
	 // OUTPUT .xyz file for visualization
	 if (xyzflag) {
	    output_xyzfile(fxyz);

	    for (moli=mol; moli<mol+NMOLS; moli++) {
	       for (i=0; i<moli->nsites; i++) {
		  fprintf(fcolor, "%f ", (moli->nuclid[i]==nmaxid ? 0.7 : 0.0));
	       }
	    }
	    fprintf(fcolor, "\n");
	 }

	 
	 // OUTPUT configuration file for further analysis
	 if (confflag) { 
	    output_conf(fconf, timestep);
	 }
	 // OUTPUT .car file for VMD visualization
	 if (carflag) {
	    output_carfile(fcar);
	 }
	 // OUTPUT .pdb file for further analysis
	 if (pdbflag) {
	    output_pdbfile(fpdb, fdat, imagen, imagex, imagey, imagez);
	 }

      } 			// finish ONE input dump file

      fclose(fin);		// close current input dump file

   }				// finish ALL input dump files

   //-------------------------------------------------------------------//
   //	Output final analysis results after ALL frames processed	//
   //-------------------------------------------------------------------//
   if (L2HDEBUG)	printf("DEBUG: ALL FRAMES PROCESSED ...\n");

   //------Frame information------//

   printf("\n# OUTPUT AFTER ALL FRAMES ARE PROCESSED\n\n");
   printf("START TIMESTEP\t%d\n", starttime);
   printf("END TIMESTEP\t%d\n", endtime);
   printf("TOTAL FRAMES\t%d\n", nframe+1);
   printf("ANALYZE EVERY\t%d\n", dnframe);
   printf("ANALYZED FRAMES\t%d\n", nframe/dnframe+1);
   
   //------Output density profile------//
   
   printf("\n# DENSITY PROFILE\n\n");
   for (i=0; i<DEN_BIN; i++) {
      printf("%ld\t%lf\n", i, density[i]);
   }
   printf("%d\t%lf\n", DEN_BIN, density[0]);
 
   //------Output velocity stat.------//

   printf("\n# VELOCITY PROFILE\n\n");
   printf("velmax= %lf\t, velmin= %lf\n", velmax, velmin);
   printf("v2= %lf\t, vf2= %lf\t, vcom2= %lf\n", v2, vf2, vcom2);
   printf("v2= %lf\t, vf2= %lf\t, vcom2= %lf\n", v2, vf2, vcom2*NSITES);
   for (i=0; i<101; i++) {
      printf("%lf\t%lf\t%lf\n", vel_distx[i], vel_disty[i], vel_distz[i]);
   }

   //------Output correlations------//

   printf("\n# TIME CORRELATION FUNCTIONS\n\n");
   corr_norm();			// normalize correlation function
   corr_print();		// print out correlation function

   //------Output distributions------//

   printf("\n# DISTRIBUTIONS\n\n");
   S_PrintAll();		// print out distributions
   //radial("print");		// print out radial distribution function
				// -- has some problem, segmentation fault
   //sq4("print");
   //sq(stdout, "print");	// print out structural factor

   /*
   for (i=0; i<180; i++) {
      printf("%ld %ld\n", i, orientdist[i]);
   }
   */

   //------Output torsion angle statistics------//

   printf("\n# TORSION STATISTICS\n");
   printf(" ## LAST FRAME\n");

   for (i=0; i<NSYSTEMS; i++) {
      for (j=0; j<3; j++) {
         for (k=0; k<3; k++) {
            Ptors2[i][j][k]	/=	((NSITES/NMOLS)-4)*NMOLS;
         }
      }
      Ptors3[i][0][0][0]	/=	((NSITES/NMOLS)-5)*NMOLS;

      printf("P\tt\t\tg+\t\tg-\n");
      printf("t\t%f\t%f\t%f\n", Ptors2[i][0][0], Ptors2[i][0][1], Ptors2[i][0][2]);
      printf("g+\t%f\t%f\t%f\n", Ptors2[i][1][0], Ptors2[i][1][1], Ptors2[i][1][2]);
      printf("g-\t%f\t%f\t%f\n", Ptors2[i][2][0], Ptors2[i][2][1], Ptors2[i][2][2]);
      printf("\nPttt = %f\n", Ptors3[i][0][0][0]);
   }

   printf("\n ## AVERAGE\n");
   avePtt	/=	((NSITES/NMOLS)-4)*NMOLS * (nframe/dnframe+1);
   avePttt	/=	((NSITES/NMOLS)-5)*NMOLS * (nframe/dnframe+1);   
   stdPtt	/=	((NSITES/NMOLS)-4)*NMOLS * (nframe/dnframe+1);
   stdPtt	/=	((NSITES/NMOLS)-4)*NMOLS;
   stdPttt	/=	((NSITES/NMOLS)-5)*NMOLS * (nframe/dnframe+1);   
   stdPttt	/=	((NSITES/NMOLS)-5)*NMOLS;   
   stdPtt	=	sqrt(stdPtt - avePtt*avePtt);
   stdPttt	=	sqrt(stdPttt - avePttt*avePttt);

   printf("Ptt = %f +- %f\t Pttt = %f +- %f\n", avePtt, stdPtt, avePttt, stdPttt);

   //------Output nucleus size distribution------//

   /*
   printf("\n# AVERAGED NUCLEUS SIZE DISTRIBUTION\n\n");
   printf(" ntot = %ld\n\n", ntot);
   printf("nsize\tN[size]\tP[size]\t\t-ln(P[size])\tRg2[size]\n");
   for (i=1; i<=maxnmax; i++) {
      printf("%ld\t%ld\t%f\t%f\t%f\n", 
		i,  Nn[i], ((double)Nn[i])/ntot, -log(((double)Nn[i])/ntot), Rg2nucl[i]/Nn[i]);
   }
   */

   //------Output xtal bead position------//

   printf("\n# XTAL BEAD POSITION STATISTICS\n\n");
   ltemp1	=	0;
   for (i=0; i<NSITES/NMOLS; i++) {
       ltemp1	+=	xbeadpos[i];
   }
   for (i=0; i<NSITES/NMOLS/10; i++) {			// every 10 beads
       ltemp2	=	0;
       for (j=0; j<10; j++) {
          ltemp2	+=	xbeadpos[i*10+j];
       }
       printf("%ld-%ld \t%ld \t%f\n", i*10, i*10+9, ltemp2, (double)ltemp2/ltemp1);
   }

   //------Output segment analysis results------//

   printf("\n# SEGMENT LENGTH STATISTICS OF LAMELLA\n");
   printf("\n ## -- for \"equilibrium\" segment length distribution\n\n");
   /* 
   printf("length\tnltail\tnlloop\tnlbrdg\tnlpbrdg\tnlxseg\n");
   for (i=0; i<NSITES/NMOLS; i++) {
      printf("segdist %ld\t %d\t %d\t %d\t %d\t %d\n", 
	i, nltail[i], nlloop[i], nlbridge[i], nlpbrdg[i], nlxseg[i]);
   }
   */
   for (i=1; i<NSITES/NMOLS; i++) {
      printf("lseg %d", i);
      for (j=2; j<8; j++) {
	 printf(" %d", lsegdist[j][i]);
      }
      printf("\n");
   }
   /*
   printf("\n### CHAIN SEGMENT INFORMATION\n\n");
   for (i=0; i<NMOLS; i++) {
      if (nsegment[i] == 1) {			// condition for output
         printf("chain %ld has %ld segments\n", i, nsegment[i]);
         for (j=0; j<(mol+i)->nsites; j++) {
            if (segment[i][j]>0 )
               printf("%ld_%ld  ", j, segment[i][j]);
         }
         printf("\n");
         for (j=0; j<nsegment[i]; j++) {
            printf("head = %ld\t id = %ld\t tail = %ld\t length = %ld\n", 
		seg_stat[i][j*3], seg_stat[i][j*3+1], seg_stat[i][j*3+2], seg_stat[i][j*3+2]-seg_stat[i][j*3]+1);
         }
         printf("\n");
      }
   }
   */

   printf("\n# BIGGEST NUCLEUS SEGMENT ANALYSIS\n\n");
   for (moli=mol; moli<mol+NMOLS; moli++) {
      k	=	moli-mol;
      
      for (m=0; m<nsegment[k]; m++) {
         if (seg_stat[k][m*3+1] == nuclid) {		// participate in the BIGGEST nucleus
            printf("chain #%ld\n", k);
            printf("type\t head\t id\t tail\t length\n\n");

            for (j=0; j<nsegment[k]; j++) {
               if (seg_stat[k][j*3+1]==-1) {
		  head	=	(j==0 ? -1 : seg_stat[k][(j-1)*3+1]);	// type of nuclei on one end
                  tail	=	(j==nsegment[k]-1 ? -1 : seg_stat[k][(j+1)*3+1]);	// on the other end

                  if ((head==-1 && tail==nuclid) || (head==nuclid && tail==-1)) {	// for nuclid
                  //if (head*tail < 0) {		// for all nuclei
                     sprintf(s, "tail\t: ");				// it's a tail
                  }
                  else if (head==nuclid && tail==nuclid) {
                  //else if (head == tail) {
		     nx  =	(int)(moli->p[seg_stat[k][j*3]-1].x / BOX[system].lx + 0.5);
		     ny	 =	(int)(moli->p[seg_stat[k][j*3]-1].y / BOX[system].ly + 0.5);
		     nz	 =	(int)(moli->p[seg_stat[k][j*3]-1].z / BOX[system].lz + 0.5);

		     nx1 =	(int)(moli->p[seg_stat[k][j*3+2]+1].x / BOX[system].lx + 0.5);
		     ny1 =	(int)(moli->p[seg_stat[k][j*3+2]+1].y / BOX[system].ly + 0.5);
		     nz1 =	(int)(moli->p[seg_stat[k][j*3+2]+1].z / BOX[system].lz + 0.5);

		     if (nx==nx1 && ny==ny1 && nz==nz1) {	// if two ends in same PBC box
			sprintf(s, "loop\t: ");			// it's a loop
		     }
		     else {					// if two ends in different PBC boxs
			sprintf(s, "pbcbridge\t:");		// it's a bridge
			printf("%f %d %f %d %f %d %f %d %f %d %f %d\n", 
			   moli->p[seg_stat[k][j*3]-1].x, nx, moli->p[seg_stat[k][j*3+2]+1].x, nx1,
			   moli->p[seg_stat[k][j*3]-1].y, ny, moli->p[seg_stat[k][j*3+2]+1].y, ny1,
			   moli->p[seg_stat[k][j*3]-1].z, nz, moli->p[seg_stat[k][j*3+2]+1].z, nz1);
		     }
                  }
                  else if (head==nuclid || tail==nuclid) {
                  //else if (head != tail) {
		     sprintf(s, "bridge\t:");		// it's a bridge
                  }   
               }
	       else {
		  sprintf(s, "xseg\t:");
	       }
               printf("%s%d\t%d\t%d\t%d\n", s,
		seg_stat[k][j*3], seg_stat[k][j*3+1], seg_stat[k][j*3+2], seg_stat[k][j*3+2]-seg_stat[k][j*3]+1);
            }
            printf("\n");
            break;
         }
      }
   }

   //------Average center of mass of chains------//
   system	=	0;
   if (pdbflag) {
      fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
      fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
 	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

      n		=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         i	=	moli-mol;
         if (ncom[i]>0) {
	 //if( fabs(com[i].x/ncom[i])<3 && fabs(com[i].y/ncom[i])<3 && fabs(com[i].z/ncom[i])<7 ) {
            n	++;
            atomname	=	'N';
            fprintf(fpdb, "ATOM  ");			// pdb command, column 1-6
            fprintf(fpdb, "%5ld ", n);			// atom number
            fprintf(fpdb, " %c  ", atomname);		// atom name
            fprintf(fpdb, " ");				// alternate location indiator
            fprintf(fpdb, "   ");			// residue name
            fprintf(fpdb, " ");				// column 21
            fprintf(fpdb, " ");				// chain identifier, column 22
            fprintf(fpdb, "    ");			// residue sequence number, 23-26
            fprintf(fpdb, " ");				// code for insertion of residues, 27
            fprintf(fpdb, "   ");			// column 28-30
            fprintf(fpdb, "%8.3f%8.3f%8.3f", com[i].x/ncom[i]*unit.LENGTH, 
			com[i].y/ncom[i]*unit.LENGTH, com[i].z/ncom[i]*unit.LENGTH);
            fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
            fprintf(fpdb, "\n"); 
         }
	 //}
      }
      fprintf(fpdb, "END\n");
      fflush(fpdb);
   }
   //______OUTPUT average center of mass of chains to .pdb file__//

   //-------------------//
   //	Closing files	//
   //-------------------//
   if (L2HDEBUG)	printf("DEBUG: Close output files ...\n");

   prog_end	=	time(NULL);
   
   printf("\n# END OF OUTPUT FILE\n");
   printf("%s", asctime(localtime(&prog_end)));
   printf("\nTotal running time was %lf seconds.\n", difftime(prog_end, prog_start));

   if (preflag)	{
      fclose(fhst);
   }
   fclose(fpre);
   if (xyzflag) {
      fclose(fxyz);
      fclose(fcolor);
   }
   if (carflag)	{
      fclose(fcar);
   }
   if (confflag) {
      fclose(fconf);
   }
   if (pdbflag) {
      fclose(fpdb); 
      fclose(fdat);
   }

   if (Zflag && MPI_myid==0) {				// ONLY process 0 do the following
      fclose(fZfile);
   }
	 
   fclose(stdout);

#ifdef myMPI
   MPI_Finalize();
#endif
   return	0;
}	
//------MAIN program ends------//


//======================//
//	SUBROUTINES	//
//======================//


//==============================================//
//	Output variables for each chain 	//
//	-- Rg, segment info. etc.		//
//	-- output format made similar		//
//	-- to .pre, so that they can		//
//	-- be easily merged			//
//==============================================//
void output_Rg()
{
   static FILE	*fPtr, *fdist;
   molstruct	*moli;
   long		i, j, k, system=0; 
   long		maxid=0, maxRgid=0, minid=0, minRgid=0;
   double	max=0, maxRg=0, min=1e6, minRg=1e6;
   long		head, tail, length, same;
   long		nloop[MAXNMOLS], 	// # of folds for each chain
		nbridge[MAXNMOLS], 
		ntail[MAXNMOLS], 
		nxseg[MAXNMOLS], 
		nfree[MAXNMOLS];
   long		nloop_tot, 		// total # of loop for the system
		nbridge_tot, ntail_tot, 
		nxseg_tot, nfree_tot;
   long		nxbead[MAXNMOLS];	// # of xtal beads on each chain
   long		nnucl[MAXNMOLS];	// # of nuclei each chain enters
   static int 	init=1;

   if (init) {
      init	=	0;
      fPtr	=	fopen("Rg.hst", "w");
      fdist	=	fopen("Rg.dist", "w");

      //=====	Print header in file Rg.hst	=====//

      fprintf(fPtr, "Timestep");
      fprintf(fPtr, " P2");
      fprintf(fPtr, " Xtal");
      fprintf(fPtr, " nmax");
      fprintf(fPtr, " 2ndnmax");
      fprintf(fPtr, " aveR02");		// average R02 of whole system
      fprintf(fPtr, " aveRg2");		// average Rg2 of whole system
      fprintf(fPtr, " avexRg2");	// average Rg2 of chains containing at least one xseg
      fprintf(fPtr, " avemRg2");	// average Rg2 of chains containing at least one xseg

      fprintf(fPtr, " <10>");		// formatting lines

      fprintf(fPtr, " smooth.y");   
      fprintf(fPtr, " smooth.z");   
      fprintf(fPtr, " ntail_tot");
      fprintf(fPtr, " nloop_tot");
      fprintf(fPtr, " nbridge_tot");
      fprintf(fPtr, " nxseg_tot");
      fprintf(fPtr, " nfree_tot");

      fprintf(fPtr, " <18>");		// formatting lines

      fprintf(fPtr, " chain1");		// ID of chain 1, the one has most xtal beads
      fprintf(fPtr, " xchain1");	// # of xtal beads on chain 1
      fprintf(fPtr, " Rg2chain1");	// Rg2 of chain 1
/*
      fprintf(fPtr, " chain2");		// ID of chain 2, the one has 2nd most xtal beads
      fprintf(fPtr, " xchain2");	// # of xtal beads on chain 2
      fprintf(fPtr, " Rg2chain2");	// Rg2 of chain 2
*/
      fprintf(fPtr, " chain3");		// ID of chain 3, the one has smallest Rg2
      fprintf(fPtr, " xchain3");	// # of xtal beads on chain 3
      fprintf(fPtr, " Rg2chain3");	// Rg2 of chain 3

      fprintf(fPtr, " chain4");		// ID of chain 4, the one has largest Rg2
      fprintf(fPtr, " xchain4");	// # of xtal beads on chain 4
      fprintf(fPtr, " Rg2chain4");	// Rg2 of chain 4

      fprintf(fPtr, " <28>");		// formatting lines

      fprintf(fPtr, " aspher");		// aphericity
      fprintf(fPtr, " dr_sys.x");
      fprintf(fPtr, " dr_sys.y");
      fprintf(fPtr, " dr_sys.z");
      fprintf(fPtr, " d2sys");

      fprintf(fPtr, " KEtot KEpec KEcom");
      fprintf(fPtr, "\n");

      //=====	Print header in Rg.dist file	=====//

      fprintf(fdist, "Timestep");
      fprintf(fdist, "\tnmax");
      fprintf(fdist, "\tchain");
      fprintf(fdist, "\tRg2");
      fprintf(fdist, "\tnxbead");
      fprintf(fdist, "\tnnucl");
      fprintf(fdist, "\tnsegment");
      fprintf(fdist, "\tnbridge");
      fprintf(fdist, "\tnxseg");
      fprintf(fdist, "\tnfree");
      fprintf(fdist, "\n");
   }

   //------Variable initialization------//

   for (i=0; i<NMOLS; i++) {
      nloop[i]	=	0;	nbridge[i]	=	0;
      ntail[i]	=	0;	nxseg[i]	=	0;
      nnucl[i]	=	0;	nxbead[i]	=	0;
      nfree[i]	=	0;
   }
   nloop_tot	=	0;	nbridge_tot	=	0;
   ntail_tot	=	0;	nxseg_tot	=	0;
   nfree_tot	=	0;

   //------Perform analysis------//

   for (moli=mol; moli<mol+NMOLS; moli++) {

      // Calculate how many xtal beads on each chain
      k	=	moli-mol;
      for (i=0; i<moli->nsites; i++) {
         nxbead[k]	+=	(moli->nuclid[i]>0 ? 1 : 0);
      }

      // Chain with max xbeads
      if (nxbead[k] > max) {
         max 	= 	nxbead[k];
         maxid	=	k;
      }

      // Chain with max/min Rg
      if (Rg2[k] > maxRg) {
         maxRg		=	Rg2[k];
         maxRgid	=	k;
      }
      if (Rg2[k] < minRg) {
         minRg		=	Rg2[k];
         minRgid	=	k;
      }
   }

   //------Find out the # of nuclei each chain participates------//
   
   for (k=0; k<NMOLS; k++) {
      for (i=0; i<nsegment[k]; i++) {

         if (seg_stat[k][3*i+1] > 0) {
            same	=	0;
            for (j=0; j<i; j++) {
               if (seg_stat[k][3*i+1]==seg_stat[k][3*j+1]) {
                  same	=	1;
                  break;
               }
            }
            if (!same) {       nnucl[k]	++; }
         }
      }
   }

   // Find out how many tails, etc, each chain has
   
   for (k=0; k<NMOLS; k++) {	
      for (j=0; j<nsegment[k]; j++) {

         length	=	seg_stat[k][j*3+2] - seg_stat[k][j*3] + 1;
	 
         if (length<3 || (length<4 && j>0 && j<nsegment[k]-1)) {
            printf("error, smooth does not work\n");
            printf("chain %ld seg %ld head %d tail %d\n", k, j, seg_stat[k][j*3], seg_stat[k][j*3+2]);
         }

         if (seg_stat[k][j*3+1]==-1) {		// amorphous segment
	    head	=	(j==0 ? -1 : seg_stat[k][(j-1)*3+1]);	
            tail	=	(j==nsegment[k]-1 ? -1 : seg_stat[k][(j+1)*3+1]);

            if (head*tail < 0) {		// tail
               ntail[k]		++;
            }
            else if (tail==-1) {		// whole chain in melt
               nfree[k]		++;
            }
            else if (head==tail) {		// fold
               nloop[k]	++;
            }
            else if (head!=tail) {		// bridge
               nbridge[k]	++;
            }
         }
         else if (seg_stat[k][j*3+1]>0) {	// nseg xtal segment
            nxseg[k]	++;
         }
      }
      ntail_tot		+=	ntail[k];
      nloop_tot		+=	nloop[k];	
      nbridge_tot	+=	nbridge[k];
      nxseg_tot		+=	nxseg[k];
      nfree_tot		+=	nfree[k];
   }

   //------Write system variables to Rg.hst------//

   fprintf(fPtr, "%-6ld", timestep);
   fprintf(fPtr, " %4.3f", P2[0]);
   fprintf(fPtr, " %ld", Xtal[0]);
   fprintf(fPtr, " %ld", nmax[0][0]);
   fprintf(fPtr, " %ld", nmax[0][1]);
   fprintf(fPtr, " %6.3f", R02);
   fprintf(fPtr, " %6.3f", aveRg2);
   fprintf(fPtr, " %6.3f", avexRg2);
   fprintf(fPtr, " %6.3f", avemRg2);

   fprintf(fPtr, "\t<10>");				// formatting lines

   fprintf(fPtr, " %ld", (long)(smoothresult.y));	// diff. in nucleus size 
   fprintf(fPtr, " %ld", (long)(smoothresult.z));	// diff. in total xtal
   fprintf(fPtr, " %ld", ntail_tot);
   fprintf(fPtr, " %ld", nloop_tot);
   fprintf(fPtr, " %ld", nbridge_tot);
   fprintf(fPtr, " %ld", nxseg_tot);
   fprintf(fPtr, " %ld", nfree_tot);

   fprintf(fPtr, "\t<18>");			// formatting lines

   fprintf(fPtr, " %ld", maxid);		// chain with max xbead
   fprintf(fPtr, " %ld", nxbead[maxid]);
   fprintf(fPtr, " %6.3f", Rg2[maxid]);

   fprintf(fPtr, " %ld", minRgid);		// chain with min Rg2
   fprintf(fPtr, " %ld", nxbead[minRgid]);
   fprintf(fPtr, " %6.3f", Rg2[minRgid]);

   fprintf(fPtr, " %ld", maxRgid);		// chain with max Rg2
   fprintf(fPtr, " %ld", nxbead[maxRgid]);
   fprintf(fPtr, " %6.3f", Rg2[maxRgid]);

   fprintf(fPtr, "\t<28>");			// formatting lines

   fprintf(fPtr, " %4.3f", asphericity);	// asphericity of the simulation box
   fprintf(fPtr, " %6.3f", dr_sys.x);
   fprintf(fPtr, " %6.3f", dr_sys.y);
   fprintf(fPtr, " %6.3f", dr_sys.z);
   //fprintf(fPtr, " %6.3f", d2initave);		// sq'ed displacement wrt initial conf.
   //fprintf(fPtr, " %6.3f", d2lastave);		// sq'ed displacement wrt last dump conf.
   fprintf(fPtr, " %6.3f", d2sys);
   //fprintf(fPtr, " %6.3f", d2syslast);

   fprintf(fPtr, " %6.3f %6.3f %6.3f", v2*16743.6, vf2*16743.6, vcom2*NSITES*16743.6);
fprintf(fPtr, " %6.3f %6.3f", lam1, amo1);
   fprintf(fPtr, "\n");

/*
   fprintf(fPtr, "TIMESTEP %ld\n", timestep);
   fprintf(fPtr, "%ld\t%ld\t%ld\n", NSYSTEMS, NMOLS, NSITES);
   fprintf(fPtr, "%f %f ", aveRg2, avexRg2);
   fprintf(fPtr, "%ld %ld %ld %ld\n", nloop_tot, nbridge_tot, ntail_tot, nxseg_tot);
*/
//average Rg2, total # of folds, total # of bridges, total # of tails, total # of xseg


   //------Write variables for each chain to Rg.dist------//

   for (i=0; i<NMOLS; i++) {
      fprintf(fdist, "%-6ld", timestep);
      fprintf(fdist, "\t%ld", nmax[0][0]);
      fprintf(fdist, "\t%ld", i);
      fprintf(fdist, "\t%f", Rg2[i]);
      fprintf(fdist, "\t%ld", nxbead[i]);	// # of xtal beads on this chain
      fprintf(fdist, "\t%ld", nnucl[i]);	// # of nuclei attached to this chain
      fprintf(fdist, "\t%d", nsegment[i]);
      fprintf(fdist, "\t%ld", ntail[i]);	// # of tails on this chain
      fprintf(fdist, "\t%ld", nloop[i]);
      fprintf(fdist, "\t%ld", nbridge[i]);
      fprintf(fdist, "\t%ld", nxseg[i]);
      fprintf(fdist, "\t%ld", nfree[i]);
      fprintf(fdist, "\n");
   }

   return;
}

//======================================================================//
//	output_pre(): write pre-analysis results for further analysis	//
//	              These pre-results should be basic but very	//
//		      time-consuming.  So by saving these results in	//
//		      a pre- file we can reuse them later without 	//
//		      calculating them every time.			//
//======================================================================//
void output_pre(FILE *fPtr) 
{
   molstruct	*moli;
   long		i, system=0;

   //------New version------//
   
   // write system variables
   Print_hst(fPtr);

   // write atomic variables
   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%ld %ld %ld", moli-mol, moli->box, moli->nsites);
      for(i=0; i<moli->nsites; i++)
         fprintf(fPtr, " %f", moli->p2[i]);
      fprintf(fPtr, "\n");
   }
   
   //------Old version------//
   /* 
   fprintf(fPtr, "TIMESTEP %ld\n", timestep);
   fprintf(fPtr, "%ld\t%ld\t%ld\n", NSYSTEMS, NMOLS, NSITES);

   // Output system variables
   fprintf(fPtr, "%8.4f ", v[0].tot);		// potential energy
   fprintf(fPtr, "%8.4f ", BOX[0].pres);	// pressure
   fprintf(fPtr, "%6.4f ", P2[0]);		// global orientation o.p.
   fprintf(fPtr, "%6.4f ", transfrac[0]);	// trans fraction
   fprintf(fPtr, "\n");

   // Output nuclei info. for the system
   fprintf(fPtr, "%ld", Xtal[0]);
   fprintf(fPtr, " %ld", realXtal[0]);
   fprintf(fPtr, " %ld", Nnucl[0]);
   for (i=0; i<10; i++)
      fprintf(fPtr, " %ld", nmax[0][i]);
   for (i=1; i<=Nnucl[0]; i++)
      fprintf(fPtr, " %ld", sizeofnucl[i]);
   for (i=1; i<=nmax[0][0]; i++) 
      fprintf(fPtr, " %ld", sizedist[i]);
   fprintf(fPtr, "\n");

   // Output atomic variables
   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%ld %ld %ld", moli-mol, moli->box, moli->nsites);
      for(i=0; i<moli->nsites; i++)
         fprintf(fPtr, " %f %ld", moli->p2[i], moli->nuclid[i]);
      fprintf(fPtr, "\n");
   }
   */

   fflush(fPtr);
   return;
}

//==============================================================//
//	read_pre(): read pre-analysis results, see output_pre()	//
//==============================================================//
void read_pre(FILE *fPtr, long timestep) 
{
   char		dummy[255];
   long		id, nsites, system, i;
   long		time	=	-999;
   float	temp;
   molstruct	*moli;

   // new version
   
   do {
      fscanf(fPtr, "%ld", &time);
      
      fscanf(fPtr, "%lf %f %lf", &(v[0].tot), &temp, &(BOX[0].pres));	// Epot, Vol, Pres
      fscanf(fPtr, "%f %f %f", &temp, &temp, &temp);			// Lx, Ly, Lz
      fscanf(fPtr, "%f %f %f", &temp, &temp, &temp);			// xy, xz, yz

      fscanf(fPtr, "%lf %lf %lf", P2+0, P2M+0, P2z+0);
      fscanf(fPtr, "%lf", transfrac+0);
      fscanf(fPtr, "%ld", Xtal+0);
      fscanf(fPtr, "%ld", realXtal+0);
      fscanf(fPtr, "%ld", Nnucl+0);
      fscanf(fPtr, "%ld %ld", &(nmax[0][0]), &(nmax[0][1]));
   
      fscanf(fPtr, "%lf", Q6+0);
      fscanf(fPtr, "%f", &rshift2);

      fscanf(fPtr, "%lf", &(v[0].stretch));
      fscanf(fPtr, "%lf", &(v[0].bending));
      fscanf(fPtr, "%lf", &(v[0].torsion));
      fscanf(fPtr, "%lf", &(v[0].lj));
      fscanf(fPtr, "%lf", &(v[0].ljcorr));
      
      fgets(dummy, sizeof(dummy), fPtr);
      
      v[0].corr		=	v[0].ljcorr;
      v[0].bonded	=	v[0].stretch + v[0].bending + v[0].torsion;
      v[0].nonbonded	=	v[0].lj + v[0].corr;
            
      for (moli=mol; moli<mol+NMOLS; moli++) {
         fscanf(fPtr, "%ld%ld%ld", &id, &system, &nsites);
         for (i=0; i<nsites; i++) {
            fscanf(fPtr, "%f", moli->p2+i);
         }
      }
   } while (time!=timestep);

   // old version
   /*
   do {
      fscanf(fPtr, "%s%ld", dummy, &time);
      fscanf(fPtr, "%ld%ld%ld", &NSYSTEMS, &NMOLS, &NSITES);
      fscanf(fPtr, "%lf%lf", &(v[0].tot), &(BOX[0].pres));
      fscanf(fPtr, "%lf%lf", &(P2[0]), &(transfrac[0]));

      fscanf(fPtr, "%ld%ld%ld", &(Xtal[0]), &(realXtal[0]), &(Nnucl[0]));
      for (i=0; i<10; i++)
         fscanf(fPtr, "%ld", nmax[0]+i);
      MAXSIZE[0]	=	nmax[0][0];

      for (i=1; i<=Nnucl[0]; i++)
         fscanf(fPtr, "%ld", sizeofnucl+i);
      for (i=1; i<=nmax[0][0]; i++) 
         fscanf(fPtr, "%ld", sizedist+i);

      for (moli=mol; moli<mol+NMOLS; moli++) {
         fscanf(fPtr, "%ld%ld%ld", &id, &system, &nsites);
         for (i=0; i<nsites; i++) {
            fscanf(fPtr, "%f%ld", moli->p2+i, moli->nuclid+i);
         }
      }
   } while (time!=timestep);
   */
   return;
}

//==============================================================//
//	output_xyzfile(): write configuration to .xyz file	//
//==============================================================//
void output_xyzfile(FILE *fPtr)
{
   molstruct	*moli;
   int		i;
   vector	r;

   fprintf(fPtr, "%ld\n", NSITES);
   fprintf(fPtr, "comments......\n");
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 // fprintf(fPtr, "O %f %f %f\n", moli->p[i].x, moli->p[i].y, moli->p[i].z);
	 r	=	MapInBox2(moli->p+i, PBC, moli->box);
	 fprintf(fPtr, "O %f %f %f\n", r.x, r.y, r.z);
      }
   }
   fflush(fPtr);
   return;
}

//==============================================================//
//	output_carfile(): write configuration to .car file	//
//==============================================================//
void output_carfile(FILE *fPtr)
{
   time_t	t	=	time(NULL);
   molstruct	*moli;
   char		s[80], ff[80];
   int		system = 0;			// only one box for now (11/18/2011)
   int		i, n;

   fprintf(fPtr, "!BIOSYM archive 3\n");
   fprintf(fPtr, "PBC=ON\n");
   fprintf(fPtr, "!TIMESTEP %ld\n", timestep);
   fprintf(fPtr, "!DATE %s", asctime(localtime(&t)));
   fprintf(fPtr, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

   n	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {

         MolInBox2(moli);
         for (i=0; i<moli->nsites; i++) {
            //moli->p[i] = MapInBox2(moli->p+i, PBC, system); //temporary

            if (moli->nuclid[i]>0)			// crystal like particle
               sprintf(s, "N%d", n++);			// N: blue color in VMD
            else
               sprintf(s, "O%d", n++);			// O: red color in VMD

	    /*
            if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system])		// note nuclid index starts from 1
               sprintf(s, "M%ld", n++);
            else if (sizeofnucl[moli->nuclid[i]] -MAXSIZE[system] > -3 && MAXSIZE[system]>=10)
               sprintf(s, "C%ld", n++);
            else if (moli->nuclid[i] >= 1)
               sprintf(s, "O%ld", n++);
            else
               sprintf(s, "H%ld", n++);
	    */

            fprintf(fPtr, "%-5.5s ", s);
            sprintf(s, "M%ld", moli-mol);
            fprintf(fPtr, "%14.8g %14.8g %14.8g ", moli->p[i].x, moli->p[i].y, moli->p[i].z);
            strcpy(ff, "O");
            fprintf(fPtr, "%-4.4s %-6ld ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[i], s));
         } 
      }   
   }
   fprintf(fPtr, "end\nend\n");
   fflush(fPtr);
   return;
}

//======================================================================//
//	output_pdbfile(): write configuration to pdb file and vmd files	//
//======================================================================//
void output_pdbfile(FILE *fPtr, FILE *fdat, long imgn, double imgx, double imgy, double imgz) 
{
   time_t	t	=	time(NULL);
   molstruct	*moli;
   char		atomname;
   int		system = 0;			// only one box for now (11/18/2011)
   int		i, m, n;
   int		drawmol;
   int		nuclid;

   fprintf(fPtr, "HEADER: file created on %s", asctime(localtime(&t)));
   fprintf(fPtr, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
        	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

   m		=	0;			// molecule sequence number
   n		=	0;			// atom sequence number
   for (moli=mol; moli<mol+NMOLS; moli++) {
      if (system==moli->box) {
         //MolInBox2(moli);

         drawmol	=	0;
         for (i=0; i<moli->nsites; i++) {
	    nuclid	=	moli->nuclid[i];
            if (sizeofnucl[nuclid] == nmax[system][0]) {	// part of the largest nucleus
               drawmol	=	1;
               break;
            }
	 }

	 //rtemp=CenterofMass(moli);
	 //if (rtemp.z < -0.25 * BOX[system].lz) {

         m	++; 
         for (i=0; i<moli->nsites; i++) {
            if (drawmol) {
	       nuclid	=	moli->nuclid[i];
               if (sizeofnucl[nuclid] == nmax[system][0]) {	// nuclid index starts from 1
                  atomname	=	'N';			// N: blue color in VMD
                  fprintf(fdat, " 10");
               }
               else {
                  atomname	=	'O';			// O: red color in VMD
                  fprintf(fdat, " 0");
               }
            }
            else {
	       atomname	=	'C';				// C: cyan color in VMD
	       fprintf(fdat," -1");
	    }

            n	++;
	    fprintf(fPtr, "ATOM  ");			// pdb command, column 1-6
            fprintf(fPtr, "%5d ", n);			// atom number
            fprintf(fPtr, " %c  ", atomname);		// atom name
            fprintf(fPtr, " ");				// alternate location indiator
  	    fprintf(fPtr, " C8");			// residue name
	    fprintf(fPtr, " ");				// column 21
            fprintf(fPtr, " ");				// chain identifier, column 22
	    fprintf(fPtr, "%4d", m);			// residue sequence number, 23-26
	    fprintf(fPtr, " ");				// code for insertion of residues, 27
            fprintf(fPtr, "   ");			// column 28-30
            fprintf(fPtr, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z);
            fprintf(fPtr, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fPtr, "%5.5s", "");
            fprintf(fPtr, "\n"); 

            if (imgn) {			// for image box
               n	++;
               fprintf(fPtr, "ATOM  ");
               fprintf(fPtr, "%5d ", n);
               fprintf(fPtr, " %c  ", atomname);
               fprintf(fPtr, " ");
               fprintf(fPtr, " C8");
               fprintf(fPtr, " ");
               fprintf(fPtr, " ");
               fprintf(fPtr, "%4d", m);
               fprintf(fPtr, " ");
               fprintf(fPtr, "   ");
               fprintf(fPtr, "%8.3f%8.3f%8.3f", 
			moli->p[i].x + imgx, moli->p[i].y+imgy, moli->p[i].z+imgz);
               fprintf(fPtr, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
               fprintf(fPtr, "%5.5s", "");
               fprintf(fPtr, "\n"); 
            }
         } 
	 //}
      }   
   }
   fprintf(fPtr, "END\n");
   fflush(fPtr);
   fflush(fdat);
   return;
}

//======================================================================//
//	output_conf(): write configuration to configuration file	//
//======================================================================//
void output_conf(FILE *fPtr, long timestep) 
{
   molstruct	*moli;
   int		i;

   fprintf(fPtr, "!TIMESTEP %ld\n", timestep);
   fprintf(fPtr, "%ld\t%ld\t%ld\n", NSYSTEMS, NMOLS, NSITES);

   for (i=0; i<NSYSTEMS; i++) {
      fprintf(fPtr, "%lf\t%lf\t%lf\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {
      fprintf(fPtr, "%ld\t%ld\t%ld\n", moli-mol, moli->box, moli->nsites);
      //fprintf(fconf, "%ld\t%ld\t%ld\t%ld\t%ld\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
      //MolInBox(moli);
      for (i=0; i<moli->nsites; i++) 
         fprintf(fPtr, "%ld\t%lf\t%lf\t%lf\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
   }
   fflush(fPtr);
   return;
}


//======================================================================//
//	output_Zinput(): write coord. as Z code input file (4/4/13)	//
//======================================================================//
void output_Zinput(FILE *fPtr)			// Output conf. as Z code input file
{
   //FILE		*fPtr;
   int		system = 0;
   int		i, j;

   //fPtr	=	fopen("Z_initconfiguration", "w");

   fprintf(fPtr, "%d  %ld  %d\n", 3, NMOLS, 0);
   fprintf(fPtr, "%lf  %lf  %lf\n", BOX[system].lx, BOX[system].ly, BOX[system].lz);
   fprintf(fPtr, "%f   %f\n", 1.0, 0.0);
   for (i=0; i<NMOLS; i++) {
      fprintf(fPtr, "%ld  ", (mol+i)->nsites);
   }
   fprintf(fPtr, "\n");

   for (i=0; i<NMOLS; i++) {
      for (j=0; j<(mol+i)->nsites; j++) {
	 fprintf(fPtr, "%lf  %lf  %lf\n", (mol+i)->p[j].x, (mol+i)->p[j].y, (mol+i)->p[j].z);
      }
   }

   //fclose(fPtr);
   return;
}
