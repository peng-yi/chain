/*
 * 
	program:	metal2hst.c
	author:		Peng Yi at JHU
	date:		Dec. 28, 2013
	purpose:	analyze LAMMPS dump file for metal simulations
	note:		require setup file
	update:
			Dec. 28, 2013	Created based on lammps2hst.c
			Mar. 04, 2014	Enable calling LAMMPS as functions
			Jan. 02, 2017	Before execute lammps module, should execute run 0 to keep
			                all computes current
	output:
*/
#define __MAIN_PROGRAM
#include "header.h"

#define M2HVERSION	"03Feb2018"
#define M2HDEBUG	0
#define MAXNINPUTFILE	256

#undef CODE

//#include "correlation.h"		// autocorrelation calculation module, can consume a lot of memory

//==============================//
//	Global variables	//
//==============================//
long		timestep;
double		pe, ke, etotal, enthalpy, vol, rho_m, atoms;	// thermo variables in LAMMPS
double		temp, pressure, presski, presske, stress[6];
double		rho_n;					// atom number density

float		R02;
float		Rg2[MAXNMOLS];		// Radius of gyration square
float		rgx2[MAXNMOLS], rgy2[MAXNMOLS], rgz2[MAXNMOLS];
float		avergx2, avergy2, avergz2;

vector		dr_sys;			// vector of system drift
float		d2initave,		// average displacement w.r.t. to the beginning
		d2lastave,		// average displacement w.r.t. to previous dump
                d2sys, d2syslast;
float		asphericity;		// asphericity of the box

//==============================================//
//	Function Declarations and Definitions	//
//==============================================//
void 		output_xyzfile(FILE *, FILE *);
void		output_carfile(FILE *);
void		output_pdbfile(FILE *);
void		output_dump(FILE *);
int 		displayatom(int n);
void		shift_coord();
void 		monitor_coord();
void 	binning(char *,char *);

//======================//
//	Main Program	//
//======================//
int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   time_t	prog_start, prog_end, start, end;
   molstruct	*moli, *molj;
   int 		i, j, k, m, n, id, maxid, size;
   int 		system;
   int		siteid, molid, type;
   int		nx, ny, nz, nx1, ny1, nz1;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo,			// box upper boundary
		xy=0.0, xz=0.0, yz=0.0;		// triclinic parameter
   long		LENGTH, accum;			// chain length variables

   // dummy variables
   
   int		itmp, itmp1, itmp2, itmp3, itmp4, itmp5, itmp6;
   int		flag;
   double	tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
   double	max, min, max1, min1, max2, min2;
   int		*itmptr1, *itmptr2;
   void		*tmptr1, *tmptr2;
   double	*tmptr;
   vector	vtmp, vtmp1, vtmp2;
   double	*data1D;
   double	**data2D;
   int		*ldata1D;
   int		*idata1D;
   int		*lmp_atomid;
   char 	line[1024];

   // file variables
   
   char		infile[MAXNINPUTFILE][80], filein[80], filename[80];
   char  	procid[16];
   char		task[80];
   char		tmpstr[80], fmtstr[80];
   char		s[80], ff[80], par[80], dummy[255];
   char		atomname;			// for visualization file output

   FILE		*fin, *fhst, *fpre;		// input and analysis output
   FILE		*fxyz;
   FILE		*fcolor;
   FILE		*fcar, *fconf, *fpdb, *fdat;	// configuration/visualization output
   FILE		*fZfile;			// Z-code input configuration file
   FILE		*fdump;
   FILE		*fsf; 				// stacking fault
   int		eof;  				// end-of-file flag
   int		file_exist;  			// file existence flag

   int		preflag		= 0;		// flags for input command line
   int		polydisperse	= 0;
   int		confflag	= 0;
   int		carflag		= 0;
   int		pdbflag		= 0;
   int		xyzflag		= 0;
   int		Zflag		= 0;
   int  	typechangeflag 	= 1; 		// if type changes during LAMMPS simulation
   int		lammpsflag	= 0;		// whether or not to invode LAMMPS
   int		nfiles=1;			// number of input files 
   int		ifile;				// index of input file
   int		nframe, dnframe = 1;		// analyze only every dnframe
   int		previous_timestep = -1;
   int		starttime = -1;
   int		endtime = -1;

   static int	firstconf = 1;			// first configuration read

   vector		con;			// center of nucleus
   static vector	rO;			// original position of nucleus

   // configuration init flag
   static int	confinit = 1;

   // cell list flag
   static int	clistinit=1;

   // hyper-distance
   static float r_MD[MAXNMOLS*MAXNMOLSITES];	// coordinates in a dynamic trajectory
   static float r_MIN[MAXNMOLS*MAXNMOLSITES];	// coordinates in a minimized trajectory
   static float	r_trj[MAXNMOLS*MAXNMOLSITES];	// coordinates in a trajectory
   static float	r_ref[MAXNMOLS*MAXNMOLSITES];	// coordinates in a reference configuration
   //static float r_trjold[MAXNMOLS*MAXNMOLSITES];	// coordinates in the previous step
   //static float r_refold[MAXNMOLS*MAXNMOLSITES];	// coordinates in the previous step

   // a group of beads

   vector	*rbead;
   rbead	=	(vector *)calloc(MAXNMOLS*MAXNMOLSITES, sizeof(vector));

   // variables to characterize nucleus, chains, etc at ONE timestep

   long		nsites;
   int 		nuclid, nuclid1, nuclid2;
   int		nmaxid;					// nuclid of the largest nucleus
   beadstruct	*nucleus;				// group beads in the same nucleus
   nucleus = (beadstruct*) calloc(MAXNMOLS*MAXNMOLSITES, sizeof(beadstruct));

   vector	com[MAXNMOLS];				// average center of mass
   long		ncom[MAXNMOLS];
   int		nchainmax;				// # of chains in the largest nucleus

   // density profile

   long		direct;
   double	density[100];
   vector	rsite[MAXNMOLS * MAXNMOLSITES];	
   //
   // composition gradient
   double	c_grad;
   double 	local_grad;
   double 	local_comp;

   // composition in the largest nucleus
   double	c_nmax;
   // radius of the largest nucleus
   double	R_nmax;
   // com of the largest nucleus (geometrical center)
   vector	com_nmax;

   // velocity statistics
 
   double	dvel=0.001;
   double	ave_vx, ave_vy, ave_vz;
   double 	velmax=0.0, velmin=0.0;
   double	vel_distx[101], vel_disty[101], vel_distz[101];
 
   // variables to AVERAGE over a number of timesteps

   int  	*Nn;				// size dist. of nucleus size
   long		ntot, maxnmax;			// total number of nuclei for normalization
   float	*Rg2nucl;			// Rg2 of nucleus as a function of size

   // dislocation location
   double	xtop, xbot;
   double	xdisl;
   double	prev_xdisl = -1e6;

   // pbc image variables
   
   long		imagen=0;
   double	imagex, imagey, imagez;	

   // LAMMPS variables

   FILE 	*flmpin;
   void		*lmp;				// pointer to LAMMPS instance
   int		lmp_comm;			// 1 if in communicator for lammps, MPI_UNDEFINED otherwise
   double	*lammps_x;
   int		nlocal, natoms;			// meaning same as in LAMMPS
   int		*atomid;			// array containing atom id, ordered by local id, then by lammps proc

   if (M2HDEBUG) {
      printf("# DEBUG: Variables declared.\n");
      fflush(stdout);
   }

   //------MPI initialization------//

#ifdef myMPI
   MPI_Init(&argc, &argv);
   MPI_Comm_size(MPI_COMM_WORLD, &MPI_numprocs);	// return number of procs
   MPI_Comm_rank(MPI_COMM_WORLD, &MPI_myid);		// return id of this proc
#else
   MPI_numprocs	=	1;				// default value
   MPI_myid	=	0;				// default value
   nprocs_lammps = 	0;
   MPI_lammpsid	=	-1;
#endif

   if (argc<2) {
      printf("metal2hst (c) 2014 by Peng Yi at JHU\n\n");
      printf("Syntax:\n");
      printf("\tmetal2hst [-option] [x= 0.1 y= 0.1 z= 0.1 n= 1] [dn= 2 nfiles= 2] lammpsdumpfile(s)\n\n");
      printf("Notes:\n");
      printf("\t* -option = \n");
      printf("\t* -pre: generate pre-analysis results\n");
      printf("\t* -poly: polydisperse system, has molid in dump file\n");
      printf("\t* -conf: configuration file output\n");
      printf("\t* -car: car file output\n");
      printf("\t* -pdb: pdb file output\n");
      printf("\t* -xyz: xyz file output\n");
      printf("\t* -lammps: invoke LAMMPS\n");
      printf("\t* -typechange: atom type change during LAMMPS\n");
      printf("\t* task=: task (dislocation, nucleation, etc)\n");
      printf("\t* x= y= z=: duplicate the system and shift unit vector\n");
      printf("\t* n=: multiple of shift vector\n");
      printf("\t* start=: the timestep that analysis stops\n");
      printf("\t* end=: the timestep that analysis stops\n");
      printf("\t* dn=: only analyze every dn frames\n");
      //printf("\t* nfiles=: number of input files if more than 1 (must be <=32)\n");
      printf("\t* files=: must be last option, followed by names of input files (no more than 32 input files)\n");
      printf("\t* \"=\" must immediately follow x or y or z or n or dn\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   char prefix[16];
   char firstfile[16];
   char lastfile[16];

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if(samestr(par, "-pre"))		preflag		=	1;
      else if (samestr(par, "-conf"))	confflag	=	1;
      else if (samestr(par, "-car"))	carflag		=	1;
      else if (samestr(par, "-pdb"))	pdbflag		=	1;
      else if (samestr(par, "-xyz"))	xyzflag		=	1;
      else if (samestr(par, "-lammps"))	lammpsflag	=	1;
      else if (samestr(par, "-typechange"))	typechangeflag	=	1;
      else if (samestr(par, "-Z"))	Zflag		=	1;
      else if (samestr(par, "x="))	imagex		=	atof(argv[i+1]);
      else if (samestr(par, "y=")) 	imagey		=	atof(argv[i+1]);
      else if (samestr(par, "z=")) 	imagez		=	atof(argv[i+1]);
      else if (samestr(par, "n=")) 	imagen		=	atol(argv[i+1]);
      else if (samestr(par, "start=")) 	starttime	=	atol(argv[i+1]);
      else if (samestr(par, "end=")) 	endtime		=	atol(argv[i+1]);
      else if (samestr(par, "dn=")) 	dnframe		=	atol(argv[i+1]);
      else if (samestr(par, "task="))	strcpy(task, argv[i+1]);
      //else if (samestr(par, "nfiles="))	nfiles		=	atoi(argv[i+1]);
      //else if (samestr(par, "files="))  nfiles		=	argc-1-i;

      else if (samestr(par, "prefix=")) {
	 strcpy(prefix, argv[i+1]); 
	 strcpy(firstfile, argv[i+2]);
	 strcpy(lastfile, argv[i+3]);
	 nfiles = atoi(lastfile) - atoi(firstfile) + 1;
      }
   }
   if (nfiles > MAXNINPUTFILE) {            				// check number of input files
      printf("Error! number of input files exceeds MAXNINPUTFILE\n");
      fflush(stdout);
      exit(0);
   }
   for (i=0; i<nfiles; i++) { 						// get input filenames
      //strcpy(infile[nfiles-1-i], argv[argc-1-i]);   
      
      strcpy(tmpstr, prefix);

      if (strcmp(firstfile, "-1")) {
	 sprintf(fmtstr, "%%0%dd", strlen(firstfile));
	 sprintf(filename, fmtstr, atoi(firstfile)+i);
	 strcat(tmpstr, filename);
      }
      strcpy(infile[i], tmpstr);
   }

   if (!strcmp(task, "dislocation") || 
       !strcmp(task, "nucleation") || 
       !strcmp(task,"profile")) { 		// check if task is valid
      // do nothing
   }
   else {
      printf("Error: task not found!\n");
      fflush(stdout);
      exit(1);
   }

   if (M2HDEBUG) {					// DEBUG output
      printf("# DEBUG: command line read. Proc id %d.\n", MPI_myid);
      fflush(stdout);
   }


   // allocate processors for lammps module
#ifdef myMPI
   if (lammpsflag) {
      if (MPI_numprocs < 2) {
	 printf("Error: MPI_numproc<2, cannot create new comm.\n");
	 fflush(stdout);
	 exit(0);
      }

      if (MPI_myid > 0) {			// all but the first proc are used by lammps
	 lmp_comm = 1;
      }
      else {
	 lmp_comm = MPI_UNDEFINED;
      }

      comm_world = MPI_COMM_WORLD;
      MPI_Comm_split(comm_world, lmp_comm, 0, &comm_lammps);	// split procs
      if (lmp_comm==1) {
	 MPI_Comm_rank(comm_lammps, &MPI_lammpsid);	// get the rank in lammps communicator
	 MPI_Comm_size(comm_lammps, &nprocs_lammps);	// get the size of lammps communicator
      }
      else {
	 MPI_lammpsid= -1;
      }
   }
   MPI_Barrier(comm_world);					// sync
   MPI_Bcast(&nprocs_lammps, 1, MPI_INT, 1, comm_world);	// pass nprocs_lammps to all procs
#endif

   //---------------------------//
   //	Open output files	//
   //---------------------------//
 
   if (nfiles==1)	strcpy(tmpstr, infile[0]);	// determine output filename
   else		strcpy(tmpstr, "multi");

   sprintf(procid, ".%d", MPI_myid);			// MPI: proc id for output files

#ifdef	myMPI
   if (MPI_myid==0) {					// ONLY process 0 do the following
#endif

      //------ Open optional output files ------//

      if (xyzflag) {
	 strcpy(filename, tmpstr);
	 strcat(filename, ".xyz");
	 strcat(filename, procid);
	 fxyz 	= 	fopen(filename, "w");
	 fcolor	=	fopen("color", "w");
      }
      if (carflag) {
	 strcpy(filename, tmpstr);
	 strcat(filename, ".car");
	 strcat(filename, procid);
	 fcar 	= 	fopen(filename, "w");
      }
      if (pdbflag) {
	 strcpy(filename, tmpstr);
	 strcat(filename, ".pdb");
	 strcat(filename, procid);
	 fpdb 	= 	fopen(filename, "w");
      }

#ifdef	myMPI
   }
#endif

   // every processor has output, so need to take care of all of them
   if (!M2HDEBUG) {			
      //strcpy(filename, tmpstr);
      //strcat(filename, ".out");	

      //strcpy(filename, "output");
      //strcat(filename, procid);		// every proc has own output file
      //freopen(filename, "w", stdout);		// redirect standard output stream to a file
      if (MPI_myid==0) {
	 freopen("output.0", "a", stdout);
	 if (!strcmp(task, "dislocation")) { 
	    fsf = fopen("stackingfault", "w");
	    fprintf(fsf, "timestep\t x_leading\t x_trailing\t width(Angstrom)\n");
	    fflush(fsf);
	 }
      }
      else {
	 freopen("output.1", "a", stdout);
      }
   }

   if (MPI_myid==0) {   
      if (M2HDEBUG) {					// DEBUG output
	 printf("# DEBUG: i/o files ready.\n");
	 fflush(stdout);
      }

      //------ Write run info. to output file ------//

      printf("# BEGINNING OF OUTPUT FILE\n");
      printf("%s\n", asctime(localtime(&t)));

      printf("# metal2hst version: %s\n", M2HVERSION);
      printf("# command: ");
      for (i=0; i<argc; i++) {	
	 printf("%s ", argv[i]);	
      }
      printf("\n");

      printf("# list of input files:\n");
      for (i=0; i<nfiles; i+=5) {
	 for (j=0; j<5; j++) {
	    printf("\t%s", infile[i+j]);
	 }
	 printf("\n");
      }
      printf("\n");

      printf("# number of processor:  %d\n", MPI_numprocs);
      if (lammpsflag) {
	 printf("# number of lammps processor:  %d\n", nprocs_lammps);
      }
      printf("# id of this processor: %d\n", MPI_myid);

      prog_start	=	time(NULL);			// to calculate running time
   }	// only proc 0

   printf("# Proc %d in lammps? %s. lammps proc id: %d\n", MPI_myid,(lmp_comm==1)?"Yes":"No",MPI_lammpsid);

   printf("\n");
   fflush(stdout);

   //-------------------//
   //	Initialization	//
   //-------------------//

   if (MPI_myid==0) {
      printf("# Initialization ...\n");
   }
   /*
   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   */
   //InitUnits();			// initialize units for calculation
   //InitForcefield();			// initialize Lennard-Jones potential mixing rule
   CalcUnits(0); 			// parameter 0 set every unit to 1.0

   NSYSTEMS	=	1;		// previously read from setup file 5/28/15

   system	=	0;		// for now 2/14/08, only one system

   //InitSample();			// initialize sampling

   nframe	=	-1;

   //------ Initialize LAMMPS module ------//
   
   if (lammpsflag) {

      /* open LAMMPS input script */

      if (MPI_lammpsid == 0) {			// FIRST lammps proc
	 flmpin = fopen("in.hybrid","r");
	 if ( flmpin ==NULL) {
	    printf("Error: analysis LAMMPS script not found.\n");
	    fflush(stdout);
	    MPI_Abort(comm_world, 1);
	    exit(0);
	 }
	 if (M2HDEBUG) {
	    printf("# DEBUG: lammps input script opened. Proc id %d.\n", MPI_myid);
	    fflush(stdout);
	 }
      }
      
      /* set lammps command line options */

      if (lmp_comm==1) {			// for ALL lammps procs
	 int	lmpargc;			// command line variable for lammps
	 char	**lmpargv;
	  
	 lmpargc 	=	5;
	 if (lmpargc>1) {
	    if ( (lmpargv = (char **)calloc(lmpargc, sizeof(char *)))==NULL) {
	       printf("Error: lmpargv allocation failed\n");
	       exit(0);
	    }
	    for (i=0; i<lmpargc; i++) {
	       if ( (lmpargv[i] = (char *)calloc(80, sizeof(char)))==NULL) {
		  printf("Error: lmpargv[%d] allocation failed\n", i);
		  exit(0);
	       }
	    }
	    strcpy(lmpargv[1], "-screen");
	    strcpy(lmpargv[2], "none");
	    strcpy(lmpargv[3], "-log");
	    strcpy(lmpargv[4], "log.analysis");
	 }

	 // lmpargc=0;
	 lammps_open(lmpargc, lmpargv, comm_lammps, &lmp);	// open a LAMMPS instance
	 if (M2HDEBUG) {
	    printf("# DEBUG: lammps instance opened\n");
	 }
      }

      /* run input lammps script one line at a time */

      while (1) {
	 if (MPI_lammpsid == 0) {				// FIRST lammps proc
	    if (fgets(line, 1024, flmpin) == NULL) n=0;
	    else n=strlen(line)+1;
	    if (n==0) {fclose(flmpin);}
	 }
	 MPI_Bcast(&n, 1, MPI_INT, 1, MPI_COMM_WORLD);		// pass command to all processors
	 //MPI_Bcast(&n, 1, MPI_INT, 0, comm_lammps);		// pass command to all lammps processors
	 if (n==0) break;
	 MPI_Bcast(line, n, MPI_CHAR, 1, MPI_COMM_WORLD);
	 //MPI_Bcast(line, n, MPI_CHAR, 0, comm_lammps);
	 if (lmp_comm==1) lammps_command(lmp, line);
      }

      //lammps_x  = (double *) calloc(3*NSITES, sizeof(double));
      if (lmp_comm==1) {
	 NSITES	= lammps_get_natoms(lmp);			// so that we don't need setup file anymore
      }
   }

   // pass variable NSITES to all processors in comm_world
   MPI_Bcast(&NSITES, 1, MPI_INT, 1, comm_world);
   MPI_Barrier(comm_world);					// sync all processors

   if (MPI_myid==0 && M2HDEBUG){ 
      printf("# DEBUG: lammps script read and executed.\n");
   }
   
   //------ Initialize system variables ------//

   NMOLS	=	NSITES;
   nsites	=	NSITES;

   if (MPI_myid==0) {
      if (M2HDEBUG) {
	 printf("NSITES = %d\n", NSITES);
	 printf("NMOLS = %d\n", NMOLS);
	 fflush(stdout);
      }

      atomid	= (int *)calloc(NSITES, sizeof(int));
      ldata1D	= (int *)calloc(NSITES, sizeof(int));
      data1D	= (double *)calloc(6*NSITES, sizeof(double));	// reason for 6 is tensor 

      if ( (mol = (molstruct*)calloc(NMOLS, sizeof(molstruct))) == NULL) {
	 printf("# Error: mol allocation failed\n");
	 fflush(stdout);
	 exit(0);
      }
      
      if ( (oldmol = (molstruct*)calloc(NMOLS, sizeof(molstruct))) == NULL) {
	 printf("# Error: oldmol allocation failed\n");
	 fflush(stdout);
	 exit(0);
      }

      if (M2HDEBUG) {
	 printf("# DEBUG: mol and oldmol allocated.\n");
      }
   }
   MPI_Barrier(comm_world);

   //------ Initialize other variables ------//

   if (MPI_myid==0) {
      Nn	=	(int *) calloc(NSITES, sizeof(int));
      ntot	=	0;
      maxnmax	=	0;
      Rg2nucl	=	(float *) calloc(NSITES, sizeof(float));
   }

   //---------------------------//
   //	Start Data Processing	//
   //---------------------------//

   for (ifile=0; ifile<nfiles; ifile++) {		// multiple input dump files

      if (MPI_myid==0) {
	 fin = fopen(infile[ifile], "r");
	 file_exist = (fin != NULL);			// file exist or not, empty file counts as exist
	 eof = 0;   					// not the end of file

	 if (file_exist) {
	    if (fgets(dummy, sizeof(dummy), fin)==NULL) {	// read one line, for the while loop next
	       eof	=	1;				// determine if end of file, like empty file
	    }
	    if (M2HDEBUG) {
	       printf("\n# DEBUG: Current dump file: %s\n", infile[ifile]);
	       printf("\n# DEBUG: eof = %d\n", eof);
	    }
	 }
	 else {
	    if (M2HDEBUG) {
	       printf("\n# DEBUG: dump file: %s does not exist.\n", infile[ifile]);
	    }
	 }
      }
      MPI_Barrier(comm_world);					// Important here !!
      MPI_Bcast(&eof, 1, MPI_INT, 0, comm_world);		// pass end-of-file flag
      MPI_Bcast(&file_exist, 1, MPI_INT, 0, comm_world);	// pass file exist flag

      while (file_exist && !eof) {

	 //------ Read one configuration from dump file	------//

	 if (MPI_myid==0) {

	    printf("\n# Reading one configuration ...\n");
	    fflush(stdout);
	 
	    /*
	    if (!fgets(dummy, sizeof(dummy), fin)) {	// read one line, if NULL then break
	       break;
	    }
	    */
	    if (strstr(dummy, "TIMESTEP")) {		// if 1st line contains string "TIMESTEP"
	       fscanf(fin, "%ld", &timestep);		// read 2nd line for the value of timestep
	       fgets(dummy, sizeof(dummy), fin);
	    }

	    //------ If first conf. read in to determine system variables used later ------//
	    //------ otherwise let LAMMPS library functions handel the configuration ------//

	    //if (firstconf) {
	    //   firstconf = 0;				// do it only once

	    fgets(dummy, sizeof(dummy), fin);						// line 3
	    fscanf(fin, "%ld", &nsites);	fgets(dummy, sizeof(dummy), fin);	// line 4
	    if (nsites!=NSITES) {
	       printf("# Error: nsites!=NSITES\n");
	       fflush(stdout);
	       exit(0);
	    }

	    fgets(dummy, sizeof(dummy), fin);						// line 5
	    if (strstr(dummy, "xy"))	PBC = 3;	// PBC=3 for triclinic box
	    fscanf(fin, "%lf%lf", &xlo, &xhi);						// line 6	
	    if (PBC==3)	fscanf(fin, "%lf", &xy);
	    else 		xy = 0.0;
	    fgets(dummy, sizeof(dummy), fin);
	    fscanf(fin, "%lf%lf", &ylo, &yhi);						// line 7
	    if (PBC==3) 	fscanf(fin, "%lf", &xz); 	
	    else		xz = 0.0;
	    fgets(dummy, sizeof(dummy), fin);
	    fscanf(fin, "%lf%lf", &zlo, &zhi);						// line 8
	    if (PBC==3) fscanf(fin, "%lf", &yz);
	    else	yz = 0.0;
	    fgets(dummy, sizeof(dummy), fin);
	    fgets(dummy, sizeof(dummy), fin);						// line 9

	    // calculate box dimension, triclinic in general

	    // xhi, xlo, etc in dump file are actually xhi_bound, xlo_bound, etc.
	    // ... need to convert to cell vectors
	    xlo -= MIN(0.0, MIN(xy, MIN(xz, xy+xz)));
	    xhi -= MAX(0.0, MAX(xy, MAX(xz, xy+xz)));
	    ylo -= MIN(0.0, yz);
	    yhi -= MAX(0.0, yz);

	    system = 0;
	    BOX[system].xlo = xlo;
	    BOX[system].xhi = xhi;
	    BOX[system].ylo = ylo;
	    BOX[system].yhi = yhi;
	    BOX[system].zlo = zlo;
	    BOX[system].zhi = zhi;

	    BOX[system].lx = xhi-xlo;
	    BOX[system].ly = yhi-ylo;
	    BOX[system].lz = zhi-zlo;
	    BOX[system].xy = xy;
	    BOX[system].xz = xz; 
	    BOX[system].yz = yz;

	    if (fabs(xy) < ZERO && fabs(xz) < ZERO && fabs(yz) < ZERO) {
	       PBC	=	1;						// for cell list to operate
	    }
       
	    // Read properties for each atom

	    LENGTH		=	NSITES/NMOLS;		// monodisperse system for now (4/26/2008)

	    for (i=0; i<nsites; i++) {			// loop over all atoms
	       fscanf(fin, "%ld", &id);			// read atom id
	       fscanf(fin, "%ld", &type);		// read atom type

	       fscanf(fin, "%lf%lf%lf", &x, &y, &z);	// read atom coordinates
	       //fscanf(fin, "%d%d%d", &nx, &ny, &nz);
	       nx = 0;					// assuming unwrapped coordinates in dump file
	       ny = 0;
	       nz = 0;
	       
	       //fscanf(fin, "%lf%lf%lf %lf%lf%lf %d%d%d", &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
	       /*
	       fscanf(fin, "%lf%lf%lf", &vx, &vy, &vz);
	       fscanf(fin, "%lf", &tmp1);
	       fscanf(fin, "%lf", &tmp2);
	       fscanf(fin, "%lf", &tmp3);
	       fscanf(fin, "%lf", &tmp4);
	       */
	       fgets(dummy, sizeof(dummy), fin);	// read the rest of the line for this atom

	       id	--; 					// -1 because lammps index starts from 1
	       molid	=	(long) (id/LENGTH);
	       siteid	=	id % LENGTH;

	       //------ Store info. only on proc 0 ------//

	       // box id and chain length
	       mol[molid].box	=	system;		// for now, only one system
	       mol[molid].nsites	=	LENGTH;		// for now, Jan/10/2010

	       // atom type 
	       mol[molid].type[siteid]	=	type - 1;	// -1 because lammps index starts from 1

	       // general triclinic box
	       mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx) 
						     + ny*(BOX[system].xy) + nz*(BOX[system].xz);
	       mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly) + nz*(BOX[system].yz);
	       mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);

	       // for density profile calculation
	       rsite[i].x	=	x;
	       rsite[i].y	=	y;
	       rsite[i].z	=	z;
       
	       /*
	       // velocity
	       mol[molid].velx[siteid]	=	vx;
	       mol[molid].vely[siteid]	=	vy;
	       mol[molid].velz[siteid]	=	vz;

	       // other atom attributes
	       mol[molid].pe[siteid]	=	tmp1;
	       mol[molid].cna[siteid]	=	(int) tmp2;
	       mol[molid].nneigh[siteid]	=	(int) tmp3;
	       */
	    }	// finish reading one configuration

	    //------ collect molecules information only on proc 0 ------//

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
	       moli->origin	=	CenterofMass(moli);
	    }

	    if (fgets(dummy, sizeof(dummy), fin)==NULL) {	// read one line
	       eof	=	1;				// determine if end of file
	    }
	 } // proc 0 read one configuration

         MPI_Barrier(comm_world);				// sync all processors
	 MPI_Bcast(&eof, 1, MPI_INT, 0, comm_world);		// pass end-of-file flag
	 MPI_Bcast(&timestep, 1, MPI_INT, 0, comm_world);	// pass timestep

	 if (M2HDEBUG) {
	    printf("proc %d ifile = %d timestep = %d previous = %d eof = %d\n", 
	            MPI_myid, ifile, timestep, previous_timestep, eof);
	    printf("proc %d start_time  = %d end_time = %d\n", 
	            MPI_myid, starttime, endtime);
	    fflush(stdout);
	 }

         //-----------------------------//
	 //	Check timestep range	//
         //-----------------------------//

	 if (timestep < starttime) {			// starttime = -1 by default
	    continue;
	 }
	 if (endtime > 0 && endtime < timestep) {	// endtime = -1 by default
	    break;
	 }

	 //------ Skip repeated frames from different input files ------//

	 if (ifile>=1 && timestep == previous_timestep) {
	    continue;
	 }
	 previous_timestep	=	timestep;	// check repeat frames

	 //------ Skip (dnframe - 1) frames ------//

	 nframe	++;
	 if (nframe%dnframe)	continue;		// analyze every dnframe frames

	 //-----------------------------------------------------//
	 //	Call LAMMPS module				//
	 //     - Feed coordinates to LAMMPS module		//
	 //	- Do calculation using LAMMPS			//
	 // 	- Extract computation results			//
	 //-----------------------------------------------------//
lammpsmodule:

         if (lammpsflag) {

	    if (MPI_myid==0) {
	       printf("\n    # Calling LAMMPS module, timestep = %d ...\n", timestep);
	       fflush(stdout);
	    }
	   
	    /* Read one configuration from dump file */
	    if (lmp_comm==1) {
	       printf("read_dump begin\n"); fflush(stdout);

	       if (typechangeflag) {  	// if type changes during LAMMPS simulation
		  sprintf(s, "read_dump %s %ld x y z wrapped no replace no purge yes add yes", infile[ifile], timestep);
	          //sprintf(s, "read_dump %s %ld x y z vx vy vz wrapped no label type v_newtype", infile[ifile], timestep);
		  printf("read_dump %s %ld x y z wrapped no replace no purge yes add yes\n", infile[ifile], timestep);
	       }
	       else {
	          sprintf(s, "read_dump %s %ld x y z type wrapped no", infile[ifile], timestep);
		  printf("read_dump %s %ld x y z type wrapped no\n", infile[ifile], timestep);
	       }

	       lammps_command(lmp, s); 

	       //lammps_command(lmp, "variable type1 atom 'type==1'");
	       //lammps_command(lmp, "group grouptype1 variable type1");
	       //lammps_command(lmp, "print count(grouptype1)");
	       //lammps_command(lmp, "compute 1 all property/atom type");
	       //lammps_command(lmp, "set atom * type c_1");
	       //lammps_command(lmp, "run 0"); 
	       //
	       //sprintf(s, "set atom * type v_newtype");
	       //lammps_command(lmp, s); 

	       //lammps_command(lmp, "change_box all ortho");

	       fflush(stdout);
	    }
	    MPI_Barrier(comm_world);		// sync because ldata1D will be reused next

	    // pass the atom types to lammps
	    /*
	    natoms = NSITES;
	    printf("proc %d natoms = %d\n", MPI_myid, natoms);
	    int *typedata = (int *)malloc(natoms*sizeof(int));
	    if (MPI_myid==0) {
	       n = 0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     typedata[n] = moli->type[i]+1;
		     n ++;
		  }
	       }
	       for (i=9292;i<9298; i++)
		  printf("type[%d] = %d\n", i, typedata[i]);
	       printf("prepare for pass types, n = %d.\n", n);
	    }
	    MPI_Barrier(comm_world);
	    MPI_Bcast(typedata, natoms, MPI_INT, 0, comm_world);	// pass atom types
	    if (lmp_comm==1){
	       if (MPI_myid==1) {
		  for (i=9292;i<9298; i++)
		     printf("type[%d] = %d\n", i, typedata[i]);
	       }
	       lammps_scatter_atoms(lmp, "type", 0, 1, typedata);

	       lammps_command(lmp, "run 1");		// important !! to keep computes current
	       
	       lammps_gather_atoms(lmp, "type", 0, 1, typedata);

	       lammps_command(lmp, "run 1");		// important !! to keep computes current
	       if (MPI_myid==1) {
		  for (i=9292;i<9298; i++)
		     printf("type[%d] = %d\n", i, typedata[i]);
	       }
	    }
	    MPI_Barrier(comm_world);
	    free(typedata);

	    printf("proc %d Atom type passed.\n", MPI_myid);
	    fflush(stdout);
	    */

	  // Test lammps gather atom-based quatities
	  /*
	  if (lmp_comm==1) {
	     printf("test lammps gather\n");
	     natoms = lammps_get_natoms(lmp);
	     int *idarray = (int *) malloc(natoms*sizeof(int));
	     double *xarray = (double *) malloc(3*natoms*sizeof(double));
	     lammps_gather_atoms(lmp,"id",0,1,idarray);
	     lammps_gather_atoms(lmp,"x",1,3,xarray);
	     if (MPI_myid==1) {
		for (i=0; i<10; i++) {
		   printf("gather id: %d %f %f %f\n", idarray[i], xarray[3*i], xarray[3*i+1], xarray[3*i+2]);
		}
		printf("\n\n");
	     }

	     xarray[0] += 0.1;
	     lammps_scatter_atoms(lmp,"x",1,3,xarray);
	     lammps_command(lmp, "run 1");

	     lammps_gather_atoms(lmp,"x",1,3,xarray);
	     if (MPI_myid==1) {
		for (i=0; i<10; i++) {
		   printf("gather id: %d %f %f %f\n", idarray[i], xarray[3*i], xarray[3*i+1], xarray[3*i+2]);
		}
	     }
	  }
	  */
	    if (lmp_comm==1) {
	       lammps_command(lmp, "run 0");		// important !! to keep computes current
	    }

	    //---------------------------------------//
	    /* Get global variables from lammps data */
	    //---------------------------------------//
	    //
	    // Get nlocal for EACH lammps proc
	    if (lmp_comm==1) {
	       itmptr1	=	lammps_extract_global(lmp, "nlocal");
	       nlocal	=	*itmptr1;

	       printf("Proc %d lammps proc %d nlocal %d\n", MPI_myid, MPI_lammpsid, nlocal);
	       fflush(stdout);

	    }
	    
	    // Get box dimensions
	    if (lmp_comm==1) {

	       if (MPI_myid==1) {				// same data for all procs, only need one

		  xlo = *((double *)lammps_extract_global(lmp,"boxxlo"));
		  xhi = *((double *)lammps_extract_global(lmp,"boxxhi"));
		  ylo = *((double *)lammps_extract_global(lmp,"boxylo"));
		  yhi = *((double *)lammps_extract_global(lmp,"boxyhi"));
		  zlo = *((double *)lammps_extract_global(lmp,"boxzlo"));
		  zhi = *((double *)lammps_extract_global(lmp,"boxzhi"));
		  xy = *((double *)lammps_extract_global(lmp,"xy"));
		  xz = *((double *)lammps_extract_global(lmp,"xz"));
		  yz = *((double *)lammps_extract_global(lmp,"yz"));

		  MPI_Send(&xlo, 1, MPI_DOUBLE, 0, 0, comm_world);		// send data to proc 0
		  MPI_Send(&xhi, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&ylo, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&yhi, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&zlo, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&zhi, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&xy, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&xz, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&yz, 1, MPI_DOUBLE, 0, 0, comm_world);
	       }
	    }
	    if (MPI_myid==0) {
	       MPI_Recv(&xlo, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);		// receive data from proc 1
	       MPI_Recv(&xhi, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&ylo, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&yhi, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&zlo, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&zhi, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&xy, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&xz, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&yz, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);

	       system = 0;
	       BOX[system].xlo = xlo;
	       BOX[system].xhi = xhi;
	       BOX[system].ylo = ylo;
	       BOX[system].yhi = yhi;
	       BOX[system].zlo = zlo;
	       BOX[system].zhi = zhi;

	       BOX[system].lx = xhi-xlo;
	       BOX[system].ly = yhi-ylo;
	       BOX[system].lz = zhi-zlo;
	       BOX[system].xy = xy;
	       BOX[system].xz = xz; 
	       BOX[system].yz = yz;
	    }

            /* Extract atom id from LAMMPS internal data structure (atom->tag) */
	    if (lmp_comm==1) {
	       ldata1D	=	lammps_extract_atom(lmp, "id");
	       
	       for (i=0; i<10; i++) {
		  //printf("# DEBUG: lammps proc %d %d\n", MPI_lammpsid, ldata1D[i]);
	       }
	       
	       MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
	       MPI_Send(ldata1D, nlocal, MPI_INT, 0, 0, comm_world);	// send atom id array to proc 0
	    }
	    if (MPI_myid==0){

	       natoms = 0;
	       for (i=1; i<MPI_numprocs; i++) {					// in the order of rank
		  MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		  MPI_Recv(ldata1D, n, MPI_INT, i, 0, comm_world, &MPIstatus);
		  //printf("n= %d\n", n);

		  for (j=0; j<n; j++) {
		     atomid[natoms+j]	=	ldata1D[j] - 1;		// atom id starts with 1 in lammps
		     //printf("id %d %d\n", natoms+j, atomid[natoms+j]);
		  }
		  natoms	+=	n;
	       }
	       //printf("\n\n");
	       if (natoms!=NSITES) {
		  printf("# Error: natoms!=NSITES.\n");
		  fflush(stdout);
		  exit(0);
	       }
	    }
	    MPI_Barrier(comm_world);		// sync because ldata1D will be reused next

            /* Extract atom type from LAMMPS internal data structure (atom->tag) */
	    /* added 10/13/2016 for nucleation project, due to atom type swapping*/
	    
	    if (lmp_comm==1) {
	       ldata1D	=	lammps_extract_atom(lmp, "type");

	       MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
	       MPI_Send(ldata1D, nlocal, MPI_INT, 0, 0, comm_world);	// send atom id array to proc 0
	    }
	    if (MPI_myid==0){

	       natoms = 0;
	       for (i=1; i<MPI_numprocs; i++) {					// in the order of rank
		  MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		  MPI_Recv(ldata1D, n, MPI_INT, i, 0, comm_world, &MPIstatus);

		  for (j=0; j<n; j++) {
		     id = atomid[natoms+j];
		     //printf("id-type %d %d %d\n", natoms+j, id, ldata1D[j]-1);

		     mol[id].type[0] = ldata1D[j] - 1;		// LAMMPS atom type starts with 1
		     
		     //if (mol[id].type[0] != ldata1D[j]-1)
	                 //printf("id-mismatch %d type %d %d\n", id, mol[id].type[0], ldata1D[j]-1);
		  }
		  natoms += n;
	       }
	       //printf("\n\n");
	    }
	    MPI_Barrier(comm_world);
	    

	    if (lmp_comm==1) {
	       lammps_command(lmp, "run 0");		// important !! to keep computes current
	    }

	    /* Get nlocal for each lammps proc */
	    /*
	    if (lmp_comm==1) {
	       itmptr1	=	lammps_extract_global(lmp, "nlocal");
	       nlocal	=	*itmptr1;

	       printf("Proc %d lammps proc %d nlocal %d\n", MPI_myid, MPI_lammpsid, nlocal);
	       fflush(stdout);
	    }
	    */
            /* Extract atom id from LAMMPS internal data structure (atom->tag) */
	    /*
	    if (lmp_comm==1) {
	       ldata1D	=	lammps_extract_atom(lmp, "id");
	       
	       //for (i=0; i<10; i++) {
		 // printf("# DEBUG: lammps proc %d %d\n", MPI_lammpsid, ldata1D[i]);
	       //}
	       
	       MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
	       MPI_Send(ldata1D, nlocal, MPI_INT, 0, 0, comm_world);	// send atom id array to proc 0
	    }
	    if (MPI_myid==0){

	       natoms = 0;
	       for (i=1; i<MPI_numprocs; i++) {					// in the order of rank
		  MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		  MPI_Recv(ldata1D, n, MPI_INT, i, 0, comm_world, &MPIstatus);

		  for (j=0; j<n; j++) {
		     atomid[natoms+j]	=	ldata1D[j] - 1;		// atom id starts with 1 in lammps
		     printf("id %d\n", atomid[natoms+j]);
		  }
		  natoms	+=	n;
	       }
	       if (natoms!=NSITES) {
		  printf("# Error: natoms!=NSITES.\n");
		  fflush(stdout);
		  exit(0);
	       }
	    }
	    MPI_Barrier(comm_world);		// sync because ldata1D will be reused next
	    */

	    // get unwrapped atom coordinates from LAMMPS
	    if (lmp_comm==1) {
	       data2D	=	lammps_extract_compute(lmp, "unwrap", 1, 2);
	       //data2D	=	lammps_extract_atom(lmp, "x");
	       if (data2D == NULL) { printf("\n# Error: compute unwrap not found.\n");  fflush(stdout); }
	       
	       data1D	=	calloc(nlocal*3, sizeof(double));
	       for (i=0; i<nlocal; i++) {
		  for (j=0; j<3; j++) {
		     data1D[3*i+j] = data2D[i][j];
		  }
	       }
	       MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
	       MPI_Send(data1D, 3*nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
	       free(data1D);
	    }
	    if (MPI_myid==0) {
	       natoms = 0;
	       for (i=1; i<MPI_numprocs; i++) {
		  MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		  MPI_Recv(data1D, 3*n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

		  for (j=0; j<n; j++) {
		     id = atomid[natoms+j];

		     // Expected to be the same ONLY if "read_dump replace yes"
		     // ... if read_dump purge and add, then not the same in general, 12/6/2018
		     /* 
		     if (fabs(mol[id].p[0].x-data1D[3*j])>1e-8)
			printf("id %d x %f %f\n", id, mol[id].p[0].x, data1D[3*j]);
		     if (fabs(mol[id].p[0].y-data1D[3*j+1])>1e-8)
			printf("id %d y %f %f\n", id, mol[id].p[0].y, data1D[3*j+1]);
		     if (fabs(mol[id].p[0].z-data1D[3*j+2])>1e-8)
			printf("id %d z %f %f\n", id, mol[id].p[0].z, data1D[3*j+2]);
	             */
		     
		     mol[id].p[0].x = data1D[3*j];
		     mol[id].p[0].y = data1D[3*j+1];
		     mol[id].p[0].z = data1D[3*j+2];
		  }
		  natoms += n;
	       }
	    }
	    MPI_Barrier(comm_world);		// sync because coordinate copying takes time


	    /* Let LAMMPS perform a run or a minimization */

	    // !!! timestep could change if a run/minimization is performed !!!

	    if (lmp_comm==1) {
	       // Note: The following LAMMPS commands will be executed for EACH snapshot read
	       //
	       //lammps_command(lmp, "fix 1 all nve");	// keep computes and thermo variables current
	       //lammps_command(lmp, "run 0");		// keep computes and thermo variables current

	       //lammps_command(lmp, "fix 3 all box/relax aniso 0.0");
	       //lammps_command(lmp, "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");
	       //lammps_command(lmp, "unfix 3");
	       //
	       if (!strcmp(task, "dislocation")) {
		  lammps_command(lmp, "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");
		  //lammps_command(lmp, "run 0");		// important !! to keep computes current
	       }
	       else if (!strcmp(task, "profile")) {
		  //lammps_command(lmp, "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");
		  lammps_command(lmp, "run 0");		// important !! to keep computes current
	       }
	       else {
		  lammps_command(lmp, "run 0");		// important !! to keep computes current
	       }
	    }

	    /* Extract computes */

	    /* Extract compute pe/atom */
	    if (lmp_comm==1) {
	       data1D	=	lammps_extract_compute(lmp, "pe", 1, 1);
	       if (data1D == NULL) {
		  printf("\n# Error: Proc. %d compute pe/atom not found.\n", MPI_myid);
		  fflush(stdout);
	       }
	       MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
	       MPI_Send(data1D, nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
	    }
	    if (MPI_myid==0) {
	       natoms = 0;
	       for (i=1; i<MPI_numprocs; i++) {
		  MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		  MPI_Recv(data1D, n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

		  for (j=0; j<n; j++) {
		     id = atomid[natoms+j];
		     /*
		     if (fabs(mol[id].pe[0]-data1D[j])>2e-4) {
			printf("id %d %f %f\n", id, mol[id].pe[0], data1D[j]);
		     }
		     */
		     mol[id].pe[0] = data1D[j];
		  }
		  natoms += n;
	       }
	    }
	    MPI_Barrier(comm_world);		// sync because copying takes time

	    /* Extract compute stress/atom, added 06/18/2017 */
	    
	    if (lmp_comm==1) {
	       data2D	=	lammps_extract_compute(lmp, "st", 1, 2);
	       if (data2D == NULL) { printf("\n# Error: compute st not found.\n");  fflush(stdout); }
	       
	       data1D	=	calloc(nlocal*6, sizeof(double));
	       if (data1D == NULL) { printf("\n# Error: compute st data1D allocation.\n");  fflush(stdout); }
	       for (i=0; i<nlocal; i++) {
		  for (j=0; j<6; j++) {
		     data1D[6*i+j] = data2D[i][j];
		  }
	       }
	       MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
	       MPI_Send(data1D, 6*nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
	       free(data1D);
	    }
	    if (MPI_myid==0) {
	       natoms = 0;
	       for (i=1; i<MPI_numprocs; i++) {
		  MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		  MPI_Recv(data1D, 6*n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

		  for (j=0; j<n; j++) {
		     id = atomid[natoms+j];

		     for (k=0; k<6; k++) {
                        mol[id].stress[0][k] = data1D[6*j+k];
		     }
		  }
		  natoms += n;
	       }
	    }
	    MPI_Barrier(comm_world);		// sync because copying takes time
	    

	    /* Extract compute orientorder/atom, added 03/13/2019 */
	    
	    if (lmp_comm==1) {
	       data2D	=	lammps_extract_compute(lmp, "bondorder", 1, 2);
	       if (data2D == NULL) { printf("\n# Error: compute bondorder not found.\n");  fflush(stdout); }
	       
	       data1D	=	calloc(nlocal*2, sizeof(double));
	       if (data1D == NULL) { printf("\n# Error: compute bondorder data1D allocation.\n");  fflush(stdout); }
	       for (i=0; i<nlocal; i++) {
		  for (j=0; j<2; j++) {
		     data1D[2*i+j] = data2D[i][j];
		  }
	       }
	       MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
	       MPI_Send(data1D, 2*nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
	       free(data1D);
	    }
	    if (MPI_myid==0) {
	       natoms = 0;
	       for (i=1; i<MPI_numprocs; i++) {
		  MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		  MPI_Recv(data1D, 2*n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

		  for (j=0; j<n; j++) {
		     id = atomid[natoms+j];

		     for (k=0; k<2; k++) {
                        mol[id].ql[0][k] = data1D[2*j+k];
		     }
		  }
		  natoms += n;
	       }
	    }
	    MPI_Barrier(comm_world);		// sync because copying takes time
	    
	    /* Extract compute cna/atom 1, cutoff 3.5A */
	    if (lmp_comm==1) {
	       data1D	=	lammps_extract_compute(lmp, "cna", 1, 1);
	       if (data1D == NULL) { printf("\n# Error: compute cna not found.\n"); fflush(stdout); }
	       /*
	       for (i=0; i<10; i++) {
		  printf("# DEBUG: lammps proc %d %f\n", MPI_lammpsid, data1D[i]);
	       }
	       */
	       MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
	       MPI_Send(data1D, nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
	    }
	    if (MPI_myid==0) {
	       natoms = 0;
	       for (i=1; i<MPI_numprocs; i++) {
		  MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		  MPI_Recv(data1D, n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

		  for (j=0; j<n; j++) {
		     id = atomid[natoms+j];
		     /*
		     if (fabs(mol[id].cna[0]-data1D[j])>0.001) {
			printf("id %d %d %f\n", id, mol[id].cna[0], data1D[j]);
		     }
		     */
		     mol[id].cna[0] = (int)(data1D[j]+0.5);
		  }
		  natoms += n;
	       }
	    }


	    /* Extract compute cna/atom 2, cutoff 3.0A */
	    //if (!strcmp(task, "nucleation")) {
	    
	    /*
	    char compute[16][32];              		// maximum 16 computes from LAMMPS
	    strcpy(compute[0], "cna_2");
	    strcpy(compute[1], "cna_3");
	    short *sptr;

	    for (i=0; i<2; i++) {
	       if (!strcmp(compute[0],"cna_2")) {
		  sptr = mol
		  */

	       if (lmp_comm==1) {
		  data1D	=	lammps_extract_compute(lmp, "cna_2", 1, 1);
		  if (data1D == NULL) { printf("\n# Error: compute cna_2 not found.\n"); fflush(stdout); }
		  MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
		  MPI_Send(data1D, nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
	       }
	       if (MPI_myid==0) {
		  natoms = 0;
		  for (i=1; i<MPI_numprocs; i++) {
		     MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		     MPI_Recv(data1D, n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

		     for (j=0; j<n; j++) {
			id = atomid[natoms+j];
			mol[id].cna2[0] = (int)(data1D[j]+0.5);
		     }
		     natoms += n;
		  }
	       }


	       /* Extract compute cna/atom 3 */
	       if (lmp_comm==1) {
		  int signal = 1;
		  data1D	=	lammps_extract_compute(lmp, "cna_3", 1, 1);
		  if (data1D == NULL) {
		     signal = 0;
		     printf("\n# Error: compute cna_3 not found.\n"); 
		     fflush(stdout); 
		  }
		  MPI_Send(&signal, 1, MPI_INT, 0, 0, comm_world);		// send signal to proc 0

		  if (signal) {
		     MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
		     MPI_Send(data1D, nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
		  }
	       }
	       if (MPI_myid==0) {
		  int signal = 0;
		  natoms = 0;
		  for (i=1; i<MPI_numprocs; i++) {
		     MPI_Recv(&signal, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive signal from proc i

		     if (signal) {
			MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
			MPI_Recv(data1D, n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

			for (j=0; j<n; j++) {
			   id = atomid[natoms+j];
			   mol[id].cna3[0] = (int)(data1D[j]+0.5);
			}
			natoms += n;
		     }
		  }
	       }

	       /* Extract compute cluster/atom */
	       if (lmp_comm==1) {
		  int signal = 1;

		  data1D	=	lammps_extract_compute(lmp, "unknownclusters", 1, 1);
		  if (data1D == NULL) {
		     signal = 0;
		     printf("\n# Error: compute beta not found.\n"); 
		     fflush(stdout); 
		  }
		  MPI_Send(&signal, 1, MPI_INT, 0, 0, comm_world);		// send signal to proc 0

		  if (signal) {
		     MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
		     MPI_Send(data1D, nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
		  }
	       }
	       if (MPI_myid==0) {
		  int signal = 0;

		  natoms = 0;
		  for (i=1; i<MPI_numprocs; i++) {
		     MPI_Recv(&signal, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive signal from proc i

		     if (signal) {
			MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);		// receive data from proc i
			MPI_Recv(data1D, n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

			for (j=0; j<n; j++) {
			   id = atomid[natoms+j];
			   mol[id].clus[0] = (int)(data1D[j]+0.5);
			}
			natoms += n;
		     }
		  }
	       }

	    //}

	    /* Extract compute nneigh/atom 1, cutoff 3.5A */
	    if (lmp_comm==1) {
	       data1D	=	lammps_extract_compute(lmp, "nneigh", 1, 1);
	       if (data1D == NULL) { printf("\n# Error: compute nneigh not found.\n"); fflush(stdout); }
	       MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
               MPI_Send(data1D, nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
	    }
	    if (MPI_myid==0) {
	       natoms = 0;
	       for (i=1; i<MPI_numprocs; i++) {
		  MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		  MPI_Recv(data1D, n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

		  for (j=0; j<n; j++) {
		     id = atomid[natoms+j];
		     mol[id].nneigh[0] = (int)(data1D[j]+0.5);
		  }
		  natoms += n;
	       }
	    }

	    /* Extract compute nneigh/atom 2, cutoff 3.0A */
	    if (!strcmp(task, "nucleation")) {
	       if (lmp_comm==1) {
		  data1D	=	lammps_extract_compute(lmp, "nneigh_2", 1, 1);
		  if (data1D == NULL) { printf("\n# Error: compute nneigh not found.\n"); fflush(stdout); }
		  MPI_Send(&nlocal, 1, MPI_INT, 0, 0, comm_world);		// send nlocal to proc 0
		  MPI_Send(data1D, nlocal, MPI_DOUBLE, 0, 0, comm_world);	// send compute array to proc 0
	       }
	       if (MPI_myid==0) {
		  natoms = 0;
		  for (i=1; i<MPI_numprocs; i++) {
		     MPI_Recv(&n, 1, MPI_INT, i, 0, comm_world, &MPIstatus);	// receive data from proc i
		     MPI_Recv(data1D, n, MPI_DOUBLE, i, 0, comm_world, &MPIstatus);

		     for (j=0; j<n; j++) {
			id = atomid[natoms+j];
			mol[id].nneigh2[0] = (int)(data1D[j]+0.5);
		     }
		     natoms += n;
		  }
	       }
	    } 

	    /* Handling atom coordinates */
	    // moved to the front

	    /* nuclei analysis to adjust atoms to display 5/18/15*/
	    /*
	    if (MPI_myid==0) {
	       system	= 0;
	       i	= 1;
	       while (sizeofnucl[i] != nmax[system][0]) {
		  i ++;
	       }
	       nmaxid	=	i;

	       //data1D[
	    }
	    if (lmp_proc) {
	       double *cna = (double *) calloc(NSITES, sizeof(double));
	    }
	    MPI_Bcast(&data1D, NSITES, MPI_DOUBLE, 0, comm_world);
	    if (lmp_proc) {
	       lammps_scatter_atom(lmp, "cna", );
	       free(cna);
	    }
	    */

	    /* 5/18/15
	    // store coordinates from the dynamics simulation trajectory
	    if (MPI_myid == 0) {
	       n = 0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     r_MD[n]	=	moli->p[i].x;
		     r_MD[n+1]	=	moli->p[i].y;
		     r_MD[n+2]	=	moli->p[i].z;
		     n	+=	3;
		  }
	       }
	    }
	    */

	    /*
	    // pass coordinates to LAMMPS
	    n = 0;
	    for (moli=mol; moli<mol+NMOLS; moli++) {
	       for (i=0; i<moli->nsites; i++) {
		  lammps_x[n]	=	moli->p[i].x;
		  lammps_x[n+1]	=	moli->p[i].y;
		  lammps_x[n+2]	=	moli->p[i].z;
		  n	=	n+3;
	       }
	    }
	    lammps_scatter_atoms(lmp, "x", 1, 3, lammps_x);
	    */


	    //--- Extract thermo variables from LAMMPS internal data structure
	    if (lmp_comm==1) {
	    //if (MPI_myid==1) {
	    
	       // Using old lammps library Feb16
	       /*
	       lammps_extract_thermo(lmp, "temp", &temp);
	       printf("extract_thermo temp = %f\n", temp);

	       lammps_extract_thermo(lmp, "pe", &pe);
	       lammps_extract_thermo(lmp, "enthalpy", &enthalpy);
	       lammps_extract_thermo(lmp, "vol", &vol);
	       lammps_extract_thermo(lmp, "density", &rho_m);
	       lammps_extract_thermo(lmp, "atoms", &atoms);
	       lammps_extract_thermo(lmp, "press", &pressure);
	       //lammps_extract_thermo(lmp, "pxx", stress);
	       //lammps_extract_thermo(lmp, "pyy", stress+1);
	       //lammps_extract_thermo(lmp, "pzz", stress+2);
	       */
	       tmptr = lammps_extract_variable(lmp, "T", "NULL");
	       temp = *tmptr;
	       //printf("extract_variable temp = %f\n", temp);

	       // lammps_get_thermo, using new lammps library Mar17
	       
	       pe 	= lammps_get_thermo(lmp, "pe");
	       enthalpy = lammps_get_thermo(lmp, "enthalpy");
	       vol 	= lammps_get_thermo(lmp, "vol");
	       rho_m 	= lammps_get_thermo(lmp, "density");
	       atoms 	= lammps_get_thermo(lmp, "atoms");
	       //atoms 	= lammps_get_natoms(lmp);
	       pressure	= lammps_get_thermo(lmp, "press");
	       stress[0]= lammps_get_thermo(lmp, "pxx");
	       stress[1]= lammps_get_thermo(lmp, "pyy");
	       stress[2]= lammps_get_thermo(lmp, "pzz");
	       
	       rho_n = atoms/vol;

               presski = 138 * temp;          		// need to multiply natoms/vol to get pressure
	       presske = atoms * presski / vol;  	// kinetic contribution to pressure, unit bar

	       printf("Proc %d step %d %f %f %f %f %f %f %f %f\n", 
	                MPI_myid, timestep, pe, enthalpy, pressure, pressure+presske, vol, rho_m, atoms, rho_n);
	       fflush(stdout);

	       if (MPI_myid==1) {				// same data for all procs, only need one
		  printf("lammps_get_thermo: pe %f\n", pe);
		  printf("lammps_get_thermo: enthalpy %f\n", enthalpy);
		  printf("lammps_get_thermo: press %f\n", pressure);
		  printf("lammps_get_thermo: vol %f\n", vol);
		  printf("lammps_get_thermo: atoms %d\n", (int)atoms);
	          printf("temperature %f\n", temp);
	          printf("presske %f\n", presske);
		  printf("pressure: %lf %lf\n", pressure, (stress[0]+stress[1]+stress[2])/3);

		  MPI_Send(&pe, 1, MPI_DOUBLE, 0, 0, comm_world);		// send data to proc 0
		  MPI_Send(&enthalpy, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&vol, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&rho_m, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&atoms, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&rho_n, 1, MPI_DOUBLE, 0, 0, comm_world);

		  MPI_Send(&pressure, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&presski, 1, MPI_DOUBLE, 0, 0, comm_world);
		  MPI_Send(&presske, 1, MPI_DOUBLE, 0, 0, comm_world);
	       }
	    }
	    if (MPI_myid==0) {
	       MPI_Recv(&pe, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);		// receive data from proc 1
	       MPI_Recv(&enthalpy, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&vol, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&rho_m, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&atoms, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&rho_n, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&pressure, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&presski, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);
	       MPI_Recv(&presske, 1, MPI_DOUBLE, 1, 0, comm_world, &MPIstatus);

	       itmp1 = 0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     if (moli->type[i]==1) {
			itmp1 ++;
		     }
		  }
	       }
	       printf("lammps_get_thermo: pe %f\n", pe);
	       printf("lammps_get_thermo: enthalpy %f\n", enthalpy);
	       printf("lammps_get_thermo: press %f\n", pressure);
	       printf("lammps_get_thermo: vol %f\n", vol);
	       printf("lammps_get_thermo: atoms %d\n", (int)atoms);
	       printf("presske %f\n", presske);
	       printf("pressure: %lf %lf\n", pressure, (stress[0]+stress[1]+stress[2])/3);
	       
	       printf("Proc %d step %d %f %f %f %f %f %f %f %f\n", 
	                MPI_myid, timestep, pe, enthalpy, pressure, pressure+presske, vol, rho_m, atoms, rho_n);
	       printf("Ntype1: %d\n", itmp1);
	       fflush(stdout);
	    }
	    MPI_Barrier(comm_world);		// sync because coordinate copying takes time

	    // Sanity check 6/18/2017
	    if (MPI_myid==0) {
	       tmp1 = 0.0;
	       tmp2 = 0.0;
	       tmp3 = 0.0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     tmp1 += moli->stress[i][0];
		     tmp2 += moli->stress[i][1];
		     tmp3 += moli->stress[i][2];
		  }
	       }
	       printf("check pressure: %f %f %f %8.3f %8.3f %8.3f\n", 
	               pressure, pressure+presske, -(tmp1+tmp2+tmp3)/(3*vol), -tmp1, -tmp2, -tmp3);
	    }

	    //--- obtain coordinates (wrapped) from LAMMPS
	    /*
	    lammps_gather_atoms(lmp, "x", 1, 3, lammps_x);
	    
	    n = 0;
	    for (moli=mol; moli<mol+NMOLS; moli++) {
	       for (i=0; i<moli->nsites; i++) {
		  moli->p[i].x	=	lammps_x[n];
		  moli->p[i].y	=	lammps_x[n+1];
		  moli->p[i].z	=	lammps_x[n+2];
		  n	=	n+3;
	       }
	    }
	    */

	    /* 5/18/15
 	    //--- Extract compute results from LAMMPS internal data structure
	    // compute rdf
	    if (lmp_comm==1) {
	       data2D	=	lammps_extract_compute(lmp, "myRDF", 0, 2);
	       if (data2D == NULL) {
		  printf("\n# Error: compute myRDF not found.\n");
		  fflush(stdout);
	       }
	       for (i=0; i<50; i++) {
		  //printf("%f %f %f\n", data2D[i][0], data2D[i][1], data2D[i][2]);
	       }
	       //fflush(stdout);
	    }
	    */

	    /* 5/18/15
	    //--- locate dislocation for velocity calculation
	    if (MPI_myid==0) {
	       n = 0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     if (moli->cna[i]==1) {
			nucleus[n].moli = moli;
			nucleus[n].site = i;
			n ++;
		     }
		  }
	       }
	       Rconn = 6.0;

	       MapInNucleus(system, nucleus, n, rbead);	// continuous coordinates in rbead

	       xdisl = 0.0;
	       for (i=0; i<n; i++) {
		  xdisl += rbead[i].x;
	       }
	       xdisl /= n;

	       if (prev_xdisl > -1e5) {		// not the first time
		  vtmp.x = xdisl-prev_xdisl;
		  vtmp.y = 0.0;
		  vtmp.z = 0.0;
		  vtmp = MapInBox2(&vtmp, PBC, system);

		  xdisl = prev_xdisl + vtmp.x;
	       }
	    }
	    */

	    /*
	    if (lmp_comm==1) {
	       lammps_extract_thermo(lmp, "pxz", &tmp1);	// pxz calc. from Virial

	       tmptr = (double *)lammps_extract_variable(lmp, "mobile_pxz", NULL);
	       tmp2 = *tmptr;
	       tmptr = (double *)lammps_extract_variable(lmp, "boundary_pxz", NULL);
	       tmp3 = *tmptr;

	       tmptr = (double *) lammps_extract_compute(lmp, "fx_top", 0, 0);
	       tmp4 = *tmptr;
	       tmptr = (double *) lammps_extract_compute(lmp, "fx_bot", 0, 0);
	       tmp5 = *tmptr;
	       tmptr = (double *) lammps_extract_variable(lmp, "stress_top", NULL);
	       tmp6 = *tmptr;
	       tmptr = (double *) lammps_extract_variable(lmp, "stress_bot", NULL);
	       tmp7 = *tmptr;
	    }
	    //... need to add code to transfer tmp1, tmp2, etc
	    if (MPI_myid==0) {
	       printf("disl: %d %f %f %f %f %f %f %f", timestep, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7);
	       printf(" %f %f\n", xdisl, xdisl-prev_xdisl);
	       fflush(stdout);
	       prev_xdisl = xdisl;
	    }
	    */

	    /* 5/18/15
	    // store coordinates of the minimized trajectory
	    if (MPI_myid==0) {
	       n = 0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     r_MIN[n]	=	moli->p[i].x;
		     r_MIN[n+1]	=	moli->p[i].y;
		     r_MIN[n+2]	=	moli->p[i].z;
		     n	+=	3;
		  }
	       }
	    }
	    */

	    /*
	    printf("lmp_atom_id\tx\ty\tz\tpe\tnneigh\tcna\tcentro\n");
	    for (moli=mol; moli<mol+NMOLS; moli++) {
	       for (i=0; i<moli->nsites; i++) {
		  printf("%d", moli-mol+1);
		  printf("\t%10.6f", moli->p[i].x);
		  printf("\t%10.6f", moli->p[i].y);
		  printf("\t%10.6f", moli->p[i].z);
		  printf("\t%10.6f", moli->pe[i]);
		  printf("\t%d", moli->nneigh[i]);
		  printf("\t%1.0f", (float)moli->cna[i]);
		  printf("\t%6.4f", moli->centro[i]);
		  printf("\n");
	       }
	    }
	    fflush(stdout);
	    */

	 }	// LAMMPS Module done
	 
	 MPI_Barrier(comm_world);
	 if (M2HDEBUG) {
	    printf("\n# DEBUG: Proc %d LAMMPS module completed.\n", MPI_myid);
	    fflush(stdout);
	 }

	 //----------------------------------------//
	 //-------- LAMMPS Module Finished --------//
	 //----------------------------------------//
	 //-------- Other Analysis Starts ---------//
	 //----------------------------------------//

	 if (MPI_myid==0) {

	    system = 0;

	    //-------- Task: nucleation -----//
	    if (!strcmp(task, "nucleation")) {

	       /* Nuclei characterization */

	       // for (j=0; j<4; j++) {		// for different crystal/clustering algorithms
	       for (j=0; j<1; j++) {		// for different crystal/clustering algorithms

		  if (j==0) {
		     find_nuclei_lmp();  		// find clusters using LAMMPS cluster/atom
		  }
		  else if (j==1){
		     find_nuclei_general("cna2", "unknown", "r", 3.0);	// Mg17Al12 in Mg/Al, at 550K
		  }
		  else if (j==2) {
		     find_nuclei_general("cna2", "bcc", "r", 3.8);	// bcc nuclei with cna cutoff 3.0A
		  }
		  else if (j==3) {
		     find_nuclei_general("cna2", "fcc/hcp", "r", 3.8);	// fcc/hcp nuclei with cna cutoff 3.0A
		  }

		  //---- Find the nuclid for the largest nucleus
		  //
		  for (i=1; i<=Nnucl[system]; i++) {
		     if (sizeofnucl[i]==nmax[system][0]) {
			nuclid1 = i;				// find nuclid of partial 1
			break;
		     }
		  }
		  for (i=1; i<=Nnucl[system]; i++) {
		     if (sizeofnucl[i]==nmax[system][1] && i!=nuclid1) {
			nuclid2 = i;				// find nuclid of partial 2
			break;
		     }
		  }

		  //---- Calcualte the radius of gyration for the largest nucleus
		  // 
		  if (nmax[system][0]>1) {
		     n=0;
		     for (moli=mol; moli<mol+NMOLS; moli++) {
			for (i=0; i<moli->nsites; i++) {
			   if (moli->nuclid[i] == nuclid1) {
			      nucleus[n].moli = moli;
			      nucleus[n].site = i;
			      n ++;
			   }
			}
		     }
		     //Rconn = 6.0;  // Rconn will change in MapInNucleus

		     MapInNucleus(system, nucleus, n, rbead, 6.0);  	// continuous coordinates in rbead
		     R_nmax = sqrt(group_Rg2(rbead, n)*5/3);   		// radius of gyration of the nucleus
		     							// 5/3 for sphere
		     vtmp1 = group_com(rbead, n);     			// geometrical center, no mass weight
		     
		     // gyration tensor of the largest nucleus
		     matrix tensor, evec;
		     vector eval;

		     tensor = grpgytensor(n, rbead); 
		     eigensol(&tensor, &eval, &evec);			// solve eigenvalue problem
		     printf("eigen: %7d %4.1f %4.1f %4.1f %4.1f\n", 
		     	timestep, R_nmax, sqrt(eval.x), sqrt(eval.y), sqrt(eval.z));

		     // Calculate # of atoms within the spherical nucleus ...
		     // ... (in case for bad xtal order parameter)
		     //
		     itmp1=0;
		     for (moli=mol; moli<mol+NMOLS; moli++) {
			for (i=0; i<moli->nsites; i++) {
			   if (DistSQ(moli->p[i], vtmp1, system) <= R_nmax*R_nmax) {
			      itmp1++; 
			   }
			}
		     }
		     while(vtmp1.x < xlo) { 			// move x coordinate within central box
			vtmp1.x += BOX[system].lx;
		     }
		     while (vtmp1.x >= xhi) {
			vtmp1.x -= BOX[system].lx;
		     }
		  }
		  else {
		     R_nmax = 0.0;
		     itmp1 = 0; 				// Number of atoms in the sphere R_nmax
		  }

		  if (j==0) {
		     tmp3 = 4*M_PI*R_nmax*R_nmax;		// Surface area of the sphere 
		     tmp4 = tmp3*R_nmax/3;			// Volume of the sphere
		     printf("analysis: %7d", timestep);
		     printf("%6d %6d %5.2f %6.3f %8.3f %8.3f %10.3f %4.3f %5.3f %5.3f %10.3f %10.3f\n",
			     nmax[system][0], itmp1, vtmp1.x, R_nmax, tmp3, tmp4, vol, rho_n,
			     itmp1/tmp4, (atoms-itmp1)/(vol-tmp4), pe, enthalpy);
		  }

		  //---- Calculate the composition in the largest nucleus
		  //
		  itmp1 = 0;
		  n=0;
		  for (moli=mol; moli<mol+NMOLS; moli++) {
		     for (i=0; i<moli->nsites; i++) {
			if (moli->nuclid[i] == nuclid1) {
			   if (moli->type[i] == 0) {			// type 0 atoms	
			      itmp1 ++;
			   }
			   n ++;
			}
		  }  }
		  c_nmax = (itmp1>0 ? 1.0*itmp1 / nmax[system][0] : 0.0);	// composition of type 0 atoms

		  //---- Realtime composition gradient in the system
		  if (j==0) {
		     //vtmp1.x = vtmp1.y = vtmp1.z = 0.0;
		     composition_profile(timestep, 'z', vtmp1, &local_comp, &local_grad, presski);
		     profile_2D(timestep, "xy", "z", presski);
		  }

		  system = 0;
		  itmp3 = 0;
		  itmp4 = 0;
		  itmp5 = 0;
		  itmp6 = 0;

		  //printf("PBC = %d system = %d box.lx = %f\n", PBC, system, BOX[system].lx);
		  for (moli=mol; moli < mol + NMOLS; moli++) {
		     for (i=0; i < moli->nsites; i++) {
	    		vtmp = MapInBox2(moli->p+i, PBC, system);
			//printf("%d %d %f %f %f %d %d\n", moli-mol, moli->type[i], 
			  //      moli->p[i].x, moli->p[i].y, moli->p[i].z,
			  //    moli->nneigh[i], moli->nneigh2[i]);

			if ( fabs(vtmp.x) <= 5.0) {
			   if (moli->type[i] == 0)
			      itmp3 ++;
			   else
			      itmp4 ++;
			}
			else if ( 0.5*BOX[system].lx - fabs(vtmp.x) <= 5.0) {
			   if (moli->type[i] == 0)
			      itmp5 ++;
			   else
			      itmp6 ++;
			}
			c_grad = (1.0*itmp5)/(itmp5+itmp6) - (1.0*itmp3)/(itmp3+itmp4);
		     }
		  }
		  
		  //---- Print out
		  if (j==0) {
		     printf("nuclei: %7d", timestep);
		  }
		  printf(" %4d", Nnucl[system]);
		  printf(" %5d", Xtal[system]);
		  printf(" %5d", nmax[system][0]);
		  printf(" %5d", nmax[system][1]);
		  printf(" %4.1f", vtmp1.x); 		// nmax position
		  printf(" %4.1f", R_nmax); 		// nmax size
		  printf(" %4.3f", c_nmax);		// nmax composition
		  if (j==0) {
		     // nmax local composition gradient, global gradient
		     printf(" %4.3f %4.3f %4.3f", local_comp, local_grad, c_grad);  
		  }
	       }
	       printf("\n");
	       fflush(stdout);
	    } // end: if task == nucleation

	    //-------- Task: profile --------//
	    //
	    else if (!strcmp(task, "profile")) {

	       //vtmp1.x = vtmp1.y = vtmp1.z = 0.0;
	       //composition_profile(timestep, 'z', vtmp1, &local_comp, &local_grad, presski);
	       profile_2D(timestep, "x", "z", presski);

	       itmp1 = 0;
	       itmp2 = 0;
	       itmp3 = 0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     /* ni-al solid-interface
		     if (moli->cna[i] == 2) {               // hcp, cutoff 3.1
			itmp1 ++;
		     }
		     if (moli->cna2[i] == 3) {              // bcc, cutoff 3.5
			itmp2 ++;
		     }
		     if (moli->cna3[i] == 3) {               // bcc, cutoff 3.9
			itmp3 ++;
		     }
		     */
		     if (moli->cna[i] == 1) {               // fcc
			itmp1 ++;
		     }
		     if (moli->cna[i] == 2) {               // fcc
			itmp2 ++;
		     }
		     if (moli->cna[i] == 3) {               // fcc
			itmp3 ++;
		     }
		  }
	       }
	       printf("nuclei: %d %d %d %d\n", timestep, itmp1, itmp2, itmp3);
	       fflush(stdout);
	    } // end: if task == profile

	    //-------- Task: dislocation --------//
	    //
	    else if (!strcmp(task, "dislocation")) {
	       printf("dislocation...\n");

	       // Exclude atoms in top and bottom slabs
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     if (moli->type[i]==1 || moli->type[i]==2) {
			moli->cna[i] = 2;
		     }
		  }
	       }

	       // Find layers in z direction, used in twin analysis
	       float zcutoff;  		// 0.5 for 77K, 0.1 for 500K
	       double dz; 		// twin layer heigh, 3.8 for 77K, 3.83 for 500K
	       int ncutoff = 160;
	       itmp1 = 500;

	       if (itmp1==77) {
		  zcutoff = 0.5;
		  dz = 3.8;
	       }
	       else if (itmp1==500) {
		  zcutoff = 0.1;
		  dz = 3.83;
	       }

	       find_nuclei_general("cna", "unknown", "z", zcutoff);

	       printf("TB atoms %6d TB layers %6d\n", Xtal[system], Nnucl[system]);

	       // Find the average z coordinate of each layer
               double zlayer[256];
	       size_t id[256];

	       for (j=1; j<=Nnucl[system]; j++) {
		  tmp1 = 0.0;
		  itmp1 = 0;
		  for (moli=mol; moli<mol+NMOLS; moli++) {
		     for (i=0; i<moli->nsites; i++) {
			if (moli->nuclid[i] == j) {
			   tmp1 += moli->p[i].z;
			   itmp1 ++;
			}
		     }
		  }
		  zlayer[j-1] = tmp1/itmp1;
		  id[j-1] = j;
	       }
	       
	       // Sort the layers based on z coordinate, in ascending order
	       gsl_sort_index(id, zlayer, 1, Nnucl[system]); 

	       // Find the right layers
	       int first = 1;
	       int firstid;

	       for (i=Nnucl[system]-1; i>=0; i--) {
		  if (sizeofnucl[id[i]+1] >= ncutoff) {
		     if (first==1) { 		// examine the first layer
			//if (max > 0) {  	// not the first configuration
			   // check max
			//   tmp1 = zlayer[id[i]] - max;
			//   if (fabs((tmp1/dz) - round(tmp1/dz)) > 0.2) { 	// height not correct
			//      continue;
			//   }
			//}
			tmp1 = zlayer[id[i]]/dz + 0.5;
			if (fabs( tmp1 - round(tmp1) ) > 0.3) {
			   continue;
			}
			first=0;
			max = zlayer[id[i]];
			firstid = id[i]+1;

			printf("layer+ %4d %6d %8.3f\n", id[i]+1, sizeofnucl[id[i]+1], zlayer[id[i]]);
		     } else {
			printf("layer %4d %6d %8.3f\n", id[i]+1, sizeofnucl[id[i]+1], zlayer[id[i]]);
		     }
		     min = zlayer[id[i]];
		  }
	       }
	       // check min
	       tmp1 = max - min;
	       if (fabs((tmp1/dz) - round(tmp1/dz)) > 0.3) {   	// TB height incorrect
		  for (i=0; i<Nnucl[system]; i++) {
		     if (zlayer[id[i]] > min + 1e-6) {  	// closest layer above
			tmp1 = max - zlayer[id[i]];
			if (fabs((tmp1/dz) - round(tmp1/dz)) < 0.3) {
			   min = zlayer[id[i]];
			   break;
			}
		     }
		  }
	       }

	       // Collect solute concentration in top layer
	       itmp1 = 0;
	       itmp2 = 0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     if ( fabs(moli->p[i].z - max) < 0.8) {
			itmp1 ++;
			if (moli->type[i] == 3) {
			   itmp2++;
			}
		     }
		  }
	       }

	       printf("\nTwin-boundary: %7d %8.3f %8.3f %8.3f %4d %4.3f", 
	                     timestep, min, max, max-min, itmp1, 1.0*itmp2/itmp1);

	       // Cluster analysis of the atoms on the top layer
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     if (moli->nuclid[i] == firstid) { moli->tmp[i] = 1; }
		     else { moli->tmp[i] = 0; }
		  }
	       }
	       find_nuclei_general("custom", "null", "r", 5.5);

               printf(" %4d %4d %4d %4d %4d %4d %4d\n",
	        		Xtal[system], Nnucl[system], realXtal[system], realNnucl[system],
				nmax[system][0], nmax[system][1], nmax[system][2]);
	       fflush(stdout);

goto dislocationdone;

	       /* Calculate stacking fault width */    // temporarily removed 2017-11-11

	       /* Calculate fraction of dislocation on the basal plane */  // temporarily removed 2017-11-11

	       /* Calculate stacking fault 2017-11-11 */

	       find_nuclei_general("cna", "hcp", "r", 3.8);

	       //printf("%d %d ", nmax[system][0], nmax[system][1]);
	       for (i=1; i<=Nnucl[system]; i++) {
		  if (sizeofnucl[i]==nmax[system][0]) {
		     nuclid1 = i;				// find nuclid of partial 1
		     break;
		  }
	       }
	       for (i=1; i<=Nnucl[system]; i++) {
		  if (sizeofnucl[i]==nmax[system][1] && i!=nuclid1) {
		     nuclid2 = i;				// find nuclid of partial 2
		     break;
		  }
	       }

	       // find the stacking fault atoms
	       n=0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     if (moli->nuclid[i] == nuclid1) {	// if in partial #1
			nucleus[n].moli = moli;
			nucleus[n].site = i;
			n ++;
		     }
		  }
	       }
	       MapInNucleus(system, nucleus, n, rbead, 15);

	       float max = -1e8; // max
	       float min = 1e8; // min

	       for (i=0; i<n; i++) {
		  if (rbead[i].x > max)  max = rbead[i].x;
		  if (rbead[i].x < min)  min = rbead[i].x;
	       }
	       //vtmp1 = group_com(rbead, n);
	       fprintf(fsf, "%d\t%f\t%f\t%f\n", timestep, max, min, max-min);
	       fflush(fsf);

	       // Profiling of bcc atoms, cna=3
	       binning("cna2", "z");


	       /* Calculate stacking fault width and partial roughening */
	       
	       /*
	       find_nuclei_general("cna", "partial");	// find atoms in dislocation
	       printf("%d %d ", nmax[system][0], nmax[system][1]);
	       for (i=1; i<=Nnucl[system]; i++) {
		  if (sizeofnucl[i]==nmax[system][0]) {
		     nuclid1 = i;				// find nuclid of partial 1
		     break;
		  }
	       }
	       for (i=1; i<=Nnucl[system]; i++) {
		  if (sizeofnucl[i]==nmax[system][1] && i!=nuclid1) {
		     nuclid2 = i;				// find nuclid of partial 2
		     break;
		  }
	       }

	       // one way to calculate average partial separation
	       n=0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     if (moli->nuclid[i] == nuclid1) {	// if in partial #1
			nucleus[n].moli = moli;
			nucleus[n].site = i;
			n ++;
		     }
		  }
	       }
	       Rconn = 4.8;
	       MapInNucleus(system, nucleus, n, rbead);
	       vtmp1 = group_com(rbead, n);

	       n=0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {
		     if (moli->nuclid[i] == nuclid2) {	// if in partial #2
			nucleus[n].moli = moli;
			nucleus[n].site = i;
			n ++;
		     }
		  }
	       }
	       Rconn = 4.8;
	       MapInNucleus(system, nucleus, n, rbead);
	       vtmp2 = group_com(rbead, n);

	       vtmp = V_Subtr(&vtmp1, &vtmp2);
	       vtmp = MapInBox2(&vtmp, PBC, system);
	       tmp3 = fabs(vtmp.y); 		// basal screw
	       tmp3 = fabs(vtmp.x);		// basal edge

	       printf("width = %f ", tmp3);
	       fprintf(fsf, "%d\t%f\n", timestep, tmp3);
	       fflush(fsf);
	       */

	       // second way to calculate partial separation // temporarily removed 2017-11-11

	       /* Calc. minimum separation between two partials */  // temporarily removed 2017-11-11

	       /* roughening */  // temporarily removed 2017-11-11

	       /* Vacancy analysis */  // temporarily removed 2017-11-11
dislocationdone:
	       
	       {} 	// do nothing

	    }// End: if task == dislocation

	    {}	 	// do nothing

	    //if (!(timestep%10000)) {
	       //radial("sample");		// sample radial distribution function
	       //printf("do sq4 sample\n");
	       //fflush(stdout);
	       //sq4("sample");
	    //}

	 } // if (MPI_myid == 0)

goto	output_conf;		// most of the analysis is done through LAMMPS functions

	 // CoorSI2System();		// Convert {x, y, z}'s, lx, ly, lz and lbox

	 // Calculate other length scales
	 for (i=0; i<NSYSTEMS; i++) {	
	    BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
	    BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
	    BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
	    BOX[i].rb		=	Rb;
	    BOX[i].rv		=	Rv;
	 }
	 // radial("sample");		// sample radial distribution function

	 //-----------------------------//
	 //	Build cell list		//
	 //-----------------------------//
#ifdef CELL_LIST
	 if (M2HDEBUG)	printf("DEBUG: build cell list ...\n");

	 if (clistinit) {
	    CL_Init();
	    CL_Build();
	    clistinit = 0;
	 }
	 /*
	 if (!clistinit) {
	    CL_Destroy();
	 }
	 clistinit = 0;

	 CL_Init();		// need to init. cell list every time due to volume change
	 CL_Build();
         */
#endif	/* CELL_LIST */

	 //-----------------------------//
	 //	Perform analysis	// 
	 //-----------------------------//

	 if (M2HDEBUG)	printf("DEBUG: perform analysis ...\n");

	 //------ Convert to Z-code input file ------//
	 
	 //------ Calculate common neighborhood parameter (CNP) ------//

	 if (M2HDEBUG)	printf("DEBUG: atom_neighbor() ...\n");
	 //atom_neighbor();					// find nearest neighbors

	 if (M2HDEBUG)	printf("DEBUG: comm_neigh_para() ...\n");
	 //comm_neigh_para();					// calculate CNP
	 
	 //------ Print out some results ------//

	 printf("lmp_atom_id\tpe\tcna\tcnp\tnneigh\n");
         for (moli=mol; moli<mol+NMOLS; moli++) {
	    for (i=0; i<moli->nsites; i++) {
	       //if (moli->type[i]==3) {
	       if (displayatom(moli-mol)) {
		  printf("%d", i+1);
		  printf("\t%f ", moli->pe[i]);
	          //printf("%f ", moli->rmin[i]);
		  //printf("%f ", moli->rmax[i]);
		  //printf("%f ", moli->rave[i]);
		  printf("\t%d", moli->cna[i]);
		  printf("\t%f", moli->cnp[i]);
		  printf("\t%d", moli->nneigh[i]);
		  printf("\n");
	       }
	    }
	 }
	 fflush(stdout);

	 //--- shift coordinates so that the com of normal atoms is at the origin
	 //shift_coord();
	
	 //--- Output the position of each column of atoms
	 //monitor_coord();

	 //------ Calculate hyper-distance ------//

	 if (!confinit ) {
	    
	    if (lammpsflag) {
	       // calculate hyper distance between MD and minimized configuration
	       m	=	0;
	       n	=	0;
	       for (moli=mol; moli<mol+NMOLS; moli++) {
		  for (i=0; i<moli->nsites; i++) {

		     if (moli->cnp[i] > 24) {
			r_ref[n]	=	r_MIN[m];
			r_ref[n+1]	=	r_MIN[m+1];
			r_ref[n+2]	=	r_MIN[m+2];

			r_trj[n]	=	r_MD[m];
			r_trj[n+1]	=	r_MD[m+1];
			r_trj[n+2]	=	r_MD[m+2];
			
			n	+=	3;
		     }
		     m	+=	3;
		  }
	       }
	       printf("MD-MIN hyper-distance = %f\n", hyper_distance(n, r_ref, r_trj));
	       printf("MD-MIN slip-distance = %f\n", slip_distance(n, r_ref, r_trj));
            }

	    
	    n = 0;
	    for (moli=oldmol; moli<oldmol+NMOLS; moli++) {	// previous conf.
	       for (i=0; i<moli->nsites; i++) {
		  if (moli->type[i]==3) {
		     r_ref[n]	=	moli->p[i].x;
		     r_ref[n+1]	=	moli->p[i].y;
		     r_ref[n+2]	=	moli->p[i].z;
		     n	+=	3;
		  }
	       }
	    }
	    

	    n = 0;
	    for (moli=mol; moli<mol+NMOLS; moli++) {		// current conf.
	       for (i=0; i<moli->nsites; i++) {
		  if (moli->type[i]==3) {
		     r_trj[n]	=	moli->p[i].x;
		     r_trj[n+1]	=	moli->p[i].y;
		     r_trj[n+2]	=	moli->p[i].z;
		     n	+=	3;
		  }
	       }
	    }

	    printf("incremental hyper-distance = %f\n", hyper_distance(n, r_ref, r_trj));
	    printf("incremental slip-distance = %f\n", slip_distance(n, r_ref, r_trj));
	    fflush(stdout);

	    /*
	    for (i=0; i<NSITES; i++) {
	       r_trjold[i]	=	r_trj[i];
	       r_refold[i]	=	r_ref[i];
	    }
	    */
	 }

	 //------ Identify and group atoms in dislocations/partials ------//

	 find_nuclei_general("cnp", "dummy", "r", 3.8);

         float	x1=0, x2=0;
	 int	n1=0, n2=0;
	 for (moli=mol; moli<mol+NMOLS; moli++) {
	    for (i=0; i<moli->nsites; i++) {
	       //if (moli->type[i]==0) {
	       if (moli->type[i]==3) {
		  if (moli->nuclid[i]==1) {
		     x1	+=	moli->p[i].x;
		     n1	++;
		  }
		  else if (moli->nuclid[i]==2) {
		     x2	+=	moli->p[i].x;
		     n2	++;
		  }
	       }
	    }
	 }
	 if (x1/n1 < x2/n2) {
	    printf("dislocation %d %d %f %d %f %d\n", timestep, n1, x1/n1, n2, x2/n2, Nnucl[0]);
	 }
	 else {
	    printf("dislocation %d %d %f %d %f %d\n", timestep, n2, x2/n2, n1, x1/n1, Nnucl[0]);
	 }
	 fflush(stdout);

	 //-------------------------------------//
	 //	Output configuration files	//
	 //-------------------------------------//
output_conf:

         if (MPI_myid==0) {
	    printf("\n    # Output configuration files ...\n");
	    fflush(stdout);

	    if (xyzflag) {
	       output_xyzfile(fxyz, fcolor);
	    }
	    if (carflag) {
	       output_carfile(fcar);
	    }
	    if (pdbflag) {
	       output_pdbfile(fpdb);
	    }
	 }

	 //-------------------------------------//
	 //	Finish analysis and clean-up	//
	 //-------------------------------------//
	 if (MPI_myid==0) {
	    for (i=0; i<NMOLS; i++) {
	       oldmol[i]	=	mol[i];		// save molecules information
	    }

	    if (confinit) {				// at least one conf. read
	       confinit	=	0;
	    }
	 }

	 //}			// if line contains string "TIMESTEP"
	 if (M2HDEBUG) 
	    printf("\n# Proc %d finished one conf. \n", MPI_myid);

      } 			// while (file_exist && !eof) for ONE input dump file
      if (M2HDEBUG) 
	 if (file_exist)
	    printf("\n# Proc %d finished file %s. \n", MPI_myid, infile[ifile]);

      if (MPI_myid==0)
	 if (file_exist) 
	    fclose(fin);	// close current input dump file
   }				// finish ALL input dump files

   //-------------------------------------------------------------------//
   //	Output final analysis results after ALL frames processed	//
   //-------------------------------------------------------------------//
   if (MPI_myid==0)	{
      //radial("print");		// print out radial distribution function
      //printf("print sq4 results\n");
      //fflush(stdout);
      //sq4("print");
   }

   if (MPI_myid==0)	{
      if (M2HDEBUG)	printf("# DEBUG: ALL FRAMES PROCESSED ...\n");

      //------ Frame information ------//

      printf("\n# OUTPUT AFTER ALL FRAMES ARE PROCESSED\n\n");
      printf("START TIMESTEP\t%d\n", starttime);
      printf("END TIMESTEP\t%d\n", endtime);
      printf("TOTAL FRAMES\t%d\n", nframe+1);
      printf("ANALYZE EVERY\t%d\n", dnframe);
      printf("ANALYZED FRAMES\t%d\n", nframe/dnframe+1);
   }

   //-------------------//
   //	Closing files	//
   //-------------------//
   if (MPI_myid==0)	{
      if (M2HDEBUG)	printf("# DEBUG: Close output files ...\n");

      prog_end	=	time(NULL);
      
      printf("\n# END OF OUTPUT FILE\n");
      printf("%s", asctime(localtime(&prog_end)));
      printf("\nTotal running time was %lf seconds.\n", difftime(prog_end, prog_start));

      if (xyzflag) {
	 fclose(fxyz);
	 fclose(fcolor);
      }
      if (carflag) {
	 fclose(fcar);
      }
      if (pdbflag) {
	 fclose(fpdb);
      }
   }
   fflush(stdout);

   if (lammpsflag) {
#ifdef myMPI
      if (lmp_comm==1) {
	 lammps_close(lmp);			// close LAMMPS instance
      }
#endif
   }

#ifdef myMPI
   MPI_Barrier(comm_world);
   /*
   if (lmp_comm==1) {
      MPI_Comm_free(&comm_lammps);
   }
   MPI_Barrier(comm_world);
   MPI_Comm_free(&comm_world);
   */
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

//======================================================================//
//	output_pre(): write pre-analysis results for further analysis	//
//	              These pre-results should be basic but very	//
//		      time-consuming.  So by saving these results in	//
//		      a pre- file we can reuse them later without 	//
//		      calculating them every time.			//
//======================================================================//

//==============================================================//
//	read_pre(): read pre-analysis results, see output_pre()	//
//==============================================================//

//==============================================================//
//	output_xyzfile(): write configuration to .xyz file	//
//			  also output a scalar value for each	//
//			  atom for color display in VMD		//
//==============================================================//
void output_xyzfile(FILE *fPtr, FILE *fclr)
{
   molstruct	*moli;
   vector	p;
   char		element;
   int		i, n, N;
   int		system=0;
   static int	init;
   static float	*color_par;
   float	tmp1, tmp2, tmp3;

   if (init) {
      init	=	0;
      color_par	=	(float *)calloc(NMOLS, sizeof(float));
   }

   //--- count the number of atoms to display
   N	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      n	=	moli-mol;
      //if (moli->type[0]>=3 && displayatom(n)) {
      if (displayatom(n) ) {
	 N	++;
      }
   }

   //--- print the header of each configuration
   fprintf(fPtr, "%d\n", N);
   fprintf(fPtr, "timestep %d.\n", timestep);

   //--- print the xyz coordinates of atoms to display
   for (moli=mol; moli<mol+NMOLS; moli++) {

      n	=	moli - mol;

      //if (moli->type[0]>=3 && displayatom(n)) {
      if (displayatom(n)) {
	 for (i=0; i<moli->nsites; i++) {
	    //p	=	MapInBox2(moli->p+i, PBC, system);
	    //fprintf(fPtr, "O %f %f %f\n", p.x, p.y, p.z);
	    //fprintf(fPtr, "O %f %f %f\n", moli->p[i].x - mol->p[0].x, moli->p[i].y-mol->p[0].y, moli->p[i].z-mol->p[0].z);
	    //if (moli->cna[i] > 2.05 || moli->cna[i] < 1.95) { 

	    if (moli->type[i] == 5) {	// solute atoms
	       element = 'C';
	    }
	    else {
	       if (moli->cna[i] > 4.5) {	// cna=5, partial
		  element = 'O';	// red
	       }
	       else if (moli->cna[i]>0.5 && moli->cna[i]<1.5) {	// cna=1, fcc
		  element = 'N';	// blue
	       }
	       else {
		  element = 'H';
	       }
	    }

	    fprintf(fPtr, "%c %f %f %f\n", element, moli->p[i].x, moli->p[i].y, moli->p[i].z);
	 }
      }
   }
   fflush(fPtr);

   //--- calculate the coloring parameter

   //--- find min and max for coloring parameter
   tmp1	=	1.0e8;
   tmp2	=	-1.0e8;

   for (moli=mol; moli<mol+NMOLS; moli++) {
      n	=	moli-mol;
      //if (moli->type[0]==0 && displayatom(moli-mol)) {
      //if (moli->type[0]>=3 && displayatom(moli-mol)) {
      if (displayatom(n)) {
	 tmp3	=	moli->nneigh[0];
	 tmp3	=	moli->pe[0];
	 tmp3	=	moli->cna[0];
	 tmp3	=	moli->rmax[0]-moli->rmin[0];
	 tmp3	=	moli->rmax[0];
	 tmp3	=	moli->rmin[0];
	 tmp3	=	moli->cnp[0];

	 if (tmp1 > tmp3) {
	    tmp1 = tmp3;		// min
	 }
	 if (tmp2 < tmp3) {
	    tmp2 = tmp3;		// max
	 }
      }
   }

   // print the coloring parameter for display in VMD

   for (moli=mol; moli<mol+NMOLS; moli++) {
      tmp3	=	moli->nneigh[0];
      tmp3	=	moli->pe[0];
      tmp3	=	moli->cna[0];
      tmp3	=	moli->rmax[0]-moli->rmin[0];
      tmp3	=	moli->rmax[0];
      tmp3	=	moli->rmin[0];
      tmp3	=	moli->cnp[0];

      n	=	moli-mol;
      //if (moli->type[0]==0 && displayatom(moli-mol)) {
      //if (moli->type[0]>=3 && displayatom(moli-mol)) {
      if (displayatom(n)) {
	 tmp3	=	0.7 * (tmp3 - tmp1)/(tmp2-tmp1);
	 //tmp3	=	0.7/35*(tmp3-5);		//cnp
      }

      //if (moli->type[0]!=0)	tmp3 = 0.8;
      //if (moli->type[0]<3)	tmp3 = 0.8;

      //if (displayatom(moli-mol)) {
      //if (moli->type[0]>=3&&displayatom(moli-mol)) {

      if (moli->cna[0] > 2.5 || moli->cna[0]<1.5) { 
	 tmp3	=	0.6;
      }
      else {
	 tmp3	=	0.1;
      }

      if (displayatom(n)) {
	 //fprintf(fclr, "%f ", tmp3);
      }
      //printf("pe %f\n", tmp3);
      //printf("cnp %f\n", tmp3);
   }
   fprintf(fclr, "\n");
   fflush(fclr);
   return;
}

//==============================================//
//	output_dump(): write LAMMPS dump file	//
//==============================================//
void output_dump(FILE *fPtr)
{
   molstruct	*moli;
   int		i, n;

   n  = 0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 n ++;
	 printf("%d %d %f %f %f %d %d\n", 
	        n, moli->type[i]+1, moli->p[i].x, moli->p[i].y, moli->p[i].z,
		moli->nuclid[i], moli->cna[i]);
      }
   }
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

         //MolInBox2(moli);
         for (i=0; i<moli->nsites; i++) {
            //moli->p[i] = MapInBox2(moli->p+i, PBC, system); //temporary

            if (moli->cna[i]>1.5)		
               sprintf(s, "O%d", n++);			// O: default red color in VMD
            else
               sprintf(s, "N%d", n++);			// N: default blue color in VMD

	    sprintf(s, "N");			// N: default blue color in VMD

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
            fprintf(fPtr, "%-4.4s%-6ld ND       C 0.000\n", ff, moli-mol);
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
void output_pdbfile(FILE *fPtr) 
{
   time_t	t	=	time(NULL);
   molstruct	*moli;
   char		atomname;
   int		system = 0;			// only one box for now (11/18/2011)
   int		i, m, n;
   int		drawmol;
   int		nuclid;

   fprintf(fPtr, "HEADER: file created on %s", asctime(localtime(&t)));
   fprintf(fPtr, "HEADER: timestep %d\n", timestep);
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

         m	++; 
         for (i=0; i<moli->nsites; i++) {

	    /*
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
	    */
	    //if ( moli->nneigh[i]>10 && (moli->cna[i]>2.5 || moli->cna[i]<1.5) ) {
	    if ( moli->nneigh[i]>10 && (moli->cna[i]>1.5) ) {
	       atomname	=	'O';
	    }
	    else {
	       atomname	=	'N';
	    }

            n	++;
	    fprintf(fPtr, "ATOM  ");			// pdb command, column 1-6
            fprintf(fPtr, "%6d" , n);			// atom number, column 7-12

            fprintf(fPtr, "    ");		// atom name, column 13-16
            //fprintf(fPtr, " %c  ", atomname);		// atom name

            fprintf(fPtr, " ");				// alternate location indiator

  	    fprintf(fPtr, "R00");			// residue name
  	    //fprintf(fPtr, "ILE");			// residue name

	    fprintf(fPtr, " ");				// column 21
            fprintf(fPtr, " ");				// chain identifier, column 22
            //fprintf(fPtr, "A");				// chain identifier, column 22
	    //fprintf(fPtr, "%4d", m);			// residue sequence number, 23-26
	    fprintf(fPtr, "%4d", 2);			// residue sequence number, 23-26
	    fprintf(fPtr, " ");				// code for insertion of residues, 27
            fprintf(fPtr, "   ");			// column 28-30
            fprintf(fPtr, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z); // col. 31-54
            fprintf(fPtr, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature	// col. 55-66
	    fprintf(fPtr, "      ");			// col. 67-72
	    fprintf(fPtr, "PROT");			// col. 73-76
            //fprintf(fPtr, "%5.5s      %c", "", atomname);				// col. 73-76
            fprintf(fPtr, "\n"); 

	    /*
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
	    */
         } 
	 //}
      }   
   }
   fprintf(fPtr, "END\n");
   fflush(fPtr);
   return;
}

//======================================================================//
//	output_conf(): write configuration to configuration file	//
//======================================================================//

//======================================================================//
//	output_Zinput(): write coord. as Z code input file (4/4/13)	//
//======================================================================//

int displayatom(int n)
{
   int	result;

   // for basal slip plane
   if ( (n<9600 && n/4/40%6==5 && (n%4==2 || n%4==3)) || 
        (n>=9600 && (n-9600)/4/41%6==0 && (n%4==0 || n%4==1)) ){ 

   // for prismatic slip plane
//   if ( (n<12000 && n/4/50%15==14 && (n%4==1 || n%4==3)) || 
//        (n>=12000 && (n-12000)/4/51%15==0 && (n%4==0 || n%4==2)) ){ 

      result	=	1;
   }
   else {
      result	=	0;
   }

   //--- show all atoms
   result	=	1;	

   return	result;
}


//==============================================================================//
//	shift_coord(): shift whole system w.r.t. the com of normal atoms	//
//			added 3/27/14						//
//==============================================================================//
void shift_coord()
{
   molstruct	*moli;
   int		i;
   int		n;
   int		system = 0;
   float	min_cnp=6.8, max_cnp=6.9;
   vector	com;

   // find the center of mass of all normal atoms
   V_Null(&com);

   /*
   n	=	0;
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 if (moli->cnp[i] > min_cnp && moli->cnp[i] < max_cnp) {	// normal atoms
	    com.x	+=	moli->p[i].x;
	    com.y	+=	moli->p[i].y;
	    com.z	+=	moli->p[i].z;
	    n	++;
	 }
      }
   }
   com	=	V_Mult(1.0/n, &com);
   */
   n	=	16000;
   for (i=0; i<4; i++) {
      moli	=	mol+n+i;
      com.x	+=	moli->p[0].x;	
      com.y	+=	moli->p[0].y;	
      com.z	+=	moli->p[0].z;	
   }
   com	=	V_Mult(1.0/4, &com);

   // shift coordinates
   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 moli->p[i].x	-=	com.x;
	 moli->p[i].y	-=	com.y;
	 moli->p[i].z	-=	com.z;
      }
   }
   printf("shift_coord(): %f %f %f\n", com.x, com.y, com.z);
   return;
}


void monitor_coord()
{
   molstruct	*moli;
   int		Ncol = 402;
   int		i, n;
   int		result;
   float	ave_cnp = 7.0;
   float	x[402];
   float	y[402];
   float	z[402];
   float	cnp[402];
   float	xdis[402];
   int		nn[402];
   int		half = 16000;
   vector	r;
   int		system = 0;

   // initialization
   for (i=0; i<Ncol; i++) {
      x[i]	=	0.0;
      cnp[i]	=	0.0;
      xdis[i]	=	0.0;
      y[i]	=	0.0;
      z[i]	=	0.0;
      nn[i]	=	0;
   }

   //shift_coord();

   for (moli=mol; moli<mol+NMOLS; moli++) {
      n	=	moli-mol;

      if (moli->type[0]==3) {			// dislocation layer

	 if (n < half && n%4==3) {
	    i	=	(n/4%100)*2;		// 0, 2, 4, ..., 198

	    x[i]	+=	moli->p[0].x;
	    y[i]	+=	moli->p[0].y;
	    z[i]	+=	moli->p[0].z;
	    cnp[i]	+=	moli->cnp[0];
	    nn[i]	++;
	 }
	 //else if (n < half && n%4==2) {	// basal
	 else if (n < half && n%4==1) {		// prismatic
	    i	=	(n/4%100)*2 + 1;	// 1, 3, 5, ..., 199

	    x[i]	+=	moli->p[0].x;
	    y[i]	+=	moli->p[0].y;
	    z[i]	+=	moli->p[0].z;
	    cnp[i]	+=	moli->cnp[0];
	    nn[i]	++;
	 }
	 else if (n >= half && n%4==0) {
	    i	=	((n-half)/4%101)*2 + 200;	// 200, 202, ..., 400

	    x[i]	+=	moli->p[0].x;
	    y[i]	+=	moli->p[0].y;
	    z[i]	+=	moli->p[0].z;
	    cnp[i]	+=	moli->cnp[0];
	    nn[i]	++;
	 }
	 //else if (n >= half && n%4==1) {	// basal
	 else if (n >= half && n%4==2) {	// prismatic
	    i	=	((n-half)/4%101)*2 + 201;	// 201, 203, ..., 401

	    x[i]	+=	moli->p[0].x;
	    y[i]	+=	moli->p[0].y;
	    z[i]	+=	moli->p[0].z;
	    cnp[i]	+=	moli->cnp[0];
	    nn[i]	++;
	 }
      }
      /*
      if (displayatom(n) && n<9600 && n%4==2) {
	 i	=	 n/4%40;

	 x[i]	+=	moli->p[0].x;
	 cnp[i]	+=	moli->cnp[0];
	 nn[i]	++;
      }
      */
   }
    
   // normalize
   for (i=0; i<Ncol; i++) {
      x[i]	/=	nn[i];
      y[i]	/=	nn[i];
      z[i]	/=	nn[i];
      cnp[i]	/=	nn[i];
   }


   printf("x position %ld", timestep);
   for (i=0; i<Ncol; i++) {
      printf(" %3d %5.3f %5.3f", i, x[i], cnp[i]);
   }
   printf("\n");
   fflush(stdout);

   for (i=0; i<Ncol; i++) {
      if (i<2 || (i>=(Ncol-2)/2 && i<(Ncol-2)/2+2)) {
	 xdis[i]	=	3.2;
      }
      else {
	 xdis[i]	=	x[i] - x[i-2];
      }
   }
   for (i=0; i<(Ncol-2)/2; i+=2){
      printf("%d %f %f %f %f %f\n", i, x[i], y[i], z[i], cnp[i], xdis[i]);
   }
   for (i=1; i<(Ncol-2)/2; i+=2){
      printf("%d %f %f %f %f %f\n", i, x[i], y[i], z[i], cnp[i], xdis[i]);
   }
   for (i=(Ncol-2)/2; i<Ncol; i+=2){
      printf("%d %f %f %f %f %f\n", i, x[i], y[i], z[i], cnp[i], xdis[i]);
   }
   for (i=(Ncol-2)/2+1; i<Ncol; i+=2){
      printf("%d %f %f %f %f %f\n", i, x[i], y[i], z[i], cnp[i], xdis[i]);
   }
   fflush(stdout);


   return;
}

void binning(char *var, char *dim)
{
   static int init = 1;
   static int *zbins, *xbins;
   static FILE *fzprof, *fxprof;
   molstruct *moli;
   int i, j, k;
   int system = 0;
   float zbinsize = (BOX[system].zhi - BOX[system].zlo)/50;
   float xbinsize = (BOX[system].xhi - BOX[system].xlo)/100;
 
   if (init) {
      zbins = (int *) calloc(50, sizeof(int));	
      xbins = (int *) calloc(100, sizeof(int));

      if (zbins == NULL || xbins == NULL)
         Exit("position", "binning", "out of memory");

      fzprof = fopen("zprofile","w");
      fxprof = fopen("xprofile","w");
      init = 0;
   }
   for (i=0; i<50; i++) {
      zbins[i] = 0;
   }
   for (i=0; i<100; i++) {
      xbins[i] = 0;
   }

   for (moli=mol; moli<mol+NMOLS; moli++) {
      for (i=0; i<moli->nsites; i++) {
	 if (moli->cna2[i] == 3) {  // bcc, cna2
	    j = (moli->p[i].z - BOX[system].zlo)/zbinsize;
	    zbins[j] ++;

	    if (j>=24 && j<=26) {
	       k = (moli->p[i].x - BOX[system].xlo)/xbinsize;
	       xbins[k] ++;
	    }
	 }
      }
   }
   fprintf(fzprof, "timestep = %d\t zbinsize = %f\n", timestep, zbinsize);
   fprintf(fzprof, "index\t z\t count(bcc atoms)\n");
   for (i=0; i<50; i++) {
      fprintf(fzprof, "%d\t %f\t %d\n", i, (i+0.5)*zbinsize, zbins[i]);
   }
   fprintf(fzprof, "\n");
   fflush(fzprof);

   fprintf(fxprof, "timestep = %d\t xbinsize = %f\n", timestep, xbinsize);
   fprintf(fxprof, "index\t x\t count(bcc atoms)\n");
   for (i=0; i<100; i++) {
      fprintf(fxprof, "%d\t %f\t %d\n", i, (i+0.5)*xbinsize, xbins[i]);
   }
   fprintf(fxprof, "\n");
   fflush(fxprof);
}
