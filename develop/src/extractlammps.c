/*
	program:	extractlammps.c
	author:		Peng Yi at MIT
	date:		May 18, 2011
	purpose:	read lammps dump file, extract the configuration of one chain
			from different frames and output one single vmd-readable file
	note:		
			require setup file

			2011/5/18 the box dimension output is not exactly correct, 
			the box size of the last snapshot is output in .car file.
*/

#define __MAIN_PROGRAM
#include "header.h"

#define CAP	5		// largest nuclei considered still the melt
#define p2threshold	0.4
#define DEBUG	0

#include "correlation.h"

long		timestep;
double		rshift2;	// shift of the biggest nucleus
long		nsitesp2;	// # of sites with p2 greater than a threshold
long		nmolsp2;
long		nmaxp2_1[10];	// nmax using p2 nucleus definition
long		nmaxp2_2[10], nmaxp2_3[10], nmaxp2_4[10];

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
} 

int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   molstruct	*moli;
   molstruct	*moln;
   long		i, j, k, system, m, n, nuclid, nsites;
   long		previous_timestep=-1;
   long		nx, ny, nz, id, siteid, molid, type;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo;			// box upper boundary
   double	temp1, temp2, temp3;		// dummy variables

   char		infile[8][255], filein[80], filename[80];
   char		s[80], ff[80], par[80], dummy[255];
   FILE		*fin, *fhst, *fout, *fconf, *fpdb, *fdat;
   long		LENGTH, accum, confflag=0, carflag=0, pdbflag=0, polydisperse=0, drawmol;
   long		nframe, dnframe=1;		// analyze only every dnframe
   long		nfiles=1, ifile;		// number of input files 
   long		nshots=1, ishots=0;		// number of snapshots
   long		chainid=0;
   char		atomname;
   static long	init=1;	

   vector		con;	// center of nucleus
   static vector	rO;	// original position of nucleus

   // chain rotational angle distribution
   vector	chainface;
   double	orient;
   long		orientdist[180];

   // segment statistics variables
   long		previous, previous_id;
   long		nlength[MAXNMOLSITES], n_sep[MAXNMOLSITES];
   long		head, tail, seg_on, seg_id, nseg, nxtal, length;
   long		segment[MAXNMOLS][MAXNMOLSITES];	// segments identification on a chain
   long		nsegment[MAXNMOLS];			// # of segments on a chain
   long		seg_stat[MAXNMOLS][MAXNMOLSITES];	// xtal segment stat for chains
   long		loose_stat[MAXNMOLS][MAXNMOLSITES];	// loose segment stat
   long		nloop, nbridge, ntail, nxseg;		// # of total loops, etc
   double	lloop, lbridge, ltail, lxseg;		// average length of loops, etc
   long		nlloop[MAXNMOLSITES], nlbridge[MAXNMOLSITES];	// segment length distribution
   long		nltail[MAXNMOLSITES], nlxseg[MAXNMOLSITES];

   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];		// group beads in the same nucleus

   vector	com[MAXNMOLS], temp;		// average center of mass
   long		ncom[MAXNMOLS];

   double	imagex, imagey, imagez;		// variables for creating pbc images
   long		imagen;

   if (argc<2) {
      printf("extractlammps (c) 2011 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\textractlammps [-option] [nshots= 1] [dn= 2] chainid= 0 lammpsdumpfile\n\n");

      printf("Notes:\n");
      printf("\t* -option = \n");
      printf("\t* -poly: polydisperse system, has molid in dump file\n");
      printf("\t* -conf: configuration file output\n");
      printf("\t* -car: car file output\n");
      printf("\t* -pdb: pdb file output\n");
      printf("\t* nshots=: number of snapshots\n");
      printf("\t* dn=: only analyze every dn frames\n");
      printf("\t* chainid: id of the chain we extract, id starts from 0\n");
      printf("\t* nfiles=: number of input files if more than 1 (must be <=8)\n");
      printf("\t* \"=\" must immediately follow x or y or z or n or dn\n");
      printf("\t* require setup file\n\n");
      exit(1);
   }

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-poly"))	polydisperse	=	1;
      else if (samestr(par, "-conf"))	confflag	=	1;
      else if (samestr(par, "-car"))	carflag		=	1;
      else if (samestr(par, "-pdb"))	pdbflag		=	1;
      else if (samestr(par, "nshots="))	nshots		=	atol(argv[i+1]);
      else if (samestr(par, "dn=")) 	dnframe		=	atol(argv[i+1]);
      else if (samestr(par, "nfiles="))	nfiles		=	atol(argv[i+1]);
      else if (samestr(par, "chainid="))chainid		=	atol(argv[i+1]);
      
   }
   for (i=0; i<nfiles; i++) {
      strcpy(infile[nfiles-1-i], argv[argc-1-i]);	// get input filenames
   }

   // Open output files
 
   if (nfiles==1)	strcpy(filein, infile[0]);
   else			strcpy(filein, "multi");

   strcpy(filename, filein);
   strcat(filename, ".car");
   if (carflag && (fout=fopen(filename, "w"))==NULL )
      Exit("lammps2hst", "main", "open car file failed.");

   strcpy(filename, filein);
   strcat(filename, ".pdb");
   if (pdbflag && ((fpdb=fopen(filename, "w"))==NULL || (fdat=fopen("vmd.dat","w"))==NULL))
      Exit("lammps2hst", "main", "open pdb file failed.");

   if (!DEBUG) {
      strcpy(filename, filein);
      strcat(filename, ".out");
      freopen(filename, "w", stdout);	// redirect standard output stream to a file
   }

   printf("nshots = %d dn = %d chainid = %d infile[0] = %s\n", 
		nshots, dnframe, chainid, infile[0]);

   ////////////////////
   // Initialization //
   ////////////////////

   if (DEBUG)	printf("Initialization: start ... \n");

   moln	=	(molstruct *) calloc(nshots, sizeof(molstruct));
  
   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   InitUnits();				// initialize units for calculation
   InitForcefield();			// initialize Lennard-Jones potential mixing rule

   system	=	0;		// for now 2/14/08, only one system
   					// initialize sampling
					// initialize chain rotation distribution
				   	// initialize segment variables
 
   nframe	=	-1;

for (ifile=0; ifile<nfiles; ifile++){		// multiple input files

   fin	=	fopen(infile[ifile], "r");

   while (!feof(fin)) {

      if (DEBUG)	printf("Read one configuration: start ...\n");

      /* Read in one configuration from dump file */

      if (!fgets(dummy, sizeof(dummy), fin))	break;			// end of file
      fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &xlo, &xhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &ylo, &yhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &zlo, &zhi);	fgets(dummy, sizeof(dummy), fin);
      fgets(dummy, sizeof(dummy), fin);

      BOX[system].lx	=	xhi-xlo;
      BOX[system].ly	=	yhi-ylo;
      BOX[system].lz	=	zhi-zlo;

      LENGTH		=	NSITES/NMOLS;	// monodisperse system for now (4/26/2008)
      accum		=	0;

      for (i=0; i<nsites; i++) {
         fscanf(fin, "%ld", &id);
         fscanf(fin, "%ld", &molid);		// Need to check the lammps.dump file format
						// because some early lammps.dump file 
						// does not have molid output
//         if (polydisperse)	fscanf(fin, "%ld", &molid);
         fscanf(fin, "%ld", &type);
         fscanf(fin, "%lf%lf%lf %lf%lf%lf %ld%ld%ld", &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
         fgets(dummy, sizeof(dummy), fin);

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
         mol[molid].box		=	system;		// for now, only one system
         mol[molid].nsites	=	LENGTH;		// for now, Jan/10/2010
         mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx);
         mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly);
         mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);
         mol[molid].type[siteid]=	type - 1;	// -1 because lammps index starts from 1
      }

      for (system=0; system<NSYSTEMS; system++) {
         NMols[system]	=	0;
         NSites[system]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;			// total # of mols in certain system
            NSites[system]	+=	moli->nsites;	// total # of sites in certain system
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

      ///////////////////////////////////////////////////
      /* Skip repeated frames in different input files */
      ///////////////////////////////////////////////////
      if (ifile>=1 && timestep==previous_timestep) {
         continue;
      }
      previous_timestep	=	timestep;	// check repeat frames in different infiles

      ///////////////////////////
      /* Skip dnframe-1 frames */
      ///////////////////////////

      nframe	++;

      if (mod(nframe, dnframe))	continue;	// analyze every dnframe frames
      
      //////////////////////////
      // Extract the molecule //
      //////////////////////////

      if (DEBUG)	printf("Extract the molecule ...\n");

      moln[ishots]	=	mol[chainid];
      MolInBox2(moln+ishots);
      ishots	++;

      if (ishots >= nshots)	break;

      ///////////////////////////////////////////////////////////////////////////
      // Start: Calculate the average position of center of mass of each chain //
      ///////////////////////////////////////////////////////////////////////////

      /////////////////////////////////////////////////////////////
      /* Convert coordinates and box size from SI to system unit */
      /////////////////////////////////////////////////////////////
      CoorSI2System();
      //_________________________________________________________//

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
      } 

      /////////////////////
      // build cell list //
      /////////////////////

      //////////////////////
      // Perform analysis //
      //////////////////////

      /////////////////////////////
      // Output analysis results //
      /////////////////////////////

      CoorSystem2SI();		// convert coordinates and box size back to SI units

   } 

   fclose(fin);		// close current input file
  }			// multiple input files

   ////////////////////////////////////////////
   // OUTPUT .car file for VMD visualization //
   ////////////////////////////////////////////
      
   if (carflag) {
      fprintf(fout, "!BIOSYM archive 3\n");
      fprintf(fout, "PBC=ON\n");
      fprintf(fout, "!TIMESTEP %d\n", timestep);
      fprintf(fout, "!DATE %s", asctime(localtime(&t)));
      fprintf(fout, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

      n	=	0;

      for (moli=moln; moli<moln+nshots; moli++) {
         if (system==moli->box) {

            //MolInBox2(moli);
            for (i=0; i<moli->nsites; i++) {

               if (moli->nuclid[i]>0)	// crystal like particle
                  sprintf(s, "N%d", n++);	// N: blue color in VMD
               else
                  sprintf(s, "O%d", n++);	// O: red color in VMD

               fprintf(fout, "%-5.5s ", s);
               sprintf(s, "M%d", moli-moln);
               fprintf(fout, "%14.8g %14.8g %14.8g ", moli->p[i].x, moli->p[i].y, moli->p[i].z);
               strcpy(ff, "O");
               fprintf(fout, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-moln, Element(moli->type[i], s));
            } 
         }   
      }
      fprintf(fout, "end\nend\n");
   }
   fflush(fout);
   //___________OUTPUT .car file_____________//
/*
   ///////////////////////////////////////////
   // OUTPUT .pdb file for further analysis //
   ///////////////////////////////////////////

   if (pdbflag) {
      fprintf(fpdb, "HEADER: file created from %s on %s", argv[1], asctime(localtime(&t)));
      fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
        	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

      m		=	0;	// molecule sequence number
      n		=	0;	// atom sequence number
#define SIZECAP	5
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if (system==moli->box) {
            //MolInBox2(moli);

            drawmol	=	0;
            for (i=0; i<moli->nsites; i++) {
               if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {
                  drawmol	=	1;		// participate the biggest nucleus
                  break;
               }
	    }
//temp=CenterofMass(moli);
//if (temp.z < -0.25 * BOX[system].lz) {
            m	++; 
            for (i=0; i<moli->nsites; i++) {
               if (drawmol) {
                  if (sizeofnucl[moli->nuclid[i]] == MAXSIZE[system]) {	// nuclid index starts from 1
                     atomname	=	'N';	// N: blue color in VMD
	             fprintf(fdat, " 10");
                  }
                  else {
                     atomname	=	'O';	// O: red color in VMD
                     fprintf(fdat, " 0");
                  }
               }
               else {
	          atomname	=	'C';		// C: cyan color in VMD
	          fprintf(fdat," -1");
	       }

               n	++;
	       fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
               fprintf(fpdb, "%5d ", n);	// atom number
               fprintf(fpdb, " %c  ", atomname);	// atom name
               fprintf(fpdb, " ");		// alternate location indiator
  	       fprintf(fpdb, " C8");		// residue name
	       fprintf(fpdb, " ");		// column 21
               fprintf(fpdb, " ");		// chain identifier, column 22
	       fprintf(fpdb, "%4d", m);	// residue sequence number, 23-26
	       fprintf(fpdb, " ");		// code for insertion of residues, 27
               fprintf(fpdb, "   ");		// column 28-30
               fprintf(fpdb, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z);
               fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
               fprintf(fpdb, "%5.5s", "");
               fprintf(fpdb, "\n"); 

               if (imagen) {			// for image box
                  n	++;
                  fprintf(fpdb, "ATOM  ");
                  fprintf(fpdb, "%5d ", n);
                  fprintf(fpdb, " %c  ", atomname);
                  fprintf(fpdb, " ");
                  fprintf(fpdb, " C8");
	          fprintf(fpdb, " ");
                  fprintf(fpdb, " ");
	          fprintf(fpdb, "%4d", m);
                  fprintf(fpdb, " ");
                  fprintf(fpdb, "   ");
                  fprintf(fpdb, "%8.3f%8.3f%8.3f", 
			moli->p[i].x + imagex, moli->p[i].y+imagey, moli->p[i].z+imagez);
                  fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                  fprintf(fpdb, "%5.5s", "");
                  fprintf(fpdb, "\n"); 
               }
            } 
//}
         }   
      }
      fprintf(fpdb, "END\n");
      fflush(fpdb);
   }	//pdbflag
   //____________OUTPUT .pdb file___________//
*/

   if (DEBUG)	printf("Closing output files ...\n");

   if (carflag)		fclose(fout);
   if (pdbflag)	     {	fclose(fpdb); fclose(fdat);}

   fflush(stdout);
   fclose(stdout);

   return	0;
}
