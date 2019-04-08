/*
	program:	period.c
	author:		Peng Yi at MIT
	date:		Aug 29, 2009
	purpose:	create periodic image 
	note:		
*/

#define __MAIN_PROGRAM
#include "header.h"

#define CAP	5		// largest nuclei considered still the melt
#define p2threshold	0.4

#include "correlation.h"

long		timestep;
double		rshift2;	// shift of the biggest nucleus
long		nsitesp2;	// # of sites with p2 greater than a threshold
long		nmolsp2;
long		nmaxp2[10];	// nmax using p2 nucleus definition

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


int main(int argc, char *argv[])
{
   time_t	t	=	time(NULL);
   long		i, j, k, system, n;
   double	x, y, z, 			// coordinates
		vx, vy, vz, 			// velocity
		fx, fy, fz,			// force
		xhi, xlo, yhi, 			// box lower boundary
		ylo, zhi, zlo;			// box upper boundary
   long		nsites, record;
   long		nx, ny, nz, id, siteid, molid, type;
   molstruct	*moli;

   char		filein[255], s[80], ff[80], par[80];
   FILE		*fin, *fhst, *fout, *fconf, *fpdb, *fdat;
   char		a[255], b[255], c[255], d[255], dummy[255];
   long		LENGTH, confflag=1, carflag=1, pdbflag=1, drawmol;
   char		atomname;
   static long	init=1;

   vector		con;	// center of nucleus
   static vector	rO;	// original position of nucleus

   vector		chainface;
   double		orient;
   long			orientdist[180];	// chain orientation distribution

   beadstruct	nucleus[MAXNMOLS*MAXNMOLSITES];

   vector	com[MAXNMOLS], temp;		// average center of mass
   long		ncom[MAXNMOLS];

   for (i=0; i<180; i++)
	orientdist[i]	=	0;

   if (argc<6) {
      printf("period (c) 2009 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tperiod x y z n dumpfile\n\n");
      printf("Notes:\n");
      printf("\t* x, y, z: the direction to create image\n");
      printf("\t* n: the number of images to create\n");
      exit(1);
   }

   shiftx	=	atof(argv[1]);
   shifty	=	atof(argv[2]);
   shiftz	=	atof(argv[3]);
   n		=	atol(argv[4]);

   for (i=1; i<argc-1; i++) {
      strcpy(par, argv[i]);
      if (samestr(par, "-noconf")) 
         confflag	=	0;
      else if (samestr(par, "-nocar"))
         carflag	=	0;
      else if (samestr(par, "-nopdb"))
         pdbflag	=	0;
   }
   strcpy(filein, argv[argc-1]);

   if ( (fin=fopen(filein, "r"))==NULL )
      Exit("lammps2hst", "main", "open input file failed.");
   if ((fhst=fopen("lammps.hst", "w"))==NULL )
      Exit("lammps2hst", "main", "open hst file failed.");
   if (confflag && (fconf=fopen("lammps.conf", "w"))==NULL )
      Exit("lammps2hst", "main", "open conf file failed.");   
   if (carflag && (fout=fopen("lammps.car", "w"))==NULL )
      Exit("lammps2hst", "main", "open car file failed.");
   if (pdbflag && ((fpdb=fopen("lammps.pdb", "w"))==NULL || (fdat=fopen("vmd.dat","w"))==NULL))
      Exit("lammps2hst", "main", "open pdb file failed.");

   InitMols(MAXNMOLS, MAXNMOLS);	// allocate memory for molecules
   GetSetup(argv);			// read in setup file
   InitUnits();				// initialize units for calculation
   InitForcefield();			// initialize Lennard-Jones potential mixing rule

   system	=	0;		// for now 2/14/08, only one system

   InitSample();			// initialize sampling

   while (!feof(fin)) {

      /* Read in one configuration from dump file */

      if (!fgets(a, sizeof(a), fin))		break;				// end of file
      fscanf(fin, "%ld", &timestep);		fgets(dummy, sizeof(dummy), fin);
      fgets(b, sizeof(b), fin);
      fscanf(fin, "%ld", &nsites);		fgets(dummy, sizeof(dummy), fin);
      fgets(c, sizeof(c), fin);
      fscanf(fin, "%lf%lf", &xlo, &xhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &ylo, &yhi); 	fgets(dummy, sizeof(dummy), fin);
      fscanf(fin, "%lf%lf", &zlo, &zhi);	fgets(dummy, sizeof(dummy), fin);
      fgets(d, sizeof(d), fin);

      BOX[system].lx	=	xhi-xlo;
      BOX[system].ly	=	yhi-ylo;
      BOX[system].lz	=	zhi-zlo;

      LENGTH		=	NSITES/NMOLS;	// only for monodisperse system (4/26/2008)

      for (i=0; i<nsites; i++) {
         fscanf(fin, "%ld%ld %lf%lf%lf %lf%lf%lf %ld%ld%ld",
			&id, &type, &x, &y, &z, &vx, &vy, &vz, &nx, &ny, &nz);
         fgets(dummy, sizeof(dummy), fin);

         id	--;				// because lammps starts from 1 not 0
         type	--;
         molid	=	(long) (id/LENGTH);
         siteid		=	id % LENGTH;
         mol[molid].box		=	system;		// for now, only one system
         mol[molid].nsites	=	LENGTH;
         mol[molid].p[siteid].x	=	x + nx*(BOX[system].lx);
         mol[molid].p[siteid].y	=	y + ny*(BOX[system].ly);
         mol[molid].p[siteid].z	=	z + nz*(BOX[system].lz);
         mol[molid].type[siteid]=	type;
      }

      // Start: Calculate the average position of center of mass
      for (moli=mol; moli<mol+NMOLS; moli++) {
      	  temp	=	CenterofMass(moli);
          if (fabs(temp.x) < 0.45 * BOX[moli->box].lx && 
		fabs(temp.y) < 0.45 * BOX[moli->box].ly &&
		fabs(temp.z) < 0.45 * BOX[moli->box].lz) {
	     com[moli-mol]	=	V_Add(&temp, com+(moli-mol));
             ncom[moli-mol]	++;
          }
      }
      // End: Calculate the average position of center of mass

      CoorSI2System();		// convert coordinates and box size from SI to system unit

      for (system=0; system<NSYSTEMS; system++) {
         NMols[system]	=	0;
         NSites[system]	=	0;
      }
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if ( (system=moli->box) >= 0) {
            NMols[system]	++;				// total # of mols in certain system
            NSites[system]	+=	moli->nsites;		// total # of sites in certain system
         }
      }

      for (i=0; i<NSYSTEMS; i++) {
         BOX[i].lbox		=	MIN(MIN(BOX[i].lx, BOX[i].ly), BOX[i].lz);
         BOX[i].vol		=	BOX[i].lx * BOX[i].ly * BOX[i].lz;
         BOX[i].rc		=	MIN(0.5*BOX[i].lbox, Rc);
         BOX[i].rb		=	Rb;
         BOX[i].rv		=	Rv;
      } 

      for (moli=mol; moli<mol+NMOLS; moli++) { 		
         for (i=0; i<moli->nsites; i++)  {
            moli->flags[i]	=	1;		// activate all the sites on this processor
            moli->parent[i]	=	i-1;		// initialize parent site
         }
         moli->flip		=	0;		// flip to the original direction
         moli->origin	=	CenterofMass(moli);
      }

#ifdef CELL_LIST					// build cell list
      if (init) {
         CL_Init();		
         init	=	0;
      } 
      else 
         CL_Destroy();

      CL_Build();
#endif	/* CELL_LIST */

      // Perform Analysis
      CalcV();				// calc energy and virial
      if (V_VIRIAL)
         Sample_Pressure();

      SampleP2All();			// sample P2 and P2m and local p2
      Dist_p2();			// put local p2 into distribution
SampleM_Q();
      SampleSpherical();		// sample spherical coordinates
      Dist_Spherical();			// put spherical coord. into distribution
      Find_Nuclei_p2(1);
      nmaxp2[0]	=	nmax[0][0];
      nmaxp2[1]	=	nmax[0][1];
      nmaxp2[2]	=	nmax[0][2];
sizeofnucl	=	sizeofnuclp2;
sizedist	=	sizedistp2;
//      Find_Nuclei(dynvar);		// find nuclei
//      Calc_Qlm(6);			// calc. Qlm for LJ system, require CELL_LIST

/*
	 for (moli=mol; moli<mol+NMOLS; moli++)		// calculate the drift of biggest nucleus
            if (sizeofnucl[moli->nuclid[0]] == nmax[0][0]) {
               con	=	CenterofNucleus(moli->nuclid[0], moli);
	       break;
	    }
         if (timestep==0) {
	    rO	=	con;
	    rshift2	=	0;
         }
         else {
	    con	=	V_Subtr(&con, &rO);
	    con	=	MapInBox2(&con, PBC, system);
	    rshift2	=	V_Dot(&con, &con);
         }
*/
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

      // chain orientation distribution, test the rotator phase

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
	 orientdist[(int)(orient/M_PI*90)]	++;
//	 printf("%f %f %f %d\n", chainface.x, chainface.y, orient, orientdist[(int)(orient/M_PI*90)]);
      }

      // calculate the # of sites with local p2 greater than a threshold

      nsitesp2	=	0;
      nmolsp2	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
//       nsitesp2	=	0;
         for (i=0; i<moli->nsites; i++) {
            if (moli->p2[i] > 0.4)
               nsitesp2	++;
         }
//            if (nsitesp2>1 && nsitesp2<6)
//	       nmolsp2	++;
      }

      Print_hst(fhst);		// print out hisgram
/*
Find_Nuclei_p2(-1);
printf("Xtal = %d\n Nmax = %d\n Nnucl = %d\n", Xtal[0], nmax[0][0], Nnucl[0]);
fprintf(fhst, " %4d\n", nmax[0][0]);
*/
      correlation();			// calculate correlation, in system units
      if (!(timestep%IRADIAL)) {
         radial("sample");		// sample radial distribution function
      }

      /* Group the segments belong to the biggest nucleus */
      n	=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         if (system == moli->box) {
            for (i=0; i<moli->nsites; i++) {
               if (sizeofnucl[moli->nuclid[i]] == nmax[system][0]) {
                  nucleus[n].moli	=	moli;
	          nucleus[n].site	=	i;
                  n	++;
	       }
            }
         }
      }

//      shiftbox(system, nucleus, n);

      /* Output */

      CoorSystem2SI();		// convert coordinates back to SI units

      // OUTPUT .car file for VMD visualization

      if (carflag) {
         fprintf(fout, "!BIOSYM archive 3\n");
         fprintf(fout, "PBC=ON\n");
         fprintf(fout, "!TIMESTEP %d\n", timestep);
         fprintf(fout, "!DATE %s", asctime(localtime(&t)));
         fprintf(fout, "PBC %9.7g %9.7g %9.7g %9.7g %9.7g %9.7g (P1)\n", 
		BOX[system].lx, BOX[system].ly, BOX[system].lz, 90.0, 90.0, 90.0);

	 n	=	0;
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {

               MolInBox2(moli);
               for (i=0; i<moli->nsites; i++) {

                  if (moli->p2[i]>critp2)
                     sprintf(s, "M%d", n++);
                  else
	             sprintf(s, "C%d", n++);

/*
                  if (sizeofnucl[moli->nuclid[i]] == nmax[system][0])		// note nuclid index starts from 1
                     sprintf(s, "M%d", n++);
                  else if (sizeofnucl[moli->nuclid[i]] -nmax[system][0] > -3 && nmax[system][0]>=10)
                     sprintf(s, "C%d", n++);
	          else if (moli->nuclid[i] >= 1)
                     sprintf(s, "O%d", n++);
                  else
                     sprintf(s, "H%d", n++);
*/
                  fprintf(fout, "%-5.5s ", s);
                  sprintf(s, "M%d", moli-mol);
                  fprintf(fout, "%14.8g %14.8g %14.8g ", moli->p[i].x, moli->p[i].y, moli->p[i].z);
                  strcpy(ff, "O");
                  fprintf(fout, "%-4.4s %-6d ND      %-3.3s 0.000\n", ff, moli-mol, Element(moli->type[i], s));
               } 
            }   
         }
         fprintf(fout, "end\nend\n");
      }

      // OUTPUT .pdb file for further analysis
      if (pdbflag) {
         fprintf(fpdb, "HEADER: file created from %s on %s", argv[1], asctime(localtime(&t)));
         fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
        	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

         n		=	0;
#define SIZECAP	5
         for (moli=mol; moli<mol+NMOLS; moli++) {
            if (system==moli->box) {
               MolInBox2(moli);

               drawmol	=	0;
               for (i=0; i<moli->nsites; i++) {
                  if (sizeofnucl[moli->nuclid[i]] == nmax[system][0]) {
	             drawmol	=	1;		// participate the biggest nucleus
		     break;
                  }
	       }

               for (i=0; i<moli->nsites; i++) {

                  if (drawmol) {
                     if (sizeofnucl[moli->nuclid[i]] == nmax[system][0]) {	// nuclid index starts from 1
                        atomname	=	'N';
		        fprintf(fdat, " 10");
                     }
                     else {
                        atomname	=	'O';
                        fprintf(fdat, " 0");
                     }
                  }
                  else {
		     atomname	=	'C';
		     fprintf(fdat," -1");
	          }

                  n	++;

	          fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
                  fprintf(fpdb, "%5d ", n);	// atom number
                  fprintf(fpdb, " %c  ", atomname);	// atom name
                  fprintf(fpdb, " ");		// alternate location indiator
  	          fprintf(fpdb, "   ");		// residue name
	          fprintf(fpdb, " ");		// column 21
                  fprintf(fpdb, " ");		// chain identifier, column 22
	          fprintf(fpdb, "    ");		// residue sequence number, 23-26
	          fprintf(fpdb, " ");		// code for insertion of residues, 27
                  fprintf(fpdb, "   ");		// column 28-30
                  fprintf(fpdb, "%8.3f%8.3f%8.3f", moli->p[i].x, moli->p[i].y, moli->p[i].z);
                  fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
                  fprintf(fpdb, "%5.5s", "");
                  fprintf(fpdb, "\n"); 
               } 
            }   
         }
         fprintf(fpdb, "END\n");
      }	//pdbflag

      // OUTPUT configuration file for further analysis

      if (confflag) { 
         fprintf(fconf, "!TIMESTEP %d\n", timestep);
         fprintf(fconf, "%d\t%d\t%d\n", NSYSTEMS, NMOLS, NSITES);
         for (i=0; i<NSYSTEMS; i++)
            fprintf(fconf, "%f\t%f\t%f\n", BOX[i].lx, BOX[i].ly, BOX[i].lz);

         for (moli=mol; moli<mol+NMOLS; moli++) {
            fprintf(fconf, "%d\t%d\t%d\n", moli-mol, moli->box, moli->nsites);
            //fprintf(fconf, "%d\t%d\t%d\t%d\t%d\n", i, moli.box, moli.nsites, moli.fix, moli.flip);
            //MolInBox(moli);
            for (i=0; i<moli->nsites; i++) 
               fprintf(fconf, "%d\t%f\t%f\t%f\n", moli->type[i], moli->p[i].x, moli->p[i].y, moli->p[i].z);
         }
      }
   } 

   // output average center of mass of chains
   if (pdbflag) {
      fprintf(fpdb, "HEADER: pdb file created from %s on %s", argv[1], asctime(localtime(&t)));
      fprintf(fpdb, "CRYST1%9.4f%9.4f%9.4f%7.2f%7.2f%7.2f P 1           1\n", 
 	BOX[system].lx, BOX[system].ly, BOX[system].lz,	90.0, 90.0, 90.0);

      n		=	0;
      for (moli=mol; moli<mol+NMOLS; moli++) {
         i	=	moli-mol;
         if (ncom[i]>0) {
            n	++;
            atomname	=	'N';
            fprintf(fpdb, "ATOM  ");		// pdb command, column 1-6
            fprintf(fpdb, "%5d ", n);		// atom number
            fprintf(fpdb, " %c  ", atomname);	// atom name
            fprintf(fpdb, " ");		// alternate location indiator
            fprintf(fpdb, "   ");	// residue name
            fprintf(fpdb, " ");		// column 21
            fprintf(fpdb, " ");		// chain identifier, column 22
            fprintf(fpdb, "    ");	// residue sequence number, 23-26
            fprintf(fpdb, " ");		// code for insertion of residues, 27
            fprintf(fpdb, "   ");	// column 28-30
            fprintf(fpdb, "%8.3f%8.3f%8.3f", com[i].x/ncom[i], com[i].y/ncom[i], com[i].z/ncom[i]);
            fprintf(fpdb, "%6.2f%6.2f", 1.0, 0.0);	// occupance and temperature
            fprintf(fpdb, "%5.5s", "");
            fprintf(fpdb, "\n"); 
         }
      }
      fprintf(fpdb, "END\n");
   }

   S_PrintAll();		// print out distributions
   corr_norm();		// normalize correlation function
   corr_print();		// print out correlation function
   radial("print");		// print out radial distribution function

   fclose(fin);		// close files
   fclose(fhst);
   if (carflag)		fclose(fout);
   if (confflag)    	fclose(fconf);
   if (pdbflag)	     {	fclose(fpdb); fclose(fdat);}

   for (i=0; i<180; i++) {
      printf("%d %d\n", i, orientdist[i]);
   }
   return	0;
}


