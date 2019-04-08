/*
    program:    builder.c
    author:     Peng Yi at MIT
    date:       November 15, 2007
    purpose:    Initialize polymer configuration for one single box.
    notes:
		July 25, 2009 enable polydispersity in Amorphous()
*/

#define __MAIN_PROGRAM
#include "header.h"

#define Rcut2	0.5*16.0801		// cutoff^2, unit Angstrom

/* Variables that are only useful in builder.c */

long		NX, NY, NZ;
double		cella, cellb, cellc, angle;
char		title1[256], title2[256], comments[256];

void BuilderReadSetup()
{
   FILE		*fPtr;
   long		i;
   char		c[256], name[256], a[256];

   if ( (fPtr=fopen("in","r"))==NULL )
      printf("Open setup file failed.\n");
   else {
      fgets(title1, sizeof(title1), fPtr);
      fgets(title2, sizeof(title2), fPtr);
      fgets(comments, sizeof(comments), fPtr);

      fscanf(fPtr, "%s", name);
      while(strcmp(name, "END_SETUP")) {

         if (!strcmp(name, "NMOLS"))			GetLVar(fPtr, 1, &NMOLS);
         else if (!strcmp(name, "NSITES"))		GetLVar(fPtr, 1, &NSITES);
         else if (!strcmp(name, "NBOX")) {		GetLVar(fPtr, 1, &NBOX);
            NSYSTEMS	=	NBOX;
	 }
         else if (!strcmp(name, "Dpoly"))		GetDVar(fPtr, 1, &Dpoly);
         else if (!strcmp(name, "Rho"))			GetDVar(fPtr, 1, &Rho);
         else if (!strcmp(name, "kT"))			GetDVar(fPtr, 1, &kT);
         else if (!strcmp(name, "PBC"))			GetLVar(fPtr, 1, &PBC);
         else if (!strcmp(name, "NX"))			GetLVar(fPtr, 1, &NX);
         else if (!strcmp(name, "NY"))			GetLVar(fPtr, 1, &NY);
         else if (!strcmp(name, "NZ"))			GetLVar(fPtr, 1, &NZ);
         else if (!strcmp(name, "a"))			GetDVar(fPtr, 1, &cella);
         else if (!strcmp(name, "b"))			GetDVar(fPtr, 1, &cellb);
         else if (!strcmp(name, "c"))			GetDVar(fPtr, 1, &cellc);
         else if (!strcmp(name, "angle"))		GetDVar(fPtr, 1, &angle);

         /* Forcefield parameters */

	 else if (!strcmp(name, "NTYPES")) 		GetLVar(fPtr, 1, &NTYPES);
         else if (!strcmp(name, "DLJ")) 		GetLVar(fPtr, 1, &DLJ);
         else if (!strcmp(name, "TYPEMASS")) {  
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].M));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "SIGMA")) {  
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].SIGMA));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "EPSILON")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].EPSILON));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "FIXBONDL")) {
	    fscanf(fPtr, "%s%ld", c, &(FIXBONDL));	fgets(comments,sizeof(comments),fPtr);
	 }
         else if (!strcmp(name, "KSTRETCH")) { 
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].KSTRETCH));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "LSTRETCH")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].LSTRETCH));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "FIXBONDA")) {
	    fscanf(fPtr, "%s%ld", c, &(FIXBONDA));	fgets(comments,sizeof(comments),fPtr);
	 }
         else if (!strcmp(name, "KBENDING")) { 
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].KBENDING));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "THETA")) { 
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].THETA));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "HS")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].HS));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION0")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[0]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION1")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[1]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION2")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[2]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION3")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[3]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION4")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[4]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else if (!strcmp(name, "TORSION5")) {
            fscanf(fPtr, "%s%s", c, a);
            for (i=0; i<MAXNTYPES; i++)
               fscanf(fPtr, "%lf", &(type[i].TORSION[5]));
            fgets(comments, sizeof(comments), fPtr);
         }
         else
            fgets(comments, sizeof(comments),fPtr);
	 
         fscanf(fPtr, "%s", name);
      }
      fgets(comments,sizeof(comments),fPtr);
      fclose(fPtr);
   }
   return;
}


void BuilderPrintSetup()
{
   FILE		*fPtr;
   long		i, j;
   char		name[256];
 
   if ( (fPtr=fopen("builder.log","a"))==NULL )
      printf("Builder log file failed.\n");
   else {
      fprintf(fPtr, "%s", title1);
      fprintf(fPtr, "%s", title2);
      fprintf(fPtr, "Program Version: \t%s\n", VERSION);
      for (i=0; i<80; i++)	
	 fprintf(fPtr, "*");
      fprintf(fPtr, "\n");

      fprintf(fPtr, "\n/* Box dimension */\n\n");

      fprintf(fPtr, "NBOX   \t=\t%ld\n", NBOX);
      fprintf(fPtr, "NMOLS  \t=\t%ld\n", NMOLS);
      fprintf(fPtr, "NSITES \t=\t%ld\n", NSITES);
      fprintf(fPtr, "Dpoly  \t=\t%f\n", Dpoly);
      fprintf(fPtr, "Rho    \t=\t%f\n", Rho);
      fprintf(fPtr, "PBC    \t=\t%ld\n", PBC);
      fprintf(fPtr, "NX     \t=\t%ld\n", NX);
      fprintf(fPtr, "NY     \t=\t%ld\n", NY);
      fprintf(fPtr, "NZ     \t=\t%ld\n", NZ);
      fprintf(fPtr, "cella  \t=\t%f\n", cella);
      fprintf(fPtr, "cellb  \t=\t%f\n", cellb);
      fprintf(fPtr, "cellc  \t=\t%f\n", cellc);
      fprintf(fPtr, "angle  \t=\t%f\n", angle);
     
      fprintf(fPtr, "\n/* Forcefield parameters */\n\n");
 
      fprintf(fPtr, "NTYPES	\t=\t%ld\n", NTYPES);
      fprintf(fPtr, "DLJ	\t=\t%ld\n", DLJ);
      fprintf(fPtr, "TYPEMASS\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].M);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "SIGMA	\t=\t{ ");
      for (i=0; i<NTYPES; i++)	        fprintf(fPtr, "%f\t", type[i].SIGMA);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "EPSILON	\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].EPSILON);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "KSTRETCH\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].KSTRETCH);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "LSTRETCH\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].LSTRETCH);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "KBENDING\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].KBENDING);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "THETA	\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].THETA);
      fprintf(fPtr, "}\n");
      fprintf(fPtr, "HS	\t=\t{ ");
      for (i=0; i<NTYPES; i++)		fprintf(fPtr, "%f\t", type[i].HS);
      fprintf(fPtr, "}\n");
      for (j=0; j<6; j++) {
         sprintf(name, "TORSION%c", '0'+j);
         fprintf(fPtr, "%s\t=\t{ ", name);
         for (i=0; i<NTYPES; i++)	fprintf(fPtr, "%f\t", type[i].TORSION[j]);
         fprintf(fPtr, "}\n");
      }
      fclose(fPtr);
   }
   return; 
}

double tors_angle()
{
   double		angle, cosb, select;
   static double	value[360], accum[360], sum;
   static long		init =1;
   long			i;
   typestruct		t = type[0];

   if (init) {
      sum	=	0;
      for (i=0; i<360; i++) {
	 value[i]	=	0;
      }
      for (i=0; i<360; i++) {
         angle	=	M_PI/180.0*i;    	// [0, 2pi]  
         cosb	=	cos(angle);
         value[i]	=	exp(-OPLS2(cosb, t.TORSION[1], t.TORSION[2], t.TORSION[3])/(kT* N_KB * N_NAV * 1e-3));
	 sum	+=	value[i];
      }
/*
      cosb	=	1.0;			// trans
      value[0]	=	exp(-OPLS2(cosb, t.TORSION[1], t.TORSION[2], t.TORSION[3])/kT);
      cosb	=	-0.5;			// g+ and g-
      value[1]	=	exp(-OPLS2(cosb, t.TORSION[1], t.TORSION[2], t.TORSION[3])/kT);
      value[2]	=	value[1];
      sum	=	value[0] + value[1] + value[2];
*/
      init	=	0; 
   }
   select	=	ran1(seed)*sum;

   i	=	0;
   while (select>value[i]) {
      select	-=	value[i];
      i	++;
   }
   return	M_PI + M_PI/180.0*i;

/*
   if (select < value[0]/sum) 
      return	0;
   else if (select < (value[0]+value[1])/sum)
      return	M_PI * 2.0/3; 
   else
      return	- M_PI * 2.0/3;
*/
}


void Amorphous(long nsystems, long nmols, long nsites)	// generate random chain conf.
{
   molstruct		*moli, *molj;
   long			i, j, n, avensites, ndelta, ibox, system, reject, nreject;
   double		r2, vol, numberdensity;
   double		LX, LY, LZ;
   vector		dp;
   sphere		s;
   typestruct		t = type[0];

   static long		Dnsites, maxnsites, minnsites;		// for polydispersity
   static long		init=1;

   if (init) {
      init	=	0;
      Dnsites	=	(int) (Dpoly*nsites/nmols);
      maxnsites	=	nsites/nmols + Dnsites;
      minnsites	=	nsites/nmols - Dnsites;
   }

   double 		l0 = type[0].LSTRETCH, 		// bond length
			theta0 = type[0].THETA,		// bond angle
                	mass = type[0].M;

   numberdensity	=	Rho*0.6022/mass;	// Rho in unit of g/cm^3
   vol	=	(nsites/nsystems)/numberdensity;	// vol in unit of Anstrom^3
   LX	=	pow(vol, 1.0/3);			// cubic
   LY	=	LX;
   LZ	=	LX;

   printf("### RANDOM NUMBER TEST (PRINT 5 RANDOM NUMBERS):\n");
   for (i=0; i<5; i++) {
      printf("%f\n", ran1(seed));
   }

   for (i=0; i<nsystems; i++) {
      BOX[i].lx		=	LX;
      BOX[i].ly		=	LY;
      BOX[i].lz		=	LZ;
      BOX[i].temp	=	kT;
   }

   for (moli=mol; moli<mol+nmols; moli++) {
      n			=	(int) (ran1(seed) * (2*Dnsites+1));	// [0, 2Dnsites]
      moli->nsites	=	minnsites + n;		// polydispersity

      moli->type[0]			=	0;	// set the type of end sites
      moli->type[moli->nsites-1]	=	0;
      for (i=1; i<moli->nsites-1; i++) {		// set the type of middle sites
         moli->type[i]			=	1;
      }
/*
      moli->type[0]			=	(moli-mol)*2+1;	// set the type of end sites
      moli->type[moli->nsites-1]	=	(moli-mol)*2+1;
      for (i=1; i<moli->nsites-1; i++) {		// set the type of middle sites
         moli->type[i]			=	(moli-mol)*2+2;
      }
*/
      for (i=0; i<moli->nsites; i++)			// set parent site
         moli->parent[i]	=	i-1;
   }

   for (moli=mol; moli<mol+nmols; moli++) {
      system	=	(moli-mol)/(nmols/nsystems);	// set box id
      moli->box	=	system;

      // Place the first atom on a chain
//printf("moli=%ld \n", moli-mol);
      do {
         reject	=	0;
         moli->p[0].x	=	BOX[system].lx * (ran1(seed)-0.5);
         moli->p[0].y	=	BOX[system].ly * (ran1(seed)-0.5);
         moli->p[0].z	=	BOX[system].lz * (ran1(seed)-0.5);
/*
         moli->p[0].x	=	0;
         moli->p[0].y	=	0;
         moli->p[0].z	=	0;
*/

         for (molj=mol; molj<moli; molj++) {		// overlap test
            for (j=0; j<molj->nsites; j++) {
               r2	=	DistSQ(moli->p[0], molj->p[j], system);
               if (r2 < Rcut2) {			// 1.25sigma ~ 5A
                  reject	=	1;
                  break;
               }
            }
            if (reject)
               break;
         }
      } while (reject);

      for (i=1; i<moli->nsites; i++) {			// subsequent atoms
         s.d		=	l0;			// atom spherical coord.
         s.alpha	=	theta0;
//printf("moli=%ld i=%ld\n", moli-mol, i);
         nreject	=	0;		// # of rejection for this site
         do {
            s.beta	=	tors_angle();
            moli->p[i]	=	SiteCartesian(moli, i, s);	// cart. coord.

            reject	=	0;

            for (molj=mol; molj<=moli; molj++) {		// overlapping test
               for (j=0; j< (molj==moli ? i-3 : molj->nsites); j++) {
                  r2	=	DistSQ(moli->p[i], molj->p[j], system);
//printf("molj=%ld j=%ld\n", molj-mol, j);
//printf("r2=%f Rcut2=%f\n", r2, Rcut2);

                  if (r2 < Rcut2) {			// 1.25 sigma ~ 5 Angstrom
                     reject	=	1;
		     nreject	++;
	             break;
                  }
               }
	       if (reject)
		  break;
            }
         } while (nreject >= 1 && nreject <= 5);

         if (nreject > 5) {		// if continuously rejected
            i	-=	2;		// back one atom
            if (i==-1) {
               do {
                  reject	=	0;
                  moli->p[0].x	=	BOX[system].lx * (ran1(seed)-0.5);
                  moli->p[0].y	=	BOX[system].ly * (ran1(seed)-0.5);
                  moli->p[0].z	=	BOX[system].lz * (ran1(seed)-0.5);

                  for (molj=mol; molj<moli; molj++) {		// overlap test
                     for (j=0; j<molj->nsites; j++) {
                        r2	=	DistSQ(moli->p[0], molj->p[j], system);
                        if (r2 < Rcut2) {			// 1.25sigma ~ 5A
                           reject	=	1;
                           break;
                        }
                     }
                     if (reject)               break;
                  }
               } while (reject);
	       i	++;
            }
         }
      }
   }
   return;
}


void Orthorhombic(long nbox, long nmols, long nsites, long nx, long ny, long nz, double a, double b, double c, double angle)
{
   // nx, ny and nz are the number of unit cells on x, y and z direction, respectively
   // a, b and c are the dimension of unit cells on x, y and z direction, respectively
   // angle is the angle between chain surface and x direction

   long		i, j, k;
   double	bondlx, bondly, bondlz, theta;
   double	cosy, siny, angletilt, temp;
   vector	com;

   theta	=	type[0].THETA / 180.0 * M_PI;	// bond angle
   angle	=	angle / 180.0 * M_PI;		// packing angle on x-y plane
//   angletilt	=	15.0 /180.0 * M_PI;		
   angletilt	=	0;

   bondlx	=	type[0].LSTRETCH * sin(0.5 * theta) * fabs(cos(angle)); 
   bondly	=	type[0].LSTRETCH * sin(0.5 * theta) * fabs(sin(angle));
   bondlz	=	type[0].LSTRETCH * cos(0.5 * theta);
//   c		=	bondlz * nsites/nmols;

   /* determine box size */

   for (i=0; i<nbox; i++) {
      BOX[i].lx		=	nx * a;
      BOX[i].ly		=	ny * b;
      BOX[i].lz		=	nz * c;
   }
   for (i=0; i<nmols; i++) {
      mol[i].nsites	=	nsites/nmols;		// monodisperse
   }

   /* Chains on the first z layer */

   // first 2 sites of the 2 chains in the first unit cell in the first z layer

   mol[0].p[0].x	=	0 - 0.5 * bondlx - 0.5*nx*a;
   mol[0].p[0].y	=	0 - 0.5 * bondly - 0.5*ny*b;
   mol[0].p[1].x	=	0 + 0.5 * bondlx - 0.5*nx*a;
   mol[0].p[1].y	=	0 + 0.5 * bondly - 0.5*ny*b;
   mol[1].p[0].x	=	0.5*a - 0.5 * bondlx - 0.5*nx*a;
   mol[1].p[0].y	=	0.5*b + 0.5 * bondly - 0.5*ny*b;
   mol[1].p[1].x	=	0.5*a + 0.5 * bondlx - 0.5*nx*a;
   mol[1].p[1].y	=	0.5*b - 0.5 * bondly - 0.5*ny*b;

   // first 2 sites of all other chains in the first z layer

   for (i=0; i<nx; i++) {
      for (j=0; j<ny; j++) {
         for (k=0; k<2; k++) {
            mol[(i*ny+j)*2].p[k].x	=	mol[0].p[k].x	+ i * a;
            mol[(i*ny+j)*2].p[k].y	=	mol[0].p[k].y	+ j * b;
            mol[(i*ny+j)*2+1].p[k].x	=	mol[1].p[k].x	+ i * a;
            mol[(i*ny+j)*2+1].p[k].y	=	mol[1].p[k].y	+ j * b;
         }
      }
   }

   // all other sites of all chains in the first z layer
 
   for (i=0; i<nx*ny*2; i++) {
      for (j=2; j<nsites/nmols; j++) {
         mol[i].p[j].x		=	mol[i].p[mod(j,2)].x;
         mol[i].p[j].y		=	mol[i].p[mod(j,2)].y;
      }
      for (j=0; j<nsites/nmols; j++)
         mol[i].p[j].z		=	j * bondlz - 0.5*NZ*c;
   }

   /* Multiple layers along z-direction */

   // chains on other z layers

   for (k=1; k<nz; k++) {
      for (i=0; i<2*nx*ny; i++) {			// one z-layer has 2*NX*NY chains
         for (j=0; j<nsites/nmols; j++) {
            if (!mod(nsites/nmols,2)) {			// even number sites 
               mol[k*2*nx*ny+i].p[j].x	=	mol[i].p[j].x;
               mol[k*2*nx*ny+i].p[j].y	=	mol[i].p[j].y;
               mol[k*2*nx*ny+i].p[j].z	=	mol[i].p[j].z + k*c;
            }
            else {
               mol[k*2*nx*ny+i].p[j].x	=	mol[i].p[(j-1)<0 ? j+1 : j-1].x;
               mol[k*2*nx*ny+i].p[j].y	=	mol[i].p[(j-1)<0 ? j+1 : j-1].y;
               mol[k*2*nx*ny+i].p[j].z	=	mol[i].p[j].z + k*c;
            }
         }
      }
   }

   /* Rotating */

   
   cosy		=	cos(angletilt);
   siny		=	sin(angletilt);

   for (i=0; i<2*nx*ny*nz; i++) {
      com	=	CenterofMass(mol+i);

      for (j=0; j<mol[i].nsites; j++) {
         mol[i].p[j].x	-=	com.x; 
         mol[i].p[j].y	-=	com.y; 
         mol[i].p[j].z	-=	com.z; 

         // rotate about y axis
         mol[i].p[j].x	=	cosy * mol[i].p[j].x + siny * mol[i].p[j].z;
         mol[i].p[j].z	=      -siny * mol[i].p[j].x + cosy * mol[i].p[j].z;

         mol[i].p[j].x	+=	com.x; 
         mol[i].p[j].y	+=	com.y; 
         mol[i].p[j].z	+=	com.z; 
      }
   } 
   
   
   /* Multiple boxes */

   for (i=1; i<nbox; i++)
      for (j=0; j<nmols/nbox; j++)
         mol[i*nmols/nbox+j]	=	mol[j];
   for (i=0; i<nmols; i++)
      mol[i].box	=	i/(nmols/nbox);

   /* molecule site type */

   for (i=0; i<nmols; i++) {
      for (j=0; j<mol[i].nsites; j++)
         mol[i].type[j]		=	0;

      if (NTYPES>=2) {
         for (j=1; j<mol[i].nsites-1; j++) 
            mol[i].type[j]	=	1;
      }
   }
   /*
   for (i=0; i<nmols; i++) {			// for melting simulation
      if ( mod(i/(nx*ny), 2)) {
         for (j=0; j<mol[i].nsites; j++)
            mol[i].type[j]	+=	2;
         mol[i].fix	=	1;		// not used yet
      }
      else
         mol[i].fix	=	0;
   }
   */
   for (i=0; i<nmols; i++)
      if (mol[i].p[0].y > 0) {		// interface along x direction
//	if (mol[i].p[0].x > 0) {		// interface along y direction
//      if ( (mol[i].p[0].y+mol[i].p[1].y) >  (mol[i].p[0].x+mol[i].p[1].x) * b/a) {	// interface along y=b/a*x direction
	 for (j=0; j<mol[i].nsites; j++)
	     mol[i].type[j]	+=	2;
      }
   /*
   for (i=0; i<nmols; i++) {			// for melting, create a cylinder
      if ( mol[i].p[0].x * mol[i].p[0].x + mol[i].p[0].y * mol[i].p[0].y < 169 ) {
         for (j=0; j<mol[i].nsites; j++) {
               mol[i].type[j]	+=	2;
         }
      }
   }
   */

   //***** Rotate about y axis *****//
   /*
   for (i=0; i<nmols; i++) {
      for (j=0; j<mol[i].nsites; j++) {
	 temp		=	mol[i].p[j].x;
	 mol[i].p[j].x	=	mol[i].p[j].z;
	 mol[i].p[j].z	=	-1 * temp;
      }
   }
   for (i=0; i<nbox; i++) {
      temp		=	BOX[i].lx;
      BOX[i].lx		=	BOX[i].lz;
      BOX[i].lz		=	temp;
   }
   */

   return;
}


void sc_lattice(long N, double rho, long PBC)		// total particle number N = NC * NC * NC
{							// box dimension L
   long		i, NC;
   double	L, cell;

   for (i=0; i<N; i++)
      mol[i].nsites	=	1;			// only one atom in one molecule

   if (PBC==1)
      NC	=	(int) rint(pow(N, 1.0/3));

   L	=	pow(N/rho, 1.0/3);			// cubic box length

   BOX[0].lx	=	L;				// only for one box now 4/18/08
   BOX[0].ly	=	L;
   BOX[0].lz	=	L;
   BOX[0].lbox	=	L;

   if (PBC==1 && N!=NC*NC*NC)
      Exit("builder", "sc_lattice", "number of particle incomp. with sc lattice.\n");

   if (PBC==1) {
      cell	=	L/NC;

      for (i=0; i<N; i++) {
         mol[i].p[0].x	=	(mod(i, NC) - NC/2) * cell;			//or (int)mod(i,M)/1
         mol[i].p[0].y	=	((int)(mod(i, NC*NC)/NC) - NC/2) * cell;	//or (int)mod(i,M*M)/M
	 mol[i].p[0].z	=	((int)(i/(NC*NC)) - NC/2) * cell;		//or (int)mod(i,M*M*M)/(M*M)
      }
   }
}


void fcc_lattice(long N, double rho, long PBC)			//reference: F.23 of Allen and Tildesley
{							// now only for one box 4/8/2008
   long		i, j, k, m, ref;
   double	L, cell, cell2;
   long		NC;

   for (i=0; i<N; i++) 
      mol[i].nsites	=	1;				// only one atom on one chain

   if (PBC==1)
      NC	=	(int) rint(pow(N/4, 1.0/3));		// rint: nearest integer value
   else if (PBC==2)
      NC	=	(int) rint(pow(N/16, 1.0/3));

   L		=	pow(N/rho, 1.0/3);			// cubic box length

   BOX[0].lx	=	L;
   BOX[0].ly	=	L;
   BOX[0].lz	=	L;
   BOX[0].lbox	=	L;

   if ( (PBC==1 && N!=4*NC*NC*NC) || (PBC==2 && N!=16*NC*NC*NC))
      Exit("builder", "fcc_lattice", "number of particles incomp. with fcc lattice.\n");

   if (PBC==1) {
      cell	=	L / NC;				// the reason to use plus sign here is to adapt with our
      cell2	=	0.5 * cell;			// p.b.c., that is if X<-LBOX/2, then X+=LBOX
							// if our pbc is if X<=-LBOX/2, then X+=LBOX, then minus sign
							// in another words, our box is [-L, +L)
      mol[0].p[0].x	=	0;	mol[0].p[0].y	=	0;	mol[0].p[0].z	=	0;
      mol[1].p[0].x	=	cell2;	mol[1].p[0].y	=	cell2;	mol[1].p[0].z	=	0;
      mol[2].p[0].x	=	0;	mol[2].p[0].y	=	cell2;	mol[2].p[0].z	=	cell2;
      mol[3].p[0].x	=	cell2;	mol[3].p[0].y	=	0;	mol[3].p[0].z	=	cell2;

      m	=	0;

      for (i=0; i<NC; i++) {
         for (j=0; j<NC; j++) {
	    for (k=0; k<NC; k++) {
	       for (ref=0; ref<4; ref++) {
		  mol[ref+m].p[0].x	=	mol[ref].p[0].x	+	cell*k;
		  mol[ref+m].p[0].y	=	mol[ref].p[0].y	+	cell*j;
		  mol[ref+m].p[0].z	=	mol[ref].p[0].z	+	cell*i;
	       }
	       m	+=	4;
	    }
	 }
      }
      for (i=0; i<N; i++) {
	  mol[i].p[0].x	-=	(double)L/2;		//same reason to use += instead of -= as mentioned above
	  mol[i].p[0].y	-=	(double)L/2;
	  mol[i].p[0].z	-=	(double)L/2;
      }
   }
   if (PBC==2) {					//truncated octahedron periodic boundary condition
      printf("Set up fcc lattice in truncated octahedron box.\n");

      cell	=	L / (2*NC);
      cell2	=	0.5 * cell;

      mol[0].p[0].x	=	0;	mol[0].p[0].y	=	0;	mol[0].p[0].z	=	0;
      mol[1].p[0].x	=	cell2;	mol[1].p[0].y	=	cell2;	mol[1].p[0].z	=	0;
      mol[2].p[0].x	=	0;	mol[2].p[0].y	=	cell2;	mol[2].p[0].z	=	cell2;
      mol[3].p[0].x	=	cell2;	mol[3].p[0].y	=	0;	mol[3].p[0].z	=	cell2;

      m	=	0;

      for (i=0; i<NC; i++) {			// (2*NC, 2*NC, NC)
         for (j=0; j<2*NC; j++) {
	    for (k=0; k<2*NC; k++) {
               for (ref=0; ref<4; ref++) {
		  mol[ref+m].p[0].x	=	mol[ref].p[0].x	+	cell*k;
		  mol[ref+m].p[0].y	=	mol[ref].p[0].y	+	cell*j;
		  mol[ref+m].p[0].z	=	mol[ref].p[0].z	+	cell*i;
	       }
	       m	+=	4;
	    }
	 }
      }
     
      for (i=0; i<N; i++) {
	 mol[i].p[0].x	-=	0.5 * L;
         mol[i].p[0].y	-=	0.5 * L;

	 if ( (fabs(mol[i].p[0].x) + fabs(mol[i].p[0].y) + fabs(mol[i].p[0].z)) > 0.75*L) {
            mol[i].p[0].x	+=	((mol[i].p[0].x >=0) ? -0.5 : 0.5) * L;
            mol[i].p[0].y	+=	((mol[i].p[0].y >=0) ? -0.5 : 0.5) * L;
            mol[i].p[0].z	+=	((mol[i].p[0].z >=0) ? -0.5 : 0.5) * L;
	 }
	 else if ( (fabs(mol[i].p[0].x) + fabs(mol[i].p[0].y) + fabs(mol[i].p[0].z) == 0.75 * L) && mol[i].p[0].z>=0) {
	    mol[i].p[0].x	+=	((mol[i].p[0].x >= 0) ? -0.5 : 0.5) * L;
	    mol[i].p[0].y	+=	((mol[i].p[0].y >= 0) ? -0.5 : 0.5) * L;
	    mol[i].p[0].z	-=	0.5 * L;
	 }
      }
   }
}


void BuilderCheckSetup(char * argv[])
{
   long	error=0;

   if (!strcmp(argv[1], "chain")) {
      if (NMOLS/NSYSTEMS != NX*NY*NZ*2) {
         printf("NMOLS mismatch!\n");
         exit(-1);
      }
   }
   if (NMOLS > MAXNMOLS) 		{printf("NMOLS>MAXNMOLS\n");	error=1;}
   if (NSITES/NMOLS > MAXNMOLSITES) 	{printf("NMOLSITES error\n");	error=1;}

   if (error) {
      printf("setup parameter error\n");
      exit(-1);
   }
}


int main(int argc, char * argv[])
{
   
   if (argc<2) {
      printf("builder (c) 2008 by Peng Yi at MIT\n\n");
      printf("Usage:\n");
      printf("\tbuilder structure\n\n");
      printf("\tstructure:	ljfcc\n");
      printf("\t	 	ljsc\n");
      printf("\t	 	chain\n");
      printf("\t	 	chainrandom\n");
      printf("Notes:\n");
      printf("\t* require in file\n\n");
      exit(1);
   }

   tim=(int *)malloc(sizeof(int));     	//random number generator
   seed=(long *)malloc(sizeof(long));
   *tim=(int)time(NULL);
   *seed= -1*(*tim);           		//seed must start out as a negative long

//   InitMols(MAXNMOLS, MAXNMOLS);	// allocate max. memory for molecules
   InitMols(MAXNMOLS, 0);

   BuilderReadSetup();
   BuilderCheckSetup(argv);
   BuilderPrintSetup();

   if (!strcmp(argv[1], "chain")) 
      Orthorhombic(NBOX, NMOLS, NSITES, NX, NY, NZ, cella, cellb, cellc, angle);	
   else if (!strcmp(argv[1], "chainrandom"))
      Amorphous(NBOX, NMOLS, NSITES);
   else if (!strcmp(argv[1], "ljfcc"))
      fcc_lattice(NMOLS, Rho, PBC);
   else if (!strcmp(argv[1], "ljsc"))
      sc_lattice(NMOLS, Rho, PBC);

   CalcUnits(0);			// unit = 1, unit needed for Write_Conf() 
   //printf("%f  %f  %f\n", mol[0].p[0].x, mol[0].p[0].y, mol[0].p[0].z);
   Write_Conf(-1);
   return 0;
}

