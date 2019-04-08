/*
    program:	varbridge.c
    author:	Peng Yi borrowed from Pieter J. in 't Veld
    date:	January 10, 2008
    purpose:	This module functions as a bridge between old style
    		system variables and new style system variables.
*/
#define __VARBRIDGE_MODULE
#include "varbridge.h"

void OrganizeMols()			// sort molecules according to their system id
{
  long			i, j, k = 0, nmols, nsites;
  molstruct		molt;
  
  for (i=0; i<NSYSTEMS; ++i)
  {
    nmols		= 0;
    nsites		= 0;
    for (j=0; j<NMOLS; ++j)
      if (mol[j].box==i)
      {
        molt		= mol[k];
	mol[k++]	= mol[j];
	mol[j]		= molt;
	nsites		+= mol[j].nsites;
        ++nmols;
      }
    NMols[i]		= nmols;
    NSites[i]		= nsites;
  }
}

bridgestruct		VB_BridgeMap;	// a collection of values and pointers pointing to system variables 
					// and molecules.  For binary history file output

bridgestruct *Bridge()
{
  return &VB_BridgeMap;
}


void BridgeSwap(long system1, long system2)
{
  molstruct		*mol_t;
  bridgestruct		*c = &VB_BridgeMap;

  mol_t			= c->mol[system1];
  c->mol[system1]	= c->mol[system2];
  c->mol[system2]	= mol_t;
}


bridgestruct *BridgeMap(long flag_mpi)
{
  long			i;
  bridgestruct		*c = &VB_BridgeMap;
  
  if (!flag_mpi)
    OrganizeMols();
    
  for (i=0; i<MAXNSYSTEMS; ++i)			// link to system and mol variables
  {
    c->system[i].n	= i;			// system.n is value, not pointer

    //GetBoxVectors(&(c->system[i].a), &(c->system[i].b), &(c->system[i].c));
    /*
    c->system[i].a.x	= BOX[i].lx;
    c->system[i].a.y	= 0;
    c->system[i].a.z	= 0;
    c->system[i].b.x	= 0;
    c->system[i].b.y	= BOX[i].ly;
    c->system[i].b.z	= 0;
    c->system[i].c.x	= 0;
    c->system[i].c.y	= 0;
    c->system[i].c.z	= BOX[i].lz;
    */
    c->system[i].lx	= &(BOX[i].lx);
    c->system[i].ly	= &(BOX[i].ly);
    c->system[i].lz	= &(BOX[i].lz);

    c->system[i].nmols	= NMols+i;		// pass pointer (address) rather than value
    c->system[i].nsites	= NSites+i;
    c->system[i].pres	= &(BOX[i].pres);
    c->system[i].vol	= &(BOX[i].vol);
    c->system[i].temp	= &(BOX[i].temp);
    c->system[i].drmax	= &(BOX[i].drmax);
    c->system[i].dlmax	= &(BOX[i].dlmax);
    c->system[i].damax	= &(BOX[i].damax);
    
    c->av[i]		= av+i;
    c->v[i]		= v+i;
    c->vir[i]		= vir+i;
    c->mol[i]		= flag_mpi ? mol : mol+i*NMols[i];
  } 
  for (i=0; i<D_NDIST; ++i)			// link to distributions
    c->dist[i]		= D_Distributions[i];

  //c->d_type		= &D_TYPE;
  c->d_type		= &NTYPES;
  c->type		= type;

  c->command.hs		= &V_HS;		// link to flags
  c->command.lj		= &V_LJ;
  c->command.ljshift	= &V_LJSHIFT;
  c->command.ljlrc	= &V_LJLRC;
  c->command.stretch	= &V_STRETCH;
  c->command.bending	= &V_BENDING;
  c->command.torsion	= &V_TORSION;
  c->command.virial	= &V_VIRIAL;
  c->command.scalecut	= &V_SCALECUTOFF;
  c->command.nvt	= &E_NVT;
  c->command.npt	= &E_NPT;
  c->command.gibbs	= &E_GIBBS;
  c->command.mpi	= &E_MPI;
  c->command.density	= &D_DENSITY;
  c->command.energy	= &D_ENERGY;
  c->command.pressure	= &D_PRESSURE;
  c->command.drift	= &D_DRIFT;
  c->command.torsional	= &D_TORSION;
  c->command.bonda	= &D_BONDA;
  c->command.bondl	= &D_BONDL;
  c->command.radial	= &D_RADIAL;
  c->command.localp2	= &D_LOCALP2;
  c->command.xtalsize	= &D_XTALSIZE;
  
/* 
  c->command.hs		= &P_HS;
  c->command.lj		= &P_LJ;
  c->command.stretch	= &P_STRETCH;
  c->command.ptorsion	= &P_TORSION;
  c->command.coulomb	= &P_COUL;
  c->command.virial	= &P_VIRIAL;
  c->command.polymer	= &E_POLYMER;
  c->command.mpi	= &E_MPI;
  c->command.async	= &S_ASYNC;
  c->command.monodisperse= &E_MONO;
  c->command.bias	= &E_BIAS;
  c->command.jacob	= &E_JACOB;
  c->command.nvt	= &E_NVT;
  c->command.npt	= &E_NPT;
  c->command.gibbs	= &E_GIBBS;
  c->command.insert	= &D_INSERT;
  c->command.widom	= &E_WIDOM;
  c->command.canonical	= &E_CANON;
  c->command.cavity	= &D_CAVITY;
  c->command.density	= &D_DENSITY;
  c->command.density3d	= &D_DENSITY3D;
  c->command.densfree	= &D_DENSITY_FREE;
  c->command.denstalobr	= &D_DENSITY_TALOBR;
  c->command.orient	= &D_ORIENTATION;
  c->command.orientfree = &D_ORIENT_FREE;
  c->command.orienttalobr= &D_ORIENT_TALOBR;
  c->command.orientcorr01 = &D_ORIENT_CORR_01;
  c->command.orientcorr02 = &D_ORIENT_CORR_02;
  c->command.orientcorr03 = &D_ORIENT_CORR_03;
  c->command.orientcorr04 = &D_ORIENT_CORR_04;
  c->command.orientcorr05 = &D_ORIENT_CORR_05;
  c->command.orientcorrCR = &D_ORIENT_CORR_CR;
  c->command.tails_etc	= &D_TAILS_ETC;
  c->command.e_n_function= &D_E_N_FUNCTION;
  c->command.radial	= &D_RADIAL;
  c->command.energy	= &D_ENERGY;
  c->command.torsion	= &D_TORSION;
  c->command.re_torsion	= &D_RE_TORSION;
  c->command.loopreentry= &D_LOOPREENTRY;
  c->command.b_length	= &D_B_LENGTH;
  c->command.b_angle	= &D_B_ANGLE;
  c->command.d_bridge	= &D_D_BRIDGE;
  c->command.e_profile	= &D_E_PROFILE;
  c->command.w_profile	= &D_W_PROFILE;
  c->command.n_profile	= &D_N_PROFILE;
  c->command.hs_dens	= &HS_DENS;
  c->command.temper	= &E_TEMPER;
  
  c->n.cycle		= &CYCLE;
  c->n.systems		= &NSYSTEMS;
  c->n.mols		= &NMOLS;
  c->n.sites		= &NSITES;
  c->n.types		= &NTYPES;
  c->n.volumes		= &NVOLUMES;
  c->n.swaps		= &NSWAPS;
  c->n.inserts		= &NINSERTS;
  c->n.cavinserts	= &NCAVINSERTS;
  c->n.blocks		= &NBLOCKS;
  c->n.cycles		= &NCYCLES;
  c->n.box		= &NBOX;
  c->n.rotation		= &NROTATION;
  c->n.reptation	= &NREPTATION;
  c->n.endbridge	= &NENDBRIDGE;
  c->n.rebridge		= &NREBRIDGE;
  c->n.bridges		= &NBRIDGES;
  c->n.fixed		= &NFIXED;
  c->n.semifixed	= &NSEMIFIXED;
  c->n.seed		= SEED;
  c->n.temper		= &NTEMPER;
  c->n.system		= &SYSTEM;
  c->n.sample		= &NSAMPLE;
  c->n.stretch		= &NSTRETCH;
  c->n.free		= &NFREE;
  c->n.ends		= &NENDS;
*/
  c->n.cycle		= &NCYCLE;		// link to numbers
  c->n.tape		= &ITAPE;
  c->n.systems		= &NSYSTEMS;
  c->n.mols		= &NMOLS;
  c->n.sites		= &NSITES;
  c->n.types		= &NTYPES;
  c->n.volchange	= &NVOLCHANGE;
  c->n.swap		= &NSWAP;
  c->n.displace		= &NDISPLACE;
  c->n.cbmc		= &NCBMC;

  return &VB_BridgeMap;
}

/*
#define VB_NLONG	2
#define VB_NDOUBLE	6
#define VB_NDIST	D_NDIST

long			**long_copy = NULL;
double			**double_copy = NULL;
diststruct		**dist_copy = NULL;
avstruct		*av_copy = NULL;
molstruct		*mol_copy = NULL;
vstruct			*v_copy = NULL;
wstruct			*w_copy = NULL;

void BridgeInitMemory(bridgestruct *bridge, long nsystems, long NMols)
{
  const char		module[32] = "varbridge",
			procedure[32] = "BridgeCreateCopy";
  long			i, nmols = 0;
  bridgestruct		*c = bridge;

  nmols			= nsystems*NMols;
  *c			= *Bridge();		// copy current bridge
  
  if (!long_copy)				// allocate copy space
  {
    if(!(long_copy = (long **) calloc(VB_NLONG, sizeof(long *))))
      Exit(module, procedure, "first long calloc error");
    for (i=0; i<VB_NLONG; ++i)
      if (!(long_copy[i] = (long *) calloc(nsystems, sizeof(long))))
        Exit(module, procedure, "second long calloc error");
  }
  if (!double_copy)
  {
    if(!(double_copy = (double **) calloc(VB_NDOUBLE, sizeof(double *))))
      Exit(module, procedure, "first double calloc error");
    for (i=0; i<VB_NDOUBLE; ++i)
      if (!(double_copy[i] = (double *) calloc(nsystems, sizeof(double))))
        Exit(module, procedure, "second double calloc error");
  }
  if (!dist_copy)
    if(!(dist_copy = (diststruct **) calloc(D_NDIST, sizeof(diststruct *))))
      Exit(module, procedure, "diststruct calloc error");
  if (!av_copy)
    if (!(av_copy = (avstruct *) calloc(nsystems, sizeof(avstruct))))
      Exit(module, procedure, "avstruct calloc error");
  if (!v_copy)
    if (!(v_copy = (vstruct *) calloc(nsystems, sizeof(vstruct))))
      Exit(module, procedure, "vstruct calloc error");
  if (!w_copy)
    if (!(w_copy = (wstruct *) calloc(nsystems, sizeof(wstruct))))
      Exit(module, procedure, "wstruct calloc error");
  if (!mol_copy)
    if (!(mol_copy = (molstruct *) calloc(nmols, sizeof(molstruct))))
      Exit(module, procedure, "molstruct calloc error");
  
  for (i=0; i<nsystems; ++i)			// set links to copy space
  {
    c->system[i].nmols	= long_copy[0]+i;
    c->system[i].nsites	= long_copy[1]+i;
    c->system[i].pres	= double_copy[0]+i;
    c->system[i].vol	= double_copy[1]+i;
    c->system[i].temp	= double_copy[2]+i;
    c->system[i].drmax	= double_copy[3]+i;
    c->system[i].dlmax	= double_copy[4]+i;
    c->system[i].damax	= double_copy[5]+i;

    c->av[i]		= av_copy+i;
    c->v[i]		= v_copy+i;
    c->vir[i]		= w_copy+i;
    c->mol[i]		= mol_copy+i*NMols;
  }
  ReinitSample(c->dist);
  for (i=0; i<D_NDIST; ++i)
    c->dist[i]		= dist_copy[i];
}


void BridgeResetVariable(bridgestruct *bridge)
{
  long			system;
  
  for (system=0; system<NSYSTEMS; ++system)
  {
    ResetAverages(bridge->av[system]);
    ResetAcceptance(bridge->av[system]);
  }
  ReinitSample(bridge->dist);
}


void BridgeCopyInstant(bridgestruct *to, bridgestruct *from, long system)
{
  long			i = system;
  
  *(to->system[i].nmols)	= *(from->system[i].nmols);
  *(to->system[i].nsites)	= *(from->system[i].nsites);
  *(to->system[i].pres)		= *(from->system[i].pres);
  *(to->system[i].vol)		= *(from->system[i].vol);
  *(to->system[i].temp)		= *(from->system[i].temp);
  *(to->system[i].drmax)	= *(from->system[i].drmax);
  *(to->system[i].dlmax)	= *(from->system[i].dlmax);
  *(to->system[i].damax)	= *(from->system[i].damax);
  to->system[i].a		= from->system[i].a;
  to->system[i].b		= from->system[i].b;
  to->system[i].c		= from->system[i].c;

  *(to->v[i])			= *(from->v[i]);
  *(to->vir[i])			= *(from->vir[i]);
  //to->av[i]->q_offset		= from->av[i]->q_offset;
  
  for (i=0; i<*(from->system[system].nmols); ++i)
    *(to->mol[system]+i)	= *(from->mol[system]+i);

  for (i=0; i<D_NDIST; ++i)
    if (from->dist[i]&&to->dist[i])
      for (system=0; system<*(from->n.systems); ++system)
        (to->dist[i]+system)->binsize = (from->dist[i]+system)->binsize;
}


void BridgeAddVariable(bridgestruct *to, bridgestruct *from)
{
  long			i, system;
  
  for (system=0; system<*(from->n.systems); ++system)
  {
    AddAvToAv(to->av[system], from->av[system]);
    for (i=0; i<D_NDIST; ++i)
      if (from->dist[i]&&to->dist[i])
        D_Add(to->dist[i]+system, from->dist[i]+system);
  }
}
*/
