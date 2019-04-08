/*
    program:	history.c
    author:	Peng Yi borrowed from Pieter J. in 't Veld
    date:	January 10, 2008
    purpose:	Binary history file module.
*/
#define __HISTORY_MODULE
#include "history.h"

int fputd(double *d, FILE *stream)
{
  fwrite(d, sizeof(double), 1, stream);
  return ferror(stream)!=0;
}


int fgetd(double *d, FILE *stream)
{
  fread(d, sizeof(double), 1, stream);
  return feof(stream)!=0;
}


int fputl(long *l, FILE *stream)
{
  fwrite(l, sizeof(long), 1, stream);
  return ferror(stream)!=0;
}


int fgetl(long *l, FILE *stream)
{
  fread(l, sizeof(long), 1, stream);
  return feof(stream)!=0;
}


int fput_vector(vector *r, FILE *stream)
{
  if (fputd(&(r->x), stream)) return 1;
  if (fputd(&(r->y), stream)) return 1;
  if (fputd(&(r->z), stream)) return 1;
  return 0;
}


int fget_vector(vector *r, FILE *stream)
{
  if (fgetd(&(r->x), stream)) return 1;
  if (fgetd(&(r->y), stream)) return 1;
  if (fgetd(&(r->z), stream)) return 1;
  return 0;
}


int fput_matrix(matrix *m, FILE *stream)
{
  if (fput_vector(&(m->x), stream)) return 1;
  if (fput_vector(&(m->y), stream)) return 1;
  if (fput_vector(&(m->z), stream)) return 1;
  return 0;
}


int fget_matrix(matrix *m, FILE *stream)
{
  if (fget_vector(&(m->x), stream)) return 1;
  if (fget_vector(&(m->y), stream)) return 1;
  if (fget_vector(&(m->z), stream)) return 1;
  return 0;
}


int fput_sphere(sphere *s, FILE *stream)
{
  if (fputd(&(s->d), stream)) return 1;
  if (fputd(&(s->alpha), stream)) return 1;
  if (fputd(&(s->beta), stream)) return 1;
  return 0;
}


int fget_sphere(sphere *s, FILE *stream)
{
  if (fgetd(&(s->d), stream)) return 1;
  if (fgetd(&(s->alpha), stream)) return 1;
  if (fgetd(&(s->beta), stream)) return 1;
  return 0;
}

/* More customized functions */

int fput_mol(molstruct *mol, FILE *stream)
{
  int			i;
  
  if (fputl(&(mol->nsites), stream)) return 1;
  if (fputl(&(mol->box), stream)) return 1;
  if (fputl(&(mol->fix), stream)) return 1;
  if (fputl(&(mol->flip), stream)) return 1;
  for (i=0; i<mol->nsites; ++i)
  {
    if (fputl(mol->type+i, stream)) return 1;
    if (fputl(mol->parent+i, stream)) return 1;
    if (fputl(mol->flags+i, stream)) return 1;
    if (fput_vector(mol->p+i, stream)) return 1;
  }
  return 0;
}
  

int fget_mol(molstruct *mol, FILE *stream, long version)
{
  int			i;
  
  if (fgetl(&(mol->nsites), stream)) return 1;
  if (fgetl(&(mol->box), stream)) return 1;
  if (fgetl(&(mol->fix), stream)) return 1;
  if (fgetl(&(mol->flip), stream)) return 1;
  if ((version<17)&&fget_vector(mol->p, stream)) return 1;
  for (i=0; i<mol->nsites; ++i)
  {
    if (fgetl(mol->type+i, stream)) return 1;
    if (fgetl(mol->parent+i, stream)) return 1;
    if (fgetl(mol->flags+i, stream)) return 1;
    if (fget_vector(mol->p+i, stream)) return 1;
  }
  return 0;
}


int fput_type(typestruct *type, FILE *stream)
{
  if (fputd(&(type->M), stream)) return 1;
  if (fputd(&(type->Q), stream)) return 1;
  if (fputd(&(type->SIGMA), stream)) return 1;
  if (fputd(&(type->EPSILON), stream)) return 1;
  if (fputd(&(type->LSTRETCH),stream)) return 1;
  if (fputd(&(type->KSTRETCH),stream)) return 1;
  if (fputd(&(type->KBENDING),stream)) return 1;
  if (fputd(&(type->THETA),stream)) return 1;
  if (fputd(type->TORSION,stream)) return 1;
  if (fputd(type->TORSION+1,stream)) return 1;
  if (fputd(type->TORSION+2,stream)) return 1;
  if (fputd(type->TORSION+3,stream)) return 1;
  if (fputd(type->TORSION+4,stream)) return 1;
  if (fputd(type->TORSION+5,stream)) return 1;
  return 0;
}


int fget_type(typestruct *type, FILE *stream, long version)
{
  if (fgetd(&(type->M), stream)) return 1;
  if (fgetd(&(type->Q), stream)) return 1;
  if (fgetd(&(type->SIGMA), stream)) return 1;
  if (fgetd(&(type->EPSILON), stream)) return 1;
  if (fgetd(&(type->LSTRETCH),stream)) return 1;
  if (fgetd(&(type->KSTRETCH),stream)) return 1;
  if (fgetd(&(type->KBENDING),stream)) return 1;
  if (fgetd(&(type->THETA),stream)) return 1;
  if (fgetd(type->TORSION,stream)) return 1;
  if (fgetd(type->TORSION+1,stream)) return 1;
  if (fgetd(type->TORSION+2,stream)) return 1;
  if (fgetd(type->TORSION+3,stream)) return 1;
  if (fgetd(type->TORSION+4,stream)) return 1;
  if (fgetd(type->TORSION+5,stream)) return 1;
  return 0;
}


int fput_system(systemstruct *system, FILE *stream)
{
  if (fputl(&(system->n), stream)) return 1;
  if (fputl(system->nmols, stream)) return 1;
  if (fputl(system->nsites, stream)) return 1;
  if (fputd(system->pres, stream)) return 1;
  if (fputd(system->vol, stream)) return 1;
  if (fputd(system->temp, stream)) return 1;
  if (fputd(system->drmax, stream)) return 1;
  if (fputd(system->dlmax, stream)) return 1;
  if (fputd(system->damax, stream)) return 1;
  /*
  if (fput_vector(&(system->a), stream)) return 1;
  if (fput_vector(&(system->b), stream)) return 1;
  if (fput_vector(&(system->c), stream)) return 1;
  */
  if (fputd(system->lx, stream)) return 1;
  if (fputd(system->ly, stream)) return 1;
  if (fputd(system->lz, stream)) return 1;
  return 0;
}


int fget_system(systemstruct *system, FILE *stream, long version)
{
  if (fgetl(&(system->n), stream)) return 1;
  if (fgetl(system->nmols, stream)) return 1;
  if (fgetl(system->nsites, stream)) return 1;
  if (fgetd(system->pres, stream)) return 1;
  if (fgetd(system->vol, stream)) return 1;
  if (fgetd(system->temp, stream)) return 1;
  if (fgetd(system->drmax, stream)) return 1;
  if (fgetd(system->dlmax, stream)) return 1;
  if (fgetd(system->damax, stream)) return 1;
  /*
  if (fget_vector(&(system->a), stream)) return 1;
  if (fget_vector(&(system->b), stream)) return 1;
  if (fget_vector(&(system->c), stream)) return 1;
  */
  if (fgetd(system->lx, stream)) return 1;
  if (fgetd(system->ly, stream)) return 1;
  if (fgetd(system->lz, stream)) return 1;

  return 0;
}


int fput_command(commandstruct *command, FILE *stream)
{
  if (fputl(command->hs, stream)) return 1;
  if (fputl(command->lj, stream)) return 1;
  if (fputl(command->ljshift, stream)) return 1;
  if (fputl(command->ljlrc, stream)) return 1;
  if (fputl(command->stretch, stream)) return 1;
  if (fputl(command->bending, stream)) return 1;
  if (fputl(command->torsion, stream)) return 1;
  if (fputl(command->virial, stream)) return 1;
  if (fputl(command->scalecut, stream)) return 1;
  if (fputl(command->nvt, stream)) return 1;
  if (fputl(command->npt, stream)) return 1;
  if (fputl(command->gibbs, stream)) return 1;
  if (fputl(command->mpi, stream)) return 1;
  if (fputl(command->density, stream)) return 1;
  if (fputl(command->energy, stream)) return 1;
  if (fputl(command->pressure, stream)) return 1;
  if (fputl(command->drift, stream)) return 1;
  if (fputl(command->torsional, stream)) return 1;
  if (fputl(command->bonda, stream)) return 1;
  if (fputl(command->bondl, stream)) return 1;
  if (fputl(command->radial, stream)) return 1;
  if (fputl(command->localp2, stream)) return 1;
  if (fputl(command->xtalsize, stream)) return 1;

/*
  if (fputl(command->hs, stream)) return 1;
  if (fputl(command->lj, stream)) return 1;
  if (fputl(command->coulomb, stream)) return 1;
  if (fputl(command->polymer, stream)) return 1;
  if (fputl(command->mpi, stream)) return 1;
  if (fputl(command->async, stream)) return 1;
  if (fputl(command->bias, stream)) return 1;
  if (fputl(command->nvt, stream)) return 1;
  if (fputl(command->npt, stream)) return 1;
  if (fputl(command->gibbs, stream)) return 1;
  if (fputl(command->insert, stream)) return 1;
  if (fputl(command->widom, stream)) return 1;
  if (fputl(command->canonical, stream)) return 1;
  if (fputl(command->cavity, stream)) return 1;
  if (fputl(command->density, stream)) return 1;
  if (fputl(command->density3d, stream)) return 1;
  if (fputl(command->densfree, stream)) return 1;
  if (fputl(command->denstalobr, stream)) return 1;
  if (fputl(command->orient, stream)) return 1;
  if (fputl(command->orientfree, stream)) return 1;
  if (fputl(command->orienttalobr, stream)) return 1;
  if (fputl(command->orientcorr01, stream)) return 1;
  if (fputl(command->orientcorr02, stream)) return 1;
  if (fputl(command->orientcorr03, stream)) return 1;
  if (fputl(command->orientcorr04, stream)) return 1;
  if (fputl(command->orientcorr05, stream)) return 1;
  if (fputl(command->orientcorrCR, stream)) return 1;
  if (fputl(command->tails_etc, stream)) return 1;
  if (fputl(command->radial, stream)) return 1;
  if (fputl(command->energy, stream)) return 1;
  if (fputl(command->hs_dens, stream)) return 1;
  if (fputl(command->temper, stream)) return 1;
  if (fputl(command->jacob, stream)) return 1;
  if (fputl(command->d_bridge, stream)) return 1;
  if (fputl(command->torsion, stream)) return 1;
  if (fputl(command->re_torsion, stream)) return 1;
  if (fputl(command->e_profile, stream)) return 1;
  if (fputl(command->n_profile, stream)) return 1;
  if (fputl(command->b_length, stream)) return 1;
  if (fputl(command->b_angle, stream)) return 1;
  if (fputl(command->ptorsion, stream)) return 1;
  if (fputl(command->loopreentry, stream)) return 1;
  if (fputl(command->virial, stream)) return 1;
  if (fputl(command->e_n_function, stream)) return 1;
  if (fputl(command->monodisperse, stream)) return 1;
  if (fputl(command->w_profile, stream)) return 1;
  if (fputl(command->stretch, stream)) return 1;
*/
  return 0;
}


int fget_command(commandstruct *command, FILE *stream, long version)
{
  long			extra_l;

  if (fgetl(command->hs, stream)) return 1;
  if (fgetl(command->lj, stream)) return 1;
  if (fgetl(command->ljshift, stream)) return 1;
  if (fgetl(command->ljlrc, stream)) return 1;
  if (fgetl(command->stretch, stream)) return 1;
  if (fgetl(command->bending, stream)) return 1;
  if (fgetl(command->torsion, stream)) return 1;
  if (fgetl(command->virial, stream)) return 1;
  if (fgetl(command->scalecut, stream)) return 1;
  if (fgetl(command->nvt, stream)) return 1;
  if (fgetl(command->npt, stream)) return 1;
  if (fgetl(command->gibbs, stream)) return 1;
  if (fgetl(command->mpi, stream)) return 1;
  if (fgetl(command->density, stream)) return 1;
  if (fgetl(command->energy, stream)) return 1;
  if (fgetl(command->pressure, stream)) return 1;
  if (fgetl(command->drift, stream)) return 1;
  if (fgetl(command->torsional, stream)) return 1;
  if (fgetl(command->bonda, stream)) return 1;
  if (fgetl(command->bondl, stream)) return 1;
  if (fgetl(command->radial, stream)) return 1;
  if (fgetl(command->localp2, stream)) return 1;
  if (fgetl(command->xtalsize, stream)) return 1;
/*  
  if (fgetl(command->hs, stream)) return 1;
  if (fgetl(command->lj, stream)) return 1;
  if (fgetl(command->coulomb, stream)) return 1;
  if (fgetl(command->polymer, stream)) return 1;
  if ((version>19)&&fgetl(command->mpi, stream)) return 1;
  if ((version>19)&&fgetl(command->async, stream)) return 1;
  if (fgetl(command->bias, stream)) return 1;
  if (fgetl(command->nvt, stream)) return 1;
  if (fgetl(command->npt, stream)) return 1;
  if (fgetl(command->gibbs, stream)) return 1;
  if (fgetl(command->insert, stream)) return 1;
  if (fgetl(command->widom, stream)) return 1;
  if (fgetl(command->canonical, stream)) return 1;
  if (fgetl(command->cavity, stream)) return 1;
  if (fgetl(command->density, stream)) return 1;
  if ((version>23)&&fgetl(command->density3d, stream)) return 1;
  if ((version>22)&&fgetl(command->densfree, stream)) return 1;
  if ((version>22)&&fgetl(command->denstalobr, stream)) return 1;
  if ((version>19)&&fgetl(command->orient, stream)) return 1;
  if ((version>22)&&fgetl(command->orientfree, stream)) return 1;
  if ((version>22)&&fgetl(command->orienttalobr, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr01, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr02, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr03, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr04, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorr05, stream)) return 1;
  if ((version>22)&&fgetl(command->orientcorrCR, stream)) return 1;
  if (fgetl(command->tails_etc, stream)) return 1;
  if (fgetl(command->radial, stream)) return 1;
  if (fgetl(command->energy, stream)) return 1;
  if (fgetl(command->hs_dens, stream)) return 1;
  if (fgetl(command->temper, stream)) return 1;
  if (version<2)
  {
    if (!fgetl(&extra_l, stream)) return 1;
    if (!fgetl(&extra_l, stream)) return 1;
    if (!fgetl(&extra_l, stream)) return 1;
    return fgetl(&extra_l, stream);
  }
  if (version<7)
    return 0;
  if (fgetl(command->jacob, stream)) return 1;
  if (fgetl(command->d_bridge, stream)) return 1;
  if (version<8)
    return 0;
  if (fgetl(command->torsion, stream)) return 1;
  if ((version>24)&&fgetl(command->re_torsion, stream)) return 1;
  if (fgetl(command->e_profile, stream)) return 1;
  if (fgetl(command->n_profile, stream)) return 1;
  if (version<9)
    return 0;
  if (fgetl(command->b_length, stream)) return 1;
  if (fgetl(command->b_angle, stream)) return 1;
  if ((version>9)&&fgetl(command->ptorsion, stream)) return 1;
  if ((version>10)&&fgetl(command->loopreentry, stream)) return 1;
  if ((version>11)&&fgetl(command->virial, stream)) return 1;
  if ((version>13)&&fgetl(command->e_n_function, stream)) return 1;
  if ((version>17)&&fgetl(command->monodisperse, stream)) return 1;
  if ((version>21)&&fgetl(command->w_profile, stream)) return 1;
  if (fgetl(command->stretch, stream)) return 1;
*/
  return 0;
}


int fput_n(nstruct *n, FILE *stream)
{
  if (fputl(n->cycle, stream)) return 1;
  if (fputl(n->tape, stream)) return 1;
  if (fputl(n->systems, stream)) return 1;
  if (fputl(n->mols, stream)) return 1;
  if (fputl(n->sites, stream)) return 1;
  if (fputl(n->types, stream)) return 1;
  if (fputl(n->volchange, stream)) return 1;
  if (fputl(n->swap, stream)) return 1;
  if (fputl(n->displace, stream)) return 1;
  if (fputl(n->cbmc, stream)) return 1;
/*
  if (fputl(n->cycle, stream)) return 1;
  if (fputl(n->systems, stream)) return 1;
  if (fputl(n->mols, stream)) return 1;
  if (fputl(n->sites, stream)) return 1;
  if (fputl(n->types, stream)) return 1;
  if (fputl(n->volumes, stream)) return 1;
  if (fputl(n->swaps, stream)) return 1;
  if (fputl(n->inserts, stream)) return 1;
  if (fputl(n->cavinserts, stream)) return 1;
  if (fputl(n->blocks, stream)) return 1;
  if (fputl(n->cycles, stream)) return 1;
  if (fputl(n->box, stream)) return 1;
  if (fputl(n->rotation, stream)) return 1;
  if (fputl(n->reptation, stream)) return 1;
  if (fputl(n->endbridge, stream)) return 1;
  if (fputl(n->rebridge, stream)) return 1;
  if (fputl(n->bridges, stream)) return 1;
  if (fputl(n->fixed, stream)) return 1;
  if (fputl(n->semifixed, stream)) return 1;
  if (fputl(n->seed, stream)) return 1;
  if (fputl(n->temper, stream)) return 1;
  if (fputl(n->sample, stream)) return 1;
  if (fputl(n->stretch, stream)) return 1;
  if (fputl(n->free, stream)) return 1;
  if (fputl(n->ends, stream)) return 1;
  if (fputl(n->system, stream)) return 1;
*/
  return 0;
}


int fget_n(nstruct *n, FILE *stream, long version)
{
  long			extra_l;

  if (fgetl(n->cycle, stream)) return 1;
  if (fgetl(n->tape, stream)) return 1;
  if (fgetl(n->systems, stream)) return 1;
  if (fgetl(n->mols, stream)) return 1;
  if (fgetl(n->sites, stream)) return 1;
  if (fgetl(n->types, stream)) return 1;
  if (fgetl(n->volchange, stream)) return 1;
  if (fgetl(n->swap, stream)) return 1;
  if (fgetl(n->displace, stream)) return 1;
  if (fgetl(n->cbmc, stream)) return 1;
/*  if (fgetl(n->cycle, stream)) return 1;
  if (fgetl(n->systems, stream)) return 1;
  if (fgetl(n->mols, stream)) return 1;
  if (fgetl(n->sites, stream)) return 1;
  if (fgetl(n->types, stream)) return 1;
  if (fgetl(n->volumes, stream)) return 1;
  if (fgetl(n->swaps, stream)) return 1;
  if (fgetl(n->inserts, stream)) return 1;
  if (fgetl(n->cavinserts, stream)) return 1;
  if (fgetl(n->blocks, stream)) return 1;
  if (fgetl(n->cycles, stream)) return 1;
  if (fgetl(n->box, stream)) return 1;
  if (fgetl(n->rotation, stream)) return 1;
  if (fgetl(n->reptation, stream)) return 1;
  if (fgetl(n->endbridge, stream)) return 1;
  if (fgetl(n->rebridge, stream)) return 1;
  if ((version<=19)&&fgetl(&extra_l, stream)) return 1; // loops
  if ((version<=19)&&fgetl(&extra_l, stream)) return 1; // tails
  if (fgetl(n->bridges, stream)) return 1;
  if (fgetl(n->fixed, stream)) return 1;
  if ((version>20)&&fgetl(n->semifixed, stream)) return 1;
  if (fgetl(n->seed, stream)) return 1;
  if (fgetl(n->temper, stream)) return 1;
  if (version<2)
  {
    if (fgetl(&extra_l, stream)) return 1;
    if (fgetl(&extra_l, stream)) return 1;
    if (fgetl(&extra_l, stream)) return 1;
    return fgetl(&extra_l, stream);
  }
  if (fgetl(n->sample, stream)) return 1;
  if (version<9) 
    return 0;
  if (fgetl(n->stretch, stream)) return 1;
  if (version<15) 
    return 0;
  if (fgetl(n->free, stream)) return 1;
  if (fgetl(n->ends, stream)) return 1;
  if ((version>19)&&fgetl(n->system, stream)) return 1;
*/
  return 0;
}


int fput_av(avstruct *av, FILE *stream)
{
  if (fputl(&(av->move), stream)) return 1;
  if (fputl(&(av->acc_move), stream)) return 1;
  if (fputl(&(av->rot), stream)) return 1;
  if (fputl(&(av->acc_rot), stream)) return 1;
  if (fputl(&(av->vol), stream)) return 1;
  if (fputl(&(av->acc_vol), stream)) return 1;
  if (fputl(&(av->cbmc), stream)) return 1;
  if (fputl(&(av->acc_cbmc), stream)) return 1;
  if (fputl(&(av->rep), stream)) return 1;
  if (fputl(&(av->acc_rep), stream)) return 1;
  if (fputl(&(av->seq), stream)) return 1;
  if (fputl(&(av->acc_seq), stream)) return 1;
  return 0;
}


int fget_av(avstruct *av, FILE *stream, long version)
{
  if (fgetl(&(av->move), stream)) return 1;
  if (fgetl(&(av->acc_move), stream)) return 1;
  if (fgetl(&(av->rot), stream)) return 1;
  if (fgetl(&(av->acc_rot), stream)) return 1;
  if (fgetl(&(av->vol), stream)) return 1;
  if (fgetl(&(av->acc_vol), stream)) return 1;
  if (fgetl(&(av->cbmc), stream)) return 1;
  if (fgetl(&(av->acc_cbmc), stream)) return 1;
  if (fgetl(&(av->rep), stream)) return 1;
  if (fgetl(&(av->acc_rep), stream)) return 1;
  if (fgetl(&(av->seq), stream)) return 1;
  if (fgetl(&(av->acc_seq), stream)) return 1;
  return 0;
}


int fput_v(vstruct *v, FILE *stream)
{
  if (fputd(&(v->tot), stream)) return 1;
  if (fputd(&(v->bonded), stream)) return 1;
  if (fputd(&(v->nonbonded), stream)) return 1;
  if (fputd(&(v->hs), stream)) return 1;
  if (fputd(&(v->lj), stream)) return 1;
  if (fputd(&(v->ljcorr), stream)) return 1;
  if (fputd(&(v->corr), stream)) return 1;
  if (fputd(&(v->stretch), stream)) return 1;
  if (fputd(&(v->bending), stream)) return 1;
  if (fputd(&(v->torsion), stream)) return 1;
  return 0;
}


int fget_v(vstruct *v, FILE *stream, long version)
{
  if (fgetd(&(v->tot), stream)) return 1;
  if (fgetd(&(v->bonded), stream)) return 1;
  if (fgetd(&(v->nonbonded), stream)) return 1;
  if (fgetd(&(v->hs), stream)) return 1;
  if (fgetd(&(v->lj), stream)) return 1;
  if (fgetd(&(v->ljcorr), stream)) return 1;
  if (fgetd(&(v->corr), stream)) return 1;
  if (fgetd(&(v->stretch), stream)) return 1;
  if (fgetd(&(v->bending), stream)) return 1;
  if (fgetd(&(v->torsion), stream)) return 1;
  return 0;
}


int fput_w(wstruct *w, FILE *stream)
{
  if (fputd(&(w->tot), stream)) return 1;
  if (fputd(&(w->lj), stream)) return 1;
  if (fputd(&(w->stretch), stream)) return 1;
  if (fputd(&(w->torsion), stream)) return 1;
  return 0;
}


int fget_w(wstruct *w, FILE *stream, long version)
{
  if (fgetd(&(w->tot), stream)) return 1;
  if (fgetd(&(w->lj), stream)) return 1;
  if (fgetd(&(w->stretch), stream)) return 1;
  if (fgetd(&(w->torsion), stream)) return 1;
  return 0;
}


int fput_dist(diststruct *dist, int top, FILE *stream)
{
  long                  i;

  if (fputl(&(dist->nbins), stream)) return 1;
  if (fputl(&(dist->n), stream)) return 1;
  if (fputl(&(dist->ncount), stream)) return 1;
  if (fputl(&(dist->startbin), stream)) return 1;
  if (fputl(&(dist->level), stream)) return 1;
  if (top)
    for (i=0; i<=dist->level; ++i)
      if (fputd(dist->binsize+i, stream)) return 1;

  if (dist->level)				// distribution fork
  {
    for (i=0; i<dist->nbins; ++i)
      if (fput_dist(dist->dist+i, FALSE, stream)) return 1;
  }
  else if (dist->nbins)				// data fork
  {
    for (i=0; i<D_NAVERAGE; ++i)
      if (fputd(dist->average+i, stream)) return 1;
    for (i=0; i<dist->nbins; ++i)
    {
      if (fputl(dist->bin+i, stream)) return 1;
      if (fputd(dist->data+i, stream)) return 1;
      if (fputd(dist->cweight+i, stream)) return 1;
    }
  }
  return 0;
}


int fget_dist(diststruct *dist, double *binsize, FILE *stream, long version)
{
  long                  i;

  if (fgetl(&(dist->nbins), stream)) return 1;
  if (fgetl(&(dist->n), stream)) return 1;
  if ((version>19)&&fgetl(&(dist->ncount), stream)) return 1;
  if (fgetl(&(dist->startbin), stream)) return 1;

  if (version<24) dist->level = 0;
  else if (fgetl(&(dist->level), stream)) return 1;
  if (!binsize)
  {
    for (i=0; i<=dist->level; ++i)
      if (fgetd(dist->binsize+i, stream)) return 1;
  }
  else dist->binsize	= binsize;

  if (dist->level)				// distribution fork
  {
    if (dist->nbins)
    {
      if (!(dist->dist = calloc(dist->nbins, sizeof(diststruct))))
	Exit("history", "fget_dist", "distribution fork calloc error");
      for (i=0; i<dist->nbins; ++i)
	if (fget_dist(dist->dist+i, dist->binsize+1, stream, version)) return 1;
    }
  }
  else						// data fork
  {
    if ((version<24)||dist->nbins)
    {
      if (!(dist->average = (double *) calloc(D_NAVERAGE, sizeof(double))))
        Exit("history", "fget_dist", "average calloc error");
      for (i=0; i<D_NAVERAGE; ++i)
        if (fgetd(dist->average+i, stream)) return 1;
    }
    if (dist->nbins)
    {
      if (!((dist->bin = (long *) calloc(dist->nbins, sizeof(long)))&&
            (dist->data = (double *) calloc(dist->nbins, sizeof(double)))&&
            (dist->cweight = (double *) calloc(dist->nbins, sizeof(double)))))
        Exit("history", "fget_dist", "data fork calloc error");
      for (i=0; i<dist->nbins; ++i)
      {
        if (fgetl(dist->bin+i, stream)) return 1;
        if (fgetd(dist->data+i, stream)) return 1;
        if (fgetd(dist->cweight+i, stream)) return 1;
      }
    }
  }
  return 0;
}


int fput_all_dist(bridgestruct *c, FILE *stream)
{
  long			i, system;

  for (i=0; i<D_NDIST; ++i) 
    for (system=0; c->dist[i]&&(system<NSYSTEMS); ++system)
      if (fput_dist(c->dist[i]+system, TRUE, stream)) return 1;
  return 0;
}


int fget_all_dist(bridgestruct *c, FILE *stream, long version)
{
  long			i, system;

  if (version<8) return 0;
  ReinitSample(c->dist);
  for (i=0; i<D_NDIST; ++i)
    for (system=0; c->dist[i]&&(system<*(c->n.systems)); ++system)
    {
      D_Reset(c->dist[i]+system);
      if (fget_dist(c->dist[i]+system, NULL, stream, version)) return 1;
    }
  return 0;
}


int fput_all(bridgestruct *c, FILE *stream)
{
  int			flag, i, j;
  
//  if (!(flag = fputl(c->d_type, stream))) {}
  if (!(flag = fput_n(&(c->n), stream)))			// 0
      
    flag		= fput_command(&(c->command), stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 1
    flag		= fput_system(c->system+i, stream);
  for (i=0; (i<*(c->n.types))&&!flag; ++i)			// 2
    flag		= fput_type(c->type+i, stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 3
    for (j=0; (j<*(c->system[i].nmols))&&!flag; ++j)
      flag		= fput_mol(c->mol[i]+j, stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 4
    flag		= fput_av(c->av[i], stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 5
    flag		= fput_v(c->v[i], stream);
  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 6
    flag		= fput_w(c->vir[i], stream);
  if (!flag) 							// 7
    flag 		= fput_all_dist(c, stream);
  return flag;
}


int fget_all(bridgestruct *c, FILE *stream, long version)
{
  int			flag = 0, i = 0;
  long			cnt = 0, debug = 0;	// cnt keeps track of where
						// read in fails
  long			data;
  
//  if (version>19) flag = fgetl(c->d_type, stream);
//  if (!flag)
  flag = fget_n(&(c->n), stream, version);			// 0
  if (!flag)
    flag		= fget_command(&(c->command), stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 1
    flag		= fget_system(c->system+i, stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.types))&&!flag; ++i)			// 2
    flag		= fget_type(c->type+i, stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.mols))&&!flag; ++i)			// 3
    flag		= fget_mol(c->mol[0]+i, stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 4
    flag		= fget_av(c->av[i], stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 5
    flag		= fget_v(c->v[i], stream, version);
  if (!flag) ++cnt;

  for (i=0; (i<*(c->n.systems))&&!flag; ++i)			// 6
    flag		= fget_w(c->vir[i], stream, version);
  if (!flag) ++cnt;

  if (!flag)							// 7
    flag		= fget_all_dist(c, stream, version);
  if (debug&&flag)
    fprintf(stderr, "history::fget_all: read error at identity #%ld, i = %d\n", cnt, i);
 
/* 
  if (!flag) {				// copy some values to corresponding system variable
    for (i=0; i<*(c->n.systems); i++) {
       BOX[i].lx	=	c->system[i].a.x;
       BOX[i].ly	=	c->system[i].b.y;
       BOX[i].lz	=	c->system[i].c.z;
    }
  }
*/
/*
  if (!flag)
  {
    U_SYSTEM		= TRUE;
    D_TYPE		= NTYPES+1;
    SetBoxVectors(c->system->a, c->system->b, c->system->c);
    CalcUnits(-1);
    if (version<16) AllCartesian();
    InitForceField();
    if (version<15)
    {
      for (i=0; i<*(c->n.systems); ++i)
      {
        av[i].v		-= v[i].corr*av[i].total/NActive[i];
        av[i].vnb	-= v[i].corr*av[i].total/NActive[i];
        v[i].total	-= v[i].corr;
        v[i].nonbonded	-= v[i].corr;
      }
      CalcVCorr();
      for (i=0; i<*(c->n.systems); ++i)
      {
        av[i].v		+= v[i].corr*av[i].total/NActive[i];
        av[i].vnb	+= v[i].corr*av[i].total/NActive[i];
        v[i].total	+= v[i].corr;
        v[i].nonbonded	+= v[i].corr;
      }
    }   
  }
*/
  return flag;
}


FILE *fcreate_hist(char *name)
{
  FILE			*fp;
  long			version	= HIST_VERSION, i;
  char			ident[5] = HIST_IDENT;
  
  if (fp = fopen(name, "w"))
  {
    for (i=0; i<4; ++i)
      fputc((int) ident[i], fp);
    fputl(&version, fp);
    return fp;
  }
  return NULL;
}


FILE *fread_hist(char *name, long *version)
{
  char			c, ident[5] = HIST_IDENT;
  long			i, flag = 1;
  FILE			*fp;
  
  if (!name[0])
    fp			= stdin;
  else 
    fp			= fopen(name, "r");
  if (fp)
  {
    for (i=0; (i<4)&&flag; ++i)
    {
      c			= fgetc(fp);
      flag		= (c==ident[i]) ? 1 : 0;
    }
    if (flag)
    {
      fgetl(version, fp);
      if (*version<=HIST_VERSION)
        return fp;
    }
  }
  return NULL;
}


bridgestruct		*VarBridge = NULL;

int StartNewHistoryFile(char *name, long flag_mpi)
{
  FILE			*fp;
  
  if (fp = fcreate_hist(name))
  {
    VarBridge		= BridgeMap(flag_mpi);
    fclose(fp);
    return 1;
  }
  return 0;
}


int H_StoreBridge(char *name, bridgestruct *bridge)
{
  int			flag;
  FILE			*fp;

  if (fp = fopen(name, "a"))
  {
    flag		= fput_all(bridge, fp);
    fclose(fp);
    return flag;
  }
  return 0;
}


int H_StoreCurrent(char *name)
{
  return H_StoreBridge(name, VarBridge);
}


int H_GetHeader(char *name, long *version)
{
  long			flag, i;
  FILE			*fp;

  if (fp = fread_hist(name, version))
  {
    VarBridge		= BridgeMap(TRUE);
    flag		= 1;
    if ((*version>19)&&fgetl(VarBridge->d_type, fp)) flag = 0;
    if (flag&&fget_n(&(VarBridge->n), fp, *version)) flag = 0;
    if (flag&&fget_command(&(VarBridge->command), fp, *version)) flag = 0;
    for (i=0; (i<*(VarBridge->n.systems))&&flag; ++i)
      flag		= !fget_system(VarBridge->system+i, fp, *version);
    fclose(fp);
    return flag ? 0 : 1;
  }
  return 1;
}
  

FILE *H_GetFirst(char *name, long *version, long flag_mpi)
{
  FILE			*fp;
  
  if (fp = fread_hist(name, version))
  {
    VarBridge		= BridgeMap(flag_mpi);
    if (!fget_all(VarBridge, fp, *version))
    {
      VarBridge		= BridgeMap(flag_mpi);
      return fp;
    }
  }
  return NULL;
}


bridgestruct *H_Bridge()		// not used
{
  return VarBridge;
}


int H_GetNext(FILE *fp, long version)
{
  if (VarBridge)
    fget_all(VarBridge, fp, version);
  return feof(fp)!=0;
}

/*
void H_InitSystem(char *argv[])		// not used
{
  InitUnits();
  InitMove();
  //InitEnsembles();
  //InitForceField();
  InitCavity(argv);
  InitSample(argv);
}
*/
