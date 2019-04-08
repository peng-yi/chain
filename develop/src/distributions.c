/*
    program:	distributions.c
    author:	Pieter J. in 't Veld for UT at Austin
    date:	January 27, June 11, 1998, April 1, 1999.
    purpose:	Definition of a multi-dimensional self-defining distribution
		formalism

    Notes:
      19980127	Creation date
      19980611	- Adaption to pointers
		- Alteration of energy functions
		- Alteration of distribution storage functions
		- Transfer to module style
      19990401	Change in finding minimum: Newton-like root finding algorithm
		through local linearization of the scaled force (see notes).
      20040306	Addition of multi-dimensionality

      20071014  Perfect multi-dimensionality, basically change InitDist() --Peng
*/
#define __DISTRIBUTIONS_MODULE
#include "distributions.h"

void D_Allocate(diststruct *d, long n)
{
  long			i, i0 = 0, dnbins = 0;
  long			*newbin;
  double		*newdata, *newcweight;
  diststruct		*newdist;

  if (!d->nbins) 
    d->startbin		= n;
  n			-= d->startbin;
  if (!(dnbins = n<0 ? i0 = -n : n>=d->nbins ? n+1-d->nbins : 0)) 
    return;
  if (d->level)
  {
    if (!(newdist = (diststruct *) calloc(d->nbins+dnbins, sizeof(diststruct))))
      Exit("distributions", "D_AllocateDist", "calloc error");
    for (i=0; i<d->nbins; ++i)
      newdist[i+i0]	= d->dist[i];
    i0			= n<0 ? 0 : d->nbins;
    for (i=i0; i<i0+dnbins; ++i)
    {
      newdist[i].binsize= d->binsize+1;
      newdist[i].level	= d->level-1;
    }
    if (d->dist) free(d->dist);	
    d->dist		= newdist;
  }
  else
  {
    if ((!d->average)&&
	  !(d->average = (double *) calloc(D_NAVERAGE, sizeof(double))))
      Exit("distributions", "D_AllocateData", "calloc error");
    if (!((newbin = (long *) calloc(d->nbins+dnbins, sizeof(long)))&&
	  (newdata = (double *) calloc(d->nbins+dnbins, sizeof(double)))&&
	  (newcweight = (double *) calloc(d->nbins+dnbins, sizeof(double)))))
      Exit("distributions", "D_AllocateData", "calloc error");
    for (i=0; i<d->nbins; ++i)
    {
      newbin[i+i0]	= d->bin[i];
      newdata[i+i0]	= d->data[i];
      newcweight[i+i0]	= d->cweight[i];
    }
    if (d->bin) free(d->bin);
    if (d->data) free(d->data);
    if (d->cweight) free(d->cweight);
    d->bin		= newbin;
    d->data		= newdata;
    d->cweight		= newcweight;
  }
  d->nbins		+= dnbins;
  if (n<0)
    d->startbin		+= n;
}


void D_Reset(diststruct *d)
{
  long			i;

  if (d->level)
  {
    if (d->dist)
    {
      for (i=0; i<d->nbins; ++i)
	D_Reset(d->dist+i);
      free(d->dist);
    }
  }
  else
  {
    if (d->bin) free(d->bin);
    if (d->data) free(d->data);
    if (d->cweight) free(d->cweight);
    if (d->average) free(d->average);
  }
  d->n			= 0;
  d->nbins		= 0;
  d->ncount		= 0;
  d->startbin		= 0;
  d->dist		= NULL;
  d->bin		= NULL;
  d->data		= NULL;
  d->cweight		= NULL;
  d->average		= NULL;
}


void D_CopyRecursive(
  diststruct *dest, diststruct *src, double *binsize, int minus)
{
  long			i;
  
  dest->nbins		= src->nbins;
  dest->startbin	= src->startbin;
  dest->ncount		= minus ? -src->ncount : src->ncount;
  dest->n		= minus ? -src->n : src->n;
  if (!binsize)
  {
    if (!(dest->binsize||
	 (dest->binsize = (double *) calloc(src->level+1, sizeof(double)))))
      Exit("distributions", "D_CopyRecursive", "binsize calloc error");
    for (i=0; i<=src->level; ++i)
      dest->binsize[i] = src->binsize[i];
  }
  if (dest->level = src->level)			// distribution fork
  {
    if (src->nbins)
    {
      if (!(dest->dist = (diststruct *) calloc(src->nbins, sizeof(diststruct))))
        Exit("distributions", "D_CopyRecursive", "dist fork calloc error");
      for (i=0; i<src->nbins; ++i)
        D_CopyRecursive(dest->dist+i, src->dist+i, dest->binsize+1, minus);
    }
  }
  else if (src->nbins)				// data fork
  {
    if (!((dest->bin = (long *) calloc(src->nbins, sizeof(long)))&&
	  (dest->data = (double *) calloc(src->nbins, sizeof(double)))&&
	  (dest->cweight = (double *) calloc(src->nbins, sizeof(double)))&&
	  (dest->average = (double *) calloc(D_NAVERAGE, sizeof(double)))))
      Exit("distributions", "D_CopyRecursive", "data fork calloc error");
    for (i=0; i<D_NAVERAGE; ++i)
      dest->average[i]	= minus ? -src->average[i] : src->average[i];
    for (i=0; i<src->nbins; ++i)
    {
      dest->bin[i]	= minus ? -src->bin[i] : src->bin[i];
      dest->data[i]	= minus ? -src->data[i] : src->data[i];
      dest->cweight[i]	= minus ? -src->cweight[i] : src->cweight[i];
    }
  }
}


void D_Copy(diststruct *dest, diststruct *src)
{
  D_Reset(dest);
  D_CopyRecursive(dest, src, NULL, FALSE);
}


void D_Add(diststruct *dest, diststruct *src)
{
  long			i, n;

  if (!dest->nbins)
  {
    D_Reset(dest);
    D_CopyRecursive(dest, src, NULL, FALSE);
    return;
  }
  for (i=0; i<src->nbins; ++i)
  {
    D_Allocate(dest, i+src->startbin);
    n			= i+src->startbin-dest->startbin;
    if (src->level)
      D_Add(dest->dist+n, src->dist+i);
    else
    {
      dest->bin[n]	+= src->bin[i];
      dest->data[n]	+= src->data[i];
      dest->cweight[n]	+= src->cweight[i];
    }
  }
  if (dest->average&&src->average)
    for (i=0; i<D_NAVERAGE; ++i)
      dest->average[i]+= src->average[i];
  dest->n		+= src->n;
  dest->ncount		+= src->ncount;
}


void D_Subtr(diststruct *dest, diststruct *src)
{
  long			i, n;

  if (!dest->nbins)
  {
    D_Reset(dest);
    D_CopyRecursive(dest, src, NULL, TRUE);
    return;
  }
  for (i=0; i<src->nbins; ++i)
  {
    D_Allocate(dest, i+src->startbin);
    n			= i+src->startbin-dest->startbin;
    if (src->level)
      D_Add(dest->dist+n, src->dist+i);
    else
    {
      dest->bin[n]	-= src->bin[i];
      dest->data[n]	-= src->data[i];
      dest->cweight[n]	-= src->cweight[i];
    }
  }
  if (dest->average&&src->average)
    for (i=0; i<D_NAVERAGE; ++i)
      dest->average[i]-= src->average[i];
  dest->n		-= src->n;
  dest->ncount		-= src->ncount;
}


void D_Submit(diststruct *d, double *x, double *y, double *weight)
{
  register long		n = floor(*x/d->binsize[0])-d->startbin;
  register double	*avg = (double *) d->level;
  double		z;

  if ((n<0)||(n>=d->nbins))
  {
    D_Allocate(d, n += d->startbin);
    n			-= d->startbin;
  }
  while (d->level)				// advance to last level
    if (((n = floor(*(++x)/(d = d->dist+n)->binsize[0])-d->startbin)<0)||
	 (n>=d->nbins))
    {
      D_Allocate(d, n += d->startbin);
      n			-= d->startbin;
    }
  d->data[n]		+= *y * *weight;
  d->cweight[n]		+= *weight;
  ++(d->bin[n]);
  ++d->n;
  if ((!avg)&&(avg = d->average))		// only average for 1D
  {
    z			= *x;
    for (n=0; n<D_NAVERAGE; ++n)
    {
      *(avg++)		+= z;
      z			*= *x;
    }
  }
}


typedef
  struct {
    double		avg[D_NAVERAGE];
    long		ntotal, nbins;
  } distheaderstruct;

//#define DISTTEST

#ifdef DISTTEST
void D_Header(
  diststruct *d, int dist_type, distheaderstruct *header, long ncount)
{
  long			i, j;
  double		x, z;

  if (!d) return;
  if (d==D_Density3D)
    fprintf(stderr, "");
  if (dist_type&&(dist_type<3))
  {
    if (d->level)
    {
      for (i=0; i<d->nbins; ++i)
	D_Header(d->dist+i, dist_type, header, ncount);
      return;
    }
    header->nbins	+= d->nbins;
    header->ntotal	+= d->n;
    for (i=0; i<d->nbins; ++i)
    {
      switch (dist_type) {
        case 1: x = z = d->cweight[i] ? d->data[i]/d->cweight[i] : 0.0; break;
        case 2: x = z = ncount ? d->data[i]/ncount : 0.0; break;
      }         
      for (j=0; j<D_NAVERAGE; ++j)
      {
	header->avg[j]	+= z;
	z		*= x;
      }
    }	 
    return;
  }
  if (!d->average) return;
  for (i=0; i<D_NAVERAGE; ++i)
    header->avg[i]	= d->average[i]*d->nbins/(d->n ? d->n : 1.0);
  header->nbins		= d->nbins;
  header->ntotal	= d->n;
}


void D_PrintHeader(diststruct *d, int dist_type, double *L)
{
  char			s[40];
  long			i, n = d->level+1;
  double		x[n], stddev = 0;
  distheaderstruct	header;

  if (!d) return;
  if (d==D_Density3D)
    fprintf(stderr, "");
  for (i=0; i<D_NAVERAGE; ++i)
    header.avg[i]	= 0.0;
  header.nbins		= 0;
  header.ntotal		= 0;
  D_Header(d, dist_type, &header, d->ncount);
  printf("/* %s */\n\n", d->header);
  for (i=0; i<n; ++i)
    x[i]		= d->binsize[i]* (L ? *L : 1.0);
  PrintDVar("BINSIZE", n, x);
  PrintLVar("NCOUNT", 1, &d->ncount);
  PrintLVar("NBINS", 1, &header.nbins);
  PrintLVar("NTOTAL", 1, &header.ntotal);
  if (header.nbins)
  {
    for (i=0; i<D_NAVERAGE; ++i)
      header.avg[i]	/= header.nbins;
    stddev		= sqrt(header.avg[1]-header.avg[0]*header.avg[0]);
  }
  PrintDVar("STDDEV", 1, &stddev);
  PrintDVar("AVERAGE", D_NAVERAGE, header.avg);
  printf("\n");
  for (i=0; i<=d->level; ++i)
  {
    sprintf(s, "x%d", i);
    if (i) printf(", %12.12s", s);
    else printf("%12.12s", s);
  }
  printf(", %12.12s\n", "y");
  PrintLine(14*(d->level+2)-2);
}
#endif /* DISTTEST */


void D_PrintData(
  diststruct *d, int dist_type, double *x, double *L, long ncount, long nlevels)
{
  long			i, j;
  double		y, binsize = *(d->binsize)*(L ? *L : 1.0);

  for (i=0; i<d->nbins; ++i)
  {
/*    printf("nbins=%d\n", d->nbins);
    printf("startbin=%d\n", d->startbin);
    printf("level=%d\n", d->level);
    printf("i+d->startbin=%d\n", i+d->startbin);
    printf("binsize=%f\n", binsize);
    printf("binsize*(i+d->startbin)=%f\n", binsize*((double)(i+d->startbin)));
    printf("%d\t%f\n", nlevels-1-d->level, binsize*(i+(d->startbin)));
*/
    x[nlevels-1-d->level] = binsize*(i+(d->startbin));

    if (d->level)
      D_PrintData(d->dist+i, dist_type, x, L, ncount, nlevels);
    else
    {
      for (j=0; j<nlevels; ++j)
	printf("%s%12.6g", j ? ", " : "", x[j]);
      if (dist_type)
      {
	switch (dist_type) {
	  case 1: y = d->cweight[i] ? d->data[i]/d->cweight[i] : 0.0; break;
          case 2: y = ncount ? d->data[i]/ncount : 0.0; break;
	  case 3: y = d->data[i]; break;
	}
        printf(", %12.6g\n", y);
      }
      else printf(", %12ld\n", d->bin[i]);
    }
  }
}


void D_Print(diststruct *d, int dist_type, double *L)
{
  double		*x;

  if (!d) return;
//  D_PrintHeader(d, dist_type, L);
  if (!(x = (double *) calloc(d->level+1, sizeof(double))))
    Exit("distributions", "D_Print", "calloc error");
  D_PrintData(d, dist_type, x, L, d->ncount, d->level+1);
  printf("\n");
  free(x);
}

#ifdef DISTTEST
void D_PrintMath(diststruct *d, int dist_type, double *L)
{
  char			s[40], mem[256];
  long			i, n = d->n ? d->n : 1;
  double		*avg = d->average;
  
  sprintf(mem, "{\nnbins -> %d, binsize -> %s, startbin -> %d, ntotal -> %d, ",
    d->nbins, DToMath(s, *(d->binsize)* (L ? *L : 1.0)), d->startbin, d->n);
  PrintToMath(mem, FALSE);
  sprintf(mem, "ncount -> %d", d->ncount); PrintToMath(mem, FALSE);
  if (d->average)
  {
    //D_Average(d, dist_type);
    for (i=0; i<D_NAVERAGE; ++i)
    {
      if (i) sprintf(mem, ", avg%d -> %s", i+1, DToMath(s, avg[i]/n));
      else sprintf(mem, ", avg -> %s", DToMath(s, avg[i]/n));
      PrintToMath(mem, FALSE);
    }
  }
  if (d->level)
  {
    PrintToMath(", dist -> {", FALSE);
    for (i=0; i<d->nbins; ++i)
      D_PrintMath(d->dist+i, dist_type, L);
    PrintToMath("}", FALSE);
  }
  else
  {
    PrintToMath(", data -> \n{\n", FALSE);
    for (i=0; i<d->nbins; ++i)
    {
      if (i) PrintToMath(", ", FALSE);
      sprintf(mem,"{%s", DToMath(s, d->data[i])); PrintToMath(mem, FALSE);
      sprintf(mem,", %s}", DToMath(s, d->cweight[i])); PrintToMath(mem, FALSE);
    }
    PrintToMath("\n}", FALSE);
  }
  PrintToMath("\n}\n", FALSE);
}
#endif /* DISTTEST */

// Only create a distribution once up initialization;  double creation leads
// to memory leaks;  distributions should only be created through D_Init.

void D_Create(diststruct *d, char *header, long nlevels, double *binsize)
{
  long			i;

  D_Reset(d);
  if (nlevels<1) 
    nlevels		= 1;
  d->level		= nlevels-1;
  if (!(d->header = (char *) calloc(strlen(header)+1, sizeof(char))))
    Exit("distributions", "D_Create", "header calloc error");
  strcpy(d->header, header);
  if (!(d->binsize = (double *) calloc(nlevels, sizeof(double))))
    Exit("distributions", "D_Create", "header/binsize calloc error");
  for (i=0; i<nlevels; ++i)
    d->binsize[i]	= binsize[i];
}


void D_Init(diststruct **dist, char *header, long nlevels, double *binsize)
{
  char			s[80];
  long			system;

  if (*dist) 					// if existed already, then reset
  { 
    for (system=0; system<NSYSTEMS; ++system) 
      D_Reset(*dist+system); 
    return; 
  }
  if (!(*dist=(diststruct *) calloc(NSYSTEMS, sizeof(diststruct))))
    Exit("distributions", "D_Init", "calloc error");
  for (system=0; system<NSYSTEMS; ++system)
  {
    sprintf(s, "%s distribution of system %ld", header, system);
    D_Create(*dist+system, s, nlevels, binsize);
  }
}


void PutInDistribution(diststruct *d, double x, double y, double weight)
{
  D_Submit(d, &x, &y, &weight);
}


void PrintDistribution(diststruct *d)
{
  D_Print(d, 0, NULL);
}


void InitDist(long flag, diststruct **dist, char *header, long nlevels, double *binsize)
{
  if (flag)
    D_Init(dist, header, nlevels, binsize);
}
