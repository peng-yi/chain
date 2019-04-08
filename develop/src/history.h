/*
  program:	history.h
  author:	Peng Yi borrowed from Pieter J. in 't Veld
  date:		January 10, 2008
  purpose:	i/o for binary history files
*/
#ifndef __HISTORY_HEADER
#define __HISTORY_HEADER

#define HIST_VERSION	25
#define HIST_IDENT	"HIST"

#include "header.h"

#ifdef __HISTORY_MODULE

#else

extern int StartNewHistoryFile(char *name, long flag_mpi);
extern int H_StoreCurrent(char *name);
extern int H_StoreBridge(char *name, bridgestruct *bridge);
extern int H_GetHeader(char *name, long *version);
extern FILE *H_GetFirst(char *name, long *version, long flag_mpi);
extern int H_GetNext(FILE *fp, long version);
extern void H_InitSystem(char *argv[]);
extern bridgestruct *H_Bridge();

extern FILE *fcreate_hist(char *);		// for hstcomb
extern FILE *fread_hist(char *, long *);	// for hstcomb
#endif

#endif

