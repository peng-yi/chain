/*
    program:	io.h
    author:     Peng Yi at MIT
    date:       October 24, 2006
    purpose:    Integral module for input/output handling
*/
#ifndef __IO_HEADER
#define __IO_HEADER


#ifdef __IO_MODULE

#include "header.h"

FILE		*foutput, *fhst, *fdump;	// output file handlers
char		file_hst[256];			// name of binary history file
long		frame;				// visualization output frame #

#else

extern FILE	*foutput, *fhst, *fdump;
extern char	file_hst[256];
extern long	frame;

extern int	samestr(char *, char *);	// compare two strings, case insensitive
extern void	GetLVar(FILE *, long, long *);
extern void	GetDVar(FILE *, long, double *);
extern void	GetSetup(char * argv[]);	// read in setup info.
extern void 	InitFile();
extern void 	CloseFile();

extern void	Exit(char *module, char *procedure, char *error); 

extern void	PrintSetup();		// print out setup for record
extern void	Print_Header(FILE *);
extern void	Print_Histogram();
extern void	Print_Nuclei();
extern void	Print_Verlet();
extern void	Print_Verlet2();
extern void	Print_Clist();
extern void	Print_q();
extern void	Print_Q();
extern void	Print_qproduct();
extern void	Printout();
extern int 	Visualize(int);
extern void	Print_gr();
extern void	Print_Summary();

extern int	Write_Conf(long timestep);
extern int	Read_Conf(char *);
extern int	Read_MultiConf(FILE *fPtr);
extern void	dump_conf(long timestep);
#endif

#endif
