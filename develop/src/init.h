/*
    program:    init.h
    author:     Peng Yi at MIT
    date:       October 23, 2006
    purpose:    Suggested initialization sequence
*/

#ifndef __INIT_HEADER
#define __INIT_HEADER

#ifdef __INIT_MODULE

#include "header.h"

char		INITCONF[256];		//initial configuration

#else

extern char	INITCONF[256];

extern void	GetCoordinates(char * filename);
extern void	InitMols(long, long);
extern void 	InitAll(char *argv[]);

#endif

#endif
