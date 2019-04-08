/*
    program:    motion.h
    author:     Peng Yi at MIT
    date:       October 26, 2007
    purpose:    move molecule, sites, functions borrowed and modified from Pieter's
*/
#ifndef __MOTION_HEADER
#define __MOTION_HEADER

#include "types.h"

#ifdef __MOTION_MODULE

#include "header.h"

#else

extern sphere	SiteSpherical(molstruct *, long);
extern vector	SiteCartesian(molstruct *, long, sphere);
extern void	MolSpherical(molstruct *);
extern void	MolCartesian(molstruct *);
extern void	AllSpherical();
extern void	AllCartesian();

// Molecule manipulators

extern void	SiteCopy(molstruct *molm, long m, molstruct *moln, long n, long p);
extern void	SiteSwap(molstruct *molm, long m, molstruct *moln, long n);
extern void	MolReverse(molstruct *mol1);	// reverse a molecule withOUT cell list update
extern void	MolFlip(molstruct *moli);	// reverse a molecule WITH cell list update
extern molstruct	MolAdd(molstruct *mol1, molstruct *mol2);
extern void	ChangeAxis(long system, vector scale);
extern void	ChangeVolume(long system, double scale);
#endif

#endif
