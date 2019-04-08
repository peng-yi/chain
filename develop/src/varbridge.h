/*
    program:	varbridge.h
    author:	Peng Yi borrowed from Pieter J. in 't Veld
    date:	January 10, 2008
    purpose:	This module functions as a bridge between old style
    		system variables and new style system variables.
*/
#ifndef __VARBRIDGE_HEADER
#define __VARBRIDGE_HEADER

#ifdef __VARBRIDGE_MODULE

#include "header.h"


#else

extern bridgestruct *Bridge();
extern bridgestruct *BridgeMap(long flag_mpi);
/*
extern void BridgeInitMemory(bridgestruct *bridge, long nsystems, long NMols);
extern void BridgeResetVariable(bridgestruct *bridge);
extern void BridgeAddVariable(bridgestruct *to, bridgestruct *from);
extern void BridgeCopyInstant(
		bridgestruct *to, bridgestruct *from, long system);
*/

#endif

#endif

