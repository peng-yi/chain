/*
    program:    header.h
    author:     Peng Yi at MIT
    date:       October 19, 2006
    purpose:    Header file for all programs.
*/
//#include </usr/mpi/intel/openmpi-1.4.3/include/mpi.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <math.h>
#include <time.h>
#include <complex.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sort.h>
#include "input.h"
#include "types.h"
#include "globals.h"
#include "distributions.h"
#include "init.h"
#include "random.h"
#include "sample.h"
#include "position.h"
#include "forcefield.h"
#include "ensembles.h"
#include "history.h"
#include "varbridge.h"
#include "lists.h"
#include "vector.h"
#include "io.h"
#include "motion.h"
#include "units.h"
#include "rebridge.h"
#include "roots.h"
#include "library.h"	// lammps library
