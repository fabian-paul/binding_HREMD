#ifndef HREMD_BOOST_H
#define HREMD_BOOST_H

#include <stdio.h>
/* TODO: use aceplug.h form CPATH */
#include "aceplug.h"

struct remd_t {
  double *E;
  double *alpha;
  int restart_freq;
  int full_energy_freq;
  unsigned short rngs[3];
  int *ff;  /* permutation ff[rep] */
  int *rep; /* permutation rep[ff] */
  FILE *f_perm;
  FILE *f_energy;
  int check_restart;
  unsigned long saved_step; /* step # stored in last restart file */
};

char* get_arg_string( int argc, char **argkey, char **argval, char *needle, char *def );
double get_arg_int( int argc, char **argkey, char **argval, char *needle, int def ); 

#define  TEST(a) { int retval = a; if( ACEPLUG_OK != retval) { printf( "ERROR %d\n", retval ); assert(0); exit(-1);  }}

#endif
