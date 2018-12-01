#ifndef BOOST_H
#define BOOST_H

#include <stdlib.h>
/* TODO: use aceplug.h form CPATH */
#include "aceplug.h"

void init_boost(struct aceplug_sim_t *s, int argc, char **argkey, char **argval);
aceplug_err_t aceplug_calcforces( struct aceplug_sim_t *s );
double calc_boost_energy(struct aceplug_sim_t *s, int i_ff, int *pairlist, size_t offset, int N);

#endif
