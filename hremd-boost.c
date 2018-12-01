#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
/* TODO: use aceplug.h form CPATH */
#include "aceplug.h"
#include "boost.h"
#include "hremd-boost.h" 

#define RT (MD_BOLTZMAN*300.0)

double get_arg_int( int argc, char **argkey, char **argval, char *needle, int def ) 
{
  int i;
  for( i=0; i<argc; i++ ) {
    if( !strcmp( argkey[i], needle ) ) return atof(argval[i]);
  }
  return def;
}

char* get_arg_string( int argc, char **argkey, char **argval, char *needle, char *def ) 
{
  int i;
  for( i=0; i<argc; i++ ) {
    if( !strcmp( argkey[i], needle ) ) return argval[i];
  }
  return def;
}

void save_restart(struct aceplug_sim_t *s) 
{
  char fname[20];
  char backup[100];
  FILE *restart;
  long f_energy_fptr, f_perm_fptr;
  struct remd_t *p = (struct remd_t*) s->privdata;

  printf("Trying to save restart file...\n");

  /* Make backup of old restart file */
  snprintf(fname, 19, "%d.restart.bin", s->ensemble_rank);
  if(p->saved_step > 0) {
    snprintf(backup, 99, "%d.restart.%lu.bak", s->ensemble_rank, p->saved_step);
    rename(fname, backup);
  }
 
  f_perm_fptr = ftell(p->f_perm);
  if(f_perm_fptr < 0) goto error;
  f_energy_fptr = ftell(p->f_energy);
  if(f_energy_fptr < 0) goto error;
  restart = fopen(fname,"w");
  if(!restart) goto error;
  assert(fwrite(&(s->step), sizeof(unsigned long), 1, restart)==1);
  assert(fwrite(&(s->ensemble_size), sizeof(int), 1, restart)==1);
  assert(fwrite(&p->rngs[0], sizeof(unsigned short), 3, restart)==3);
  assert(fwrite(p->rep, sizeof(int), s->ensemble_size, restart)==s->ensemble_size);
  assert(fwrite(p->ff, sizeof(int), s->ensemble_size, restart)==s->ensemble_size);
  assert(fwrite(&f_perm_fptr, sizeof(long), 1,  restart)==1);
  assert(fwrite(&f_energy_fptr, sizeof(long), 1, restart)==1);
  assert(fclose(restart)==0);
  
  printf("Wrote restart information at step %lu.\n", s->step);
  p->saved_step = s->step;
  return;
 
error: 
  printf("Failed to write restart information at step %lu.\n", s->step);
}

void load_restart(struct aceplug_sim_t *s)
{
  int n;
  char fname[20];
  FILE *restart;
  long f_energy_fptr, f_perm_fptr;
  struct remd_t *p = (struct remd_t*) s->privdata;

  snprintf(fname, 19, "%d.restart.bin", s->ensemble_rank);
  restart = fopen(fname,"r");
  assert(restart);
  assert(fread(&(p->saved_step), sizeof(unsigned long), 1, restart)==1);
  assert(fread(&n, sizeof(int), 1, restart)==1);
  assert(n==s->ensemble_size);
  assert(fread(&p->rngs[0], sizeof(unsigned short), 3, restart)==3);
  assert(fread(p->rep, sizeof(int), s->ensemble_size, restart)==s->ensemble_size);
  assert(fread(p->ff, sizeof(int), s->ensemble_size, restart)==s->ensemble_size);
  assert(fread(&f_perm_fptr, sizeof(long), 1,  restart)==1);
  assert(fread(&f_energy_fptr, sizeof(long), 1, restart)==1);
  assert(fclose(restart)==0);
  assert(fseek(p->f_perm,f_perm_fptr,SEEK_SET)==0);
  assert(fseek(p->f_energy,f_energy_fptr,SEEK_SET)==0);
  
  printf("HREMD plugin loaded restart information from step %lu.\n", p->saved_step);
}

aceplug_err_t aceplug_init(struct aceplug_sim_t *s, int argc, char **argkey, char **argval) 
{
  struct remd_t *p;
  int i;

  p = (struct remd_t*) malloc( sizeof( struct remd_t ) );
  s->privdata = p;
  if( s->ensemble_size <= 1 ) {
    fprintf( stderr, "# REMD : Need to run in a parallel context.\n" );
    return ACEPLUG_ERR;
  }

  init_boost(s, argc, argkey, argval);

  p->ff = (int*) malloc( sizeof(int) * s->ensemble_size );
  assert(p->ff);
  p->rep = (int*) malloc( sizeof(int) * s->ensemble_size );
  assert(p->rep);


  /* get frequency for writing restart info  */
  p->restart_freq = get_arg_int(argc, argkey, argval, "restartfreq", 250000);
  
  /* get frequency for calculating the full matrix of exchange energies */
  p->full_energy_freq = get_arg_int(argc, argkey, argval, "fullenergyfreq", 250);
  
  //if(s->ensemble_rank==0) {
  char fname_perm[20];
  snprintf(fname_perm, 19, "%d.perm.log", s->ensemble_rank);
  char fname_energy[20];
  snprintf(fname_energy, 19, "%d.energy.log", s->ensemble_rank);
  //}
  
  
  if(strcmp(get_arg_string(argc, argkey, argval, "restart", "off"),"on")==0) {
    /* restore from restart file */
    p->f_perm = fopen(fname_perm,"r+");
    assert(p->f_perm);
    p->f_energy = fopen(fname_energy,"r+");
    assert(p->f_energy);
    load_restart(s);
    p->check_restart = 1;
  } else {
    /* normal initialization */
    /* open log files */
    p->f_perm = fopen(fname_perm,"w");
    assert(p->f_perm);
    fprintf(p->f_perm,"# cumulative permutation, column=Hamiltonian, value=replica\n");
    p->f_energy = fopen(fname_energy,"w");
    assert(p->f_energy);
    fprintf(p->f_energy,"# self & foreign energies: column=Hamiltonian, row=replica (not slot) \n");
    
    /* initialize assignement of replicas to Hamiltonians */
    for(i=0; i<s->ensemble_size; i++) { p->ff[i]=i; p->rep[i]=i; }
    
    /* initialize random number generator */
    p->rngs[0]= 0x330E;
    p->rngs[1]= 666;
    p->rngs[2]= 0;
    
    p->saved_step = -1;
    p->check_restart = 0;
  }
  
  printf("Using Hamiltonian %d in replica(=rank) %d.\n", p->ff[s->ensemble_rank], s->ensemble_rank); // debug
  
  return ACEPLUG_OK;
}


aceplug_err_t aceplug_terminate( struct aceplug_sim_t *s ) {
  return ACEPLUG_OK;
}

/* Given a Hamiltonian index i and the swapping scheme definition (first,last)
 * return the peer's Hamitonian index. The index must have a peer, i.e.
 * has_peer(i,first,last) should return true. */
int get_peer(int i, int first, int last)
{
  if((i+first)%2==0) { /* i-first ? */
    assert(i>=first && i+1<last);
    return i + 1;
  } else { 
    assert(i-1>=first && i<last);
    return i - 1;
  }
}

/* For Hamiltonian index i, return true if this Hamiltonian participates
 * in the current exchange scheme (first,last). */
int has_peer(int i, int first, int last) 
{
   if((i+first)%2==0 && i>=first && i+1<last) return 1;
   if((i+first)%2==1 && i-1>=first && i<last) return 1; 
   return 0;
}


aceplug_err_t aceplug_endstep( struct aceplug_sim_t *s) 
{
  struct remd_t *r = (struct remd_t*) s->privdata;
  double dU, rn;
  double energies[ENERGY_MAX];
  double temp_energies[s->ensemble_size];
  double all_energies[s->ensemble_size * s->ensemble_size];
  int *pairlist; size_t offset; int Npair;
  int N;
  int i,j, k, l, temp, other;
  int first, last;
  
  /* display debug info about restart file */
  if(r->check_restart) {
    printf("Restart info was from step %lu, current step is %lu\n.", r->saved_step, s->step);
    if(r->saved_step+1 != s->step) {
      printf("They don\'t agree. Terminating.\n");
      return ACEPLUG_ERR;
    }
    r->check_restart = 0;
  }

  /* not an energy step */
  if( ACEPLUG_OK !=  s->plugin_get_energy_temp( energies ) ) 
    return ACEPLUG_OK; 
  
  N = s->ensemble_size;

  /* decide on swapping attempt scheme */
  if(erand48(r->rngs)>0.5) 
    first = 0; /* 0-1 2-3 4-5 ... */
  else 
    first = 1; /* 0 1-2 3-4 5-6 ... */
  
  last = ((N-first)/2)*2+first;
  /* first and last uniquely define the scheme 
   * Hamiltonians with indices i: first<=i<last participate in the exchange. */
  
  /* failsafe: set all energies to sentinel, so we can later check if we missed some calculations */
  for ( i=0; i < s->ensemble_size; i++ ) temp_energies[i] = NAN;

  /* prepare ACEMD energy calculation */
  TEST(s->plugin_load_positions()); 
  TEST(s->plugin_get_pairlist(&pairlist,&offset,&Npair)); 
  
  /* TODO: refactor here:
   * (1) calculate the physical energy V once per replica and store it in a local variable V_phys
   * (2) use this value in all calls to calc_boost_energy
   * (3) we could go a step further and store V_phys in r, then we can use it in the force evaluation of the next step.
   *     I would be important to know, at what time the pairlist is updated.
   * */
  
  /* get self energy */
  temp_energies[r->ff[s->ensemble_rank]] = calc_boost_energy(s,r->ff[s->ensemble_rank],pairlist,offset,Npair);

  /* evaluate this systems's energies under the  */
  /* Hamiltonian of all other peers */
  if(s->step % r->full_energy_freq==0) {
    /* calculate all combinations of energies */
    for ( i=0; i < s->ensemble_size; i++ ) {
      if(i!=r->ff[s->ensemble_rank]) { /* we already know the self energy */
        temp_energies[i] = calc_boost_energy(s,i,pairlist,offset,Npair);
      }
    } 
  } else { 
    if(has_peer(r->ff[s->ensemble_rank],first,last)) {
      other = get_peer(r->ff[s->ensemble_rank],first,last);
      temp_energies[other] = calc_boost_energy(s,other,pairlist,offset,Npair);
    }
  }

  /* communinate energies to the swarm */
  TEST( s->plugin_ensemble_allgather(temp_energies, N*sizeof(double), all_energies) );
  /* all_energies[i*N+j] is the i'th replica's energy evaluated with Hamiltonian j */
  
  /* i,j replica indices; k,l: Hamiltonian indices */
  for(k=first; k<last; k+=2) {
    l = get_peer(k,first,last);
    rn = erand48(r->rngs);
    i = r->rep[k];
    j = r->rep[l];
    assert(!isnan(all_energies[i*N+k]));
    assert(!isnan(all_energies[i*N+l])); 
    assert(!isnan(all_energies[j*N+l])); 
    assert(!isnan(all_energies[j*N+k]));
    dU = all_energies[i*N+k] - all_energies[i*N+l] + all_energies[j*N+l] - all_energies[j*N+k];
    //printf("Attempting %d %d %f %f\n",k,l,rn,dU/RT);
    if(rn < exp(dU/RT))  {
      //printf("Exchanged replicas %d and %d (or Hamiltonians %d and %d)\n",i,j,k,l);
      /* update ff permuation... */
      temp = r->ff[i];
      r->ff[i] = r->ff[j];
      r->ff[j] = temp;
      /* ...and its inverse, the replica permutation */
      temp = r->rep[k];
      r->rep[k] = r->rep[l];
      r->rep[l] = temp;
    }
  }
  /* Potentials are automatically swapped, because calc_forces respects r->ff */
  
  /* debug: check that the permuations are each others inverses and that they contain no duplicate numbers */
  for(i=0;i<N;i++) assert(i==r->rep[r->ff[i]] && i==r->ff[r->rep[i]]);
  
  /* restrict ourselves to full_energy_freq when writing permutations and energies */
  /* note that the simulation starts with step=1 (not 0) */
  if(s->step % r->full_energy_freq==0) {
  //if(s->ensemble_rank==0) { // debug
    /* write permuation log file */
    for(i=0; i<N; i++) {
      fprintf(r->f_perm, "%d", r->rep[i]);
      if(i==N-1) fprintf(r->f_perm, "\n");
      else fprintf(r->f_perm, " ");
    }
    fflush(r->f_perm);
    
    /* write energy log file */
    for(i=0; i<N; i++) {
      for(j=0; j<N; j++) { 
        fprintf(r->f_energy, "%f", all_energies[i*N+j]);
        if(j==N-1) fprintf(r->f_energy, "\n");
        else fprintf(r->f_energy, " ");
      }
    }
    fprintf(r->f_energy, "&\n");
    fflush(r->f_energy);
  //}
  }
  
  if(s->step % r->restart_freq==0) save_restart(s);
  
  return 0;
}
