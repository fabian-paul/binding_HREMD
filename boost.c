#include <stdio.h> 
#include <assert.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
/* TODO: use aceplug.h form CPATH */
#include "aceplug.h"
#include "hremd-boost.h" 
#include "boost.h" 

#define N_pmi 195
#define atom_in_mdm2(x) ((x)<1449)
/* in a more general version pmi and atom_in_mdm2 should come from plugin parameters as lists of atoms */


/* debug function */
double slow_pbc_diff(double xi, double xj, double xsize) {
  double dx;
  xi=xi-floor(xi/xsize)*xsize; 
  xj=xj-floor(xj/xsize)*xsize; 
  dx=xi-xj;
  dx=dx-(int)round((dx/xsize))*xsize;
  return dx;
}
/* end debug function */

/* calculate xi-xj taking into account periodic boundary conditions */
inline double pbc_diff(double xi, double xj, double xsize)
{
  double v = xi-xj, half = xsize / 2.0;
  while(v >= half) v -= xsize;
  while(v < -half) v += xsize;
  //double v2 = slow_pbc_diff(xi,xj,xsize); /* debug */
  //if(v!=v2) fprintf(stderr, "@ pbc mismatch %e\n, ",v-v2); /* debug */
  return v;
}

#define NUMBER 0
#define SPACE 1
static int parse_string(char *str, double *dest, int size)
{
  char buffer[0x100];
  char *bufptr;
  int state;
  int ocount;
  double d;
  
  ocount = 0;
  bufptr = &buffer[0];
  if(isspace(*str)) state = SPACE; else state = NUMBER;

  while(*str!='\0') {
    switch(state) {
      case NUMBER:
        if(isdigit(*str) || *str=='.' || *str=='-' || *str=='+' || *str=='e' || *str=='E') {
          *bufptr = *str; bufptr++;
          if(bufptr-buffer >= 0x100) return -1;
        }
        else if(isspace(*str) || *str==',') { 
          *bufptr = '\0';
          d = atof(buffer);
          if(size>0 && ocount>size) return -1;
          if(size>0) dest[ocount] = d;
          //printf("Converted %lf\n",d);
          ocount ++;
          bufptr = &buffer[0];
          state = SPACE; 
        } 
        else return -1;
        break;
      case SPACE:
        if(!isspace(*str)) { state = NUMBER; str--; }
        break;       
    }
    str++;
  }
  
  if(bufptr!=buffer) { /* last number was not converted */
    *bufptr = '\0';
    d = atof(buffer);
    if(size>0 && ocount>size) return -1;
    if(size>0) dest[ocount] = d;
    ocount ++;
  }  
  //printf("Ocount %d\n",ocount);
  return ocount;
}

void init_boost(struct aceplug_sim_t *s, int argc, char **argkey, char **argval) {
  int i, j;
  int pmi[N_pmi];
  struct remd_t *r = (struct remd_t*) s->privdata;
  char *E_string;
  char *alpha_string;
  
  for(i=1449, j=0; i < 1644; i++, j++) pmi[j]=i;
  TEST(s->plugin_set_pairlist(pmi, N_pmi));  
  
  r->E = (double*)malloc(sizeof(double)*s->ensemble_size);
  assert(r->E);
  r->alpha = (double*)malloc(sizeof(double)*s->ensemble_size);
  assert(r->alpha);

  E_string = get_arg_string(argc, argkey, argval, "E", NULL);
  assert(E_string);
  //printf("E_string [%s]\n",E_string);
  alpha_string = get_arg_string(argc, argkey, argval, "alpha", NULL);
  assert(alpha_string);
  //printf("alpha_string [%s]\n",alpha_string);
  
  assert( parse_string(E_string, r->E, s->ensemble_size) == s->ensemble_size);
  assert( parse_string(alpha_string, r->alpha, s->ensemble_size) == s->ensemble_size);
  
  for(i=0; i<s->ensemble_size; i++) printf("E[%d] = %lf\n", i, r->E[i]);
  for(i=0; i<s->ensemble_size; i++) printf("alpha[%d] = %lf\n", i, r->alpha[i]);

  /*E = 0.0;*/
  /*alpha = 927.574/4.184;*/
  /*sigma = 20 / 4.184 ;
  V_min = -250  / 4.184 ;  approximate minimal LJ interaction 
  V_min_new = -50 / 4.184 ;
  E = V_min_new + sigma;
  alpha = (E-V_min)*(E-V_min)/(V_min_new-V_min) + V_min - E;  
	return 0;*/
}

/* Calculate (unmodified) LJ interaction energy between ligand and protein */
double calc_LJ_energy(const struct aceplug_sim_t *const s, const int* const pairlist, const size_t offset, const int N) {
  double dx, dy, dz;
  double r2, r6, r12;
  double a, b;
  int k, l, i, j;
  double V;
  const double4* const pos = s->pos;
  const float4 box=s->box;

  V = 0.0;
  /* l and k are just loop indices, i and j refer to atom numbers */
//#pragma omp parallel for default(none) private(k,i,j,dx,dy,dz,r2,r6,r12,a,b) shared(N,pairlist,offset,box,pos,s) reduction(+:V)
#pragma omp parallel for private(k,i,j,dx,dy,dz,r2,r6,r12,a,b) reduction(+:V)
  /* access to V is the only concurrent write operation, use OpenMP summation reduction */
  for(l=0; l<N; l++) {
    k = 1; /* zeroth element is PMI atom index i  */
    i = pairlist[l*offset];
    while((j=pairlist[l*offset+k])!=-1) { 
      if(atom_in_mdm2(j)) {
        dx = pbc_diff(pos[i].x,pos[j].x,box.x);
        dy = pbc_diff(pos[i].y,pos[j].y,box.y);
        dz = pbc_diff(pos[i].z,pos[j].z,box.z);
        r2 = dx*dx+dy*dy+dz*dz;
        if(r2 < 9*9) {
          r6 = r2*r2*r2;
          r12 = r6*r6;
          TEST(s->plugin_vdw_param(i,j,0,&a,&b));
          V += a/r12 - b/r6;
        }
      } 
      k++; 
    }  
  }
  return V;
}

/*  Apply a modified Lennard-Jones force to protein and ligand atoms.
 *  The modification consists of a prefactor that multiplies the
 *  force. (factor=-1 would mean undoing the force that was added on the GPU).
 */ 
void apply_forces(const struct aceplug_sim_t *const s, const double factor, const int *const pairlist, const size_t offset, const int N) {
  double r2, r6, r12;
  double dx, dy, dz;
  double frc_scalar;
  double a,b;
  int k, l, i, j;
  const double4* const pos = s->pos;
  float4* const frc = s->frc;
  const float4 box = s->box;

  /* l and k are just loop indices, i and j refer to atom number */
#pragma omp parallel for private(k,i,j,dx,dy,dz,r2,r6,r12,a,b,frc_scalar)
  for(l=0; l<N; l++) {
    k = 1; /* zeroth element is PMI atom index i  */
    i = pairlist[l*offset];
    while((j=pairlist[l*offset+k])!=-1) { 
      if(atom_in_mdm2(j)) {
        dx = pbc_diff(pos[i].x,pos[j].x,box.x);
        dy = pbc_diff(pos[i].y,pos[j].y,box.y);
        dz = pbc_diff(pos[i].z,pos[j].z,box.z);
        r2 = dx*dx+dy*dy+dz*dz;
        if(r2 < 9*9) {
          r6 = r2*r2*r2;
          r12 = r6*r6;
          TEST(s->plugin_vdw_param(i,j,0,&a,&b)); 
          frc_scalar = factor * (12*a/r12 - 6*b/r6) / r2;
          frc[i].x += dx*frc_scalar;
          frc[i].y += dy*frc_scalar;
          frc[i].z += dz*frc_scalar;
#pragma omp atomic
          frc[j].x -= dx*frc_scalar;
#pragma omp atomic
          frc[j].y -= dy*frc_scalar;
#pragma omp atomic
          frc[j].z -= dz*frc_scalar;
          /* The frc[i] operations don't need to be atomic, because there is a */
          /* bijection from i to l and l is the loop variable that is managed */
          /* by OpenMP. Therefore every thread gets a non-overlapping chunk of i's. */
          /* However nothing can be said about the j's and they certainly overlap. */
          /* Therefore write acces to frc[j] needs to be atomic. */
        }
      } 
      k++; 
    }  
  }
}

/* Correct forces according to current boosting level. */ 
aceplug_err_t aceplug_calcforces( struct aceplug_sim_t *s ) {
  int *pairlist;
  size_t offset;
  int N;
  double V;
  double factor;
  double E, alpha;
  struct remd_t *r = (struct remd_t*) s->privdata;
  
  E = r->E[r->ff[s->ensemble_rank]];
  alpha = r->alpha[r->ff[s->ensemble_rank]];

  TEST(s->plugin_load_positions()); 
  TEST(s->plugin_get_pairlist(&pairlist,&offset,&N)); 
  assert(N<=N_pmi);
  
  V=calc_LJ_energy(s,pairlist,offset,N); 
  if(V<E) {
    factor = (alpha/(alpha+E-V))*(alpha/(alpha+E-V)) - 1.0;
    TEST(s->plugin_load_forces()); 
    apply_forces(s,factor,pairlist,offset,N);
    TEST(s->plugin_update_forces()); 
  }

  /*if((s->step % 250)==0) printf("@@@@@ LJ energy %f kcal/mol = %f kJ/mol\n", V, V*4.184);*/
  
  return ACEPLUG_OK;
}

/* Calculate modified energy of current molecular conformation with Hamiltonian number i_ff */
double calc_boost_energy(struct aceplug_sim_t *s, int i_ff, int *pairlist, size_t offset, int N) {
  double V, E, alpha;
  struct remd_t *r = (struct remd_t*) s->privdata;
  
  V = calc_LJ_energy(s, pairlist, offset, N);
  E = r->E[i_ff];
  alpha = r->alpha[i_ff];
 
  if(V>=E) 
    return V;
  else 
    return V + (E-V)*(E-V)/(alpha+E-V);
}
