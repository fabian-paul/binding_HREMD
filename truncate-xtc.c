/* compile: gcc -Wall truncate-xtc.c -lxdrfile -o truncate-xtc */
/* run:     truncate-xtc Nframes input.xtc output.xtc */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <xdrfile/xdrfile.h>
#include <xdrfile/xdrfile_xtc.h>
#include <xdrfile/xdrfile_trr.h>

void error(const char *str) 
{
  fputs("truncate-xtc: ", stderr);
  fputs(str, stderr);
  fputs("\n", stderr);
  exit(1);
}

int main(int argc, char* argv[])
{
  XDRFILE *input;
  XDRFILE *output;
  int step;
  float time, lambda;
  matrix box;
  rvec *frame, *v, *f;
  float prec;
  int N_atoms;
  int N_frames;
  int res, i;
  int xtc = 0;

  if(argc<4) error("Argument error");
 
  N_frames = atoi(argv[1]);
 
  /* determine input file type */ 
  i = strlen(argv[2])-4;
  if(i>0 && strcmp(argv[2]+i,".xtc")==0) xtc=1;
  else xtc=0;
 
  if(xtc) read_xtc_natoms(argv[2], &N_atoms);
  else    read_trr_natoms(argv[2], &N_atoms);
  input = xdrfile_open(argv[2], "r");
  if(!input) error("Open file (read) error");
  
  frame = malloc(N_atoms*sizeof(rvec));
  if(!frame) error("Memory error");
  if(!xtc) {
    v = malloc(N_atoms*sizeof(rvec));
    f = malloc(N_atoms*sizeof(rvec));
    if(!v || !f) error("Memory error");
  }
  
  output = xdrfile_open(argv[3], "w");
  if(!output) error("Open file (write) error");
  
  for(i=0; i<N_frames; i++) {
    if(xtc) res = read_xtc(input, N_atoms, &step, &time, box, frame, &prec);
    else    res = read_trr(input, N_atoms, &step, &time, &lambda, box, frame, v, f);
    if(res==exdrENDOFFILE) break;
    if(res!=exdrOK) error("Read error");
    res = write_xtc(output, N_atoms, step, time, box, frame, prec);
    if(res!=exdrOK) error("Write error");
  }
  if(i<N_frames) fprintf(stderr, "WARNING: file %s contains less than %d frames!\n", argv[2], N_frames);
  xdrfile_close(input);
  xdrfile_close(output);
  
  return 0;
}
