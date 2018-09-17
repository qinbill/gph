#include <getopt.h>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <vector>
#include "sys/types.h"


#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <pthread.h>
#include <time.h>
#include <string.h>

#include "types.h"
#include "linscan.h"
#include "histgram.h"
#include "bitops.h"

using namespace std;

const char* toBinary(uint64_t code, int size)
{
  static char bs[128];
  bs[size] = '\0';
  for (int i = 0; i < size; i++) {
    if (code % 2 == 0) {
      bs[size-i-1] = '0';
    }else {
      bs[size-i-1] = '1';
    }
    code = code / 2;
  }
  return bs;
}

void print_usage(){
  fprintf(stderr, "usage: -o <histgram file prefix>\n");
  fprintf(stderr, "       -m number of parts for mih\n");
  fprintf(stderr, "       -t the distance\n");
  fprintf(stderr, "\n");
  exit(0);
}

int main(int argc, char* argv[])
{
  char c;    
  char *outputprefix = NULL;
  int M = 2;
  int pid = 0;
  while ((c = getopt(argc,argv, "ho:m:i:")) != -1)
    switch (c){
    case 'o':
      outputprefix = optarg;
      break;
    case 'm':
      M = atoi(optarg);
      break;
    case 'i':
      pid = atoi(optarg);
      break;
    case 'h':
      print_usage();
      break;
    case '?':
      if ( optopt == 't' || optopt == 'd' || optopt == 'D' || optopt == 'l' || optopt == 'q' || optopt == 't' )
	cerr << "Error: Option -" << optopt << "requires an argument." << endl;
      else if ( isprint(optopt))
	cerr << "Error: Unknown Option -" << optopt << endl;
      else
	cerr << "Error: Unknown Option character" <<endl;
      return 1;
    default:
      print_usage();
    }
  
  if (outputprefix == NULL){
    print_usage();
  }

  histgram* Hist = new histgram[M];
  
  char outfile[256];
  for (int k = 0; k < M; k++) {
    sprintf(outfile, "%s_histgram_m=%d_p=%d.bin", outputprefix, M, k);
    Hist[k].LoadFromFile(outfile);
  }

  for (uint64_t i = 0; i < Hist[pid].bucketsize; i ++) {
    fprintf(stdout, "%s", toBinary(i, Hist[pid].bucketbits));
    for (int err = 0; err < Hist[pid].max_error; err++) {
      fprintf(stdout, " %.9f", (1.0*Hist[pid].mHistgram[i*Hist[pid].max_error +err])/Hist[pid].totalcount);
    }
    fprintf(stdout, "\n");
  }
  delete [] Hist;
  return 0;
}
