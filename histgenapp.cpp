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
#include "hmdata.h"
#include "histgram.h"
#include "bitops.h"

using namespace std;

void print_usage(){
  fprintf(stderr, "usage: -d <Location of Data file>\n");
  fprintf(stderr, "       -o <histgram file prefix>\n");
  fprintf(stderr, "       -m number of parts for mih\n");
  fprintf(stderr, "       -N number of data documents\n");
  fprintf(stderr, "       -r rate of sample\n");
  fprintf(stderr, "\n");
  exit(0);
}

int main(int argc, char* argv[])
{
  struct timeval begin, end;
  char c;    
  char *dataFile = NULL;
  char *outputprefix = NULL;
  int M = 2;
  int64_t N = -1;
  double rateofsample = 0.1;
  while ((c = getopt(argc,argv, "hd:o:m:N:r:")) != -1)
  {
    switch (c)
    {
      case 'd':
        dataFile = optarg;
        break;
      case 'o':
        outputprefix = optarg;
        break;
      case 'm':
        M = atoi(optarg);
        break;
      case 'N':
        N = atoi(optarg);
        break;
      case 'r':
        rateofsample = atof(optarg);
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
  }
  
  if (dataFile == NULL){
    cerr << "Need data input and query input file name" <<endl;
    print_usage();
  }
  if (outputprefix == NULL){
    outputprefix = dataFile;
  }
  
  HmData *indata = new HmData(dataFile);
  
  if (N < 0 || N > indata->mNumData) {
    N = indata->mNumData;
  } 

  uint8_t* codes = indata->mData;
  uint64_t * chunks = new uint64_t[M];
  int b = ceil((double)indata->mDim/M);
  int mplus = indata->mDim - M * (b-1);

  int histgramrate = ceil(1.0 / rateofsample);

  histgram *Hist = new histgram[M];  
  for (int i=0; i<mplus; i++) {
    Hist[i].init(b);
  }
  for (int i=mplus; i<M; i++) {
    Hist[i].init(b-1);
  }

  gettimeofday(&begin, NULL);
  for (int64_t i=0; i<N; i++) {
    if (i % histgramrate == 0) {
      // printf("process doc %lld of %lld\n", i, N);
      split(chunks, codes+i*indata->mDataBytes, M, mplus, b);
      for (int k=0; k<M; k++) {
	Hist[k].insert(chunks[k], 1);
      }
      
      if (i % (int)ceil(N/1000) == 0) {
	printf("%.2f%%\r", (double)i/N * 100);
	fflush(stdout);
      }
    }
  }
  gettimeofday(&end, NULL);
  float indextime = end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1.0 / CLOCKS_PER_SEC;
  
  char outfile[256];
  for (int k = 0; k < M; k++) {
    sprintf(outfile, "%s_histgram_m=%d_p=%d.bin", outputprefix, M, k);
    Hist[k].WriteToFile(outfile);
  }
  fprintf(stdout, "Histgenapp for data %s M= %d time= %.3f N= %lld histgramrate= %d\n",
          dataFile, M, indextime, N, histgramrate);
  
  delete [] Hist;
  delete [] chunks;
  delete indata;
  return 0;
}
