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
#include "histallo.h"
#include "bitops.h"

using namespace std;

void print_usage(){
  fprintf(stderr, "usage: -d <Location of Data file>\n");
  fprintf(stderr, "       -o <histgram file prefix>\n");
  fprintf(stderr, "       -m number of parts for mih\n");
  fprintf(stderr, "       -t the distance\n");
  fprintf(stderr, "       -N number of data documents\n");
  fprintf(stderr, "\n");
  exit(0);
}

int32_t *slots;
histgram *Hist;
int m;
#define kMaxBuckets 64

int main(int argc, char* argv[])
{

  char c;    
  char *dataFile = NULL;
  char *outputprefix = NULL;
  int M = 2;
  int64_t N = -1;
  int tau = 33;
  double rateofsample = 0.1;
  while ((c = getopt(argc,argv, "hd:o:m:N:r:t:")) != -1)
    switch (c){
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
    case 't':
      tau = atoi(optarg);
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
  m = M;

  Hist = new histgram[M];
  slots = new int32_t[M];

  uint64_t greedytotal=0;
  uint64_t redgreedytotal=0;
  uint64_t dptotal=0;
  uint64_t reddptotal=0;
  uint64_t floortotal=0;
  uint64_t robintotal=0;
  uint64_t robinreducetotal=0;
      

  char outfile[256];
  for (int k = 0; k < M; k++) {
    sprintf(outfile, "%s_histgram_m=%d_p=%d.bin", outputprefix, M, k);
    Hist[k].LoadFromFile(outfile);
  }
  
  for (uint64_t i=0; i<N; i++) {
    split(chunks, codes+i*indata->mDataBytes, M, mplus, b);

    uint64_t totalcost = floorevenallocator(Hist, m, slots, chunks, tau);
    floortotal += totalcost;
    fprintf(stderr, "Floor[");
    for (int k=0; k<M; k++) {
      fprintf(stderr, "%d,", slots[k]);
    }
    fprintf(stderr, "]  %lld  ", totalcost);

    totalcost = roundrobinallocator(Hist, m, slots, chunks, tau);
    robintotal += totalcost;
    fprintf(stderr, "Robin[");
    for (int k=0; k<M; k++) {
      fprintf(stderr, "%d,", slots[k]);
    }
    fprintf(stderr, "]  %lld  ", totalcost);

    totalcost = roundrobinreduceallocator(Hist, m, slots, chunks, tau);
    robinreducetotal += totalcost;
    fprintf(stderr, "RedRobin[");
    for (int k=0; k<M; k++) {
      fprintf(stderr, "%d,", slots[k]);
    }
    fprintf(stderr, "]  %lld  ", totalcost);
    
    totalcost = greedyallocator(Hist, m, slots, chunks, tau);
    greedytotal += totalcost;
    fprintf(stderr, "Greedy[");
    for (int k=0; k<M; k++) {
      fprintf(stderr, "%d,", slots[k]);
    }
    fprintf(stderr, "]  %lld  ", totalcost);
    
    totalcost = greedyallocatorreduce(Hist, m, slots, chunks, tau);
    redgreedytotal += totalcost;
    fprintf(stderr, "RedGreedy[");
    for (int k=0; k<M; k++) {
      fprintf(stderr, "%d,", slots[k]);
    }
    fprintf(stderr, "]  %lld  ", totalcost);
    
    totalcost = DP_optimal_allocator(Hist, m, slots, chunks, tau);
    dptotal += totalcost;
    fprintf(stderr, "DP[");
    for (int k=0; k<M; k++) {
      fprintf(stderr, "%d,", slots[k]);
    }
    fprintf(stderr, "]  %lld  ", totalcost);
    
    totalcost = Reduced_DP_optimal_allocator(Hist, m, slots, chunks, tau);
    reddptotal += totalcost;
    fprintf(stderr, "RedDP[");
    for (int k=0; k<M; k++) {
      fprintf(stderr, "%d,", slots[k]);
    }
    fprintf(stderr, "]  %lld  ", totalcost);   
    
    for (int k=0; k<M; k++) {
      fprintf(stderr, "k=%d ->(", k);
      for (int error = 0; error < Hist[k].max_error; error ++) {
	fprintf(stderr, "%d, ", Hist[k].search(chunks[k], error));
      }
      fprintf(stderr, ") ");
    }
    fprintf(stderr, "\n");
  }

  printf("Allocation of %lld quries of %lld Dim of tau %d IN Floor %.3f Robin %.3f RedRobin %.3f Greedy %.3f ReducedGreedy %.3f DP %.3f ReducedDP %.3f\n",
	 N, indata->mDim, tau, floortotal*1.0/N, robintotal*1.0/N, robinreducetotal*1.0/N, greedytotal*1.0/N, redgreedytotal*1.0/N, dptotal*1.0/N, reddptotal*1.0/N);

  delete [] Hist;
  delete [] chunks;
  delete indata;
  return 0;
}
