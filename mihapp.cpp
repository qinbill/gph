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
#include "mihasher.h"


using namespace std;

void print_usage(){
  fprintf(stderr, "usage: -t <Hamming Distance>\n");
  fprintf(stderr, "       -d <Location of Data file>\n");
  fprintf(stderr, "       -q <Location of query file>\n");
  fprintf(stderr, "       -p <projected dimension>\n");
  fprintf(stderr, "       -k calculate top k instead of ranger search\n");
  fprintf(stderr, "       -l use linear scan\n");
  fprintf(stderr, "       -n use number of queries\n");
  fprintf(stderr, "       -g threshold gap\n");
  fprintf(stderr, "       -m number of parts for mih\n");
  fprintf(stderr, "       -N number of data strings, default all.\n");
  fprintf(stderr, "\n");
  exit(0);
}

int main(int argc, char* argv[])
{

  struct timeval begin, end;
  char c;    
  int tau = 0; 
  char *dataFile = NULL;
  char *queryFile = NULL;
  int K = 100;
  int NQ = 100000000;
  int M = 2;
  int gap = 1;
  int64_t N = -1;
  bool uselinscan = false;
  bool dotopk = false;
  int projdim = 0;
  while ((c = getopt(argc,argv, "ht:d:q:lk:n:m:N:p:g:")) != -1)
    switch (c){
      case 't':
        tau = atoi(optarg);
        break;        
      case 'd':
        dataFile = optarg;
        break;
      case 'q':
        queryFile = optarg;
        break;
      case 'n':
        NQ = atoi(optarg);
        break;
      case 'g':
	gap = atoi(optarg);
        break;
      case 'm':
        M = atoi(optarg);
        break;
      case 'N':
        N = atoi(optarg);
        break;
      case 'p':
        projdim = atoi(optarg);
        break;	
      case 'l':
        uselinscan = true;
        break;
      case 'k':
        K = atoi(optarg);
        dotopk = true;
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
  
  if (dataFile == NULL || queryFile == NULL){
    cerr << "Need data input and query input file name" <<endl;
    print_usage();
  }


  HmData *indata = new HmData(dataFile);
  HmData *inquery = new HmData(queryFile);
  HmData *projdata, *projquery;
  
  if (projdim != 0) {
    char newfilename[1024];
    sprintf(newfilename, "%s.p%d.bin", dataFile, projdim);
    //fprintf(stderr, "Load File %s\n", newfilename);
    projdata = new HmData(newfilename);
    sprintf(newfilename, "%s.p%d.bin", queryFile, projdim);
    //fprintf(stderr, "Load File %s\n", newfilename);
    projquery = new HmData(newfilename);
  }else {
    projdata = indata;
    projquery = inquery;
  }
  
  if (N < 0 || N > indata->mNumData) {
    N = indata->mNumData;
  }

  uint32_t B = indata->mDim;
  NQ=min((uint32_t)NQ, (uint32_t)inquery->mNumData);
  UINT32 *res = new UINT32[K*N];
  UINT32 *counter = new UINT32[(B+1)*NQ];

  mihasher* MIH = NULL;
  qstat *stats = NULL;
  uint32_t *results = NULL;
  uint32_t *numres = NULL;

  if (dotopk) {
    K = min((uint32_t)K, (uint32_t)N);
  } else {
    K = N;
  }
    
  // Preprocessing. 
  gettimeofday(&begin, NULL);
  if (uselinscan) {
    // Do nothing.
  } else {
    // Do the him indexing. 
    MIH = new mihasher(projdata->mDim, M);
    MIH->populate(projdata->mData, N, projdata->mDataBytes);
    stats = (qstat*) new qstat[NQ];
    results = new uint32_t[(K)*NQ];
    numres = new uint32_t[(indata->mDim+1)*NQ];
    MIH->setK(K);
  }
  gettimeofday(&end, NULL);
  float indextime = end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1.0 / CLOCKS_PER_SEC;

  UINT32 numresults = 0;
  for (int t = 0; t <=tau; t += gap) {
    gettimeofday(&begin, NULL);
    // Perform the query processing. 
    if (dotopk) {
      // Perform top k search here.
      if (uselinscan) {
        linscan_query(counter, res, indata->mData, inquery->mData, indata->mNumData, NQ, B, K, \
                      indata->mDataBytes, inquery->mDataBytes);
      } else {
        MIH->batchquery (results, numres, stats, inquery->mData, NQ, inquery->mDataBytes);
      }
    } else {
      // Perform range search here.
      if (uselinscan) {
        results = new uint32_t[indata->mNumData];
        UINT32 numres = 0;
        linscan_rangequery(results, &numres, indata->mData, inquery->mData, N, NQ, B,
                           indata->mDataBytes, inquery->mDataBytes, t);
        printf("numResults %d\n", numres);
        delete [] results;
      } else {
        numresults = MIH->batchrangequery (results, stats, projquery->mData, NQ, projquery->mDataBytes, t, inquery->mData, indata->mData, indata->mDataBytes);
      }    
    }
    gettimeofday(&end, NULL);
    float querytime = end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1.0 / CLOCKS_PER_SEC;
    fprintf(stdout, "Alg: MIHRING #dim: %lld Tau: %d qrytime %.3f idxtime %.f #data %lld #qry %d K: %d #rst %d\n", indata->mDim, t, querytime, indextime, N, NQ, K, numresults);
  }

  if (!uselinscan) {
    delete MIH;
    delete [] stats;
    delete [] results;
    delete [] numres;
  }
  delete counter;
  delete indata;
  delete inquery;
  return 0;
}
