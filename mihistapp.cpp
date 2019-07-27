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
#include "histallo.h"


// Add dlib regressor 
#include "mlregressor.h"

using namespace std;



void print_usage(){
  fprintf(stderr, "usage: -t <Hamming Distance>\n");
  fprintf(stderr, "       -d <Location of Data file>\n");
  fprintf(stderr, "       -o <Location of Histgram file>\n");
  fprintf(stderr, "       -q <Location of Query file>\n");
  fprintf(stderr, "       -n use number of queries\n");
  fprintf(stderr, "       -p print out single query candidate.\n");
  fprintf(stderr, "       -m number of parts for mih\n");
  fprintf(stderr, "       -N number of data strings, default all.\n");
  fprintf(stderr, "       -g threshold gap default 1\n");
  fprintf(stderr, "       -s threshold start default 0\n");
  fprintf(stderr, "       -e threshold extension, default 0\n");
  fprintf(stderr, "       -a <Algorithm for allocation>\n");
  fprintf(stderr, "          Floor    -> floorevenallocator>\n");
  fprintf(stderr, "          Greedy   -> greedyallocator>\n");
  fprintf(stderr, "          Rgreedy  -> greedyallocatorreduce>\n");
  fprintf(stderr, "          Robin    -> roundrobinallocator>\n");
  fprintf(stderr, "          Rrobin   -> roundrobinreduceallocator>\n");
  fprintf(stderr, "          Dp       -> DP_optimal_allocator>\n");
  fprintf(stderr, "          Rdp      -> Reduced_DP_optimal_allocator>\n");
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
  char *outputprefix = NULL;
  int K = 100;
  int NQ = 100000000;
  int ext = 0;
  int M = 2;
  bool printcand=false;
  int64_t N = -1;
  bool uselinscan = false;
  bool dotopk = false;
  char* algorcode = "rdp";
  int taustart=0;
  int taugap = 1;
  histgram *Hist = NULL;

  // Prefix of regressor file
  char* regfile = NULL;
  while ((c = getopt(argc,argv, "ht:d:q:a:o:n:m:s:g:N:e:pr:")) != -1)
    switch (c){
      case 't':
        tau = atoi(optarg);
        break;
      case 's':
        taustart = atoi(optarg);
        break;
      case 'g':
        taugap = atoi(optarg);
        break;
      case 'e':
        ext = atoi(optarg);
        break;                	
      case 'd':
        dataFile = optarg;
        break;
      case 'o':
        outputprefix = optarg;
        break;	
      case 'q':
        queryFile = optarg;
        break;
      case 'n':
        NQ = atoi(optarg);
        break;
      case 'm':
        M = atoi(optarg);
        break;
      case 'N':
        N = atoi(optarg);
        break;
      case 'a':
        algorcode = optarg;
        break;
      case 'p':
	printcand = true;
	break;
      case 'h':
        print_usage();
        break;
      case 'r':
 	regfile = optarg;	
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

  if (outputprefix == NULL) {
    outputprefix = dataFile;
  }

  HmData *indata = new HmData(dataFile);
  HmData *inquery = new HmData(queryFile);
  if (N < 0 || N > indata->mNumData) {
    N = indata->mNumData;
  }

  Hist = new histgram[M];
  char outfile[256];
  for (int k = 0; k < M; k++) {
    sprintf(outfile, "%s_histgram_m=%d_p=%d.bin", outputprefix, M, k);
    Hist[k].LoadFromFile(outfile);
  }


  uint32_t B = indata->mDim;
  NQ=min((uint32_t)NQ, (uint32_t)inquery->mNumData);
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
    MIH = new mihasher(B, M);
    //MIH->populate(indata->mData, N, indata->mDataBytes);
    stats = (qstat*) new qstat[NQ];
    results = new uint32_t[(K)*NQ];
    numres = new uint32_t[(indata->mDim+1)*NQ];
    MIH->setK(K);
  }
  gettimeofday(&end, NULL);
  float indextime = end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1.0 / CLOCKS_PER_SEC;

  UINT32 numresults = 0;
  uint64_t * chunks = new uint64_t[M];
  allocator_Func allor = NULL;
  int32_t* slots = new int32_t[M];
  
  if (strcmp(algorcode, "floor")==0) {
    allor = &floorevenallocator;
  } else if (strcmp(algorcode, "greedy")==0) {
    allor = &greedyallocator;
  } else if (strcmp(algorcode, "rgreedy")==0) {
    allor = &greedyallocatorreduce;
  } else if (strcmp(algorcode, "robin")==0) {
    allor = &roundrobinallocator;
  } else if (strcmp(algorcode, "rrobin")==0) {
    allor = &roundrobinreduceallocator;
  } else if (strcmp(algorcode, "dp")==0) {
    allor = &DP_optimal_allocator;
  } else if (strcmp(algorcode, "rdp")==0) {
    allor = &Reduced_DP_optimal_allocator;
  }

  const char* alloName[] = {"Floor", "Greedy", "Rgreedy", "Robin", "Rrobin", "Dp", "Rdp", NULL};
  for (int al=0; alloName[al]!=NULL;al++) {
    if (strstr(algorcode, alloName[al])==NULL) continue;    
    if (strcmp(alloName[al], "Floor")==0) {
      allor = &floorevenallocator;
    } else if (strcmp(alloName[al], "Greedy")==0) {
      allor = &greedyallocator;
    } else if (strcmp(alloName[al], "Rgreedy")==0) {
      allor = &greedyallocatorreduce;
    } else if (strcmp(alloName[al], "Robin")==0) {
      allor = &roundrobinallocator;
    } else if (strcmp(alloName[al], "Rrobin")==0) {
      allor = &roundrobinreduceallocator;
    } else if (strcmp(alloName[al], "Dp")==0) {
      allor = &DP_optimal_allocator;
    } else if (strcmp(alloName[al], "Rdp")==0) {
      allor = &Reduced_DP_optimal_allocator;
    }

    mlregressor* mlr = NULL;
    // Create Regressor
    if(regfile != NULL)
    {
       mlr = new mlregressor(tau, regfile, MIH->b, M, 2);

       // Generate feature vectors for queries
       for(int i = 0; i < NQ; i++) {
       	split(chunks, inquery->mData + i * inquery->mDataBytes, M, MIH->mplus, MIH->b); 
	mlr->gen_feature_vecs(chunks, i);

       }
    }

    for (int ee = 0; ee <= ext; ee++) {
      for (int t = taustart; t <=tau; t +=taugap) {
	gettimeofday(&begin, NULL);
	// Perform the query processing.
	//fprintf(stderr, "Process %d error\n", t);
	// numresults = MIH->batchrangequery (results, stats, inquery->mData, NQ, inquery->mDataBytes, t);
	qstat *pstats = stats;
	// memset(stats, 0, sizeof(qstat)*NQ);
	UINT32 totalnumres = 0;
	uint64_t allocost = 0;
	uint64_t totalookups = 0;
	uint64_t totalcand = 0;
	uint64_t totaldup = 0;
	
	for (int i=0; i<NQ; i++) {
	  split(chunks, inquery->mData + i * inquery->mDataBytes, M, MIH->mplus, MIH->b);      

          if(regfile != NULL) {
	         // Generate feature vectors for the current query i
	 	 //mlr->gen_feature_vecs(chunks);
		 allocost = Reduced_dpml_allocator(Hist, M, slots, chunks, t+ee, mlr, i);
	  }
	  else
	  	allocost = (*allor)(Hist, M, slots, chunks, t+ee);
	  // if (ext == 0) 
	  //   MIH->rangequerywithallo(results, pstats, inquery->mData + i * inquery->mDataBytes, chunks, slots, t);
	  // else 
	  MIH->rangequerywithalloandext(results, pstats, inquery->mData + i * inquery->mDataBytes, chunks, slots, t, ee);
	  totalnumres += pstats->numres;
	  totalookups += pstats->numlookups;
	  totaldup += pstats->numdups;
	  totalcand += pstats->numcand;
	  if (printcand) {
 	    fprintf(stdout, "INDIV CAND Alg: MIHRING Allor: %s #dim: %lld Tau: %d Ext: %d #Rst %d #Lkup %d #Cnd %d #Dup %d\n",
		    alloName[al], indata->mDim, t, ee, pstats->numres, pstats->numlookups, pstats->numcand, pstats->numdups);
	  }
	  pstats ++;
	}
	
	gettimeofday(&end, NULL);       
	float querytime = end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1.0 / CLOCKS_PER_SEC;
	// fprintf(stderr, "Alg: MIHRING #dim: %lld Tau: %d qrytime %.3f idxtime %.f #data %lld #qry %d K: %d #rst %lld #lkup %lld #cnd %lld #dup %lld\n", indata->mDim, t, querytime, indextime, N, NQ, K, totalnumres, totalookups, totalcand, totaldup);
	fprintf(stderr, "Alg: MIHRING Allor: %s #dim: %lld Tau: %d Ext: %d qrytime %.3f idxtime %.f #data %lld #qry %d K: %d #avgrst %.2f #avglkup %.2f #avgcnd %.2f #avgdup %.2f\n", alloName[al], indata->mDim, t, ee, querytime, indextime, N, NQ, K, totalnumres/(NQ*1.0), totalookups/(NQ*1.0), totalcand/(NQ*1.0), totaldup/(NQ*1.0));
      }
    }
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
