#include <getopt.h>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <vector>
#include "sys/types.h"

#include <stdio.h>
#include <algorithm>
#include <time.h>
#include <string.h>

#include "types.h"
#include "hmdata.h"
#include "bitops.h"

using namespace std;

void print_usage(){
  fprintf(stderr, "usage: -d <Location of Data file>\n");
  fprintf(stderr, "       -q <Location of Query file>\n");
  fprintf(stderr, "       -m number of parts for mih\n");
  fprintf(stderr, "       -M number of parts as a block\n");
  fprintf(stderr, "       -N number of data documents\n");
  fprintf(stderr, "       -Q number of query documents\n");
  fprintf(stderr, "\n");
  exit(0);
}

const int MAXERR = 32;
const int BLOCKMAXERR = 32;
const int kMaxBuckets = 64;
const int INFINIT=10000000000;

// #define myabs(x) ((x)<0 ? -(x) : (x))

int64_t myabs(int64_t x) {
  if (x < 0) {
    return -x;
  }else {
    return x;
  }
}

int64_t greedyallo(int32_t *allo, int64_t** Hist, int tau, int MM) {
  int64_t nextCost[kMaxBuckets];
  int32_t available_dists = tau;
  int64_t totalCost = 0;
  
  for (int32_t pid = 0; pid < MM; pid ++) {
    allo[pid] = -1;
  }

  for (int32_t pid = 0; pid < MM; pid ++) {
    nextCost[pid] = Hist[pid][0];
  }
  
  while (available_dists > 0) {
    int minpid = 0;
    for (int pid = 1; pid < MM; pid ++) {
      if (nextCost[pid] < nextCost[minpid]) {
	minpid = pid;
      }
    }
    // Then update the corrent find pid.
    totalCost += nextCost[minpid];
    allo[minpid] ++;
    available_dists --;
    nextCost[minpid] = Hist[minpid][allo[minpid]+1>BLOCKMAXERR?BLOCKMAXERR:allo[minpid]+1];
  }
  return totalCost;

}


#define kMaxBuckets 256
int64_t histsearch(int64_t** Hist, int pid, int err){
  if (err > BLOCKMAXERR) {
    return INFINIT;
  } else if (err < 0) {
    return 0;
  } else {
    return Hist[pid][err];
  }
}

int64_t Reduced_DP_optimal_allocator(int32_t *allo, int64_t** Hist, int tau, int MM) {
  // uint64_t Reduced_DP_optimal_allocator(int32_t *allo, uint64_t * chunks, int tau) {
  // int32_t available_dists = tau;
  int32_t pid = 0;
  int32_t ttau=tau;
  int64_t DistCost[kMaxBuckets][kMaxBuckets*2];
  int16_t DistPath[kMaxBuckets][kMaxBuckets*2];
  int64_t cumulativecost= 0; 

  for (int32_t dst = 0; dst <= ttau; dst++) {
    cumulativecost += histsearch(Hist, 0, dst -1);
    DistCost[0][dst] = cumulativecost;
    DistPath[0][dst] = dst;
  }

  // Process Intermediate buckets.
  for (pid = 1; pid < MM; ++pid) {
    for (int32_t dst = 0; dst <= ttau; ++dst) {
      cumulativecost = histsearch(Hist, pid, -1);
      DistCost[pid][dst] = DistCost[pid-1][dst] + cumulativecost;
      DistPath[pid][dst] = 0;
      for (int32_t err = 1; err <= ttau && dst - err >=0; err ++) {
	cumulativecost += histsearch(Hist, pid, err-1);
        int64_t current_cost = DistCost[pid-1][dst-err] + cumulativecost;
        if (DistCost[pid][dst] > current_cost) {
          DistCost[pid][dst] = current_cost;
          DistPath[pid][dst] = err;
        }
      }
    }
  }
  
  // Backtrace the error path.
  int32_t err = ttau;
  for (pid = MM - 1; pid >= 0; --pid) {
    allo[pid] = DistPath[pid][err] -1;
    err = err - DistPath[pid][err];
  }

  // for (err = 0; err <= ttau; err ++) {
  //   fprintf(stderr, "\nRederr %d-> ", err);
  //   for (pid = 0; pid < MM; pid ++) {
  //     fprintf(stderr, "(%lld|%d) ", DistCost[pid][err], DistPath[pid][err]);
  //   }
  //   fprintf(stderr, "\n");
  // }  
  return DistCost[MM-1][ttau];  
}




void recursiveEstimate(int i, int64_t**Hist, int error, double cost, int MM, int64_t* totalcost, int N) {
  if (i >= MM) {
    if (error <= BLOCKMAXERR) {
      totalcost[error] += cost * N;
    } else {
      totalcost[BLOCKMAXERR] += cost * N;
    }
  } else {
    for (int err = 0; err <= MAXERR; err++) {
      if (Hist[i][err] > 0) {
        double localcost = cost * ((double)(Hist[i][err])/ N);
        recursiveEstimate(i+1, Hist, error+err, localcost, MM, totalcost, N);
      }
    }
  }
}

int main(int argc, char* argv[])
{
  struct timeval begin, end;
  char c;    
  char *dataFile = NULL;
  char *queryFile = NULL;
  int M = 2;
  int MM = 1;
  int64_t dataN = -1;
  int queryN = -1;
  while ((c = getopt(argc,argv, "hd:m:N:q:M:Q:")) != -1)
  {
    switch (c)
    {
      case 'd':
        dataFile = optarg;
        break;
      case 'q':
        queryFile = optarg;
        break;
      case 'Q':
        queryN = atoi(optarg);
        break;
        // case 'o':
        //   outputprefix = optarg;
        //   break;
      case 'm':
        M = atoi(optarg);
        break;
      case 'M':
        MM = atoi(optarg);
        break;        
      case 'N':
        dataN = atoi(optarg);
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
  
  if (dataFile == NULL || queryFile == NULL){
    cerr << "Need data input and query input file name" <<endl;
    print_usage();
  }
  
  HmData *indata = new HmData(dataFile);
  HmData *querydata = new HmData(queryFile);

  if (dataN < 0 || dataN > indata->mNumData) {
    dataN = indata->mNumData;
  }

  if (queryN < 0 || queryN > querydata->mNumData) {
    queryN = querydata->mNumData;
  }
  
  uint8_t* codes = indata->mData;
  uint64_t * chunks = new uint64_t[M];
  int b = ceil((double)indata->mDim/M);
  int mplus = indata->mDim - M * (b-1);
  int M_d_MM = M / MM;

  uint8_t* querycodes = querydata->mData;
  uint64_t* querychunks = new uint64_t[queryN * M];

  int64_t** queryHistM = new int64_t*[queryN * M];
  int64_t** queryHistMM = new int64_t*[queryN * M_d_MM];
  int64_t** queryHistMMEsti = new int64_t*[queryN * M_d_MM];
  
  for (int i = 0; i < queryN; i++) {
    split(querychunks+i*M, querycodes+i*querydata->mDataBytes, M, mplus, b);
    for (int j = 0; j < M; j++) {
      queryHistM[i*M + j] = new int64_t[MAXERR+1];
    }
    for (int j = 0; j < M_d_MM; j++) {
      queryHistMM[i*M_d_MM + j] = new int64_t[BLOCKMAXERR+1];
      queryHistMMEsti[i*M_d_MM + j] = new int64_t[BLOCKMAXERR+1];
    }
  }
  
  gettimeofday(&begin, NULL);
  for (int64_t i=0; i<dataN; i++) {
    // printf("process doc %lld of %lld\n", i, N);
    split(chunks, codes +i*indata->mDataBytes, M, mplus, b);
    for (int q = 0; q < queryN; q ++) {
      int blockerr = 0;
      for (int k=0; k<M; k++) {
        int err = match((UINT8*)(chunks+k), (UINT8*)(querychunks+q*M+k), 8);
        blockerr += err;
        if (err<=MAXERR) {
          queryHistM[q*M+k][err] ++;
        }else{
          // printf("largeerror %d\n", err);
          queryHistM[q*M+k][MAXERR]++;
        }
        if ((k+1)%MM==0){
          int blockid = ((k+1)/MM)-1;
          if (blockerr<=BLOCKMAXERR){
            queryHistMM[q*M_d_MM+blockid][blockerr] ++;
          } else{
            // printf("largeerror %d\n", blockerr);
            queryHistMM[q*M_d_MM+blockid][BLOCKMAXERR] ++;
          }
          blockerr = 0;
        }
      }
    }
    if (i % (int)ceil(dataN/1000) == 0) {
      printf("%.2f%%\r", (double)i/dataN * 100);
      fflush(stdout);
    }
  }
  
  for (int i  = 0; i < queryN; i++) {
    for (int k = 0; k < M_d_MM; k ++) {
      recursiveEstimate(0, (queryHistM + i*M + k*MM), 0, 1, MM, queryHistMMEsti[i*M_d_MM+k], dataN);
    }
  }  
  
  gettimeofday(&end, NULL);
  float indextime = end.tv_sec - begin.tv_sec + (end.tv_usec - begin.tv_usec) * 1.0 / CLOCKS_PER_SEC;
  
  fprintf(stdout, "Histgenapp for data %s M= %d time= %.3f N= %lld histgramrate= %d\n",
          dataFile, M, indextime, dataN, 1);
  
  for (int i  = 0; i < queryN; i++) {
    printf("Query %d -> [", i);
    for (int k = 0; k < M; k ++) {
      printf("M%d->(", k);
      for (int err = 0; err <= MAXERR; err ++) {
        printf("%lld,", queryHistM[i*M+k][err]);        
      }
      printf("),");
    }
    for (int k = 0; k < M_d_MM; k ++) {
      printf("MM%d->(", k);
      for (int err = 0; err <= BLOCKMAXERR; err ++) {
        printf("%lld,", queryHistMM[i*M_d_MM+k][err]);        
      }
      printf("),");
    }

    for (int k = 0; k < M_d_MM; k ++) {
      printf("MMEsti%d->(", k);
      for (int err = 0; err <= BLOCKMAXERR; err ++) {
        printf("%lld,", queryHistMMEsti[i*M_d_MM+k][err]);        
      }
      printf("),");
    }
    printf(" ]\n");
  }

  for (int err = 0; err < BLOCKMAXERR; err+=4) {
    int64_t totalGrountCost = 0;
    int64_t totalEstiCost = 0;
    int64_t totalPartCost = 0;
    
    for (int i = 0; i < queryN; i++) {
      int32_t allo[kMaxBuckets];
      int32_t allo2[kMaxBuckets];
      // Get greedy ground truth.
      //uint64_t grandcost = greedyallo(allo, queryHistMM + i*M_d_MM, err + 1, M_d_MM);
      int64_t grandcost = Reduced_DP_optimal_allocator(allo, queryHistMM + i*M_d_MM, err + 1, M_d_MM);
      grandcost = 0;
      for (int k = 0; k < M_d_MM; k++) {
        if (allo[k]>0 && allo[k]<=BLOCKMAXERR){
	  for (int err=0;err<=allo[k]; err++) 
	    grandcost+=queryHistMM[i*M_d_MM+k][err];
        }
      }
      totalGrountCost += grandcost;

      printf("Error: %d, GrandAllo->(", err);
      for (int k = 0; k < M_d_MM; k++) {
        printf("%d, ", allo[k]);
      }
      printf(") ");
      
      int64_t costEsti = Reduced_DP_optimal_allocator(allo, queryHistMMEsti + i*M_d_MM, err + 1, M_d_MM);
      //uint64_t costEsti = greedyallo(allo, queryHistMMEsti + i*M_d_MM, err + 1, M_d_MM);
      costEsti = 0;
      for (int k = 0; k < M_d_MM; k++) {
        if (allo[k]>0 && allo[k]<=BLOCKMAXERR){
	  for (int err=0;err<=allo[k]; err++) 
	    costEsti+=queryHistMM[i*M_d_MM+k][err];
        }
      }
      totalEstiCost += costEsti;
      printf("Error: %d, Estiallo->(", err);
      for (int k = 0; k < M_d_MM; k++) {
        printf("%d, ", allo[k]);
      }
      printf(") ");

      
      //uint64_t costPart = greedyallo(allo, queryHistM+i*M, err + 1, M);
      int64_t costPart = Reduced_DP_optimal_allocator(allo, queryHistM+i*M, err + 1, M);
      
      for (int k = 0; k < M_d_MM; k++) {
        allo2[k] = -1;
      }
      for (int k = 0; k < M; k++) {
        allo2[k/MM] += allo[k]-(-1);
      }      
      costPart = 0;
      for (int k = 0; k < M_d_MM; k++) {
        if (allo2[k]>0 && allo2[k]<=BLOCKMAXERR){
	  for (int err=0;err<=allo2[k]; err++) 
          costPart+=queryHistMM[i*M_d_MM+k][err];
        }
      }
      totalPartCost += costPart;
      printf("Error: %d, Partallo->(", err);
      for (int k = 0; k < M_d_MM; k++) {
        printf("%d, ", allo2[k]);
      }
      printf(") ");

      if (grandcost > costEsti) {
        printf(" Grand>Esti");
      }
      printf("\n");
    }     
    double estiErrate = (totalGrountCost > 0)? (myabs(totalEstiCost-totalGrountCost)*1.0/totalGrountCost):0;
    double partErrate = (totalGrountCost > 0)? (myabs(totalPartCost-totalGrountCost)*1.0/totalGrountCost):0;
    
    
    fprintf(stderr, "Cost for Error %d: Ground: %lld Esti: %lld RelErr: %.4f, Part: %lld RelErr: %.4f\n", err, totalGrountCost, totalEstiCost, estiErrate, totalPartCost, partErrate);    
  }

  delete [] chunks;
  delete indata;
  delete querydata;
  return 0;
    print(mnist_train.shape, mnist_test.shape)
}
