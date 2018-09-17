#include <iostream>

using namespace std;

const int MAXERR = 32;
const int BLOCKMAXERR = 32;
const int kMaxBuckets = 64;
const int INFINIT=10000000000;

/* PartEst is to estimate candidate 
 * according to partial bits
 */
class PartEsti {
private:
    // estimate candidate size
    void recursiveEstimate(int i, int32_t** Hist, int error, double cost, int MM, int32_t* totalcost, int N) {
        if(i >= MM) {
	    if(error <= BLOCKMAXERR) {
	        totalcost[error] += cost * N;
	    } else {
	        totalcost[BLICKMAXERR] += cost * N;
	    }
	} else {
	    for(it err = 0; err <= MAXERR; err++) {
	        if(Hist[i][err] > 0) {
		   double localcost = cost * ((double)(Hist[i][err]) / N);
		   recursiveEstimate(i + 1, Hist, error + err, localcost, MM, totalcost, N);
		}
	    }
	}
    }
public:
    int64_t** queryHistMMEsti;
    int64_t** queryHistM;
    int b;
    int mplus;
    // larger partNum
    int M_d_MM;
    // smaller partNum
    int M;
    int MM;
    int32_t* bitsPart;
    HmData* indata;
    HmData* querydata;
    int dataN = -1;
    int queryN = -1;

    PartEsti(int _M_d_MM, int _MM, HmData* indata, HmData* querydata) {
    	this->queryN = _queryN;
	this->b = _b;
	this->M_d_MM = _M_d_MM;
	this->MM = _MM;
	this->M = this->M_d_MM * MM;
        if(dataN < 0 || data > indata->mNumData) {
	    dataN = indata->mNumdata;
	}
	if(queryN < 0 || queryN > querydata->mNumData) {
	    queryN = querydata->mNumData;
	}
	uint8_t* codes = indata->mData;
	uint64_t* chunks = new uint64_t[M];
	b = ceil((double)(indata->mDim) / M);
	mplus = indata->mDim - M * (b - 1);
	uint8_t* querycodes = querydata->mData;
	uint64_t* querychunks = new uint64_t[queryN * M];

	// malloc
	this->queryHistM = new int64_t*[queryN * M];
	this->queryHistMMEsti = new int32_t*[queryN * M_d_MM];
        for (int i = 0; i < queryN; i++) {
            split(querychunks+i*M, querycodes+i*querydata->mDataBytes, M, mplus, b);
            for (int j = 0; j < M; j++) {
                queryHistM[i*M + j] = new int64_t[MAXERR+1];
            }
            for (int j = 0; j < M_d_MM; j++) {
                queryHistMMEsti[i*M_d_MM + j] = new int64_t[BLOCKMAXERR+1];
            }
        }      
 
        // fill into histgrams
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
             } 
          }
        }  
 
        for (int i  = 0; i < queryN; i++) {
            for (int k = 0; k < M_d_MM; k ++) {
            recursiveEstimate(0, (queryHistM + i*M + k*MM), 0, 1, MM, queryHistMMEsti[i*M_d_MM+k], dataN);
            }
        }  
    }    

    /*
    // fill out estimated candidate number
    void fillEstimateCands(histgram* Hist, uint64_t* chunk, int MM) {
	int32_t* candsMM = new int32_t[MM];
	// fill into partial candidate number
        for(int i = 0; i < queryN; i++) {
	    for(int j = 0; j < M_d_MM; j++) {
		// fill in candidate numbers
		for(int z = 0; z < MM; z++) candsMM[z] = Hist[j * MM + z].
		recursiveEstimate(0, ());
	    }
	}
    }
    */

    // allocate thresholds of queries
    int64_t Reduced_DP_optimal_allocator(int qid, int32_t *allo, int tau, int MM) {
	    int64_t** Hist = this->queryHistMMEsti + qid*M_d_MM;
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


};
