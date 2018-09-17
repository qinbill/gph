#include "histgram.h"
#include "histallo.h"

#include "mlregressor.h"

#define kMaxBuckets 64

uint32_t combine_taus[50];

// Calcualte the value of selecting K from N
int cal_combination(int n, int k)
{
	if(k == 0)
	       	return 1;
	return (n * cal_combination(n - 1, k - 1)) / k;
}

void assign_combine_taus(int N)
{
	int maxtau = N / 2 + 1;
	for(int tau = 0; tau <= maxtau; tau++)
	{
		combine_taus[tau] = cal_combination(N, tau);
	}
}

uint64_t greedyallocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau) {
  uint64_t nextCost[kMaxBuckets];
  uint32_t available_dists = tau;
  uint64_t totalCost = 0;
  
  for (int32_t pid = 0; pid < m; pid ++) {
    allo[pid] = 0;
  }

  for (int32_t pid = 0; pid < m; pid ++) {
    totalCost +=  Hist[pid].search(chunks[pid], 0);
    nextCost[pid] = Hist[pid].search(chunks[pid], 1);    
  }
  
  while (available_dists > 0) {
    int minpid = 0;
    for (int pid = 1; pid < m; pid ++) {
      if (nextCost[pid] < nextCost[minpid]) {
	minpid = pid;
      }
    }
    // Then update the corrent find pid.
    totalCost += nextCost[minpid];
    allo[minpid] ++;
    available_dists --;
    nextCost[minpid] = Hist[minpid].search(chunks[minpid], allo[minpid]+1);
  }
  return totalCost;
}


uint64_t floorevenallocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau) {
  uint64_t totalCost = 0;  
  for (int32_t pid = 0; pid < m; pid ++) {
    allo[pid] = floor(tau*1.0/m);
    for (int err = 0; err <= allo[pid]; err++) 
      totalCost += Hist[pid].search(chunks[pid], err);
  }
  return totalCost;
}


uint64_t roundrobinallocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau) {
  uint64_t totalCost = 0;
  for (int32_t pid =0; pid < m; pid ++) {
    if (pid < tau % m) {
      allo[pid] = (tau / m) + 1;
    } else {
      allo[pid] = tau / m;
    }
    for (int err = 0; err <= allo[pid]; err++) 
      totalCost += Hist[pid].search(chunks[pid], err);
  }
  return totalCost;
}



uint64_t roundrobinreduceallocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau) {
  uint64_t totalCost = 0;
  uint64_t mincost = 0xffffffff;
  int minpid = -1;
  for (int32_t pid =0; pid < m; pid ++) {
    if (pid < tau % m) {
      allo[pid] = (tau / m) + 1;
    } else {
      allo[pid] = tau / m;
    }
    if (Hist[pid].search(chunks[pid], allo[pid])<mincost) {
      mincost = Hist[pid].search(chunks[pid], allo[pid]);
      minpid = pid;
    }
  }
  for (int32_t pid = 0; pid < m; pid ++) {
    if (pid != minpid) {
      allo[pid] = allo[pid]-1;
    }
    for (int err = 0; err <= allo[pid]; err++) 
      totalCost += Hist[pid].search(chunks[pid], err);    
  }
  return totalCost;
}


uint64_t greedyallocatorreduce(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau) {
  uint64_t nextCost[kMaxBuckets];
  uint32_t available_dists = tau-(m-1)+m;
  uint64_t totalCost = 0;
  
  for (int32_t pid = 0; pid < m; pid ++) {
    allo[pid] = -1;
  }

  for (int32_t pid = 0; pid < m; pid ++) {
    nextCost[pid] = Hist[pid].search(chunks[pid], 0); // + 0.01 * combine_taus[0];
  }
  
  while (available_dists > 0) {
    int minpid = 0;
    for (int pid = 1; pid < m; pid ++) {
      if (nextCost[pid] < nextCost[minpid]) {
	minpid = pid;
      }
    }
    // Then update the corrent find pid.
    totalCost += nextCost[minpid];
    allo[minpid] ++;
    available_dists --;
    nextCost[minpid] = Hist[minpid].search(chunks[minpid], allo[minpid]+1); // + 0.01 * combine_taus[allo[minpid] + 1];
  }

/*
  // display the allo threshold
  fprintf(stderr, "%d\t", tau);
  for(int i = 0; i < m; i ++)
  {
	  fprintf(stderr, "%d\t", allo[i]);
  }
  fprintf(stderr, "\n");
*/
  return totalCost;
}

// Written by Yaoshu
// Naive DP of Machine Learning
uint64_t Reduced_dpml_allocator(histgram *Hist, int m, int32_t *allo, uint64_t *chunk, int tau, mlregressor* mlr, int queryid)
{
	int32_t pid = 0;
	int32_t ttau = tau + 1;
	double DistCost[kMaxBuckets][kMaxBuckets * 2];
	int64_t DistPath[kMaxBuckets][kMaxBuckets * 2];
	double cumulativecost = 0.0;

	for(int32_t dst = 0; dst <= ttau; dst++) {
		cumulativecost += mlr->search(0, dst - 1, queryid);
		DistCost[0][dst] = cumulativecost;
		//printf("1. %f ", cumulativecost);
		DistPath[0][dst] = dst;	
	}

	// Process Intermediate buckets.
	for(pid = 1; pid < m; pid++)
	{
		for(int32_t dst = 0; dst <= ttau; ++dst) {
			cumulativecost = mlr->search(pid, -1, queryid);
			//printf("2. %f ", cumulativecost);
			DistCost[pid][dst] = DistCost[pid - 1][dst] + cumulativecost;
			DistPath[pid][dst] = 0;
			for(int32_t err = 1; err <= ttau && dst - err >= 0; err ++) {
				cumulativecost += mlr->search(pid, err - 1, queryid);
				//printf("3. %f ", cumulativecost);
				double current_cost = DistCost[pid-1][dst-err] + cumulativecost;
				if(DistCost[pid][dst] > current_cost) {
					DistCost[pid][dst] = current_cost;
					DistPath[pid][dst] = err;
				}
			}
		}
	}

	// Backrace the error path.
	int32_t err = ttau;
	for(pid = m - 1; pid >= 0; --pid) {
		allo[pid] = DistPath[pid][err] - 1;
       		//fprintf(stderr, "allo: %d \n", allo[pid]);
		err = err - DistPath[pid][err];
	}
/*
  for (err = 0; err <= ttau; err ++) {
    fprintf(stderr, "\nRederr %d-> ", err);
    for (pid = 0; pid < m; pid ++) {
       fprintf(stderr, "(%lf|%d) ", DistCost[pid][err], DistPath[pid][err]);
     }
     fprintf(stderr, "\n");
   } 
  
	return (uint64_t)DistCost[m-1][ttau];
*/

/*
	// Malloc a DP table
	uint64_t** DP = new uint64_t*[m + 1];
	int **path = new int*[m + 1];
	for(int i = 0; i < m; i++)
	{
		DP[i] = new uint64_t[tau + 1];
		path[i] = new int[tau + 1];
	}	
	// Initialize the table
	for(int i = 0; i < tau + 1; i++) DP[0][i] = 0;
	for(int i = 1; i < m + 1; i++)
	{
		for(int j = 0; j < tau + 1; j++)
		{
			uint64_t cur_cost = 0;
			int min_tau = 0;
			for(int z = 0; z <= j; z++)
			{
				uint64_t T = DP[i - 1][z] + (mlr->df[j - z])(mlr->ft[i - 1]);
				if(T <= cur_cost)
				{
					cur_cost = T;
					min_tau = z;
				}
			}
			DP[i][j] = cur_cost;
			path[i][j] = min_tau;
		}
	}
	
	// Backtracking to find the results
	int t = tau;
	for(int i = m; i > 0; i++)
	{
		int c = t - path[i][t];
		allo[i - 1] = c;
		t = path[i][t];
	}
	
	return DP[m][tau];
*/
}


uint64_t DP_optimal_allocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau) {
  // int32_t available_dists = tau;
  int32_t pid = 0;
  int64_t DistCost[kMaxBuckets][kMaxBuckets*2];
  int16_t DistPath[kMaxBuckets][kMaxBuckets*2];
  int64_t cumulativecost= 0; 
  
  for (int32_t dst = 0; dst <= tau; dst++) {
    cumulativecost += Hist[0].search(chunks[0], dst);
    DistCost[0][dst] = cumulativecost;
    DistPath[0][dst] = dst;
  }

  // Process Intermediate buckets.
  for (pid = 1; pid < m; ++pid) {
    for (int32_t dst = 0; dst <= tau; ++dst) {
      cumulativecost = Hist[pid].search(chunks[pid], 0);
      DistCost[pid][dst] = DistCost[pid-1][dst] + Hist[pid].search(chunks[pid], 0);
      DistPath[pid][dst] = 0;
      for (int32_t err = 1; err <= tau && dst - err >=0; err ++) {
	cumulativecost +=  Hist[pid].search(chunks[pid], err);
        int64_t current_cost = DistCost[pid-1][dst-err] + cumulativecost;
        if (DistCost[pid][dst] > current_cost) {
          DistCost[pid][dst] = current_cost;
          DistPath[pid][dst] = err;
        }
      }
    }
  }
  
  // Backtrace the error path.
  int32_t err = tau ;
  for (pid = m - 1; pid >= 0; --pid) {
    allo[pid] = DistPath[pid][err];
    err = err - DistPath[pid][err];
  }

  // for (err = 0; err <= tau; err ++) {
  //   fprintf(stderr, "\nDP err %d-> ", err);
  //   for (pid = 0; pid < m; pid ++) {
  //     fprintf(stderr, "(%lld|%d) ", DistCost[pid][err], DistPath[pid][err]);
  //   }
  //   fprintf(stderr, "\n");
  // }
  
  return DistCost[m-1][tau];  
}



uint64_t Reduced_DP_optimal_allocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau) {
  // int32_t available_dists = tau;
  int32_t pid = 0;
  int32_t ttau=tau+1;
  int64_t DistCost[kMaxBuckets][kMaxBuckets*2];
  int16_t DistPath[kMaxBuckets][kMaxBuckets*2];
  int64_t cumulativecost= 0; 

  for (int32_t dst = 0; dst <= ttau; dst++) {
    cumulativecost += Hist[0].search(chunks[0], dst-1);
    DistCost[0][dst] = cumulativecost;
    DistPath[0][dst] = dst;
  }

  // Process Intermediate buckets.
  for (pid = 1; pid < m; ++pid) {
    for (int32_t dst = 0; dst <= ttau; ++dst) {
      cumulativecost = Hist[pid].search(chunks[pid], -1);
      DistCost[pid][dst] = DistCost[pid-1][dst] + cumulativecost;
      DistPath[pid][dst] = 0;
      for (int32_t err = 1; err <= ttau && dst - err >=0; err ++) {
	cumulativecost +=  Hist[pid].search(chunks[pid], err-1);
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
  for (pid = m - 1; pid >= 0; --pid) {
    allo[pid] = DistPath[pid][err] -1;
    err = err - DistPath[pid][err];
  }

/*
  // display histogram
  fprintf(stderr, "histogram Info: \n");
  for(pid = 0; pid < m; ++ pid) { 
	 fprintf(stderr, "Partition %d-> ", pid);
	 for(int k = 0; k <= tau; k++) {
		fprintf(stderr, "(%d, %d)\t",Hist[pid].search(chunks[pid], k), k);
	 }
	 fprintf(stderr, "\n");
  }  

  // display the allo threshold
  fprintf(stderr, "threshold: %d\t", tau);
  for(int i = 0; i < m; i ++)
  {
	  fprintf(stderr, "%d\t", allo[i]);
  }
  fprintf(stderr, "\n");

   for (err = 0; err <= ttau; err ++) {
     fprintf(stderr, "\nRederr %d-> ", err);
     for (pid = 0; pid < m; pid ++) {
       fprintf(stderr, "(%lld|%d) ", DistCost[pid][err], DistPath[pid][err]);
     }
     fprintf(stderr, "\n");
   }
  */
  return DistCost[m-1][ttau];  
}

