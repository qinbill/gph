#include <algorithm>
#include "mihasher.h"

using namespace std;

#define kMaxBuckets 64


#ifdef GREEDYALLO
#define histgramrate 100
#else
#define histgramrate 1000
#endif

/*
 * Inputs: query, numq, dim1queries
 */

/*
 * Outputs: results, numres, stats
 *
 *   results: an array of indices (1-based) of the K-nearest neighbors
 *   for each query. So the array includes K*numq uint32 integers.
 *
 *   numres: includes the number of database entries that fall at any
 *   specific Hamming distance from the query until the K nearest
 *   neighbors are reached. So from this array you can figure out the
 *   Hamming distances of the K-nearest neighbors.
 */

void mihasher::batchquery(UINT32 *results, UINT32 *numres, qstat *stats, UINT8 *queries, UINT32 numq, int dim1queries){
  UINT32 *res  = new UINT32[K*(D+1)];
  uint64_t *chunks = new uint64_t[m];

  UINT32 *presults = results;
  UINT32 *pnumres = numres;
  qstat *pstats = stats;
  UINT8 *pq = queries;
  UINT32 totalnumres = 0;

  for (int i=0; i<numq; i++) {
    query(presults, pnumres, pstats, pq, chunks, res);
    totalnumres += pstats->numres;
    presults += K;
    pnumres += B+1;
    pstats ++;
    pq += dim1queries;
  }
  printf("The Total number of results %d\n", totalnumres);
	
  delete [] res;
  delete [] chunks;

}

UINT32 mihasher::batchrangequery(UINT32 *results, qstat *stats, UINT8 *queries, UINT32 numq, int dim1queries, int tau, UINT8 *origqueries, UINT8 *origcodes, int origdimbytes)
{
    // UINT32 *res  = new UINT32[K*(D+1)];
  uint64_t *chunks = new uint64_t[m];

  UINT32 *presults = new UINT32[N];
  qstat *pstats = stats;
  UINT8 *pq = queries;
  UINT8 *origpq = origqueries;
  UINT32 totalnumres = 0;
  UINT32 totalnumlookup = 0;
  UINT32 totalnumcand = 0;
  UINT32 totalnumdups = 0;

  for (int i=0; i<numq; i++) {
    rangequery(presults, pstats, pq, chunks, tau, origpq, origcodes, origdimbytes);
    totalnumres += pstats->numres;
    totalnumcand += pstats->numcand;
    totalnumlookup += stats->numlookups;
    totalnumdups += stats->numdups;
    pstats ++;
    pq += dim1queries;
    origpq += origdimbytes;
  }

  printf("The Total number of results %d of numquery %d  lookups %d  Cand %d  dups %d\n", totalnumres, numq, totalnumlookup, totalnumcand, totalnumdups);

  delete [] presults;
  delete [] chunks;
  // delete [] res;
  return totalnumres;
}




// Temp variables: chunks, res -- we did not want to malloc inside
// query, so these arrays are passed from outside
void mihasher::rangequery(UINT32 *results, qstat *stats, UINT8 *query, uint64_t *chunks, int tau, UINT8 *origquery, UINT8 *origcodes, int origdimbytes)
{
  
  //  UINT32 n = 0; 		// number of results so far obtained (up to a distance of s per chunk)
  UINT32 nc = 0;       		// number of candidates tested with full codes (not counting duplicates)
  UINT32 nd = 0;              	// counting everything retrieved (duplicates are counted multiple times)
  UINT32 nl = 0;	       	// number of lookups (and xors)
  UINT32 nr = 0;                // number of results so far obtained (up to a distance of s per chunk)


  UINT32 index;
  int hammd;
  clock_t start, end;
  start = clock();
  
  counter->erase(N);
  split(chunks, query, m, mplus, b);
    
  int s;			// the growing search radius per substringo
  int curb = b;	        	// current b: for the first mplus substrings it is b, for the rest it is (b-1)
  //  int t = tau;
  
  int32_t slots[kMaxBuckets];

#ifdef GREEDYALLO
  greedyallocator (slots, chunks, tau);
#else
  for (int k =0; k < m; k ++) {
    if (k < tau % m) {
      slots[k] = (tau / m) + 1;
    } else {
      slots[k] = tau / m;
    }
  }
#endif

#ifdef REDUCTION  
  for (int k =0; k < m-1; k ++) {
    slots[k] = slots[k]-1;
  }
#endif
  
  for (int k = 0; k < m; k ++) {
    if (slots[k]<0) continue;
    
    if (k < mplus)
      curb = b;
    else
      curb = b-1;
    
    //uint64_t chunksk = chunks[k];
    for (s = 0; s <= slots[k]; s++) {
      if (k < mplus)
        curb = b;
      else
        curb = b-1;
      
      uint64_t chunksk = chunks[k];
      nl += xornum[s+1] - xornum[s];	// number of bit-strings with s number of 1s
      
      uint64_t bitstr = 0; 	       	// the bit-string with s number of 1s
      for (int i=0; i<s; i++)
        power[i] = i;			// power[i] stores the location of the i'th 1
      power[s] = curb+1;	       	// used for stopping criterion (location of (s+1)th 1)
      
      int bit = s-1;			// bit determines the 1 that should be moving to the left
      // we start from the left-most 1, and move it to the left until it touches another one
      
      while (true) {			// the loop for changing bitstr
        if (bit != -1) {
          bitstr ^= (power[bit] == bit) ? (uint64_t)1 << power[bit] : (uint64_t)3 << (power[bit]-1);
          power[bit]++;
          bit--;
        } else { // bit == -1
          /* the binary code bitstr is available for processing */
#ifdef USEDENSE
          int size = 0;
          uint64_t *arr;
          arr = DH[k].query(chunksk ^ bitstr, &size); // lookup
#else
          int size = 0;
          uint32_t *arr;
          arr = H[k].query(chunksk ^ bitstr, &size); // lookup
#endif          
          // arr = H[k].query(chunksk ^ bitstr, &size); // lookup
          if (size) {			// the corresponding bucket is not empty
            nd += size;
            for (int c = 0; c < size; c++) {
              index = arr[c];
              if (!counter->get(index)) { // if it is not a duplicate
                counter->set(index);
                // hammd = match(codes + (uint64_t)index*(B_over_8), query, B_over_8);
		hammd = match(origcodes + index * origdimbytes, origquery, origdimbytes);
                nc++;
                if (hammd <= tau) {
                  results[nr++] = index + 1;
                }
              }
            }
          }
          /* end of processing */
          
          while (++bit < s && power[bit] == power[bit+1]-1) {
            bitstr ^= (uint64_t)1 << (power[bit]-1);
            power[bit] = bit;
          }
          if (bit == s)
            break;
        }
      }
    }
  }
  
  end = clock();
  
  // printf("%d\n", nr);
  stats->ticks = end-start;
  stats->numcand = nc;
  stats->numdups = nd;
  stats->numlookups = nl;  
  stats->numres = nr;
}





// UINT32 mihasher::batchrangequeryprojected(UINT32 *results, qstat *stats, UINT8 *queries, UINT32 numq, int dim1queries, int tau)
// {
//     // UINT32 *res  = new UINT32[K*(D+1)];
//   uint64_t *chunks = new uint64_t[m];

//   UINT32 *presults = new UINT32[N];
//   qstat *pstats = stats;
//   UINT8 *pq = queries;
//   UINT32 totalnumres = 0;
//   UINT32 totalnumlookup = 0;
//   UINT32 totalnumcand = 0;
//   UINT32 totalnumdups = 0;

//   for (int i=0; i<numq; i++) {
//     rangequery(presults, pstats, pq, chunks, tau);
//     totalnumres += pstats->numres;
//     totalnumcand += pstats->numcand;
//     totalnumlookup += stats->numlookups;
//     totalnumdups += stats->numdups;
//     pstats ++;
//     pq += dim1queries;
//   }

//   printf("The Total number of results %d of numquery %d  lookups %d  Cand %d  dups %d\n", totalnumres, numq, totalnumlookup, totalnumcand, totalnumdups);

//   delete [] presults;
//   delete [] chunks;
//   // delete [] res;
//   return totalnumres;
// }





// // Temp variables: chunks, res -- we did not want to malloc inside
// // query, so these arrays are passed from outside
// void mihasher::rangequeryprojected(UINT32 *results, qstat *stats, UINT8 *query, uint64_t *chunks, int tau)
// {
  
//   //  UINT32 n = 0; 		// number of results so far obtained (up to a distance of s per chunk)
//   UINT32 nc = 0;       		// number of candidates tested with full codes (not counting duplicates)
//   UINT32 nd = 0;              	// counting everything retrieved (duplicates are counted multiple times)
//   UINT32 nl = 0;	       	// number of lookups (and xors)
//   UINT32 nr = 0;                // number of results so far obtained (up to a distance of s per chunk)


//   UINT32 index;
//   int hammd;
//   clock_t start, end;
//   start = clock();
  
//   counter->erase(N);
//   split(chunks, query, m, mplus, b);
    
//   int s;			// the growing search radius per substringo
//   int curb = b;	        	// current b: for the first mplus substrings it is b, for the rest it is (b-1)
//   //  int t = tau;
  
//   int32_t slots[kMaxBuckets];

// #ifdef GREEDYALLO
//   greedyallocator (slots, chunks, tau);
// #else
//   for (int k =0; k < m; k ++) {
//     if (k < tau % m) {
//       slots[k] = (tau / m) + 1;
//     } else {
//       slots[k] = tau / m;
//     }
//   }
// #endif

// #ifdef REDUCTION  
//   for (int k =0; k < m-1; k ++) {
//     slots[k] = slots[k]-1;
//   }
// #endif
  
//   for (int k = 0; k < m; k ++) {
//     if (slots[k]<0) continue;
    
//     if (k < mplus)
//       curb = b;
//     else
//       curb = b-1;
    
//     //uint64_t chunksk = chunks[k];
//     for (s = 0; s <= slots[k]; s++) {
//       if (k < mplus)
//         curb = b;
//       else
//         curb = b-1;
      
//       uint64_t chunksk = chunks[k];
//       nl += xornum[s+1] - xornum[s];	// number of bit-strings with s number of 1s
      
//       uint64_t bitstr = 0; 	       	// the bit-string with s number of 1s
//       for (int i=0; i<s; i++)
//         power[i] = i;			// power[i] stores the location of the i'th 1
//       power[s] = curb+1;	       	// used for stopping criterion (location of (s+1)th 1)
      
//       int bit = s-1;			// bit determines the 1 that should be moving to the left
//       // we start from the left-most 1, and move it to the left until it touches another one
      
//       while (true) {			// the loop for changing bitstr
//         if (bit != -1) {
//           bitstr ^= (power[bit] == bit) ? (uint64_t)1 << power[bit] : (uint64_t)3 << (power[bit]-1);
//           power[bit]++;
//           bit--;
//         } else { // bit == -1
//           /* the binary code bitstr is available for processing */
// #ifdef USEDENSE
//           int size = 0;
//           uint64_t *arr;
//           arr = DH[k].query(chunksk ^ bitstr, &size); // lookup
// #else
//           int size = 0;
//           uint32_t *arr;
//           arr = H[k].query(chunksk ^ bitstr, &size); // lookup
// #endif          
//           // arr = H[k].query(chunksk ^ bitstr, &size); // lookup
//           if (size) {			// the corresponding bucket is not empty
//             nd += size;
//             for (int c = 0; c < size; c++) {
//               index = arr[c];
//               if (!counter->get(index)) { // if it is not a duplicate
//                 counter->set(index);
//                 hammd = match(codes + (uint64_t)index*(B_over_8), query, B_over_8);
//                 nc++;
//                 if (hammd <= tau) {
//                   results[nr++] = index + 1;
//                 }
//               }
//             }
//           }
//           /* end of processing */
          
//           while (++bit < s && power[bit] == power[bit+1]-1) {
//             bitstr ^= (uint64_t)1 << (power[bit]-1);
//             power[bit] = bit;
//           }
//           if (bit == s)
//             break;
//         }
//       }
//     }
//   }
  
//   end = clock();
  
//   // printf("%d\n", nr);
//   stats->ticks = end-start;
//   stats->numcand = nc;
//   stats->numdups = nd;
//   stats->numlookups = nl;  
//   stats->numres = nr;
// }




// Temp variables: chunks, res -- we did not want to malloc inside
// query, so these arrays are passed from outside
void mihasher::rangequerywithallo(UINT32 *results, qstat *stats, UINT8 *query, uint64_t *chunks, int* slots, int tau)
{
  
  //  UINT32 n = 0; 		// number of results so far obtained (up to a distance of s per chunk)
  UINT32 nc = 0;       		// number of candidates tested with full codes (not counting duplicates)
  UINT32 nd = 0;              	// counting everything retrieved (duplicates are counted multiple times)
  UINT32 nl = 0;	       	// number of lookups (and xors)
  UINT32 nr = 0;                // number of results so far obtained (up to a distance of s per chunk)

  UINT32 index;
  int hammd;
  clock_t start, end;
  start = clock();
  
  counter->erase(N);
  
  int s;			// the growing search radius per substringo
  int curb = b;	        	// current b: for the first mplus substrings it is b, for the rest it is (b-1)
  //  int t = tau;
  
  for (int k = 0; k < m; k ++) {
    if (slots[k]<0) continue;
    
    if (k < mplus)
      curb = b;
    else
      curb = b-1;
    
    //uint64_t chunksk = chunks[k];
    for (s = 0; s <= slots[k]; s++) {
      if (k < mplus)
        curb = b;
      else
        curb = b-1;
      
      uint64_t chunksk = chunks[k];
      nl += xornum[s+1] - xornum[s];	// number of bit-strings with s number of 1s
      
      uint64_t bitstr = 0; 	       	// the bit-string with s number of 1s
      for (int i=0; i<s; i++)
        power[i] = i;			// power[i] stores the location of the i'th 1
      power[s] = curb+1;	       	// used for stopping criterion (location of (s+1)th 1)
      
      int bit = s-1;			// bit determines the 1 that should be moving to the left
      // we start from the left-most 1, and move it to the left until it touches another one
      
      while (true) {			// the loop for changing bitstr
        if (bit != -1) {
          bitstr ^= (power[bit] == bit) ? (uint64_t)1 << power[bit] : (uint64_t)3 << (power[bit]-1);
          power[bit]++;
          bit--;
        } else { // bit == -1
          /* the binary code bitstr is available for processing */
#ifdef USEDENSE
          int size = 0;
          uint64_t *arr;
          arr = DH[k].query(chunksk ^ bitstr, &size); // lookup
#else
          int size = 0;
          uint32_t *arr;
          arr = H[k].query(chunksk ^ bitstr, &size); // lookup
#endif          
          // arr = H[k].query(chunksk ^ bitstr, &size); // lookup
          if (size) {			// the corresponding bucket is not empty
            nd += size;
            for (int c = 0; c < size; c++) {
              index = arr[c];
              if (!counter->get(index)) { // if it is not a duplicate
                counter->set(index);
                hammd = match(codes + (uint64_t)index*(B_over_8), query, B_over_8);
                nc++;
                if (hammd <= tau) {
                  results[nr++] = index + 1;
                }
              }
            }
          }
          /* end of processing */
          
          while (++bit < s && power[bit] == power[bit+1]-1) {
            bitstr ^= (uint64_t)1 << (power[bit]-1);
            power[bit] = bit;
          }
          if (bit == s)
            break;
        }
      }
    }
  }
  
  end = clock();
  
  // printf("%d\n", nr);
  
  stats->ticks = end-start;
  stats->numcand = nc;
  stats->numdups = nd;
  stats->numlookups = nl;  
  stats->numres = nr;
}

void mihasher::rangequerywithalloandext(UINT32 *results, qstat *stats, UINT8 *query, uint64_t *chunks, int* slots, int tau, int ext, int _id)
{
  
  //  UINT32 n = 0; 		// number of results so far obtained (up to a distance of s per chunk)
  UINT32 nc = 0;       		// number of candidates tested with full codes (not counting duplicates)
  UINT32 nd = 0;              	// counting everything retrieved (duplicates are counted multiple times)
  UINT32 nl = 0;	       	// number of lookups (and xors)
  UINT32 nr = 0;                // number of results so far obtained (up to a distance of s per chunk)

  UINT32 index;
  int hammd;
  clock_t start, end;
  start = clock();
  
  bytecounter->erase(N);
  //bytecounter->erase();
  
  int s;			// the growing search radius per substringo
  int curb = b;	        	// current b: for the first mplus substrings it is b, for the rest it is (b-1)
  //  int t = tau;
  
  for (int k = 0; k < m; k ++) {
    if (slots[k]<0) continue;
    
    if (k < mplus)
      curb = b;
    else
      curb = b-1;
    
    //uint64_t chunksk = chunks[k];
    for (s = 0; s <= slots[k]; s++) {
      if (k < mplus)
        curb = b;
      else
        curb = b-1;
      
      uint64_t chunksk = chunks[k];
      nl += xornum[s+1] - xornum[s];	// number of bit-strings with s number of 1s
      
      uint64_t bitstr = 0; 	       	// the bit-string with s number of 1s
      for (int i=0; i<s; i++)
        power[i] = i;			// power[i] stores the location of the i'th 1
      power[s] = curb+1;	       	// used for stopping criterion (location of (s+1)th 1)
      
      int bit = s-1;			// bit determines the 1 that should be moving to the left
      // we start from the left-most 1, and move it to the left until it touches another one
      
      while (true) {			// the loop for changing bitstr
        if (bit != -1) {
          bitstr ^= (power[bit] == bit) ? (uint64_t)1 << power[bit] : (uint64_t)3 << (power[bit]-1);
          power[bit]++;
          bit--;
        } else { // bit == -1
          /* the binary code bitstr is available for processing */
#ifdef USEDENSE
          int size = 0;
          uint64_t *arr;
          arr = DH[k].query(chunksk ^ bitstr, &size); // lookup
#else
          int size = 0;
          uint32_t *arr;
          arr = H[k].query(chunksk ^ bitstr, &size); // lookup
#endif          
          // arr = H[k].query(chunksk ^ bitstr, &size); // lookup
          if (size) {			// the corresponding bucket is not empty
            //printf("SIZE : %d, %d\n", size, _id);
            nd += size;
            for (int c = 0; c < size; c++) {
              index = arr[c];
              if(index >= _id) break;
	      if (bytecounter->checkcount(index, slots[k] - s, ext)) {
		hammd = match(codes + (uint64_t)index*(B_over_8), query, B_over_8);
                nc++;
                if (hammd <= tau) {
                  //results[nr++] = index + 1;
                  nr++;
                }		
	      }
              // if (!counter->get(index)) { // if it is not a duplicate
              //   counter->set(index);
              // }
            }
          }
          /* end of processing */
          
          while (++bit < s && power[bit] == power[bit+1]-1) {
            bitstr ^= (uint64_t)1 << (power[bit]-1);
            power[bit] = bit;
          }
          if (bit == s)
            break;
        }
      }
    }
  }
  end = clock();
  
  // printf("%d\n", nr);
  
  stats->ticks = end-start;
  stats->numcand = nc;
  stats->numdups = nd;
  stats->numlookups = nl;  
  stats->numres = nr;
}





// Temp variables: chunks, res -- we did not want to malloc inside
// query, so these arrays are passed from outside
void mihasher::rangequerywithalloandext(UINT32 *results, qstat *stats, UINT8 *query, uint64_t *chunks, int* slots, int tau, int ext)
{
  
  //  UINT32 n = 0; 		// number of results so far obtained (up to a distance of s per chunk)
  UINT32 nc = 0;       		// number of candidates tested with full codes (not counting duplicates)
  UINT32 nd = 0;              	// counting everything retrieved (duplicates are counted multiple times)
  UINT32 nl = 0;	       	// number of lookups (and xors)
  UINT32 nr = 0;                // number of results so far obtained (up to a distance of s per chunk)

  UINT32 index;
  int hammd;
  clock_t start, end;
  start = clock();
  
  //bytecounter->erase(N);
  bytecounter->erase();
  
  int s;			// the growing search radius per substringo
  int curb = b;	        	// current b: for the first mplus substrings it is b, for the rest it is (b-1)
  //  int t = tau;
  
  for (int k = 0; k < m; k ++) {
    if (slots[k]<0) continue;
    
    if (k < mplus)
      curb = b;
    else
      curb = b-1;
    
    //uint64_t chunksk = chunks[k];
    for (s = 0; s <= slots[k]; s++) {
      if (k < mplus)
        curb = b;
      else
        curb = b-1;
      
      uint64_t chunksk = chunks[k];
      nl += xornum[s+1] - xornum[s];	// number of bit-strings with s number of 1s
      
      uint64_t bitstr = 0; 	       	// the bit-string with s number of 1s
      for (int i=0; i<s; i++)
        power[i] = i;			// power[i] stores the location of the i'th 1
      power[s] = curb+1;	       	// used for stopping criterion (location of (s+1)th 1)
      
      int bit = s-1;			// bit determines the 1 that should be moving to the left
      // we start from the left-most 1, and move it to the left until it touches another one
      
      while (true) {			// the loop for changing bitstr
        if (bit != -1) {
          bitstr ^= (power[bit] == bit) ? (uint64_t)1 << power[bit] : (uint64_t)3 << (power[bit]-1);
          power[bit]++;
          bit--;
        } else { // bit == -1
          /* the binary code bitstr is available for processing */
#ifdef USEDENSE
          int size = 0;
          uint64_t *arr;
          arr = DH[k].query(chunksk ^ bitstr, &size); // lookup
#else
          int size = 0;
          uint32_t *arr;
          arr = H[k].query(chunksk ^ bitstr, &size); // lookup
#endif          
          // arr = H[k].query(chunksk ^ bitstr, &size); // lookup
          if (size) {			// the corresponding bucket is not empty
            nd += size;
            for (int c = 0; c < size; c++) {
              index = arr[c];
	      if (bytecounter->checkcount(index, slots[k] - s, ext)) {
		hammd = match(codes + (uint64_t)index*(B_over_8), query, B_over_8);
                nc++;
                if (hammd <= tau) {
                  //results[nr++] = index + 1;
                }		
	      }
              // if (!counter->get(index)) { // if it is not a duplicate
              //   counter->set(index);
              // }
            }
          }
          /* end of processing */
          
          while (++bit < s && power[bit] == power[bit+1]-1) {
            bitstr ^= (uint64_t)1 << (power[bit]-1);
            power[bit] = bit;
          }
          if (bit == s)
            break;
        }
      }
    }
  }
  end = clock();
  
  // printf("%d\n", nr);
  
  stats->ticks = end-start;
  stats->numcand = nc;
  stats->numdups = nd;
  stats->numlookups = nl;  
  stats->numres = nr;
}





// Temp variables: chunks, res -- we did not want to malloc inside
// query, so these arrays are passed from outside

void mihasher::query(UINT32 *results, UINT32* numres, qstat *stats, UINT8 *query, uint64_t *chunks, UINT32 *res)
{
  UINT32 maxres = K ? K : N;	// if K == 0 that means we want everything to be processed.
                                // So maxres = N in that case. Otherwise K limits the results processed.

  UINT32 n = 0; 		// number of results so far obtained (up to a distance of s per chunk)
  UINT32 nc = 0;       		// number of candidates tested with full codes (not counting duplicates)
  UINT32 nd = 0;              	// counting everything retrieved (duplicates are counted multiple times)
  UINT32 nl = 0;	       	// number of lookups (and xors)
  UINT32 *arr;
  int size = 0;
  UINT32 index;
  int hammd;
  clock_t start, end;

  start = clock();

  counter->erase(N);
  memset(numres, 0, (B+1)*sizeof(*numres));

  split(chunks, query, m, mplus, b);
    
  int s;			// the growing search radius per substringo
  int curb = b;	        	// current b: for the first mplus substrings it is b, for the rest it is (b-1)

  for (s = 0; s <= d && n < maxres; s++) {
    for (int k=0; k<m; k++) {
      if (k < mplus)
        curb = b;
      else
        curb = b-1;
      uint64_t chunksk = chunks[k];
      nl += xornum[s+1] - xornum[s];	// number of bit-strings with s number of 1s

      uint64_t bitstr = 0; 	       	// the bit-string with s number of 1s
      for (int i=0; i<s; i++)
        power[i] = i;			// power[i] stores the location of the i'th 1
      power[s] = curb+1;	       	// used for stopping criterion (location of (s+1)th 1)

      int bit = s-1;			// bit determines the 1 that should be moving to the left
      // we start from the left-most 1, and move it to the left until it touches another one

      while (true) {			// the loop for changing bitstr
        if (bit != -1) {
          bitstr ^= (power[bit] == bit) ? (uint64_t)1 << power[bit] : (uint64_t)3 << (power[bit]-1);
          power[bit]++;
          bit--;
        } else { // bit == -1
          /* the binary code bitstr is available for processing */
          arr = H[k].query(chunksk ^ bitstr, &size); // lookup
          if (size) {			// the corresponding bucket is not empty
            nd += size;
            for (int c = 0; c < size; c++) {
              index = arr[c];
              if (!counter->get(index)) { // if it is not a duplicate
                counter->set(index);
                hammd = match(codes + (uint64_t)index*(B_over_8), query, B_over_8);
                nc++;
                if (hammd <= D && numres[hammd] < maxres) {
                  res[hammd * K + numres[hammd]] = index + 1;
                }
                numres[hammd]++;
              }
            }
          }
          /* end of processing */

          while (++bit < s && power[bit] == power[bit+1]-1) {
            bitstr ^= (uint64_t)1 << (power[bit]-1);
            power[bit] = bit;
          }
          if (bit == s)
            break;
        }
      }

      n = n + numres[s*m+k]; // this line is very tricky ;)
      // The k'th substring (0 based) is the last chance of an
      // item at a Hamming distance of s*m+k to be
      // found. Because if until the k'th substring, an item
      // with distance of s*m+k is not found, then it means that
      // all of the substrings so far have a distance of (s+1)
      // or more, and the remaining substrings have a distance
      // of s or more (total > s*m+k).
	    
      if (n >= maxres)
        break;
    }
  }
  
  end = clock();
  
  stats->ticks = end-start;
  stats->numcand = nc;
  stats->numdups = nd;
  stats->numlookups = nl;
  
  n = 0;
  for (s = 0; s <= D && n < K; s++ ) {
    for (int c = 0; c < numres[s] && n < K; c++)
      results[n++] = res[s*K + c];
  }

  UINT32 total = 0;
  stats->maxrho = -1;
  for (int i=0; i<=B; i++) {
    total += numres[i];
    if (total >= K && stats->maxrho == -1)
      stats->maxrho = i;
  }
  stats->numres = n;
}

mihasher::mihasher(int _B, int _m)
{
  B = _B;
  B_over_8 = B/8;
  m = _m;
  b = ceil((double)B/m);
 
  D = ceil(B/2.0);		// assuming that B/2 is large enough radius to include all of the k nearest neighbors
  d = ceil((double)D/m);
   
  mplus = B - m * (b-1);
  // mplus     is the number of chunks with b bits
  // (m-mplus) is the number of chunks with (b-1) bits

  counter = new bitarray;
  counter->init(1024);

  bytecounter = new ByteCounter;
  bytecounter->init(1024);
    
  xornum = new UINT32 [d+2];
  xornum[0] = 0;
  for (int i=0; i<=d; i++)
    xornum[i+1] = xornum[i] + choose(b, i);
  
#ifdef USEDENSE
  DH = new DenseHashtable[m];
  for (int i=0; i<mplus; i++)
    DH[i].init(b);
  for (int i=mplus; i<m; i++)
    DH[i].init(b-1);  
#else
  H = new SparseHashtable[m];
  // H[i].init might fail
  for (int i=0; i<mplus; i++)
    H[i].init(b);
  for (int i=mplus; i<m; i++)
    H[i].init(b-1);
#endif

  Hist = new histgram[m];
  for (int i=0; i<mplus; i++) {
    Hist[i].init(b);
  }
  for (int i=mplus; i<m; i++) {
    Hist[i].init(b-1);
  }
}

void mihasher::setK(int _K)
{
  K = _K;
}

mihasher::~mihasher()
{
  delete[] xornum;
  delete[] Hist;
  delete counter;
  delete bytecounter;
#ifdef USEDENSE
  delete [] DH;
#else
  delete[] H;  
#endif
}

#ifdef USEDENSE
void mihasher::populate(UINT8 *_codes, UINT32 _N, int dim1codes)
{
  N = _N;
  codes = _codes;
  uint64_t * chunks = new uint64_t[m];

  for (uint64_t i=0; i<N; i++) {
    split(chunks, codes + i * dim1codes, m, mplus, b);
    for (int k=0; k<m; k++) {
      uint64_t entry = chunks[k];
      DH[k].insertcount(entry);
    }
    if (i % (int)ceil(N/1000) == 0) {
      //printf("%.2f%%\r", (double)i/N * 100);
      //fflush(stdout);
    }
  }

  // Allocate space for further insertion. 
  for (int k=0; k<m; k++) {
    DH[k].insertstart();
  }

  for (uint64_t i=0; i<N; i++) {
    split(chunks, codes+ i * dim1codes, m, mplus, b);
    
    for (int k=0; k<m; k++) {
      // Now, we need to conbine those index to form a way.
      uint64_t entry = chunks[k];
      uint64_t padedi = i;
      DH[k].insert(entry, padedi);
      if (i % histgramrate == 0) {
	Hist[k].insert(chunks[k], 1);
      }
    }
    
    if (i % (int)ceil(N/1000) == 0) {
      //printf("%.2f%%\r", (double)i/N * 100);
      //fflush(stdout);
    }
  }

  for (int k=0; k<m; k++) {
    DH[k].insertfinish();
  }
    
  delete [] chunks;
}

void mihasher::populate_one(UINT8 *_codes, UINT32 _id, UINT32 _N, int dim1codes, bool init)
{
  N = _id;
  codes = _codes;
  uint64_t * chunks = new uint64_t[m];

  //for (uint64_t i=0; i<N; i++) {
    uint64_t i = (uint64_t)(_id);
    split(chunks, codes + i * dim1codes, m, mplus, b);
    for (int k=0; k<m; k++) {
      uint64_t entry = chunks[k];
      DH[k].insertcount(entry);
    }
    if (i % (int)ceil(N/1000) == 0) {
      //printf("%.2f%%\r", (double)i/N * 100);
      //fflush(stdout);
    }
  //}

  // Allocate space for further insertion. 
 if(init == true) { 
    for (int k=0; k<m; k++) {
       DH[k].insertstart();
    }
  }

  //for (uint64_t i=0; i<N; i++) {
    split(chunks, codes+ i * dim1codes, m, mplus, b);
    
    for (int k=0; k<m; k++) {
      // Now, we need to conbine those index to form a way.
      uint64_t entry = chunks[k];
      uint64_t padedi = i;
      DH[k].insert(entry, padedi);
      if (i % histgramrate == 0) {
	Hist[k].insert(chunks[k], 1);
      }
    }
    
    if (i % (int)ceil(N/1000) == 0) {
      //printf("%.2f%%\r", (double)i/N * 100);
      //fflush(stdout);
    }
  //}

  for (int k=0; k<m; k++) {
    DH[k].insertfinish();
  }
    
  delete [] chunks;
}

#else

void mihasher::populate_one(UINT8 *_codes, UINT32 _id, UINT32 _N, int dim1codes, bool init)
{
  N = _N;
  codes = _codes;
  uint64_t * chunks = new uint64_t[m];
  //for (uint64_t i=0; i<N; i++) {
  uint64_t i = (uint64_t)(_id);
    split(chunks, codes+i*dim1codes, m, mplus, b);

    for (int k=0; k<m; k++) {
      H[k].insert(chunks[k], i);
      if (i % histgramrate == 0) {
	Hist[k].insert(chunks[k], 1);
      }
    }
    
    if (i % (int)ceil(N/1000) == 0) {
      //printf("%.2f%%", (double)i/N * 100);
      //fflush(stdout);
    }
  //}
  
  delete [] chunks;
}

void mihasher::populate(UINT8 *_codes, UINT32 _N, int dim1codes)
{
  N = _N;
  codes = _codes;
  uint64_t * chunks = new uint64_t[m];
  for (uint64_t i=0; i<N; i++) {
    split(chunks, codes+i*dim1codes, m, mplus, b);

    for (int k=0; k<m; k++) {
      H[k].insert(chunks[k], i);
      if (i % histgramrate == 0) {
	Hist[k].insert(chunks[k], 1);
      }
    }
    
    if (i % (int)ceil(N/1000) == 0) {
      printf("%.2f%%\r", (double)i/N * 100);
      fflush(stdout);
    }
  }
  
  delete [] chunks;
}


#endif

void mihasher::greedyallocator(int32_t *allo, uint64_t * chunks, int tau) {
  uint64_t nextCost[kMaxBuckets];
  uint32_t available_dists = tau;
  uint64_t totalCost = 0;
  
  for (int32_t pid = 0; pid < m; pid ++) {
    allo[pid] = 0;
  }

  for (int32_t pid = 0; pid < m; pid ++) {
    nextCost[pid] = Hist[pid].search(chunks[pid], 0);
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
    nextCost[minpid] =  Hist[minpid].search(chunks[minpid], allo[minpid]+1);
  }
}



