#ifndef __MIHASHER_H
#define __MIHASHER_H

#include <stdio.h>
#include <pthread.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include "types.h"
#include "bitops.h"
#include "histgram.h"

#include "sparse_hashtable.h"
#include "dense_hashtable.h"
#include "bitarray.h"
#include "bytecounter.h"

struct qstat {

    UINT32 numres;		// Total number of returned results
    UINT32 numcand;		// Number of hamming distance computations executed
    UINT32 numdups;		// Number of candidates skipped because they were duplicates
    UINT32 numlookups;
    UINT32 maxrho;		// Largest distance that was searched exhaustively
    clock_t ticks;		// Number of clock ticks spent on each query
	
};

class mihasher {
 public:

    int B;			// Bits per code

    int B_over_8;

    int b;			// Bits per chunk (must be less than 64)

    int m;			// Number of chunks



    int D;			// Maximum hamming search radius (we use B/2 by default)

    int d;			// Maximum hamming search radius per substring

    int K;			// Maximum results to return

    uint64_t N;			// Number of codes
	
    UINT8 *codes;		// Table of original full-length codes

    /* is not thread safe */
    bitarray *counter;		// Counter for eliminating duplicate results

    ByteCounter *bytecounter;   // byte counter for candidate counting. 
	
    SparseHashtable *H;		// Array of m hashtables;

    DenseHashtable *DH;		// Array of m hashtables;

    histgram *Hist;		// Array of m hashtables;
		
    UINT32 *xornum;		// Volume of a b-bit Hamming ball with radius s (for s = 0 to d)

    int power[100];		// Used within generation of binary codes at a certain Hamming distance

    // public:
    int mplus;			// Number of chunks with b bits (have 1 bit more than others)
    
    mihasher();

    ~mihasher();

    mihasher(int B, int m);

    void setK(int K);

    void populate(UINT8 *codes, UINT32 N, int dim1codes);

    void batchquery (UINT32 *results, UINT32 *numres, qstat *stats, UINT8 * q, UINT32 numq, int dim1queries);



    void query(UINT32 *results, UINT32* numres, qstat *stats, UINT8 *q, uint64_t * chunks, UINT32 * res);

    //    UINT32 batchrangequery (UINT32 *results, qstat *stats, UINT8 * q, UINT32 numq, int dim1queries, int tau);
    UINT32 batchrangequery (UINT32 *results, qstat *stats, UINT8 *queries, UINT32 numq, int dim1queries, int tau, UINT8 *origqueries, UINT8 *origcodes, int origdimbytes);
    void rangequery(UINT32 *results, qstat *stats, UINT8 *query, uint64_t *chunks, int tau, UINT8 *origquery, UINT8 *origcodes, int origdimbytes);
    // void rangequery(UINT32 *results, qstat *stats, UINT8 *q, uint64_t * chunks, int tau);

    void rangequerywithallo(UINT32 *results, qstat *stats, UINT8 *query, uint64_t *chunks, int* slots, int tau);
    void rangequerywithalloandext(UINT32 *results, qstat *stats, UINT8 *query, uint64_t *chunks, int* slots, int tau, int ext);

    void greedyallocator(int32_t *allo, uint64_t* chunks, int tau);
};

#endif
