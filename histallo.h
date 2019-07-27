#ifndef __HISTALLO_H__
#define __HISTALLO_H__

#include "histgram.h"

#include "mlregressor.h"

typedef uint64_t(*allocator_Func)(histgram *, int, int32_t*, uint64_t*, int);

/* Use greedy allocation algorithm to allocate */
uint64_t greedyallocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau);

uint64_t floorevenallocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau);

uint64_t roundrobinallocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau);

uint64_t roundrobinreduceallocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau);

uint64_t greedyallocatorreduce(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau);

uint64_t DP_optimal_allocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau);

uint64_t Reduced_DP_optimal_allocator(histgram *Hist, int m, int32_t *allo, uint64_t * chunks, int tau);

uint64_t Reduced_dpml_allocator(histgram *Hist, int m, int32_t *allo, uint64_t *chunk, int tau, mlregressor* mlr, int queryid);

#endif
