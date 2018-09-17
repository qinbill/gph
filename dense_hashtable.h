#ifndef DENSEHTABLE_H__
#define DENSEHTABLE_H__

#include <stdio.h>
#include <math.h>
#include "types.h"
#include "bucket_group.h"

class DenseHashtable {

 private:
    static const int MAX_B;	// Maximum bits per key before folding the table	

 public:

    int b;			// Bits per index
	
    uint64_t size;		// Number of bins

    uint64_t tablesize;

    UINT32 *indexstart;

    UINT32 *indexlen;

    uint64_t *table;

    DenseHashtable();

    ~DenseHashtable();

    int init(int _b);

    void insertcount(uint64_t index);

    void insertstart();
    
    void insertfinish();
	
    void insert(uint64_t index, uint64_t padeddata);

    uint64_t* query(uint64_t index, int* size);

};

#endif
