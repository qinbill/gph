#ifndef SPHASHTABLE_H__
#define SPHASHTABLE_H__

#include <stdio.h>
#include <math.h>
#include "types.h"
#include "bucket_group.h"

class SparseHashtable {

 private:
    static const int MAX_B;	// Maximum bits per key before folding the table	

    BucketGroup *table;		// Bins (each bin is an Array object for duplicates of the same key)

 public:

    int b;			// Bits per index
	
    uint64_t size;		// Number of bins

    SparseHashtable();

    ~SparseHashtable();

    int init(int _b);
	
    void insert(uint64_t index, UINT32 data);

    UINT32* query(uint64_t index, int* size);

};

#endif
