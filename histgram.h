#ifndef HISTGRAM_H__
#define HISTGRAM_H__

#include <stdio.h>
#include <math.h>
#include "types.h"
#include "bucket_group.h"

class histgram {
public:
  static const int MAX_BUCKET;
  static const int MAX_ERROR;
  static const int MAX_COUNT;
  
  //uint32_t **mHistgram;
  uint64_t bucketsize;
  uint64_t bucketbits;
  uint64_t totalcount;
  int64_t  max_error;
  uint64_t valuemask;
  uint32_t *mHistgram;
  // public:
  
  int b;
  
  uint64_t size;

  histgram() {
  }
  
  ~histgram() {
    if (mHistgram !=NULL) {
      /* for (uint32_t i = 0; i < bucketsize; i++ ) { */
      /* 	delete [] mHistgram[i]; */
      /* } */
      delete [] mHistgram;
    }
  }

  void insert(uint64_t index, UINT32 count);
  uint32_t search(uint64_t index, int error);
  void init(int _bucketbits);
  void LoadFromFile(const char *filename);
  void WriteToFile(const char *filename);
};

#endif
