#ifndef __BYTECOUNTER_H
#define __BYTECOUNTER_H

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "types.h"


class ByteCounter {

public:

  UINT8 *arr;
  UINT32 length;
	
  ByteCounter() {
    arr = NULL;
    length = 0;
  }
	
  ByteCounter(uint64_t _bits) {
    init(_bits);
  }
	
  void init(uint64_t _bytes) {
    length =  _bytes * 2; // (UINT32)ceil(_bits/2);
    arr = new UINT8[length];
    erase();
  }
  
  ~ByteCounter() {
    if (arr)
      delete[] arr;
  }

  inline int addcount(uint64_t index, int count) {
    arr[(index<<1)]++;
    arr[((index<<1)) + 1] = arr[(index<<1) + 1] + (UINT8)count;
    return arr[(index<<1)] + arr[(index<<1)+1];
  }

  inline bool checkcount(uint64_t index, int count, int ext) {
    bool precheck = false;
    if (arr[(index<<1)] + arr[(index<<1)+1] > ext) {
      precheck = true; // Already satisfied.
    }
    arr[(index<<1)]++;
    arr[(index<<1) + 1] = arr[(index<<1) + 1] + (UINT8)count;
    if ((arr[(index<<1)] + arr[(index<<1)+1] > ext) && !precheck) {
      return true;
    } else {
      return false;
    }
  }
  
  
  inline UINT8 get(uint64_t index) {
    return arr[(index<<1)] + arr[(index<<1)+1];
  }
	
  inline void erase() {
    memset(arr, 0, sizeof(UINT8) * length);
  }

  inline void erase(uint64_t _bytes) {
    if ((length > 1) < _bytes){
      delete[] arr;
      init(_bytes);
    } else {
      memset(arr, 0, sizeof(UINT8) * length);
    }
  }  
};

#endif
