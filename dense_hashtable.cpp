#include "dense_hashtable.h"

const int DenseHashtable::MAX_B = 32;

DenseHashtable::DenseHashtable() {
  size = 0;
  tablesize = 0;
  b = 0;
  indexstart = NULL;
  indexlen = NULL;
  table = NULL;
}

int DenseHashtable::init(int _b) {
  b = _b;    
  if ( b > MAX_B || b > sizeof(uint64_t)*8)
    return 1;
    
  size = (((uint64_t)1) << b) + 1;
  indexstart = new UINT32[size];
  indexlen = new UINT32[size];
  return 0;
}

DenseHashtable::~DenseHashtable () {
  delete [] table;
  delete [] indexstart;
  delete [] indexlen;
}

void DenseHashtable::insertcount(uint64_t index) {
  if (index < size)
    indexlen[index] ++;
}
void DenseHashtable::insertstart() {
  uint64_t totalcount = 0;
  for (uint64_t i = 0; i < size; i++) {
    indexstart[i] = totalcount;
    totalcount += indexlen[i];
    indexlen[i] = 0;
  }
  tablesize = totalcount;
  table = new uint64_t[totalcount];
}

void DenseHashtable::insert(uint64_t index, uint64_t padeddata) {
  if (index < size) {
    table[indexstart[index]+indexlen[index]] = padeddata;
    indexlen[index] ++;
  }
}

void DenseHashtable::insertfinish() {
  delete [] indexlen;
  indexlen = NULL;
}

uint64_t* DenseHashtable::query(uint64_t index, int *sizep) {
  if (index < size) {
    *sizep = indexstart[index + 1] - indexstart[index];
    return table + indexstart[index];
  } else {
    *sizep = 0;
    return NULL;
  }
}
