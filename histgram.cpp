#include "histgram.h"


const int histgram::MAX_BUCKET = 20;
const int histgram::MAX_ERROR = 10;
const int histgram::MAX_COUNT = 1024*1024*128;

void histgram::insert(uint64_t index, UINT32 count) {    
  for (uint64_t i = 0; i < bucketsize; i++) {
    uint32_t err = __builtin_popcount((index ^ i) & valuemask);
    if (err < max_error)
      mHistgram[i*max_error + err] += count;
  }
  totalcount += count;
}


uint32_t histgram::search(uint64_t index, int error) {
  if (error < 0) {
    return 0;
  }    
  if (error>=max_error) {
    return MAX_COUNT;
  }
  return mHistgram[((index&valuemask) * max_error) + error];
}


void histgram::init(int _bucketbits) {
  totalcount = 0;
  bucketbits = _bucketbits;
  if (_bucketbits > MAX_BUCKET) {
    bucketbits = MAX_BUCKET;
  }
  bucketsize = pow(2, bucketbits);
  if (_bucketbits > MAX_ERROR) {
    max_error = MAX_ERROR;
  } else {
    max_error = _bucketbits;
  }
  mHistgram = new uint32_t[bucketsize * max_error];
  valuemask=bucketsize-1;
}

void histgram::LoadFromFile(const char *filename) {
  FILE *fp = fopen(filename, "rb");
  
  // Read the dimentionality informations.
  fread(&bucketsize, sizeof(uint64_t), 1, fp);
  fread(&bucketbits, sizeof(uint64_t), 1, fp);
  fread(&totalcount, sizeof(uint64_t), 1, fp);
  fread(&max_error, sizeof(uint64_t), 1, fp);
  fread(&valuemask, sizeof(uint64_t), 1, fp);
  mHistgram = new uint32_t[bucketsize*max_error];
  fread(mHistgram, sizeof(uint32_t), bucketsize * max_error, fp);
  fclose(fp);
}

void histgram::WriteToFile(const char *filename) {
  FILE *fp = fopen(filename, "wb");
  
  // Read the dimentionality informations.
  fwrite(&bucketsize, sizeof(uint64_t), 1, fp);
  fwrite(&bucketbits, sizeof(uint64_t), 1, fp);
  fwrite(&totalcount, sizeof(uint64_t), 1, fp);
  fwrite(&max_error, sizeof(uint64_t), 1, fp);
  fwrite(&valuemask, sizeof(uint64_t), 1, fp);
  fwrite(mHistgram, sizeof(uint32_t), bucketsize * max_error, fp);
  fclose(fp);
}
