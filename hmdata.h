#ifndef __HMDATA_H__
#define __HMDATA_H__

#include <stdint.h>

#define kMaxDim  32*1024

class HmData {
 public:
  uint64_t mDim; // The member dimension.
  uint8_t *mData;
  uint64_t mDataBytes;
  uint64_t mNumData;

  HmData() {
    mDim = 0;
    mData = NULL;
    mDataBytes = 0;
    mNumData = 0;
  }

  HmData(const char* filename);  

  ~HmData();
  
  const char* DataToString(uint32_t did);

  HmData * RandProject(int dim, int *projector);

  void WriteToFile(const char* filename);

};

#endif // __HMDATA_H__
