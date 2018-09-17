#include <stdio.h>
#include <string.h>
#include <string>       // std::string
#include <iostream>     // std::cout
#include <sstream>      // std::stringstream
#include <time.h>
#include <math.h>


#include "hmdata.h"

using namespace std;


HmData::HmData(const char *filename) {
  FILE *fp = fopen(filename, "rb");
  uint64_t N = 0;
  int B, D;

  // Read the dimentionality informations.
  fread(&N, sizeof(uint64_t), 1, fp);
  fread(&D, sizeof(int), 1, fp);
  fread(&B, sizeof(int), 1, fp);

  mDim = D;
  mDataBytes = B;
  mNumData = N;

  mData = new uint8_t[mNumData*mDataBytes];
  fread(mData, sizeof(uint8_t), mNumData*mDataBytes, fp);
  fclose(fp);
}


void HmData::WriteToFile(const char* filename) {
  FILE *fp = fopen(filename, "wb");
  
  // Write meta information.
  fwrite(&mNumData, sizeof(uint64_t), 1, fp);
  fwrite(&mDim, sizeof(uint32_t), 1, fp);
  fwrite(&mDataBytes, sizeof(uint32_t), 1, fp);

  fwrite(mData, 1, mNumData * mDataBytes, fp);
  fclose(fp);
}


HmData* HmData::RandProject(int dim, int *projector) {
  int projected[kMaxDim];
  
  HmData *hmdataProj = new HmData();
  hmdataProj->mDim = dim;
  hmdataProj->mNumData = mNumData;
  hmdataProj->mDataBytes = dim / 8;
  hmdataProj->mData = 
    new uint8_t[hmdataProj->mNumData * hmdataProj->mDataBytes];
  for (int i = 0; i < mNumData; i++) {
    memset(projected, 0, sizeof(int) * kMaxDim);
    for (int b = 0; b < mDataBytes; b++) {
      uint8_t orig = mData[i * mDataBytes + b];
      for (int c = 0; c < 8; c ++) {
	int cval = 0;
	if (((orig >> c) % 2) == 1){
	  cval = 1;
	}
	projected[projector[b * 8 + c]] += cval;
      }
    }
    for (int b = 0; b < hmdataProj->mDataBytes; b++) {
      uint8_t val = 0;
      for (int c = 0; c < 8; c++) {
	val = (val << 1) + (projected[b*8+c]==0?0:1);
      }
      hmdataProj->mData[i * hmdataProj->mDataBytes + b] = val;
    }
  }
  return hmdataProj;
}






HmData::~HmData() {
  delete [] mData;
  mData = NULL;
}
