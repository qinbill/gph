#ifndef BITOPTS_H__
#define BITOPTS_H__

#define popcntll __builtin_popcountll
#define popcnt __builtin_popcount

#include <stdio.h>
#include <math.h>
#include "types.h"


const int lookup [] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8};

inline int match(UINT8*P, UINT8*Q, int codelb) {
  switch(codelb) {
    case 4: // 32 bit
      return popcnt(*(UINT32*)P ^ *(UINT32*)Q);
      break;
    case 8: // 64 bit
      return popcntll(((uint64_t*)P)[0] ^ ((uint64_t*)Q)[0]);
      break;
    case 16: // 128 bit
      return popcntll(((uint64_t*)P)[0] ^ ((uint64_t*)Q)[0]) \
          + popcntll(((uint64_t*)P)[1] ^ ((uint64_t*)Q)[1]);
      break;
    case 32: // 256 bit
      return popcntll(((uint64_t*)P)[0] ^ ((uint64_t*)Q)[0]) \
          + popcntll(((uint64_t*)P)[1] ^ ((uint64_t*)Q)[1]) \
          + popcntll(((uint64_t*)P)[2] ^ ((uint64_t*)Q)[2]) \
          + popcntll(((uint64_t*)P)[3] ^ ((uint64_t*)Q)[3]);
      break;
    case 64: // 512 bit
      return popcntll(((uint64_t*)P)[0] ^ ((uint64_t*)Q)[0]) \
          + popcntll(((uint64_t*)P)[1] ^ ((uint64_t*)Q)[1]) \
          + popcntll(((uint64_t*)P)[2] ^ ((uint64_t*)Q)[2]) \
          + popcntll(((uint64_t*)P)[3] ^ ((uint64_t*)Q)[3]) \
          + popcntll(((uint64_t*)P)[4] ^ ((uint64_t*)Q)[4]) \
          + popcntll(((uint64_t*)P)[5] ^ ((uint64_t*)Q)[5]) \
          + popcntll(((uint64_t*)P)[6] ^ ((uint64_t*)Q)[6]) \
          + popcntll(((uint64_t*)P)[7] ^ ((uint64_t*)Q)[7]);
      break;
    default:
      int output = 0;
      for (int i=0; i < codelb; i++) 
        output+= lookup[P[i] ^ Q[i]];
      return output;
      break;
  }

  return -1;
}

/* b <= 64 */
inline void split (uint64_t *chunks, UINT8 *code, int m, int mplus, int b) {
  uint64_t temp = 0x0;
  int nbits = 0;
  int nbyte = 0;
  uint64_t mask = b==64 ? 0xFFFFFFFFFFFFFFFFLLU : ((uint64_t_1 << b) - uint64_t_1);

  for (int i=0; i<m; i++) {
    while (nbits < b) {
      temp |= ((uint64_t)code[nbyte++] << nbits);
      nbits += 8;
    }
    chunks[i] = temp & mask;
    temp = b==64 ? 0x0 : temp >> b; 
    nbits -= b;
    if (i == mplus-1) {
      b--;		/* b <= 63 */
      mask = ((uint64_t_1 << b) - uint64_t_1);
    }
  }
}

/* generates the next binary code (in alphabetical order) with the
   same number of ones as the input x. Taken from
   http://www.geeksforgeeks.org/archives/10375 */
inline uint64_t next_set_of_n_elements(uint64_t x) {
  uint64_t smallest, ripple, new_smallest;

  smallest     = x & -x;
  ripple       = x + smallest;
  new_smallest = x ^ ripple;
  new_smallest = new_smallest / smallest;
  new_smallest >>= 2;
  return ripple | new_smallest;
}

inline void print_code(uint64_t tmp, int b) {
  for (int j=(b-1); j>=0; j--) {
    printf("%llu", (long long int) tmp/(1 << j));
    tmp = tmp - (tmp/(1 << j)) * (1 << j);
  }
  printf("\n");
}

inline uint64_t choose(int n, int r) {
  uint64_t nchooser = 1;
  for (int k=0; k < r; k++) {
    nchooser *= n-k;
    nchooser /= k+1;
  }
  return nchooser;
}

#endif
