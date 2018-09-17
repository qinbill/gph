#include <getopt.h>
#include <iostream>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#include <vector>
#include "sys/types.h"

#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <pthread.h>
#include <time.h>
#include <string.h>

#include "types.h"
#include "hmdata.h"
#include "bitops.h"

using namespace std;

void print_usage(){
  fprintf(stderr, "usage: -d <Location of Data file>\n");
  fprintf(stderr, "       `-q <Location of Data file>\n");
  fprintf(stderr, "\n");
  exit(0);
}

int main(int argc, char* argv[])
{
  char c;    
  char *dataFile = NULL;
  char *queryFile = NULL;
  int64_t N = -1;
  int64_t NQ = -1;
  //double rateofsample = 0.1;
  while ((c = getopt(argc,argv, "hd:q:")) != -1)
  {
    switch (c)
    {
      case 'd':
        dataFile = optarg;
        break;
      case 'q':
        queryFile = optarg;
        break;
      case 'h':
        print_usage();
        break;
      case '?':
        if ( optopt == 't' || optopt == 'd' || optopt == 'D' || optopt == 'l' || optopt == 'q' || optopt == 't' )
          cerr << "Error: Option -" << optopt << "requires an argument." << endl;
        else if ( isprint(optopt))
          cerr << "Error: Unknown Option -" << optopt << endl;
        else
          cerr << "Error: Unknown Option character" <<endl;
        return 1;
      default:
        print_usage();
    }
  }
  
  if (dataFile == NULL){
    cerr << "Need data input and query input file name" <<endl;
    print_usage();
  }

  
  HmData *indata = new HmData(dataFile);
  HmData *inquery = new HmData(queryFile);
  
  if (N < 0 || N > indata->mNumData) {
    N = indata->mNumData;
  } 

  if (NQ < 0 || NQ > inquery->mNumData) {
    NQ = inquery->mNumData;
  } 

  uint8_t* datacodes = indata->mData;
  uint8_t* querycodes = inquery->mData;

  int numbytes=inquery->mDataBytes;
  
  int **dimposstats = new int*[indata->mDim];
  int **dimnegstats = new int*[indata->mDim];
  
  double **simposscore = new double*[indata->mDim];
  double **simnegscore = new double*[indata->mDim];
  
  double *dimposrate = new double[indata->mDim];
  double *dimnegrate = new double[indata->mDim];

  for (int i = 0; i < indata->mDim; i++) {
    dimposstats[i] = new int[indata->mNumData];
    dimnegstats[i] = new int[indata->mNumData];    
    simposscore[i] = new double[indata->mDim];
    simnegscore[i] = new double[indata->mDim];
  }

  for (int i=0; i< indata->mNumData; i++) {
    for (int b=0; b< numbytes; b++) {
      uint8_t abyte = indata->mData[i * indata->mDataBytes + b];
      for (int t = 0; t < 8; t++) {
	int dim = b*8+t;
	if ((abyte >> (7-t))%2==1) {
	  dimposstats[dim][i] = i;
	  dimnegstats[dim][i] = -i;
	} else {
	  dimposstats[dim][i] = -i;
	  dimnegstats[dim][i] = i;
	}
      }
    }
  }

  for (int i = 0; i < indata->mDim; i++) {
    dimposrate[i] = 0;
    dimnegrate[i] = 0;
    
    for (int k = 0; k < indata->mNumData; k++) {
      if (dimposstats[i][k] > 0) {
	dimposrate[i] +=1;
      } else {
	dimnegrate[i] +=1;
      }
    }

    for (int j = 0; j < indata->mDim; j++) {
      int posdiff = 0;
      int poscomm = 0;
      int negdiff = 0;
      int negcomm = 0;
      
      for (int k = 0; k < indata->mNumData; k++) {
	if (dimposstats[i][k]!=dimposstats[j][k]) {
	  posdiff ++;
	} else {
	  poscomm ++;
	}
	if (dimposstats[i][k]!=dimnegstats[j][k]) {
	  negdiff ++;
	} else {
	  negcomm ++;
	}
      }
      simposscore[i][j] = (poscomm * 1.0)/(indata->mNumData + posdiff);
      simnegscore[i][j] = (negcomm * 1.0)/(indata->mNumData + negdiff);
    }
  }
  
  
//   for (int i = 0; i < indata->mDims; i++) { 
//   }

  for (int i = 0; i < indata->mDim; i++) {
    int mostpossim = -1;
    // int secdpossim = -1;
    int mostnegsim = -1;
    // int secdnegsim = -1;

    printf("dim: %d posrate: %.6f negrate: %.6f minrate: %.6f ", i,
	   dimposrate[i]/indata->mNumData, dimnegrate[i]/indata->mNumData, 
	   (dimposrate[i]<dimnegrate[i]?dimposrate[i]:dimnegrate[i])/indata->mNumData);
    
    for (int j = 0; j < indata->mDim; j++) { 
      if (i==j) continue;
      if (mostpossim==-1 || simposscore[i][mostpossim] < simposscore[i][j]) {
	mostpossim = j;
      }
      if (mostnegsim==-1 || simnegscore[i][mostnegsim] < simnegscore[i][j]) {
	mostnegsim = j;
      }
    }
    printf(" mostsimpos: %d %.6f mostnegpos: %d %.6f\n", 
	   mostpossim, simposscore[i][mostpossim], 
	   mostnegsim, simnegscore[i][mostnegsim]);
  }

  delete indata;
  delete inquery;
  return 0;
}
