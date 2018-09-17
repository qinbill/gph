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
#include "linscan.h"
#include "hmdata.h"
#include "mihasher.h"


using namespace std;

void print_usage(){
  fprintf(stderr, "usage: -d <Location of Data file>\n");
  fprintf(stderr, "       -q <Location of query file>\n");
  fprintf(stderr, "       -p <projected dimension>\n");
  fprintf(stderr, "\n");
  exit(0);
}

int main(int argc, char* argv[])
{

  char c;    
  char *dataFile = NULL;
  char *queryFile = NULL;
  int projdim = 0;
  while ((c = getopt(argc,argv, "hd:q:p:")) != -1)
    switch (c){
    case 'd':
      dataFile = optarg;
      break;
    case 'q':
      queryFile = optarg;
      break;
    case 'p':
      projdim = atoi(optarg);
      break;	
    case 'h':
      print_usage();
      break;
    case '?':
      if ( optopt == 'd' || optopt == 'q' || optopt == 'D' || optopt == 'Q' || optopt == 'p' )
	cerr << "Error: Option -" << optopt << "requires an argument." << endl;
      else if ( isprint(optopt))
	cerr << "Error: Unknown Option -" << optopt << endl;
      else
	cerr << "Error: Unknown Option character" <<endl;
      return 1;
    default:
      print_usage();
    }
  
  if (dataFile == NULL || queryFile == NULL){
    cerr << "Need data input and query input file name" <<endl;
    print_usage();
  }


  HmData *indata = new HmData(dataFile);
  HmData *inquery = new HmData(queryFile);
  HmData *projdata, *projquery;
  
  if (projdim != 0) {
    srand (time(NULL));
    int projector[kMaxDim];
    for (int i = 0; i < indata->mDim; i++) {
      projector[i] = rand() % projdim;
    }
    projdata = indata->RandProject(projdim, projector);
    projquery = inquery->RandProject(projdim, projector);
  }else {
    projdata = indata;
    projquery = inquery;
  }

  char newfilename[1024];
  sprintf(newfilename, "%s.p%d.bin", dataFile, projdim);
  projdata->WriteToFile(newfilename);
  sprintf(newfilename, "%s.p%d.bin", queryFile, projdim);
  projquery->WriteToFile(newfilename);

  if (projdata != indata)
    delete projdata;
  delete indata;
  if (projquery != inquery)
    delete projquery;
  delete inquery;

  return 0;
}
