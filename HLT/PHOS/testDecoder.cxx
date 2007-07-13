#include "AliHLTDDLDecoder.h"
#include "AliHLTAltroData.h"
#include "Rtypes.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <time.h>
#include "stdio.h"

using namespace std;

///#define N_LOOPS 1000

int main(int argc, const char** argv)
{
  int n_loops = 1000000;

  clock_t  start;
  clock_t  end;

  AliHLTAltroData altrodata;
  AliHLTAltroBunch altrobunch;

  ifstream fin;
  int length;
  
  AliHLTDDLDecoder *decoder = new AliHLTDDLDecoder();
  
  fin.open(argv[1], ios::binary);
  
  fin.seekg (0, ios::end);
  length = fin.tellg();
  fin.seekg (0, ios::beg);

  char *dataPtr = new char[length];

  fin.read (dataPtr,length);
  fin.close();

  start =clock();
 
  for(int i=0; i < n_loops; i++)
    {
      decoder->SetMemory((UChar_t*)dataPtr, length);
      decoder->Decode();
      while(decoder->NextChannel(&altrodata) == true)
	{
	  //	  altrodata.Reset();
	  //      decoder->PrintInfo(altrodata, altrodata.fDataSize, 4);
 
	  //     cout << endl;


	  	  
	    while( altrodata.NextBunch(altrobunch) == true)
	    {
	    //	  cout << "altrobunch.fDataSize = "    << altrobunch.fBunchSize  << endl;
	    //	  cout << "altrobunch. fEndTimeBin = " << altrobunch. fEndTimeBin  << endl;
	    }
	  
	}

      //     end = clock();

    }

  end = clock();

  float mikro = (float)(((float)end -(float)start)/((float)n_loops));

  

  printf("\nProcessing time per event is %f  us\n", mikro);

  decoder->GetFailureRate();

  //  cnt ++;
  return 0;
}  
