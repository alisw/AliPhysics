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



int main(int argc, const char** argv)
{
  int n_loops = 1;

  clock_t  start;
  clock_t  end;

  AliHLTAltroData altrodata;
  AliHLTAltroBunch *altrobunchPtr = new AliHLTAltroBunch;

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
	  //	  printf("\n\n");
	  altrodata.Reset();
	 
	  if(  altrodata.fDataSize != 0 )
	    {
	      
	      //          printf("\n\n");
	      Double_t tmpMax = 0;
	      
	      for(int i = 0; i < altrodata.fDataSize -3 ; i++)
		{
		  if(altrodata.fData[i] > tmpMax)
		    {
		      tmpMax = altrodata.fData[i];
		    }
		}
	   
	      
	      //      cout <<"tmpMax = "<< tmpMax<<endl;

	      if (tmpMax > 100)
		{ 
		  cout <<"tmpMax = "<< tmpMax<<endl;
		  decoder->PrintInfo(altrodata, altrodata.fDataSize, 4);
		}
	  	  
	      while( altrodata.NextBunch(altrobunchPtr) == true)
		{
		  //		  cout << "altrobunch.fDataSize = "    << altrobunchPtr->fBunchSize   << endl;
		  //	  cout << "altrobunch.fEndTimeBin = "  << altrobunchPtr->fEndTimeBin  << endl;
		}
	      //        printf("\n\n");
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
