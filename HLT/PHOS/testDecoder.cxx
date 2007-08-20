
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
  
  //  int n_loops = 200000;
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

  int cnt = 0;
  int channelCnt = 0;

  start =clock();
 
  for(int i=0; i < n_loops; i++)
    {
      decoder->SetMemory((UChar_t*)dataPtr, length);
      decoder->Decode();
 
      
      while(decoder->NextChannel(&altrodata) == true && channelCnt < 10)
	{
	  channelCnt ++;
  decoder->PrintInfo(altrodata, altrodata.fDataSize +4);	   
	 
	  if(  altrodata.fDataSize != 0 )
	    {
	      cnt = 0;
	      
	      while( altrodata.NextBunch(altrobunchPtr) == true)
		{
		  printf("\n");

		  if(cnt < 5)
		    { 
		      printf("\n");
		      cout <<"cnt = "<< cnt <<endl;
		      cout << "altrobunch.fDataSize = "    << altrobunchPtr->fBunchSize   << endl;
		      cout << "altrobunch.fEndTimeBin = "  << altrobunchPtr->fEndTimeBin  << endl;
		     
		      for(int i=0; i<altrobunchPtr->fBunchSize+20; i++)
			{
			  if(i != 0 && i%4==0)
			    {
			      printf("\n");
			    }
			  printf("%d\t", altrobunchPtr->fData[i]);
			}
			
		      printf("\n\n");
		    }
		  cnt ++; 
		}

	    }
	}
      
    }

  end = clock();

  float mikro = (float)(((float)end -(float)start)/((float)n_loops));

  printf("\nProcessing time per event is %f  us\n", mikro);

  decoder->GetFailureRate();

  //  cnt ++;
  
 return 0;
 
}  
