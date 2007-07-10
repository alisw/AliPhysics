#include "AliHLTDDLDecoder.h"
#include "AliHLTAltroData.h"
#include "Rtypes.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>

using namespace std;

int main(int argc, const char** argv)
{
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
 
  decoder->SetMemory((UChar_t*)dataPtr, length);
  decoder->Decode();

  int cnt = 0; 

  while(decoder->NextChannel(&altrodata) == true)
    {

      altrodata.Reset();

   
      //      decoder->PrintInfo(altrodata, altrodata.fDataSize, 4);
 
      //     cout << endl;

      while( altrodata.NextBunch(altrobunch) == true)
	{
	  //	  cout << "altrobunch.fDataSize = "    << altrobunch.fBunchSize  << endl;
	  //	  cout << "altrobunch. fEndTimeBin = " << altrobunch. fEndTimeBin  << endl;
	}
	   
    }

  decoder->GetFailureRate();

  cnt ++;
  return 0;
}  
