// $Id$

// Author: Constantin Loizides <loizides@ikf.uni-frankfurt.de
//*-- Copyright &copy ALICE HLT Group


#include "AliHLTStandardIncludes.h"

#include "AliHLTRootTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTLogger.h"
#include "AliHLTMemHandler.h"
#include "AliHLTAltroMemHandler.h"
#include "AliHLTDigitData.h"
#include "AliHLTTransform.h"

#if __GNUC__ == 3
using namespace std;
#else
#include <stream.h>
#endif

#include <libgen.h>


/**
   Example program how to open and read a raw datafile.
   In addition it shows, how to store (and read) the
   digits in an Altro like data format. 
*/

int main(Int_t argc,Char_t **argv)
{
  Int_t slice=0;
  Int_t patch=0;
  Bool_t altroout=kFALSE;
  FILE *afile=0;
  
  AliHLTLogger l;
  l.Set(AliHLTLogger::kAll);
  //l.UseStderr();
  //l.UseStdout();
  //l.UseStream();

  if(argc<2)
    {
      cout<<"Usage: read datafile [slice] [patch] [altrodatfile]"<<endl;
      exit(1);
    }
  if (argc>2) {
    slice=atoi(argv[2]);
  }
  if (argc>3) {
    patch=atoi(argv[3]);
  }
  if (argc>4) {
    altroout=kTRUE;
    afile=fopen(argv[4],"w");
  }
  
  //Loading all specific aliroot version quantities, needed.
  Char_t fname[1024];
  strcpy(fname,argv[1]);
  AliHLTTransform::Init(dirname(fname));
  strcpy(fname,argv[1]);

  //Filehandler object:
  AliHLTMemHandler file;

  //Give slice and patch information (see filename convention)
  if((patch>=0)&&(patch<6)) file.Init(slice,patch);
  else {
    Int_t srows[2]={0,175};
    patch=0;
    file.Init(slice,patch,srows);
  }

  //Open the data file:
  if(!file.SetBinaryInput(argv[1]))
    {
      cerr<<"Error opening file "<<argv[1]<<endl;
      return -1;
    }

  //Create an RowData object to access the data
  AliHLTDigitRowData *digits=0;
  UInt_t nrows=0;
  
  //Read the file, and store the data in memory. Return value is a pointer to the data.
  digits = file.CompBinary2Memory(nrows);

  //Create an ALtroMemHandler object
  AliHLTAltroMemHandler altromem;
  if(altroout) altroout=altromem.SetASCIIOutput(afile);

  UShort_t time,charge;
  UChar_t pad;
  Int_t row=file.GetRowMin()-1,crows=0,lrow=row;

  for(UInt_t r=0; r<nrows; r++) //Loop over padrows
    {
      //Get the data on this padrow:
      AliHLTDigitData *dataPt = (AliHLTDigitData*)digits->fDigitData;
      row++;
      if(lrow+1==row) crows++;
      
      //Loop over all digits on this padrow:
      for(UInt_t ndig=0; ndig<digits->fNDigit; ndig++)
	{
	  lrow=row;
	  pad = dataPt[ndig].fPad;
	  time = dataPt[ndig].fTime;
	  charge = dataPt[ndig].fCharge;
	  cout << "Padrow " << r << " pad " << (int)pad << " time " <<(int) time << " charge " << (int)charge << endl;
	  //cout << "Padrow " << row << " pad " << (int)pad << " time " <<(int) time << " charge " << (int)charge << endl;
	  if(altroout) altromem.Write(r,pad,time,charge);
	}
      
      //Move the pointer to the next padrow:
      file.UpdateRowPointer(digits);
    }
  
  if(afile) {
    altromem.WriteFinal();
    fclose(afile);
    
#if 0
    //test Altro read
    UShort_t rrow=0,rtime=0,rcharge=0;
    UChar_t rpad=0;
    afile=fopen(argv[4],"r");
    altromem.SetASCIIInput(afile);
    while(altromem.Read(rrow,rpad,rtime,rcharge)){
      cout << "Padrow " << (int)rrow << " pad " << (int)rpad << " time " <<(int)rtime << " charge " << (int)rcharge << endl;
    }
    fclose(afile);  
#endif
#if 0
    //test Altro read sequence
    UShort_t rrow=0,rtime=0;
    UChar_t rpad=0,n=100,i=100;
    UShort_t *charges=new UShort_t[100];
    afile=fopen(argv[4],"r");
    altromem.SetASCIIInput(afile);
    while(altromem.ReadSequence(rrow,rpad,rtime,i,&charges)){
      cout << "Padrow " << (int)rrow << " pad " << (int)rpad << " time " <<(int)rtime << " charges ";
      for(UChar_t ii=0;ii<i;ii++) cout << (int)charges[ii] << " ";
      cout << endl;
      i=n;
    }
    fclose(afile);  
#endif
  }

  //cerr << "Rows: " << (file.GetRowMax()-file.GetRowMin()+1) << " Consecutive: " << crows << endl;
  return 0;
}


