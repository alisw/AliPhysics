/* $Id$
Author: Constantin Loizides <mailto: loizides@ikf.physik.uni-frankfurt.de
*/

#include <stream.h>
#include <libgen.h>
#include "AliL3RootTypes.h"
#include "AliL3Logger.h"
#include "AliL3MemHandler.h"
#include "AliL3AltroMemHandler.h"
#include "AliL3DigitData.h"
#include "AliL3Transform.h"

/**
Example program how to open and read a raw datafile.
In addition it shows, how to store results in an Altro like 
data format. 
*/

int main(int argc,char **argv)
{
  Int_t slice=0;
  Int_t patch=0;
  Bool_t altroout=kFALSE;
  FILE *afile=0;
  
  /*
    AliL3Logger l;
    l.Set(AliL3Logger::kAll);
    //l.UseStdout();
    l.UseStream();
  */

  if(argc<2)
    {
      cout<<"Usage: read datafile [slice] [patch] [altrodatfile]"<<endl;
      return -1;
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
  AliL3Transform::Init(dirname(fname));
  strcpy(fname,argv[1]);

  //Filehandler object:
  AliL3MemHandler file;

  //Give slice and patch information (see filename convention)
  if((patch>=0)&&(patch<6)) file.Init(slice,patch);
  else {
    Int_t srows[2]={0,175};
    file.Init(slice,0,srows);
  }

  //Open the data file:
  if(!file.SetBinaryInput(argv[1]))
    {
      cerr<<"Error opening file "<<argv[1]<<endl;
      return -1;
    }

  //Create an RowData object to access the data
  AliL3DigitRowData *digits=0;
  UInt_t nrows=0;
  
  //Read the file, and store the data in memory. Return value is a pointer to the data.
  digits = file.CompBinary2Memory(nrows);

  //Create an ALtroMemHandler object
  AliL3AltroMemHandler altromem;
  if(altroout) altroout=altromem.SetASCIIOutput(afile);

  UShort_t time,charge;
  UChar_t pad;
  Int_t row=file.GetRowMin()-1,crows=0,lrow=row;

  for(UInt_t r=0; r<nrows; r++) //Loop over padrows
    {
      //Get the data on this padrow:
      AliL3DigitData *dataPt = (AliL3DigitData*)digits->fDigitData;
      row++;
      if(lrow+1==row) crows++;
      
      //Loop over all digits on this padrow:
      for(UInt_t ndig=0; ndig<digits->fNDigit; ndig++)
      //for(UInt_t ndig=digits->fNDigit;ndig>0;ndig--)
	{
	  lrow=row;
	  pad = dataPt[ndig].fPad;
	  time = dataPt[ndig].fTime;
	  charge = dataPt[ndig].fCharge;
	  cout << "Padrow " << row << " pad " << (int)pad << " time " <<(int) time << " charge " << (int)charge << endl;
	  if(altroout) altromem.Write(row,pad,charge,time);
	}
      
      //Move the pointer to the next padrow:
      file.UpdateRowPointer(digits);
    }
  
  if(afile) {
    altromem.WriteFinal();
    fclose(afile);
  }

  //cerr << "Rows: " << (file.GetRowMax()-file.GetRowMin()+1) << " Consecutive: " << crows << endl;
  return 0;
}


