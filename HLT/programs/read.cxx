/* $id$
Author: Constantin Loizides <mailto: loizides@ikf.physik.uni-frankfurt.de
*/

#include <stream.h>
#include "AliL3MemHandler.h"
#include "AliL3AltroMemHandler.h"
#include "AliL3DigitData.h"

/**
Example program how to open and read a raw datafile.
And to store results in an Altro like data format. 
*/

int main(int argc,char **argv)
{
  UInt_t nrows=175;
  Bool_t altroout=kFALSE;
  FILE *afile=0;
  if(argc<2)
    {
      cout<<"Usage: read datafile [padrows]"<<endl;
      return -1;
    }
  if (argc>2) {
    nrows=atoi(argv[2]);
  }
  if (argc>3) {
    altroout=kTRUE;
    afile=fopen(argv[3],"w");
  }

  //Filehandler object:
  AliL3MemHandler file;
  
  //Open the data file:
  if(!file.SetBinaryInput(argv[1]))
    {
      cerr<<"Error opening file "<<argv[1]<<endl;
      return -1;
    }
  
  //Create an RowData object to access the data
  AliL3DigitRowData *digits=0;
  UInt_t ndigits=0;
  
  //Read the file, and store the data in memory. Return value is a pointer to the data.
  digits = file.CompBinary2Memory(ndigits);

  //Create an ALtroMemHandler object
  AliL3AltroMemHandler altromem;
  if(altroout) altroout=altromem.SetBinaryOutput(afile);
  
  UShort_t time,charge;
  UChar_t pad;
  for(UInt_t row=0; row<nrows; row++) //Loop over padrows
    {
      //Get the data on this padrow:
      AliL3DigitData *dataPt = (AliL3DigitData*)digits->fDigitData;
      
      //Loop over all digits on this padrow:
      for(UInt_t ndig=0; ndig<=digits->fNDigit; ndig++)
	{
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
  return 0;
}
