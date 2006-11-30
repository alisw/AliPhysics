// @(#) $Id$

// Author: Anders Vestbo <mailto:vestbo@fi.uib.no>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTStandardIncludes.h"

#include "AliHLTRootTypes.h"
#include "AliHLTLogging.h"
#include "AliHLTMemHandler.h"
#include "AliHLTLogging.h"
#include "AliHLTTransBit.h"
#include "AliHLTTransform.h"
#include "AliHLTDataHandler.h"
#include "AliHLTDigitData.h"

#if __GNUC__ >= 3
using namespace std;
#endif

/** \class AliHLTDataHandler
<pre>
//_____________________________________________________________
// AliHLTDataHandler
//
// HLT Binary file handler.
//
// This class have more or less the same functionality as AliHLTMemHandler,
// except that it handles 8 bit ADC-values. Reading and writing is done in the same way
// as illustrated in example 1) and 2) in AliHLTMemHandler.
//
// For converting 10 bit data files to 8 bit data files, do:
//
// AliHLTMemHandler *file = new AliHLTDataHandler();
// file->Init(slice,patch);
// file->SetBinaryInput(inputfile);    //10 bit data file
// file->SetBinaryOutput(outputfile);  //8 bit data file
// file->Convert10to8Bit();
// file->CloseBinaryInput();
// file->CloseBinaryOutput();
// delete file;
//
// Compress data format
// --------------------
//
// The data is RLE encoded, using _8_bit representation of the ADC-values.
// Conversion is done in the class AliHLTTransBit.
//
// In the beginning of every row, the row number if written and the number of pads
// containing data on that row. For every pad with data the pad number is written,
// and then comes the ADC-values on that pad. When a serie of zeros occure, a zero
// is written followed by the number of zeros. If the number of zeros is more than
// 255 (8 bit), another 8 bit word is written for the remaining. At the end of one 
// pad, 2 zeros are written. Example:
//
// ROW NPADSWITHDATA PAD 0 NZEROS ADC ADC ADC ADC 0 NZEROS ADC ADC 0 0
//
// Everything is written using 8 bit;
// (ROW < 176, PAD < 200, ADC < 255, if(NZEROS > 255) write 2 words;)
</pre>
*/

ClassImp(AliHLTDataHandler)
  
AliHLTDataHandler::AliHLTDataHandler()
{
  // default constructor
  fBitTransformer = 0;
  LOG(AliHLTLog::kInformational,"AliHLTDataHandler::AliHLTDataHandler","Data format")
    <<"8 bit data handler initialized"<<ENDLOG;
}

AliHLTDataHandler::~AliHLTDataHandler()
{
  // destructor
  if(fBitTransformer)
    delete fBitTransformer;
}

void AliHLTDataHandler::Convert10to8Bit()
{
  //Convert from 10 bit data in inputfile, to 8 bit data written to outputfile.
  
  if(!fInBinary)
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::Convert10to8Bit","File")
	<<AliHLTLog::kHex<<"Pointer to input file : "<<(void*)fInBinary<<ENDLOG;
      return;
    }
  if(!fOutBinary)
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::Convert10to8Bit","File")
	<<AliHLTLog::kHex<<"Pointer to output file : "<<(void*)fOutBinary<<ENDLOG;
      return;
    }
  
  
  //Initialize the bit transformation class:
  fBitTransformer = new AliHLTTransBitV1();
  Int_t b0=10;  // original number of bits
  Int_t b1=8;   // compressed
  fBitTransformer->SetBits(b0,b1);
  fBitTransformer->FindOptimumX0();
  fBitTransformer->Update();
  
  AliHLTMemHandler *memory = new AliHLTMemHandler();
  memory->Init(fSlice,fPatch);
  memory->SetBinaryInput(fInBinary);
  UInt_t nrow;
  AliHLTDigitRowData *data = (AliHLTDigitRowData*)memory->CompBinary2Memory(nrow);
  
  Memory2CompBinary(nrow,data);
  
  delete memory;
}

Bool_t AliHLTDataHandler::Memory2CompBinary(UInt_t nrow,AliHLTDigitRowData *data)
{
  //Compress data by RLE, and write to a binary file.
  
  UInt_t size = GetCompMemorySize(nrow,data);
  Byte_t *comp = Allocate(size);
  Memory2CompMemory(nrow,data,comp);
  if(!CompMemory2CompBinary(nrow,comp,size))
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::Memory2CompBinary","File")
	<<"Error writing to file "<<ENDLOG;
      return 0;
    }
  Free();
  return kTRUE;
}

AliHLTDigitRowData *AliHLTDataHandler::CompBinary2Memory(UInt_t &nrow)
{
  //Read RLE compressed binary file, unpack it and return pointer to it.
  
  AliHLTMemHandler *memory = new AliHLTMemHandler();
  memory->SetBinaryInput(fInBinary);
  Byte_t *comp = memory->Allocate();
    
  if(!CompBinary2CompMemory(nrow,comp))
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::CompBinary2Memory","File")
	<<"Error reading from file "<<ENDLOG;
      return 0;
    }

  UInt_t size = GetMemorySize(nrow,comp);
  AliHLTDigitRowData *data = (AliHLTDigitRowData*)Allocate(size);
  CompMemory2Memory(nrow,data,comp);
  delete memory;
  return data;
}

void AliHLTDataHandler::Write(Byte_t *comp,UInt_t &index,UShort_t value)
{
  //Write one value (=1 byte) to array comp.

  if(value > 255)
    {
      LOG(AliHLTLog::kFatal,"AliHLTDataHandler::Write","Bitnumbers")
	<<"Value too big for storing in 1 byte, something is wrong: "<<value<<" "<<index<<ENDLOG;
    }
  comp[index] = (Byte_t)value;
  index++;
}

Short_t AliHLTDataHandler::Read(Byte_t *comp,UInt_t &index)
{
  //Read one value (=1 byte) from array comp

  Short_t value = (Short_t)comp[index];
  index++;
  return value;
}

Short_t AliHLTDataHandler::Test(Byte_t *comp,UInt_t index)
{
  //Check the value (=1 byte) in array comp, but not read.

  Short_t value = (Short_t)comp[index];
  return value;
}

Bool_t AliHLTDataHandler::Memory2CompMemory(UInt_t nrow,AliHLTDigitRowData *data,Byte_t *comp)
{
  //Perform RLE.
  
  if(!data)
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::Memory2CompMemory","Data")
	<<AliHLTLog::kHex<<" Pointer to data = "<<(void*)data<<ENDLOG;
      return 0;  
    }
  if(!comp)
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::Memory2CompMemory","Data")
	<<AliHLTLog::kHex<<" Pointer to compressed data = "<<(void*)comp<<ENDLOG;
      return 0;  
    }

  AliHLTDigitRowData *rowPt = data;
  
  UInt_t index = 0;
  Int_t npads[200];      
  
  for(UInt_t i=0; i<nrow; i++)
    {
      //Write the row number:
      UShort_t value = rowPt->fRow;
      Write(comp,index,value);
      
      UShort_t numberOfPads=0;
      UShort_t maxPad = 0;
      
      for(Int_t j=0; j<200; j++)
	npads[j]=0;
      for(UInt_t dig=0; dig<rowPt->fNDigit; dig++)
	{
	  if(rowPt->fDigitData[dig].fPad < 200)
	    npads[rowPt->fDigitData[dig].fPad]++;
	}
      for(Int_t j=0; j<200; j++)
	{
	  if(npads[j])
	    {
	      numberOfPads++;
	      maxPad = j;
	    }
	}
      
      //Write the number of pads on this row:
      Write(comp,index,numberOfPads);
      UInt_t digit=0;
      
      for(UShort_t pad=0; pad <= maxPad; pad++)
	{
	  
	  if(digit >= rowPt->fNDigit || rowPt->fDigitData[digit].fPad !=  pad)
	    continue;
	
	  //Write the current pad:
	  Write(comp,index,pad);
	  
	  if(digit < rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad)
	    {
	      if(rowPt->fDigitData[digit].fTime > 0)
		{
		  //If first time!=0, write the number of following zeros, 
		  //and then the first timebin:
		  Write(comp,index,0);
		  
		  //Check if we have to use more than 1 byte to write the zeros:
		  Int_t numberOfZeroIntervals=0;
		  if(rowPt->fDigitData[digit].fTime >= 255)
		    {
		      numberOfZeroIntervals++;
		      Write(comp,index,255);
		      if(rowPt->fDigitData[digit].fTime >= 2*255)
			{
			  cerr<<"AliHLTDataHandler::Memory2CompMemory : Should not happen "<<(Int_t)rowPt->fDigitData[digit].fTime<<endl;
			  Write(comp,index,255);
			  numberOfZeroIntervals++;
			}
		    }
		  Write(comp,index,(rowPt->fDigitData[digit].fTime - numberOfZeroIntervals*255));
		}
	    }
	  
	  while(digit < rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad)
	    {
	      UShort_t charge = rowPt->fDigitData[digit].fCharge;
	      
	      if(fBitTransformer)
		charge = fBitTransformer->Get0to1(charge); //Transform 10 to 8 bit.
	      
	      //Check for saturation:
	      if(charge>255)
		{
		  LOG(AliHLTLog::kWarning,"AliHLTDataHandler::Memory2CompMemory","Digit")
		    <<"ADC-value saturated : "<<charge<<ENDLOG;
		  charge=255;
		}
	      
	      //Write the charge:
	      Write(comp,index,charge);
	      
	      //Check if the next digit is zero:
	      if(digit+1 < rowPt->fNDigit && rowPt->fDigitData[digit+1].fPad == pad)
		{
		  if(rowPt->fDigitData[digit].fTime + 1 != rowPt->fDigitData[digit+1].fTime)
		    {
		      Write(comp,index,0);
		      UShort_t nzero = rowPt->fDigitData[digit+1].fTime - (rowPt->fDigitData[digit].fTime + 1);
		      
		      //Check if we have to use more than one byte to write the zeros:
		      Int_t numberOfZeroIntervals=0;
		      if(nzero >= 255)
			{
			  numberOfZeroIntervals++;
			  Write(comp,index,255);
			  if(nzero >= 2*255)
			    {
			      cerr<<"AliHLTDataHandler::Memory2CompMemory : Should not happen "<<(Int_t)rowPt->fDigitData[digit].fTime<<endl;
			      Write(comp,index,255);
			      numberOfZeroIntervals++;
			    }
			}
		      Write(comp,index,(nzero - numberOfZeroIntervals*255));
		    }  
		}
	      digit++;
	    }
	  
	  //This is the end of the pad, state it with 2 zeros:
	  Write(comp,index,0);
	  Write(comp,index,0);
	}
      
      UpdateRowPointer(rowPt);
      
    }
  
  return index * sizeof(Byte_t);
    
}

UInt_t AliHLTDataHandler::GetCompMemorySize(UInt_t nrow,AliHLTDigitRowData *data)
{
  //Calculate the size (in bytes) of RLE data.
  
  if(!data)
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::GetCompMemorySize","Data")
	<<AliHLTLog::kHex<<" Data pointer = "<<(void*)data<<ENDLOG;
      return 0;
    }
  
  AliHLTDigitRowData *rowPt = data;
  
  UInt_t index = 0;
  Int_t npads[200];
  
  for(UInt_t i=0;i<nrow;i++)
    {
      //Write the row number:
      index++;
      
      UShort_t maxPad=0; 
      UShort_t numberOfPads = 0;
      
      for(Int_t j=0; j<200; j++) 
	npads[j]=0;
      
      for(UInt_t dig=0; dig<rowPt->fNDigit; dig++)
	{
	  if(rowPt->fDigitData[dig].fPad <200)
	    npads[rowPt->fDigitData[dig].fPad]++;
	}
      for(Int_t j=0; j<200; j++)
	{ 
	  if(npads[j])
	    {
	      numberOfPads++;
	      maxPad = j;
	    }
	}
      
      //Write the number of pads on this row:
      index++;
      
      UInt_t digit=0;
      for(UShort_t pad=0; pad <= maxPad; pad++)
	{
	  if(digit>=rowPt->fNDigit || rowPt->fDigitData[digit].fPad !=  pad)
	    continue;
	  
	  //Write the current pad:
	  index++;
	  
	  
	  if(digit<rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad)
	    {
	      if(rowPt->fDigitData[digit].fTime > 0)
		{
		  //If first time!=0, write the number of following zeros, 
		  //and then the first timebin:
		  
		  index++;
		  index++;
		  
		  //Check if we have to use more than 1 byte to write the zeros:
		  if(rowPt->fDigitData[digit].fTime >= 255)
		    index++;
		  if(rowPt->fDigitData[digit].fTime >= 2*255)
		    index++;
		}
	    }
	  
	  while(digit < rowPt->fNDigit && rowPt->fDigitData[digit].fPad == pad)
	    {
	      //Write the charge:
	      index++;
	      
	      //Check if the next digit is zero:
	      if(digit+1 < rowPt->fNDigit && rowPt->fDigitData[digit+1].fPad == pad)
		{
		  if(rowPt->fDigitData[digit].fTime +1 != rowPt->fDigitData[digit+1].fTime)
		    {
		      index++;
		      index++;
		      
		      //Check if we have to use more than 1 byte to write the zeros:
		      UInt_t nzeros = rowPt->fDigitData[digit+1].fTime - rowPt->fDigitData[digit].fTime + 1;
		      if(nzeros >= 255)
			index++;
		      if(nzeros >= 2*255)
			index++;
		    }  
		}
	      digit++;
	    }
	  
	  //Mark the end of the pad with 2 zeros:
	  index++;
	  index++;
	}
      
      UpdateRowPointer(rowPt);
    }
  
  return index * sizeof(Byte_t);
  
}

UInt_t AliHLTDataHandler::CompMemory2Memory(UInt_t nrow,AliHLTDigitRowData *data,Byte_t *comp)
{
  //Uncompress RLE data.
  
  if(!data)
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::CompMemory2Memory","Array")
	<<AliHLTLog::kHex<<"Pointer to data: "<<(void*)data<<ENDLOG;
      return 0;
    }
  if(!comp)
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::CompMemory2Memory","Array")
	<<AliHLTLog::kHex<<"Pointer to compressed data: "<<(void*)comp<<ENDLOG;
      return 0;
    }
  
  Int_t outsize=0;
  
  AliHLTDigitRowData *rowPt = data;
  UInt_t index=0;

  UShort_t pad,time,charge;
  for(UInt_t i=0; i<nrow; i++)
    {
      UInt_t ndigit=0;
      
      //Read the row:
      rowPt->fRow = Read(comp,index);

      //Read the number of pads:
      UShort_t npads = Read(comp,index);
      
      for(UShort_t p=0; p<npads; p++)
	{
	  //Read the current pad:
	  pad = Read(comp,index);
	  
	  time = 0;
	  
	  //Check for zeros:
	  if(Test(comp,index) == 0) //Zeros
	    {
	      //Read the first zero
	      Read(comp,index);
	      
		
	      if(Test(comp,index) == 0)//end of pad.
		{
		  time = Read(comp,index); 
		  continue;
		}
	      if( (time = Read(comp,index)) == 255 )
		if( (time += Read(comp,index)) == 2*255)
		  time += Read(comp,index);
	    }

	  while(1)
	    {
	      while( (charge = Read(comp,index)) != 0)
		{
		  if(time >= AliHLTTransform::GetNTimeBins())
		    cerr<<"AliHLTDataHandler::CompMemory2Memory : Time out of range "<<time<<endl;
		  rowPt->fDigitData[ndigit].fPad = pad;
		  rowPt->fDigitData[ndigit].fTime = time;
		  rowPt->fDigitData[ndigit].fCharge = charge;
		  ndigit++;
		  if(Test(comp,index) != 0)
		    time++;
		}
	      if(Test(comp,index) == 0)
		{
		  Read(comp,index); //end of pad
		  break;
		}
	      UShort_t timeShift;
	      if( (timeShift = Read(comp,index)) == 255)
		if( (timeShift += Read(comp,index)) == 2*255)
		  timeShift += Read(comp,index);
	      time += timeShift;
	      
	    }
	}
      rowPt->fNDigit = ndigit;
      UpdateRowPointer(rowPt);
      outsize += sizeof(AliHLTDigitData)*ndigit + sizeof(AliHLTDigitRowData);
    }
  
  return outsize;
}

UInt_t AliHLTDataHandler::GetMemorySize(UInt_t nrow,Byte_t *comp)
{
  //Calculate size (in bytes) of unpacked data.

  UInt_t index=0;
  Int_t outsize=0;
  
  for(UInt_t i=0; i<nrow; i++)
    {
      UInt_t ndigit=0;//Digits on this row.

      //Row number:
      Read(comp,index);
      
      UShort_t npad = Read(comp,index);
      
      for(UShort_t pad=0; pad<npad; pad++)
	{
	  //Read the pad number:
	  Read(comp,index);
	  
	  //Check for zeros:
	  if(Test(comp,index)==0) //Zeros are coming
	    {
	      Read(comp,index);
	      if(Test(comp,index) == 0) 
		{
		  Read(comp,index); //This was the end of pad.
		  continue; 
		}
	      if(Read(comp,index) == 255) //There can be up to 3 bytes with zero coding.
		if(Read(comp,index) == 255)
		  Read(comp,index);
	    }
	  
	  while(1)
	    {
	      while(Read(comp,index) != 0) ndigit++;
	      
	      if(Test(comp,index) == 0)
		{
		  Read(comp,index); //2 zeros = end of pad.
		  break;
		}
	      if(Read(comp,index) == 255) //There can be up to 3 bytes with zero coding.
		if(Read(comp,index) == 255)
		  Read(comp,index);
	      
	    }
	  
	}
      Int_t size = sizeof(AliHLTDigitData)*ndigit + sizeof(AliHLTDigitRowData);
      outsize += size;
    }
  return outsize;
}

Bool_t AliHLTDataHandler::CompBinary2CompMemory(UInt_t &nrow,Byte_t *comp)
{
  //Read RLE data from binary file into array comp.
  rewind(fInBinary);
  UInt_t size = GetFileSize() - 2;
  Byte_t type;
  if(fread(&type,1,1,fInBinary)!=1) return kFALSE;
  if(type > 0)
    {
      LOG(AliHLTLog::kError,"AliHLTDataHandler::CompBinary2CompMemory","Filetype")
	<<"Inputfile does not seem to contain 8 bit data : "<<type<<ENDLOG;
      return kFALSE;
    }
  if(fread(&nrow,1,1,fInBinary)!=1) return kFALSE;
  if(fread(comp,size,1,fInBinary)!=1) return kFALSE;

  return kTRUE;
}

Bool_t AliHLTDataHandler::CompMemory2CompBinary(UInt_t nrow,Byte_t *comp,UInt_t size)
{
  //Write RLE data in comp to binary file.
  //In order to distinguish these files from 10 bit data, 
  //a zero is written to the beginning of the file.

  Byte_t length = (Byte_t)nrow;
  Byte_t type = 0;
  if(fwrite(&type,1,1,fOutBinary)!=1) return kFALSE; //Write a zero, to mark that this file contains 8 bit data.
  if(fwrite(&length,1,1,fOutBinary)!=1) return kFALSE;
  if(fwrite(comp,size,1,fOutBinary)!=1) return kFALSE;
  return kTRUE;
}
