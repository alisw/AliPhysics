//-*- Mode: C++ -*-
// $Id$

///**************************************************************************
///* This file is property of and copyright by the ALICE HLT Project        * 
///* All rights reserved.                                                   *
///*                                                                        *
///* Primary Author: Jenny Wagner  (jwagner@cern.ch)                        *
///*                                                                        *
///* Permission to use, copy, modify and distribute this software and its   *
///* documentation strictly for non-commercial purposes is hereby granted   *
///* without fee, provided that the above copyright notice appears in all   *
///* copies and that both the copyright notice and this permission notice   *
///* appear in the supporting documentation. The authors make no claims     *
///* about the suitability of this software for any purpose. It is          * 
///* provided "as is" without express or implied warranty.                  *
///**************************************************************************

/// @file   AliHLTCOMPHuffmanAltro.cxx
/// @author Jenny Wagner
/// @date   29-08-2007
/// @brief  The Huffman compressor
///

#include "AliHLTCOMPHuffmanAltro.h"
#include "AliRawReaderMemory.h"
#include "AliAltroRawStreamV3.h"
#include <memory>

#if __GNUC__ >= 3
using namespace std;
#endif

#include <numeric>
using std::accumulate;

namespace
{
  // Helper class for std::accumulate algorithm.
  class AliHLTCOMPHuffmanOccurrenceSum {
  public:
    typedef int first_argument_type;
    typedef AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct second_argument_type;
    typedef bool result_type;
    int operator() (int a, AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct b) {
      return a+b.fabundance;
    }
  };
} // end of namespace

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTCOMPHuffmanAltro)

AliHLTCOMPHuffmanAltro::AliHLTCOMPHuffmanAltro()
  : AliHLTLogging()
  , fpRawReader(NULL)
  , fpAltroRawStream(NULL)
  , fTrainingMode(0)
  , fCompressionSwitch(0)
  , fPointer2InData(NULL)
  , fPointer2OutData(NULL)
  , fInputDataSize(0)
  , fOutputDataSize(0)
  , fNrcuTrailerwords(0)
  , fEntropy(0.0)
  , fVerbosity(0)
  , fTrainingTable(NULL)
  , fTranslationTable(NULL)
{
  // standard constructor
}

AliHLTCOMPHuffmanAltro::AliHLTCOMPHuffmanAltro(Bool_t compressionswitch, Bool_t trainingmode, AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* translationtable, Int_t nrcutrailerwords)
  : AliHLTLogging()
  , fpRawReader(NULL)
  , fpAltroRawStream(NULL)
  , fTrainingMode(trainingmode)
  , fCompressionSwitch(compressionswitch)
  , fPointer2InData(NULL)
  , fPointer2OutData(NULL)
  , fInputDataSize(0)
  , fOutputDataSize(0)
  , fNrcuTrailerwords(nrcutrailerwords)
  , fEntropy(0.0)
  , fVerbosity(0)
  , fTrainingTable(NULL)
  , fTranslationTable(translationtable)
{
  // constructor
}

AliHLTCOMPHuffmanAltro::~AliHLTCOMPHuffmanAltro()
{
  /// destructor
  if (fpAltroRawStream) delete fpAltroRawStream;
  fpAltroRawStream=NULL;
  if (fpRawReader) delete fpRawReader;
  fpRawReader=NULL;

  if (fTrainingTable) delete [] fTrainingTable;
  fTrainingTable=NULL;
}

int AliHLTCOMPHuffmanAltro::AddInputData(UChar_t* memory, ULong_t datasize, Int_t equipmentId)
{
  /// SetInputData takes input data pointer and input data size from component class
  if (memory == NULL) {
    return -EINVAL;
  };
  
  // check datasize < header + trailer
  if (datasize <= 32 + fNrcuTrailerwords) { // 8*32 (headerwords) = 8* 4*8 (headerbytes)
    HLTError("Error! Size of input data less than header and trailer");
    return -ENODATA;
  }

  // FIXME: fPointer2InData is obsolete
  fPointer2InData = (AliHLTUInt8_t*)memory;
  fInputDataSize = datasize; // size in bytes

  if (!fpRawReader) {
    std::auto_ptr<AliRawReaderMemory> rawReader(new AliRawReaderMemory);
    std::auto_ptr<AliAltroRawStreamV3> rawStream(new AliAltroRawStreamV3(rawReader.get()));

    if (!rawReader.get() || !rawStream.get()) {
      return -ENOMEM;
    }

    fpRawReader=rawReader.release();
    fpAltroRawStream=rawStream.release();
  }

  fpRawReader->AddBuffer(memory, datasize, equipmentId);

  return 0;
}

int AliHLTCOMPHuffmanAltro::Reset()
{
  /// Reset and prepare for new event
  if (fpRawReader) {
    fpRawReader->ClearBuffers();
  }
  if (fpAltroRawStream) {
    fpAltroRawStream->Reset();
  }

  InitNewTrainingTable();

  return 0;
}

int AliHLTCOMPHuffmanAltro::SetOutputData(AliHLTUInt8_t* outputdata, unsigned long outputsize)
{
  /// SetOuptputData takes output data pointer and outputsize from component class
  if (outputdata == NULL) {
    return -EINVAL;
  };

  fPointer2OutData = (AliHLTUInt64_t*) outputdata;
  fOutputDataSize = outputsize;

  // errors: compression: when outputdatasize < inputdatasize
  //         decompress.: when outputdatasize < inputMultiplier * inputdatasize
  // this should not happen with current config in component (inputMultiplier = 4.0 for decomp and 1.0 for comp)
  if (fCompressionSwitch) { // compression
    if(outputsize < fInputDataSize) {
      HLTError("Error! Size of output data not equal to size of input data for compression: reserved memory space might be too small");
      return -ENOSPC;
    }
  } else { // decompression
    if(outputsize < (fInputDataSize << 1)) { // smaller than inputdatasize * 2 not possible
      HLTError("Error! Size of output data too small for decompression");
      return -ENOSPC;
    }
  }

  return 0;
}

unsigned long AliHLTCOMPHuffmanAltro::GetOutputDataSize()
{
  // GetOutputDataSize delivers output data size for component class to arrange blocks for output writing
  // error when fOuputDataSize = 0; (not by training as there is no output when training)
  if((fOutputDataSize == 0) && (fTrainingMode == kFALSE)) {
    HLTError("Error! Calculated output data size is zero");
  };

  return fOutputDataSize;
}

int AliHLTCOMPHuffmanAltro::InitNewTrainingTable()
{
  // Initialise training table
  if (!fTrainingTable) fTrainingTable =  new AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct[TIMEBINS];
  if (!fTrainingTable) return -ENOMEM;
  
  //initialise this array:
  for(Int_t ii = 0; ii < TIMEBINS; ii++) {
    fTrainingTable[ii].famplitude = ii;   
    fTrainingTable[ii].fabundance = 0;
    fTrainingTable[ii].fcode = 2; //to assure that each leaf is assigned either 0 or 1!!!
  }

  return 0;
}

int AliHLTCOMPHuffmanAltro::SetHuffmanData(AliHLTCOMPHuffmanData* huffmandata) const
{
  // write produced code and occurrence table to Huffman Data
  return huffmandata->InitHuffmanData(fTrainingTable, fTranslationTable);
}

int AliHLTCOMPHuffmanAltro::GetTrainingTable(const AliHLTCOMPHuffmanData* huffmandata)
{
  // take translation table from HuffmanData
  huffmandata->GetOccurrenceTable(fTrainingTable);
  if (fVerbosity>0) Print("trainingtable");

  return 0;
}

int AliHLTCOMPHuffmanAltro::GetTranslationTable(const AliHLTCOMPHuffmanData* huffmandata)
{
  // initialise fTranslationTable
  if (fTranslationTable) delete fTranslationTable;
  fTranslationTable = new AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct[TIMEBINS];
  if (!fTranslationTable) return -ENOMEM;

  for(Int_t kk = 0; kk < TIMEBINS; kk++) {
    fTranslationTable[kk].famplitude = 0;
    fTranslationTable[kk].fhuffmancode = 0;
    fTranslationTable[kk].fvalidcodelength = 0;
  }

  // take translation table from HuffmanData
  huffmandata->GetCodeTable(fTranslationTable);
  if (fVerbosity>0) Print("translationtable");

  return 0;
}

/** ProcessData function to process data compression or decompression */
void AliHLTCOMPHuffmanAltro::ProcessData()
{
  // see header file for class documentation
  // set mode to process:
  if(fCompressionSwitch == kTRUE && fTrainingMode == kFALSE) // encoding without training
    {
      EntropyCompression();

      // information about compression ratio
       HLTDebug("original filesize (bytes) = %d", fInputDataSize);
       HLTDebug("compressed filesize (bytes) = %d", fOutputDataSize);
       HLTDebug("compression ratio after Huffman encoding is %f", ((Double_t) fOutputDataSize)/fInputDataSize);
    }
  else
    {
      if(fCompressionSwitch == kTRUE && fTrainingMode == kTRUE) // encoding with training
	{
	  TrainingData();

	  // information about calculated entropy (contained in GetEntropy()) and compression ratio
	   CalcEntropy();
	   GetEntropy();
	}
      else // decoding (no training possible, which is checked for errors by component)
	{
	 EntropyDecompression();

	 // information about compression ratio
	 HLTDebug("compressed filesize (bytes) = %d", fInputDataSize);
	 HLTDebug("filesize after decoding (bytes) = %d", fOutputDataSize);
	 HLTDebug("compression ratio calculated from Huffman decompressing has been %f",((Double_t) fInputDataSize)/fOutputDataSize);
	}
    }
}

/** GetEntropy returns entropy and prints it to screen */
Double_t AliHLTCOMPHuffmanAltro::GetEntropy()
{
  // see header file for class documentation
  // Information about calculated entropy
  HLTInfo("Calculated entropy =  %f",fEntropy);

  return fEntropy;
}

/** Compression: take raw input data, code table and put out entropy encoded data */
Int_t AliHLTCOMPHuffmanAltro::EntropyCompression()
{
  // see header file for class documentation
  // initialise input and output pointers
  AliHLTUInt32_t* pointer2InData = (AliHLTUInt32_t*) fPointer2InData;
  AliHLTUInt64_t* pointeroutputprop = fPointer2OutData;

  //first output word has to be initialised
  pointeroutputprop[0] = 0;

  // initialise output size
  fOutputDataSize = 0;

  // initialise required variables for input reading
  AliHLTUInt32_t data10 = 0; // temporary number to read out 10bit word
  UInt_t idxwd = 8;           // counts lines in input array (start after 8 x 32 bits header words)
  UInt_t idx10 = 0;           // counts 10-bits words
  UInt_t bitpswd = 0;         // counts current bitposition in idxwd

  // number of 32-bits trailer and 32-bits common data header words:
  // (taken out of compression!)
  // UInt_t fNrcuTrailerwords = 1; // determined by user input as argument
  UInt_t headerwords = 8;

  // initialise required variables for ouput creation:
  UInt_t idxoutwd = headerwords >> 1;  // counts lines in output array (start after header words)
  UInt_t bitpsoutwd = 0;               // counts current bitposition in idxoutwd

  // validation test variables:
  //  Int_t bitcounter = 0;           // additional test for number of read in 10 bit words
  // Int_t sumlength = 0;             // additional test for length of encoded file

  // error and abort when total number of bits cannot be divided by 32:
  // remember that fInputDataSize is given in BYTES!
  if (((fInputDataSize << 3) & 0x1F) != 0)
    {
      HLTError("Error! Input data size (bytes) %i cannot be divided into 32 bit words", fInputDataSize);

      return 1;
    };
  
  // error and abort when input data without trailer words cannot be divided into 10 bit words:
  if ((( (fInputDataSize << 3) - (fNrcuTrailerwords << 5) - (headerwords << 5)) % 10) != 0)
    {
      HLTError("Error! Input data (bytes) %i without trailer (bytes) %i and header words (bytes) %i cannot be divided into 10 bit words", fInputDataSize, (fNrcuTrailerwords << 2),(headerwords << 2));

      return 1;
    };
  
  // last line of 32-bits which is to be encoded
  UInt_t endNdx = (fInputDataSize>>2)-fNrcuTrailerwords;

  // write header to output (two 32-bits words input header in one 64-bits word output header):
  UInt_t headcntr = 0;
  UInt_t headwordcntr = 0;

  while(headcntr != (headerwords >> 1))
    {
      // write first header word to the left
      pointeroutputprop[headcntr] = ((AliHLTUInt64_t) pointer2InData[headwordcntr]);

      // write second header word to the right
      pointeroutputprop[headcntr] |= ( ((AliHLTUInt64_t) pointer2InData[headwordcntr+1]) << 32);

      // goto next header line and next two header words
       ++headcntr;
       headwordcntr += 2;
    }

    // EntropyCompression() for 10-bit values of ALTRO-Blocks: 
    while (idxwd < endNdx)
      {
	// write input data in data10 and cut out 10 bits with mask 0x3FFU = 0011 1111 1111
	data10 = (pointer2InData[idxwd] >> bitpswd ) & 0x3FFU;

	// validation test
	// cout << "10 bit data cut out from data10 =  " << data10 << endl;
	
	// if number of bits left in read input less than 10, get remaining bits from next 32-bits word
	if(bitpswd > 21)
	  {
	    // validation test
	    // cout << "adding bits from next word at bit position  " << bitpswd << endl;

	    data10 |= (pointer2InData[idxwd + 1]) << (32 - bitpswd);
	    data10 &= 0x3FFU;

	    // validation test
	    // cout << "complete 10-bits word = << data10 << endl;
	  };
	    
	// validation test
	// print first 100 10-bits words with resp.code and validcodelength to screen...
	// if (idx10 < 101)
	//  {
	//    cout << "read 10-bits word (HEX) =  " << hex << data10 << dec <<" with code (DEC) =  " << fTranslationTable[data10].fhuffmancode << "  and validcodelength (DEC) =  " << fTranslationTable[data10].fvalidcodelength << endl;
	// };

	  // start encoding of codes with less than 65 bits codelength
	  if(fTranslationTable[data10].fvalidcodelength < 65)
	    {
	      // error and abortion when validcodelength = 0
	      if (fTranslationTable[data10].fvalidcodelength == 0)
		{
		  HLTError("Error! Valid codelength of read 10 bit word (DEC) %i is zero", data10);

		  return 1;
		}

	      // validation test
	      // bitcounter += 1;
	      // sumlength += fTranslationTable[data10].fvalidcodelength;
	            
	      if(fTranslationTable[data10].fvalidcodelength + bitpsoutwd < 64) //if new code fits in current output line
		{
		  // validation test
		  // cout << "writing code (HEX) =  " << hex << fTranslationTable[data10].fhuffmancode << dec << " to line (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;

		  // shift line to right
		  pointeroutputprop[idxoutwd] <<= fTranslationTable[data10].fvalidcodelength;
	      
		  // append code at left side of line
		  pointeroutputprop[idxoutwd] |= fTranslationTable[data10].fhuffmancode;

		  // validation test
		  // cout << "code written and current line updated to (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;
		}
	      else //if not, write lowest bits to next line
		{
		  // validation test
		   // cout << "writing bits of code (HEX) =  " << hex << fTranslationTable[data10].fhuffmancode << dec << " to line (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;

		  // create temporary variable for performance gain:
		  Int_t restcode = fTranslationTable[data10].fvalidcodelength - 64 + bitpsoutwd;
		  
		  // shift line completely to the right
		  pointeroutputprop[idxoutwd] <<= (64 - bitpsoutwd);
		  
		  // append upper bits of code to this line
		  pointeroutputprop[idxoutwd] |= (fTranslationTable[data10].fhuffmancode >> (restcode));
		  
		  // validation test
		  // cout << "code bits written and current line updated to (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;

		  // write lowest code bits to next line
		  
		  ++idxoutwd;
	      
		  // trick to cut out already written bits from code		  
		  AliHLTUInt64_t tempbuffer = fTranslationTable[data10].fhuffmancode;
		  tempbuffer <<= 64-bitpsoutwd; // shift away written code
		  tempbuffer >>= 64-bitpsoutwd; // shift back to get lowest bits
		  
		  pointeroutputprop[idxoutwd] = tempbuffer;

		  // validation test
		  // cout << "code bits written and next line started as (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;
		}
	    }
	  else   // if code is longer than 64 bits: use *morecode (not tested therefore left away)
	    {
	      HLTError("Error! Huffman Code Table does not contain codes with valid codelength longer than 64 bits");

	      return 1;

	    }
      
	  // validation test
	  // cout << " current positions:" << endl;
	  // cout << "input :  idxwd (DEC) =  " << idxwd << "   bitpswd (DEC) =  " << bitpswd << "   idx10 (DEC) =  " << idx10 << endl;
	  // cout << "output :  idxoutwd (DEC) =  " << idxoutwd << "   bitpsoutwd (DEC) =  " << bitpsoutwd << endl;
	 
	  // calculate new positions
	  bitpsoutwd = (fTranslationTable[data10].fvalidcodelength + bitpsoutwd) & 0x3F;
	  
	  ++idx10; // go on reading 10 bit word
	  
	  // bit position word reads position in 32bit word when 10 bits are read out  
	  bitpswd = (idx10 * 10) & 0x1F;

	  // index word reads position of 32bit word in input array
	  idxwd = ((idx10 * 10)>>5)+8; 

	  // validation test
	  //if(idxoutwd > 2635)
	  // {
	  //   cout << " encoding value (HEX) = " << hex << data10 << dec << endl;
	  //   cout << " code (HEX) = " << hex << fTranslationTable[data10].fhuffmancode << dec << endl;
	  //   cout << " valid code length (DEC) = " << fTranslationTable[data10].fvalidcodelength << endl;
	  //   cout << " new positions:" << endl;
	  //   cout << "input :  idxwd (DEC) =  " << idxwd << "   bitpswd (DEC) =  " << bitpswd << "   idx10 (DEC) =  " << idx10 << endl;
	  //   cout << "output :  idxoutwd (DEC) =  " << idxoutwd << "   bitpsoutwd (DEC) =  " << bitpsoutwd << endl; 
	  // }

      }// end of encoding (while loop)
    
    // if new line is already started
    if(bitpsoutwd == 0)
      {
	--idxoutwd;
	bitpsoutwd = 64;
      }
    else
      {
	// fill rest of line with zeroes
	pointeroutputprop[idxoutwd] <<= (64 - bitpsoutwd);
      }

    // next 32 bits are to memorise first bitposition after the code
    UInt_t lastvalidbitps = bitpsoutwd;     
    
    // validation test
    // cout << "last bitpsoutwd = lastvalidbitps (DEC) =  " << bitpsoutwd << endl;
    
    // validation test
    // cout << "second last code line (HEX) = " <<  hex << pointeroutputprop[idxoutwd-1] << dec << endl;
    // cout << "last code line (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;
    // cout << "idxoutwd (DEC) =  " << idxoutwd << endl;
    
    // start new line and write last valid bit position and trailer words there:
    bitpsoutwd = 32;
    ++idxoutwd;
    
    pointeroutputprop[idxoutwd] = lastvalidbitps;
    
    // validation test
    // cout << "first trailer line with lastvalidbitps (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;
    // cout << "idxoutwd (DEC) =  " << idxoutwd << endl;
    
    // write first trailerword
    bitpsoutwd = 64;
    pointeroutputprop[idxoutwd] <<= 32;
    pointeroutputprop[idxoutwd] |= pointer2InData[endNdx];
    
    // validation test
    // cout << "first trailer line with lastvalidbitpos and first trailerword (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;
    
    if(fNrcuTrailerwords > 1) // if there is more than one trailerword
      { 
	++idxoutwd;
	bitpsoutwd = 32;
	
	pointeroutputprop[idxoutwd] = pointer2InData[endNdx + 1];  // write second
	
	// validation test
	// cout << "last trailer line with second trailerword (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;
	// cout << "idxoutwd (DEC) =  " << idxoutwd << endl;
	
	if (fNrcuTrailerwords == 3) // write third
	  {
	    bitpsoutwd = 64;
	    
	    pointeroutputprop[idxoutwd] <<= 32;
	    pointeroutputprop[idxoutwd] |= pointer2InData[endNdx +2];
	    
	    // validation test
	    // cout << "last trailer line with last trailerword (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;
	  }
      }
    
    // validation tests
    // cout << "file ended at current positions:" << endl;
    // cout << "input :  idxwd (DEC) =  " << idxwd << "   bitpswd (DEC) =  " << bitpswd << "   idx10 (DEC) =  " << idx10 << endl;
    // cout << "output :  idxoutwd (DEC) =  " << idxoutwd << "   bitpsoutwd (DEC) =  " << bitpsoutwd << endl; 
    // cout << "current output line (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;
    // cout << "sum of codelengths (DEC) =  " << sumlength << endl;
    // cout << "number of read 10 bit data (DEC) =  " << bitcounter << endl;
    // cout << "first input trailer (HEX) =  " << hex << pointer2InData[endNdx] << dec << endl;

    // set fOutputDataSize in byte:
    // idxoutwd*8 = number of whole 64 bits words in byte
    // bitpsoutwd/8 = number of wohle 8 bits words within last idxoutwd
     fOutputDataSize = (Int_t) ((idxoutwd << 3) + (bitpsoutwd >> 3));

    // error and abortion when output data size for encoding is larger than input data size
    if(fOutputDataSize > fInputDataSize)
      {
	HLTError("Error! Calculated data size for compressed output stream is larger than input data size");

	return 1;
      };
    
    // validation test
    // cout << "output data size (DEC) =  " << fOutputDataSize << endl;

    return 0;
}

/** Mergesort used by TrainingData to sort an array with n entries from high abundance to low abundance */
AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* AliHLTCOMPHuffmanAltro::Mergesort(AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct *unsortedarray, Int_t n)
{
  // see header file for class documentation
  //divide array into two halfs: left and right (fleft.size = divlength)
  Int_t divlength = n >> 1;
      
  if (n==1)
    {
      // error when unsorted array = NULL
      if (unsortedarray == NULL)
	{
	  HLTError("Error! Pointer to final merge sorted abundance array = NULL");
	};

      return unsortedarray;
    }
  else
    {
      // sort both halfs recursively:
      Mergesort(unsortedarray, divlength);  //left half
      Mergesort(&unsortedarray[divlength], n-divlength); //right half

      // merge together:
      // create temporary result array:   
      AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* temp = new AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct[n];

      // counters:
      Int_t ii, jj, kk;

      // if left and right halves both have elements: chose the smaller one:
      for (ii = 0, jj = divlength, kk = 0; ii < divlength && jj < n;)
        {
	  if (unsortedarray[ii].fabundance > unsortedarray[jj].fabundance)
	    { 
	      temp[kk] = unsortedarray[ii];
	      ++ii;
	    }
          else 
	    {
	      temp[kk] = unsortedarray[jj];
	      ++jj;
	    }
	  
          // increase kk
	  ++kk;
	}     
      
      // if one half is empty take elements of the other one:
      while (ii < divlength)
        {
          temp[kk] = unsortedarray[ii];
	  ++kk;
	  ++ii;
        }
      
      while (jj < n) 
        {
          temp[kk] = unsortedarray[jj];
	 ++kk;
	 ++jj;
        }
      
      // copy sorted temp array back into original data array
      for (Int_t ll = 0; ll < n; ll++)
        {
          unsortedarray[ll] = temp[ll];
        }
   
      // free space
      delete[] temp;
         
      // return pointer to original data array (which is sorted now)
      return unsortedarray;
    }
}

/** CreateHuffmanTree used by TrainingData to create the binary Huffman tree */
AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* AliHLTCOMPHuffmanAltro::CreateHuffmanTree(AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* listroot,AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* lastelement, Int_t n)
{
  // see header file for class documentation
  // initialise pointer to go through list
  //AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* nextpointer = listroot; // pointer for validation test below
 AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* lastpointer = lastelement;

  // build Huffman tree while there are still elements in the list
  while(n > 2)
    {
      // validation test
      // cout << "current number of remaining list elements n (DEC) =  " << n << endl;
      // cout << "current abundance of last list element (DEC) =  " << lastpointer->fleafcontents.fabundance << endl;
      // cout << "print all list elements from high to low to screen:" << endl;
      // while(nextpointer != NULL)
      //    { 
      //      cout << nextpointer->fleafcontents.fabundance << endl; 
      //      nextpointer = nextpointer->fnext; 
      //    }
      // cout << " end of list " << endl;

      // create new tree element from last two entries in orderedarray
     AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* temptree = new AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct;
     // error if temptree = NULL
     if (temptree == NULL)
       {
	 HLTError("Error! Allocation of Huffman binary tree failed");
	 return NULL;
       };

      // initialise the new element:
      temptree->fleafcontents.famplitude = TIMEBINS; //internal nodes mustn't be confused with real amplitudes (from 0 to bitsize)!

      // validation test
      // cout <<"current abundance of last list element (DEC) =  " << lastpointer->fleafcontents.fabundance << endl;

      // initialise abundance starting with last element
      temptree->fleafcontents.fabundance = lastpointer->fleafcontents.fabundance;

      // code assignment small ones get 0, large ones get 1:
      lastpointer->fleafcontents.fcode = 0;

      // append smallest element on the left and set pointer temptree = parent
      temptree->fleft = lastpointer;
      lastpointer->fparent = temptree;

      // go on to second last element and set abundance, code and append on right side of the new element
      lastpointer = lastpointer->fprevious;

      // validation test
      //  cout << "current abundance of second last list element (DEC) =  " << lastpointer->fleafcontents.fabundance << endl;

      cout.flush();

      temptree->fleafcontents.fabundance += lastpointer->fleafcontents.fabundance;
      lastpointer->fleafcontents.fcode = 1;
      temptree->fright = lastpointer;
      lastpointer->fparent = temptree;

      temptree->fnext = NULL;
      temptree->fprevious = NULL; 
      temptree->fparent = NULL;

      // validation test
      // cout <<   "current abundance of temptree element = sum of leaf abundances (DEC) =  " << temptree->fleafcontents.fabundance << endl;

      // write temptree at the suitable place according to its abundance in the list
      // initialise searchpointer to go through list
     AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* search = listroot;
     AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* sorting = listroot; // pointer from previous element
    
      // check if listroot[0].fleafcontents.fabundance < temptree.fleafcontents.fabundance,
      // if so, make temptree new root of list
      if (temptree->fleafcontents.fabundance > listroot[0].fleafcontents.fabundance)
	{
	  // validation test
	  // cout << "abundance of new element (DEC) =  " << temptree->fleafcontents.fabundance << "   is larger than root abundance (DEC) =  " << listroot[0].fleafcontents.fabundance << endl;

	  // temptree = root
	  temptree->fnext = search; 
	  search->fprevious = temptree;
	  //temptree->fprevious = NULL // for first list element is already done!
	  listroot = temptree;
	}
      else //sort in at the suitable point in the list
	{
	  search = listroot->fnext; // go one further than listroot!

	  while(search != NULL)
	    {
	      // validation test
	      // cout << "current abundance of searchpointer (DEC) =  " << searchpointer->fleafcontents.fabundance << endl;

	      // if abundance of list entry > abundance of new element, go to next element,
	      // else sort in at that point
	      if(search->fleafcontents.fabundance > temptree->fleafcontents.fabundance)
		{
		  // validation test		  
		  // cout << "abundance of searchpointer (DEC) =  " << search->fleafcontents.fabundance << "   is larger than temptree abundance (DEC) =  " << temptree->fleafcontents.fabundance << endl;

		  sorting = search;        // remember previous element
		  search = sorting->fnext;  // goto next element
		}
	      else 
		{
		  sorting->fnext = temptree;     //insert temptree to previous
		  search->fprevious = temptree;  //insert temptree to next
		  temptree->fprevious = sorting; // insert temptree to previous
		  temptree->fnext = search;      //insert temptree to next 
		 
		  search = NULL;                //stop sorting in
		}
	    }	  
	} //end of sorting in

      // cut the two elements out of the list:
  
      lastpointer = lastpointer->fprevious;
      lastpointer->fnext = NULL;
      
      // validation test
      // cout << "cutting out last two list elements, abundance of current last list element is (DEC) =  " << lastpointer->fleafcontents.fabundance << endl;
      
      cout.flush();
      
      // validation test
      // if(n <= 15)
      // {
      //   cout << "current list with new (temptree) element sorted in:" << endl;
	  
      //  AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct * showpointer = listroot;
      //   while((showpointer->fnext != NULL) && (showpointer->fleafcontents.fabundance != 0))	  
      //    {
      //      cout << "amplitude (DEC) =  " << showpointer->fleafcontents.famplitude << "   abundance (DEC) =  " <<  showpointer->fleafcontents.fabundance << endl; 
	      
      //      showpointer = showpointer->fnext;
      //    }
	  
      //   cout << "amplitude (DEC) =  " << showpointer->fleafcontents.famplitude << "   abundance (DEC) =  " <<  showpointer->fleafcontents.fabundance << endl; 
      // };
      
      // perform createHuffmanTree with n-1 elements:
      --n;
    }
      
  // termination for n = 2:
 
  // create new tree element from last two entries
 AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* temptree = new AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct;
  // error if temptree = NULL
  if (temptree == NULL)
    {
      HLTError("Error! Pointer to root of Huffman binary tree = NULL");
      return NULL;
    };
  
  // initialise the new element:
  temptree->fleafcontents.famplitude = TIMEBINS; //internal nodes mustn't be confused with real event names (from 0 to TIMEBINS)!
 
  // validation test 
  // cout << "assemble last two elements to final tree:" << endl;
  // cout << "abundance of last list element (DEC) =  " << lastpointer->fleafcontents.fabundance << endl;
  
  // initialise abundance starting with last element
  temptree->fleafcontents.fabundance = lastpointer->fleafcontents.fabundance;
  
  // code assignment small ones get 0, large ones get 1:
  lastpointer->fleafcontents.fcode = 0;
  
  // append smallest element on the left:
  temptree->fleft = lastpointer;
  lastpointer->fparent = temptree;
 
  // validation test
  // cout << "abundance of first list element (DEC) =  " << listroot->fleafcontents.fabundance << endl;

  cout.flush();
 
  // go on to second last element and set abundance, code and append on right side of the new element 
  temptree->fleafcontents.fabundance +=listroot->fleafcontents.fabundance;
  listroot->fleafcontents.fcode = 1;
  temptree->fright = listroot;
  listroot->fparent = temptree;
  
  temptree->fnext = NULL;
  temptree->fprevious = NULL; 
  temptree->fparent = NULL;
  
  // validation test
  // cout << "  final abundance of tree root (DEC) =  " << temptree->fleafcontents.fabundance << endl;
  
  return temptree;
}

/** HuffmanCode used by TrainingData to create the Huffman code for all leaves of the HuffmanTree */
Int_t AliHLTCOMPHuffmanAltro::HuffmanCode(AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* treeroot, AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* HuffmanCodearray)
{
  // see header file for class documentation
  while ((treeroot->fleft != NULL) || (treeroot->fright != NULL))
    {
      // define pointer to go through tree
     AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* finder = treeroot;
      
      // define temporary HuffmanCode struct to sort into the Huffmanarray later
      // array necessary for codelengths > 64 bits
      Int_t tempstructsize = (Int_t) ceil((double) TIMEBINS/64); // maximal codelength: TIMEBINS
      
      AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* tempstruct = new AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct[tempstructsize];
      if (!tempstruct) return -1;

      for(Int_t jj = 0; jj < tempstructsize; jj++)
	{
	  tempstruct[jj].fhuffmancode = 0;
	}

      Int_t passednodes = 0;

      // start searching for leaves
      while (1)
	{
	  // if finder points at a leaf:
	  if (finder->fright == NULL && finder->fleft == NULL)
	    { 
	      //if it is an internal node
	      if(finder->fleafcontents.famplitude == TIMEBINS) 
		{
		  // validation test  
		  // cout << "found internal node with following abundance to delete (DEC) =  " << finder->fleafcontents.fabundance <<  endl; 

		  //go back to parent and delete node with 256
		  if (finder->fleafcontents.fcode == 0) // node was left one
		    {
		     AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* parent = finder->fparent;
		      parent->fleft = NULL;
		    }
		  else // node was right one
		    {
		     AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* parent = finder->fparent;
		      parent->fright = NULL;
		    }	 
		  
		  break;
		}
	      else // if it is a leaf
		{
		  Int_t checkright = 0; // so leaf is left children
		  
		  // leaf can also be a right leaf:
		  if (finder->fleafcontents.fcode == 1)
		    {
		      checkright = 1; // if code == 1, leaf is right children
		    };
		  
		  // write collected data in the respective array:

		  // validation test if amplitudes are correctly found (after being initialised to TIMEBINS+1 in TrainingData()
		  // HuffmanCodearray[finder->fleafcontents.famplitude].famplitude = finder->fleafcontents.famplitude;
		  
		  // validation test
		  // cout << "found leaf with amplitude (DEC) =  " << finder->fleafcontents.famplitude << endl;
		  
		  HuffmanCodearray[finder->fleafcontents.famplitude].fvalidcodelength = passednodes;

		  // validation test
		  // cout << "found leaf with validlength (DEC) =  " << HuffmanCodearray[foundleaves].fvalidcodelength << endl;
		 
		  //write code:
		  HuffmanCodearray[finder->fleafcontents.famplitude].fhuffmancode = tempstruct[0].fhuffmancode;

		  // validation test
		  //cout << "found leaf with code (HEX) =  " << hex << HuffmanCodearray[finder->fleafcontents.famplitude].fhuffmancode << dec << endl;

		  if(passednodes > 64) // if there is more code to write (not usual in normal zero-suppressed data streams!) 
		    {

		      HLTError("Error! Valid codelength for datum (DEC) %d is larger than 64 bits, which is not usual for normal input data", finder->fleafcontents.famplitude);

		      delete[] tempstruct;  
		      return 1;
		    }
		  
		  // validation test  
		  // cout << "found leaf written after passed nodes (DEC) =  " << passednodes << endl;
		  
		  // deleting found leaf from tree
		 AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* parent = finder->fparent; // go one up and set pointer to leaf to zero:
		  
		   // validation test 		  
		   // cout << "parent abundance from found leaf (DEC) =  " << finder->fleafcontents.fabundance << endl;

		  if (checkright == 1) // if leaf is right children -> delete right
		    {
		      parent->fright = NULL;   
		    }
		  else 	//else: delete fleft
		    {
		      parent->fleft = NULL;
		    }
		  
		  break;		  
		} 
	      // validation test
	      // cout << " finder pointer at end of deleting the found leaf (DEC) =  " << finder << endl; 
	    }
	  else // if node is not a leaf     
	    {
	      // validation test
	      //cout << "finder pointer to left child (DEC) =  " << finder->fleft << endl;
	      //cout << "finder pointer to right child (DEC) =  " << finder->fright << endl;
	      
	      if(finder->fleft != NULL)
		{
		  finder = finder->fleft;  // pointer to children

		  // validation test
		  //cout << "left child has abundance (DEC) =  " << finder->fleafcontents.fabundance << endl;
		 
		  // shift left to arrange space for new code element which is zero (->left!)
		  // i.e. nothing further needed to create code

		  // validation test
		  //cout << "code for left child is created as (DEC) =  " << tempstruct->fhuffmancode << endl;
		  
		  //problem if passednodes > 63 --> inserting another bit becomes tricky:
		  if(passednodes < 64)
		    {
		      tempstruct[0].fhuffmancode = tempstruct[0].fhuffmancode << 1; 
		    }
		  else // error as valid codelength should not exceed 64 bits
		    {

		      HLTError("Error! Valid codelength for current encoding process becomes larger than 64 bits, which is not usual for normal input data");

		      delete[] tempstruct;  
		      return 1;
		    }

		  // validation test
		  //cout << "code for left child has been written as (DEC) =  " << tempstruct->fhuffmancode << endl;
		  
		 ++passednodes;
		}
	      else
		{
		  finder = finder->fright; // goto right children

		  // validation test
		  //cout << "right child has abundance (DEC) =  " << finder->fleafcontents.fabundance << endl;

		  // shift left to arrange space for new code element which is one (->right!)
		  // i.e. append 1 at the right end of the code via OR-function

		  // validation test
		  //cout << "code for right child is created as (DEC) =  " << tempstruct->fhuffmancode << endl;

		  //problem if passednodes > 63 --> inserting another bit becomes tricky:
		  if(passednodes < 64)
		    {
		      tempstruct[0].fhuffmancode = tempstruct[0].fhuffmancode << 1; 
		      tempstruct[0].fhuffmancode = (tempstruct[0].fhuffmancode) | 1;
		    }
		  else
		    {

		      HLTError("Error! Valid codelength for current encoding process becomes larger than 64 bits, which is not usual for normal input data");

		    }
				  
		  // validation test
		  //cout << "code for right children has been written as (DEC) =  " << tempstruct->fhuffmancode << endl;
		  
		  ++passednodes;
		} // end of right children search
	    } // end of no-leaf-search
	} // end of while-loop for leaf-search

      //freespace
      delete[] tempstruct;  
    } // end of while-loop for tree != root only
    
    // error when there is an amplitude = TIMEBINS in the HuffmanArray
    for(int jj = 0; jj < TIMEBINS; jj++)
      {
    	if(HuffmanCodearray[jj].famplitude == TIMEBINS)
	  {
	    HLTError("Error! Internal node from Huffman tree in Huffmanarray for 10 bit word (DEC) %i", jj);
	  };
      }
    
    return 0;
}

Int_t AliHLTCOMPHuffmanAltro::TrainingData()
{ 
  // fill training data to create the HuffmanCode from a binary tree
  int iResult=0;

  // total number of events counted
  Int_t totalnumber = 0;

  if (!fpAltroRawStream) return -ENODATA;
  fpAltroRawStream->Reset();
 
  // initialise required variables for input reading:
  AliHLTUInt32_t data10 = 0; // temporary number to read out 10bit word with

  if (!fTrainingTable && (iResult=InitNewTrainingTable())<=0)
    return -iResult;

  // read altro data
  while (fpAltroRawStream->NextDDL()) {
    while (fpAltroRawStream->NextChannel()) {
      int iNofBunches=0;
      // include payload size in training table
      data10=fpAltroRawStream->GetChannelPayloadSize();
      //fTrainingTable[(Int_t)(data10)].fabundance += 1;
      if (fpAltroRawStream->NextBunch()) {
	do {
	  // include bunch length in training table
	  int iBunchLen=fpAltroRawStream->GetBunchLength();
	  data10=iBunchLen;
	  fTrainingTable[(Int_t)(data10)].fabundance += 1;

	  // include start time in training table
	  int iStartTime=fpAltroRawStream->GetStartTimeBin();
	  data10=iStartTime;
	  //fTrainingTable[(Int_t)(data10)].fabundance += 1;

	  const UShort_t* signals=fpAltroRawStream->GetSignals();
	  for (int sigpos=0; sigpos<iBunchLen-2; sigpos++) {
	    data10=signals[sigpos];
	    // include signals in training table
	    fTrainingTable[(Int_t)(data10)].fabundance += 1;
	  }
	  iNofBunches++;
	} while (fpAltroRawStream->NextBunch());
      }
      
      // increase total number of events
      totalnumber += 1;
    }
  }

  fOutputDataSize = 0;

  return 0;
}

/** Create Huffman code table out of abundance table filled by training data in DoEvent function of component */
Int_t AliHLTCOMPHuffmanAltro::CreateCodeTable()
{
  // see header file for class documentation
  // using merge sort to sort the list: (stable algorithm, quick even in worst case O(nlogn))
  AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* HuffmanArraySorted = Mergesort(fTrainingTable, TIMEBINS);

  // abort when there is no pointer to sorted array (error produced in Mergesort function)
  if(HuffmanArraySorted == NULL)
    {
      return 1;
    }

  // error and abort if list not properly sorted:
  // go through list and watch if current abundance of HuffmanData is always larger than next element
  for(int kk = 0; kk < TIMEBINS - 1; kk++)
    {
      if(HuffmanArraySorted[kk].fabundance < HuffmanArraySorted[kk+1].fabundance)
	{
	  HLTError("Error! List of 10 bit word abundances not properly sorted by merge sort", "Element with 10 bit word (DEC) %i has abundance %i which is smaller than element with 10 bit word (DEC) %i with abundance %i",HuffmanArraySorted[kk].famplitude, HuffmanArraySorted[kk].fabundance, HuffmanArraySorted[kk+1].famplitude, HuffmanArraySorted[kk+1].fabundance);

	  return 1;
	};
    }

  // validation test
  //cout << "merge sort for training table is now done " << endl;
  
  // after mergesorting, zero arrays are last ones in array with abundance = 0
  //  -> count filled array entries
  UInt_t filled = TIMEBINS; // to ensure that every amplitude (no matter if it appears in training or not) gets a code
  UInt_t zeroentries = 0;  // number of non used 10 bit words (with abundance = 0)

  // check how many 10 bit words are used ("missing channel problem)
  for(UInt_t kk = 0; kk < TIMEBINS; kk++)
    {
      if((HuffmanArraySorted[kk].fabundance == 0))
  	{
	  ++zeroentries;
  	};
    }
  
  // -> new array size is "filled"
  // validation test
  // cout << "number of filled array entries (DEC) =  " << filled << endl;
  
  // warning when there are a lot of non used 10 bit words -> this training table cannot be used in general case and taken to encode other files
  if(zeroentries > 50)
    {
      HLTWarning("Warning! Only %i different 10 bit words out of %i are used, i.e. created Huffman Code table might not be optimal to encode other data", TIMEBINS-zeroentries, TIMEBINS);
    }
  
  // initialise leaves of the tree as list (= queue) ofAliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct,
 AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct * HuffmanTreeList = new AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct[filled];
  
  // initialise first element
  HuffmanTreeList[0].fleafcontents = HuffmanArraySorted[0];
  HuffmanTreeList[0].fleft = NULL;
  HuffmanTreeList[0].fright = NULL;
  HuffmanTreeList[0].fnext =  &HuffmanTreeList[1];
  HuffmanTreeList[0].fprevious = NULL;
  
  // validation test
  // write list to screen
  // cout <<"Amplitude  " << HuffmanTreeList[0].fleafcontents.famplitude << " | " << "Abundance   "<<  HuffmanTreeList[0].fleafcontents.fabundance << " | " << "Left " << HuffmanTreeList[0].fleft << " | " << "Right " << HuffmanTreeList[0].fright << " | " << "Next " << HuffmanTreeList[0].fnext << endl; 
  
  // initialise 1 to filled-2 elements
  UInt_t kk = 1;
  while(kk < filled-1)
    {
      HuffmanTreeList[kk].fleafcontents = HuffmanArraySorted[kk];
      HuffmanTreeList[kk].fleft = NULL;
      HuffmanTreeList[kk].fright = NULL;
      HuffmanTreeList[kk].fnext =  &HuffmanTreeList[kk+1];
      HuffmanTreeList[kk].fprevious = &HuffmanTreeList[kk-1];
      
      // validation test
      // write list to screen
      // cout <<"Amplitude  " << HuffmanTreeList[kk].fleafcontents.famplitude << " | " << "Abundance   "<<  HuffmanTreeList[kk].fleafcontents.fabundance << " | " << "Left " << HuffmanTreeList[kk].fleft << " | " << "Right " << HuffmanTreeList[kk].fright << " | " << "Next " << HuffmanTreeList[kk].fnext << endl; 
      
      kk += 1;
    }
  
  // initialise last list entry with next pointer to NULL
  HuffmanTreeList[filled-1].fleafcontents = HuffmanArraySorted[filled-1];
  HuffmanTreeList[filled-1].fleft = NULL;
  HuffmanTreeList[filled-1].fright = NULL;
  HuffmanTreeList[filled-1].fnext = NULL;
  HuffmanTreeList[filled-1].fprevious = &HuffmanTreeList[filled-2];
  
  // initialise pointer to last list element:
 AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* lastelement = &HuffmanTreeList[filled-1];
  
  // validation test
  // write last element to screen
  //  cout <<"Amplitude  " << HuffmanTreeList[filled-1].fleafcontents.famplitude << " | " << "Abundance   "<<  HuffmanTreeList[filled-1].fleafcontents.fabundance << " | " << "Left " << HuffmanTreeList[filled-1].fleft << " | " << "Right " << HuffmanTreeList[filled-1].fright << " | " << "Next " << HuffmanTreeList[filled-1].fnext << endl; 
  
  //use this sorted list to build up the binary tree with root pointer to the beginning:
 AliHLTCOMPHuffmanData::AliHLTCOMPHuffmanTreeDataStruct* root = CreateHuffmanTree(HuffmanTreeList, lastelement, filled);

  // abort if root = NULL (error already produced in CreateHuffmanTree function)
 delete [] HuffmanTreeList;

  if(root == NULL)
    {
      return 1;
    };

  // validation test
  // cout << "binary tree is now created " << endl;

  // create an array for the HuffmanCode data:
  fTranslationTable = new AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct[TIMEBINS];

  // initialise the array with validcodelengths 0 and codes = 0
  for(Int_t kk1 = 0; kk1 < TIMEBINS; kk1++)
    {
      // validation test for correct HuffmanCode()
      // fTranslationTable[kk1].famplitude = TIMEBINS+1; // check if eventnames are written

      fTranslationTable[kk1].famplitude = kk1;
      fTranslationTable[kk1].fhuffmancode = 0;
      fTranslationTable[kk1].fvalidcodelength = 0;
      // fTranslationTable[kk1].morecode = NULL;
    }

  // create HuffmanCode and abort when HuffmanCode produces errors
  if(HuffmanCode(root,fTranslationTable) == 1)
    {
      return 1;
    };

  // validation test
  // cout << "Huffman coding is now done " << endl;

  // validation test
  // print the Huffman code table to screen
  // for(Int_t kk = 0; kk < TIMEBINS; kk++)
  //  {
  //    if (fTranslationTable[kk].fvalidcodelength != 0)
  //     cout << fTranslationTable[kk].famplitude << "   |   " << fTranslationTable[kk].fhuffmancode << "   |   " << fTranslationTable[kk].fvalidcodelength << endl;
  //  }
  
  // findout maximal and minimal codelength and print them out
    UInt_t maxcodelength = fTranslationTable[0].fvalidcodelength;
    UInt_t mincodelength = TIMEBINS;

    for (Int_t kk2 = 0; kk2 < TIMEBINS; kk2++)
      {
	//	if(fTranslationTable[kk2].fvalidcodelength != 0) {cout << kk2 << " |  " << fTranslationTable[kk2].fhuffmancode << " |  " << fTranslation "Table[kk2].fvalidcodelength << endl;}

	if(fTranslationTable[kk2].fvalidcodelength > maxcodelength)
	  { maxcodelength = fTranslationTable[kk2].fvalidcodelength;};
	if( (fTranslationTable[kk2].fvalidcodelength != 0) && (fTranslationTable[kk2].fvalidcodelength < mincodelength) )
	  { mincodelength = fTranslationTable[kk2].fvalidcodelength;};
      }

    // print results to screen
    HLTInfo("maximal codelength (DEC) = %i ", maxcodelength);
    HLTInfo("minimal codelength (DEC) = %i ", mincodelength);
  
    return 0;
}

/** TTMergesort used by EntropyDecoding to sort translation table according to validcodelength from low to high */
AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* AliHLTCOMPHuffmanAltro::TTMergesort(AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct *unsortedarray, Int_t n)
{
  // see heaeder file for class documentation
  //divide array into two halfs: left and right (fleft.size = divlength)
  Int_t divlength = n >> 1;
      
  if (n==1)
    {
      // error when unsorted array = NULL
      if (unsortedarray == NULL)
	{
	  HLTError("Error! Pointer to final merge sorted code table = NULL");
	};

      return unsortedarray;
    }
  else
    {
      // sort both halfs recursively:
      TTMergesort(unsortedarray, divlength);  //left half
      TTMergesort(&unsortedarray[divlength], n-divlength); //right half

      // merge together:
      // create temporary result array:   
      AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct* temp = new AliHLTCOMPHuffmanCodeData::AliHLTCOMPHuffmanCodeStruct[n];

      // counters:
      Int_t ii, jj, kk;

      // if left and right halves both have elements: chose the smaller one:
      for (ii = 0, jj = divlength, kk = 0; ii < divlength && jj < n;)
        {
	  if (unsortedarray[ii].fvalidcodelength < unsortedarray[jj].fvalidcodelength)
	    { 
	      temp[kk] = unsortedarray[ii];
	      ++ii;
	    }
          else 
	    {
	      temp[kk] = unsortedarray[jj];
	      ++jj;
	    }
	  
          // increase kk
	  ++kk;
	}     
      
      // if one half is empty take elements of the other one:
      while (ii < divlength)
        {
          temp[kk] = unsortedarray[ii];
	  ++kk;
	  ++ii;
        }
      
      while (jj < n) 
        {
          temp[kk] = unsortedarray[jj];
	 ++kk;
	 ++jj;
        }
      
      // copy sorted temp array back into original data array
      for (Int_t ll = 0; ll < n; ll++)
        {
          unsortedarray[ll] = temp[ll];
        }
   
      // free space
      delete[] temp;
         
      // return pointer to original data array (which is sorted now)
      return unsortedarray;
    }
}

/** function for the entropy decoding, needed input: code table, input data pointer and size and pointer where to write decoded data */
Int_t AliHLTCOMPHuffmanAltro::EntropyDecompression()
{
  // see heaeder file for class documentation
  // validation test
  // print translation table to screen
  // for(int kk= 0; kk < TIMEBINS; kk++)
  // {
  //   cout << fTranslationTable[kk].famplitude << "  |  " << fTranslationTable[kk].fhuffmancode << "  |  " << fTranslationTable[kk].fvalidcodelength << endl;
  // }
  
  // sort translation table according to validcodelengths from low to high (codes that occur often are found quickly)
   fTranslationTable = TTMergesort(fTranslationTable, TIMEBINS);

  // validation test
  // print merge sorted translation table to screen
  // for(int kk= 0; kk < TIMEBINS; kk++)
  // {
  //     cout << fTranslationTable[kk].famplitude << "  |  " << fTranslationTable[kk].fhuffmancode << "  |  " << fTranslationTable[kk].fvalidcodelength << endl;
  // }

  // do set output data first in order to know the size of the decompressed array!
  // take care of trailers and headers when decoding!

  // initialise reading pointer on compressed data 
  AliHLTUInt64_t* pointerprop = (AliHLTUInt64_t*)fPointer2InData;

  // initialise writing pointer for decompressed output array
  AliHLTUInt64_t* pointeroutputprop = (AliHLTUInt64_t*) fPointer2OutData;

  // outputdata size initialised
  fOutputDataSize = 0;
  
  // validation test
  // cout << "input data size of encoded data (DEC) =  " << fInputDataSize << endl;
  // cout << "output data size of decoded data (DEC) =  " << fOutputDataSize << endl;

  // create temporary word to read out 1 bit from, maximal length: 64 bits */
  AliHLTUInt64_t codedata = 0;
  UInt_t codedatalength = 0;

  // number of 32-bit trailer and common data header words:
  UInt_t headerwords = 8;
  UInt_t header64 = headerwords >> 1;

  // initialise counter variables for input array:
  UInt_t idxwd = header64; // due to 8*32 = 4*64 bit header words
  UInt_t bitpswd = 0;

  // initialise counter variables for output array:
  UInt_t idxoutwd = header64; // due to 8*32 = 4*64 bit header words
  // UInt_t idx10out = 0; // only used for validation test
  UInt_t bitpsoutwd = 0;

  // write header to output (two words input header in one word output header):
  UInt_t head64cntr = 0;
  
  while(head64cntr < header64)
    {
      pointeroutputprop[head64cntr] = pointerprop[head64cntr];

      // validation test
      //cout << "header line in output stream (HEX) =  " << hex << pointeroutputprop[head64cntr] << dec << endl;
    
      ++head64cntr;
    }
  
  // end of encoded ALTRO-blocks [Byte] = finputdatasize - (number of 32-bit-trailer-words + 1) / 4 
  // match this with idxwd as number of 64-bit-words (endNdx[#64-bit-wds] = endNdx[#8-bit-wds]*8/64)
  UInt_t lastfullline = (fInputDataSize) >> 3;
  UInt_t lastbit = (fInputDataSize << 3) & 0x3F;
  UInt_t lastcodeline = lastfullline - ((fNrcuTrailerwords +1) >> 1);

  --lastfullline;
  --lastcodeline;

  if (lastbit == 0)
    {
      lastbit = 64;
    }

  // initialise last valid bit position and trailer array
  AliHLTUInt32_t lastvalidbitpos = 0;
  AliHLTUInt64_t* trailer = new AliHLTUInt64_t[fNrcuTrailerwords];
  if (!trailer) return -1;

  // validation test
  // cout << "last bit (DEC) =  " << lastbit << endl;
  // cout << "lastfullline (DEC) =  " << lastfullline << endl;
  // cout << "lastcodeline(DEC) =  " << lastcodeline << endl;
  // cout << "second last code line data (HEX) = " <<  hex << pointerprop[lastcodeline-1] << dec << endl;
  //  cout << "last code line data (HEX) =  " << hex << pointerprop[lastcodeline] << dec << endl;
  // cout << "last full line data (HEX) =  " << hex << pointerprop[lastfullline] << dec << endl;
 
  // get last valid bit position (first 32 bit word in first trailer line after code
  // and trailerwords
  if (fNrcuTrailerwords < 3) // then lastfullline = only trailerline
    {
      lastvalidbitpos = (pointerprop[lastfullline] >> 32);

      // first trailerword
      trailer[0] = ((pointerprop[lastfullline] << 32 ) >> 32);

      // second trailerword
      if(fNrcuTrailerwords == 2)
	{
	  trailer[1] = ((pointerprop[lastfullline+1] << 32) >> 32);
	};
    }
  else
    {
      lastvalidbitpos = (pointerprop[lastfullline-1] >> 32);

      // first trailerword
      trailer[0] = ((pointerprop[lastfullline-1] << 32 ) >> 32);

      //second trailerword
      trailer[1] = (pointerprop[lastfullline] >> 32);

      // third trailerword
      trailer[2] =((pointerprop[lastfullline] << 32 ) >> 32);
    }

  // validation test
  // cout << "last valid bit position (DEC) =  " << lastvalidbitpos << endl;

  // warning if one of the trailer words is zero:
  for (UInt_t kk = 0; kk < fNrcuTrailerwords; kk++)
    {
      if(trailer[kk] == 0)
	{
	  HLTWarning("Warning! Trailer word %i is zero",kk+1);
	};
    }

  // validation test
  // print trailer array to screen
  // for(Int_t ii=0; ii < fNrcuTrailerwords; ii++)
  // {
  // cout << "trailerword " << ii+1 << " after getting whole trailer (HEX) =  " << hex << trailer[ii] << dec << endl;
  //  }

  // find minmal validcodelength from first entry of translation table and
  // read in minimal validcodelength bits if codedata is zero...
  UInt_t minimalcodelength = fTranslationTable[0].fvalidcodelength;

  // validation test
  // cout << "minimal codelength (DEC) =  " << minimalcodelength << endl;

  Int_t readnewcode = 1; // read new code = 1; go on with old code = 0

  // start decoding
  // NEW  (bitpswd >= lastvalidbitpos) instead of  (bitpswd > lastvalidbitpos)
  while (!((idxwd >= lastcodeline) && (bitpswd >= lastvalidbitpos)) )
    {
      if((idxwd == lastcodeline+1) && (bitpswd == 0)) break; // abortion when lastvalidbitpos = 64

      if (codedatalength < 64) // if codedata can still get one further bit
	{
	  // if new code word is read before end position
	  if((readnewcode == 1) && !((idxwd >= lastcodeline) && (bitpswd + minimalcodelength > lastvalidbitpos)))
	    {
	      // codedata gets next bits from encoded input file:
	      codedata <<= minimalcodelength; //shift bits left
	      
	      if (bitpswd + minimalcodelength < 65) // all bits in the same line
		{
		  codedata |= ((pointerprop[idxwd] << bitpswd)) >> (64 - minimalcodelength); //append bits from input file to the right
		}
	      else // some of this bits in the next line
		{
		  // append bits of current line to the right end of codedata
		  codedata |= (((pointerprop[idxwd] << bitpswd)) >> (bitpswd)) << (minimalcodelength - 64 + bitpswd);
		  
		  // append bits of next line to the right end of codedata
		  codedata |= (pointerprop[idxwd+1] >> (128 - minimalcodelength - bitpswd)); // 128 - mcl - bitpswd = 64 - (mcl - (64 - bitpswd))
		}
	      
	      codedatalength += minimalcodelength;
	      
	      // new code is read, set readnewcode back to 0 again
	      readnewcode = 0;
	      
	      // go and get next input bit
	      bitpswd += minimalcodelength;
	    }
	  else
	    {
	      // codedata gets on bit after another from encoded input file:
	      codedata <<= 1; //shift one bit left
	      codedata |= ((pointerprop[idxwd] << bitpswd)) >> 63; //append next bit from input file to the right
	      
	      ++codedatalength;
	      
	      // go and get next input bit
	      ++bitpswd;
	    }
	  
	  // go and get next input bit
	  if(bitpswd > 63) 
	    {
	      ++idxwd;
	      bitpswd = bitpswd & 0x3F;
	    };
	  
	  // compare if new codedata is in translation table:
	  for(UInt_t kk = 0; kk < TIMEBINS; kk++)
	    {
	      // stopping when current codedatalength smaller than lookup validcodelength (i.e. current code not in table, read next bit)
	      if(fTranslationTable[kk].fvalidcodelength > codedatalength)
		{
		  break;
		};

	      if(fTranslationTable[kk].fhuffmancode == codedata) //lookup bit pattern
		{
		  if(fTranslationTable[kk].fvalidcodelength == codedatalength) //lookup correct codelength
		    {
		      // validation test
		      // if( idxoutwd >= 2636)
		      //{
		      //  cout << "write 10 bit word to decoded output (DEC) =  " << fTranslationTable[kk].famplitude << endl;
		      //  cout << "current idxwd (DEC) =  " << idxwd << endl;
		      //  cout << "current idxoutwd (DEC) =  " << idxoutwd << endl;
		      //}
		      if(bitpsoutwd + 10 < 65) //fits in old line
			{
			  // valdidation test for end of decoded data
			  // print all relevant quantities to screen
			  //if(idxwd > 2635)
			  // {
			  //   cout << "current input in codedata (HEX) =  " << hex << (pointerprop[idxwd] << bitpswd) << dec << endl;
			  //  cout << "current idxwd (DEC) =  " << idxwd << endl;
			  //  cout << "current bitpswd (DEC) =  " << bitpswd << endl;
			  //  cout << "current idxoutwd (DEC) =  " << idxoutwd << endl;
			  //  cout << "current bitpsoutwd (DEC) =  " << bitpsoutwd << endl;
			  //  cout << "decoding value (DEC) =  " << codedata << endl;
			  //  cout << "value length (DEC) =  " << codedatalength << endl;
			  //  cout << "10 bit value (HEX) =  " << hex << fTranslationTable[kk].famplitude << dec << endl;
			  // }
			  
			  pointeroutputprop[idxoutwd] |= ((AliHLTUInt64_t)(fTranslationTable[kk].famplitude)) << bitpsoutwd;
			  
			  // validation test
			  // if (idxoutwd > 362880)
			  // {
			  //   cout << "updated output line (HEX) =  " <<hex << pointeroutputprop[idxoutwd] <<dec << endl;
			  // }
			  
			  // goto next 10-bits output position
			  bitpsoutwd += 10;
			  
			  if(bitpsoutwd == 64) {bitpsoutwd = 0; ++idxoutwd;};
			  
			}
		      else //needs start of new line
			{
			  
			  pointeroutputprop[idxoutwd] |= ((AliHLTUInt64_t)(fTranslationTable[kk].famplitude)) << (bitpsoutwd); 
			  
			  ++idxoutwd; //start new line
			  
			  // validation test
			  // if(idxwd > 2635)
			  //{
			  //  cout << "current input in codedata (HEX) =  " << hex << (pointerprop[idxwd] << bitpswd) << dec << endl;
			  //  cout << "current idxwd (DEC) =  " << idxwd << endl;
			  //  cout << "current bitpswd (DEC) =  " << bitpswd << endl;
			  //  cout << "current idxoutwd (DEC) =  " << idxoutwd << endl;
			  //  cout << "current bitpsoutwd (DEC) =  " << bitpsoutwd << endl;
			  //  cout << "decoding value (DEC) =  " << codedata << endl;
			  //  cout << "value length (DEC) =  " << codedatalength << endl;
			  //  cout << "10 bit value (HEX) =  " << hex << fTranslationTable[kk].famplitude << dec << endl;
			  //}
			  
			  // write next line
			  AliHLTUInt64_t buffer = fTranslationTable[kk].famplitude;	  
			  buffer >>= (64 - bitpsoutwd);
			  
			  pointeroutputprop[idxoutwd] = buffer;
		      
			  // validation test
			  // if (idxoutwd > 362880)
			  // {
			  //   cout << "new line (HEX) =  " << hex << pointeroutputprop[idxoutwd] << dec << endl;
			  // }

			  // go to next bit position
			  bitpsoutwd -= 54;
			}
		      
		      // set buffers to zero again
		      codedata = 0;
		      codedatalength = 0;
		      
		      // prepare codedata to read new code word
		      readnewcode = 1;

		      // validation test
		      // ++idx10out;

		      // go on with next code bits from input
		      break;
		    };
		};
	    }	
	}
      else // if *morecode is used
	{

	  HLTError("Error! Valid codelength for current decoding is larger than 64 bits, error in recognising the correct code");

	  delete [] trailer;
	  return 1;

	}
    }

  // last (incomplete) line from input encoded data
  // after all code is read, bitpsoutwd at lastvalidbitpos
  // -> go to next byte (i.e. start of trailer)

  // validation test
  // cout << "bit position after decoding (DEC) =  " << bitpsoutwd << endl;
  // cout << "line position after decoding (DEC) =  " << idxoutwd << endl;

  // warning if byte position is not 8-bits aligned
  if ((bitpsoutwd & 0x07) != 0 )
    {
      HLTWarning("Warning! Bit position after decoding %i is not aligned to 8 bits", bitpsoutwd);
    };

  if(bitpsoutwd + 32 < 65) // trailer fits in current line, append on the right
    {
      pointeroutputprop[idxoutwd] <<= 32;

      pointeroutputprop[idxoutwd] |= trailer[0]; 

      bitpsoutwd += 32;
    }
  else
    {
      pointeroutputprop[idxoutwd] <<= 64 - bitpsoutwd;
      pointeroutputprop[idxoutwd]|= (trailer[0] >> (bitpsoutwd - 32));

      // go to next line
      ++idxoutwd;

      // get rid of upper, already written bits
      trailer[0] <<= 96 - bitpsoutwd;
      trailer[0] >>= 96 - bitpsoutwd;

      pointeroutputprop[idxoutwd] = trailer[0]; // write lower bits to next line

      bitpsoutwd -= 32;
    }

  if(fNrcuTrailerwords > 1)
    {

      if(bitpsoutwd == 64)
	{
	  bitpsoutwd = 0;
	  ++idxoutwd;
	}

      for(UInt_t ii = 1; ii < fNrcuTrailerwords; ii++)
	{
	  // write second trailer to output data
	  if(bitpsoutwd + 32 < 65) // trailer fits in current line, append on the right
	    {
	      pointeroutputprop[idxoutwd] <<= 32;
	      
	      pointeroutputprop[idxoutwd] |= trailer[ii]; 

	      bitpsoutwd += 32;

	      if(bitpsoutwd == 64)
		{
		  bitpsoutwd = 0;
		  ++idxoutwd;
		}
	    }
	  else
	    {
	      pointeroutputprop[idxoutwd] <<= 64 - bitpsoutwd;
	      pointeroutputprop[idxoutwd]|= (trailer[ii] >> (bitpsoutwd - 32));
	      
	      // go to next line
	      ++idxoutwd;
	      
	      // get rid of upper, already written bits
	      trailer[ii] <<= 96 - bitpsoutwd;
	      trailer[ii] >>= 96 - bitpsoutwd;
	      
	      pointeroutputprop[idxoutwd] = trailer[ii]; // write lower bits to next line
	      
	      bitpsoutwd -= 32;
	    }
	}
    }
      // validation test
  // cout << "bitpswd after decoding (DEC) =  " << bitpswd << endl;
  // cout << "bitpsoutwd after decoding (DEC) =  " << bitpsoutwd << endl;
  // cout << "idxwd after decoding (DEC) =  " << idxwd << endl;
  // cout << "idxoutwd after decoding (DEC) =  " << idxoutwd << endl;

  // write trailer to decoded output

  // set fOutputDataSize in byte:
  // idxoutwd*8 = number of whole 64 bits words in byte
  // bitpsoutwd/8 = number of wohle 8 bits words within last idxoutwd
  // (bitpsoutwd%8)?1:0 = determine whether one more byte has to be used for "left over" bits
  fOutputDataSize = (Int_t) ((idxoutwd << 3) + (bitpsoutwd >> 3) + ((bitpsoutwd & 0x7) ? 1 : 0) );
  
  delete [] trailer;
  // error and abort when output data size smaller than input data size (impossible when decompressing)
  if(fOutputDataSize < fInputDataSize)
    {
      HLTError("Error! Data size for decompressed output stream (bytes) %i is smaller than compressed input data size (bytes)  %i", fOutputDataSize, fInputDataSize);

      return 1;
    };

  // validation test
  // cout << "output data size (DEC) =  " << fOutputDataSize << endl;
  
  return 0;
}

Int_t AliHLTCOMPHuffmanAltro::CalcEntropy(const AliHLTCOMPHuffmanOccurrenceData::AliHLTCOMPHuffmanDataStruct* occurrencetable)
{
  // calculate entropy of 10 bit values
  if (!occurrencetable) occurrencetable=fTrainingTable;
  if (!occurrencetable) return -ENODATA;

  int totalnumber=accumulate(occurrencetable, occurrencetable+TIMEBINS, int(0), AliHLTCOMPHuffmanOccurrenceSum());
  if (totalnumber<1) {
    HLTWarning("can not calculate entropy from empty table");
    fEntropy=0.0;
    return 0;
  }

  double entropy=0.0;  
  const double l2 = log(2.0);
  for(UInt_t jj = 0; jj < TIMEBINS; jj++) {
    double value=occurrencetable[jj].fabundance;
    if (value<1.0) continue;
    entropy += (- (Double_t) value / (Double_t) totalnumber ) * log( ( (Double_t) value / (Double_t) totalnumber )) / (l2);
  }
 
  fEntropy=entropy;
  return 0;
}

/* CopyData is just a function for performance testing */
Int_t AliHLTCOMPHuffmanAltro::CopyData()
{
  // see heaeder file for class documentation
  // output data: fPointer2OutData (to 8-bit words)
  // initialise propagation pointer for output array
  AliHLTUInt32_t* pointer2InData = (AliHLTUInt32_t*) fPointer2InData;
  AliHLTUInt32_t* pointeroutprop = (AliHLTUInt32_t*) fPointer2OutData;

  //first output word has to be initialised
  pointeroutprop[0] = 0;

  fOutputDataSize = fInputDataSize;

  // initialise required variables for input reading:
  AliHLTUInt32_t data10 = 0; // temporary number to read out 10bit word with
  UInt_t idxwd = 8;
  UInt_t idx10 = 0;
  UInt_t bitpswd = 0;

  // number of 32-bit trailer and common data header words:
  // (taken out of compression!)
  // UInt_t headerwords = 8;

  // initialise required variables for ouput creation:
  UInt_t idxoutwd = 0;
  //  Int_t bitcounter = 0; // additional test for number of read in 10 bit words
  // Int_t sumlength = 0;   // additional test for length of encoded file
  
  UInt_t endNdx = (fInputDataSize>>2)-fNrcuTrailerwords;

  // start reading and encoding input data (which are 10 bit words):  
  while (idxwd < endNdx)
    {
      // write input data in temp and cut out 10 bits with mask 0x3FF
      data10 = (pointer2InData[idxwd] >> bitpswd ) & 0x3FFU;
      
      // validation test
      //cout << "current 10 bits word (DEC) =  " << data10 << endl;
	  
	  // if number of bits left in read input less than ten, get remaining bits from next 32bit word
	  if(bitpswd > 21)
	    {
	      data10 |= (pointer2InData[idxwd + 1]) << (32 - bitpswd);
	      data10 &= 0x3FFU;
	    };

	  // validation test
	  // print first 10 bits words to screen
	  // if (idx10 < 100)
	  // {
	  //   cout << "read data:  " << hex << data10 << dec << endl;
	  // };

	  // write data to 32 bit output:
	  pointeroutprop[idxoutwd] |= (data10 << bitpswd);

	  if(bitpswd > 21)
	    {
	      pointeroutprop[idxoutwd + 1] = (data10 >> (32 - bitpswd));
	    }

	  // validation test
	  // cout << "next bitpswd (DEC) =  " << bitpswd << endl;

	  ++idx10; // go on reading 10 bit word
	  
      	  // index word reads position of 32bit word in input array
	  idxoutwd = (idx10 * 10);
	  bitpswd = idxoutwd & 0x1F;
	  idxoutwd >>= 5;
	  idxwd = idxoutwd+8;
	  // bit position word reads position in 32bit word when 10 bits are read out  
    }

  // validation test
  // cout << "output data size (DEC) =  " << fOutputDataSize <<  endl;
  // cout << "sum of codelengths (DEC) =  " << sumlength << endl;
  // cout << "number of 10 bits words (DEC) =  " << bitcounter << endl;
   
  return 0;
}

void AliHLTCOMPHuffmanAltro::Print(Option_t* option) const
{
  // print content

  if (strcmp(option, "trainingtable")==0) {
    // print fTrainingTable to screen (non-zero entries only):
    if (!fTrainingTable) {
      cout << "no training table" << endl;
      return;
    }
    for(Int_t jj = 0; jj < TIMEBINS; jj++) {
      if(fTrainingTable[jj].fabundance != 0)
	cout << jj << "  |  " << fTrainingTable[jj].famplitude << "  |  " << fTrainingTable[jj].fabundance << endl;
    }
    return;
  }

  if (strcmp(option, "translationtable")==0) {
    // print fTranslationTable to screen (non-zero entries only):
    if (!fTranslationTable) {
      cout << "no translation table" << endl;
      return;
    }
    for(Int_t jj = 0; jj < TIMEBINS; jj++) {
      if(fTranslationTable[jj].fvalidcodelength != 0)
	cout << jj << "  |  " << fTranslationTable[jj].fhuffmancode << "  |  " << fTranslationTable[jj].fvalidcodelength << endl;
    }
    return;
  }

}
