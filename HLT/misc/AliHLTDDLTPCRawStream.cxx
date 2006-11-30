// @(#) $Id$

// Author: Constantin Loizides <mailto:loizides@ikf.uni-frankfurt.de>
//*-- Copyright &copy ALICE HLT Group

#include "AliHLTRootTypes.h"
#include "AliHLTStandardIncludes.h"
#include "AliHLTLogging.h"
#include "AliHLTDDLRawReader.h"

#include "AliHLTDDLTPCRawStream.h"
//#include "AliTPCHuffman.h"


/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** \class AliHLTDDLTPCRawReaderStream
<pre>
//_____________________________________________________________
// AliHLTDDLTPCRawReaderStream (taken from the offline AliROOT code,
// original authors: D.Favretto and A.K.Mohanty)
//
// This is a base class for reading TPC raw data 
// and providing information about digits
</pre>
*/

#if __GNUC__ >= 3
using namespace std;
#endif

ClassImp(AliHLTDDLTPCRawStream)

AliHLTDDLTPCRawStream::AliHLTDDLTPCRawStream(AliHLTDDLRawReader* rawReader)
{
  // create an object to read TPC raw digits

  fRawReader = rawReader;
  fRawReader->Select(3);
  fData = new UShort_t[fgkDataMax];
  fDataSize = fPosition = 0;
  fCount = fBunchLength = 0;

  fSector = fPrevSector = fRow = -1;
  fPrevRow = fPad = fPrevPad = -1;
  fTime = fSignal = -1;
}

AliHLTDDLTPCRawStream::~AliHLTDDLTPCRawStream()
{
  // clean up
  delete[] fData;
}

Bool_t AliHLTDDLTPCRawStream::SetDDLID(Int_t d)
{
  // sets DDL ID
  if((d<0)||(d>216)){
    LOG(AliHLTLog::kFatal,"AliHLTDDLTPCRawStream::SetDDLID","DDL")
      <<AliHLTLog::kDec<<"DDL number out of range "<<d<<ENDLOG;
    return kFALSE;
  }

  //partial clean
  fDataSize = fPosition = 0;
  fCount = fBunchLength = 0;

  fSector = fPrevSector = fRow = fPrevRow = fPad = fPrevPad = fTime = fSignal = -1;

  fRawReader->Reset();
  fRawReader->Select(3,d,d+1);
	
  return kTRUE;
}

Bool_t AliHLTDDLTPCRawStream::Next()
{
  // read the next raw digit
  // returns kFALSE if there is no digit left

  fPrevSector = fSector;
  fPrevRow = fRow;
  fPrevPad = fPad;

  while (fCount == 0) {  // next trailer
    if (fPosition >= fDataSize) {  // next payload
      UChar_t* data;
      do {
	if (!fRawReader->ReadNextData(data)) return kFALSE;
      } while (fRawReader->GetDataSize() == 0);

      if (fRawReader->IsCompressed()) {  // compressed data
 	LOG(AliHLTLog::kFatal,"AliHLTDDLTPCRawStream::Next","Compression")
      <<"Compression is not implemented (yet)!"<<ENDLOG;
	return kFALSE;
      } else {                           // uncompressed data
	fDataSize = 0;
	Int_t pos = (fRawReader->GetDataSize() * 8) / 10;
	while (Get10BitWord(data, pos-1) == 0x2AA) pos--;
	while (pos > 0) {
	  for (Int_t i = 0; i < 4; i++) {  // copy trailer
	    fData[fDataSize++] = Get10BitWord(data, pos-4+i);
	  }
	  pos -= 4;
	  Int_t count = fData[fDataSize-4];
	  pos -= (4 - (count % 4)) % 4;  // skip fill words

	  while (count > 0) {
	    UShort_t bunchLength = Get10BitWord(data, pos-1);
	    fData[fDataSize++] = bunchLength;
	    fData[fDataSize++] = Get10BitWord(data, pos-2);  // time bin

	    // copy signal amplitudes in increasing order on time
	    for (Int_t i = 0; i < bunchLength-2; i++) {
	      fData[fDataSize++] = Get10BitWord(data, pos - bunchLength + i);
	    }
	    pos -= bunchLength;
	    count -= bunchLength;
	  }
	}
      }

      fPosition = 0;
    }
    if (fPosition + 4 >= fDataSize) {
      LOG(AliHLTLog::kError,"AliHLTDDLTPCRawStream::Next","Data")
      <<"Could not read trailer"<<ENDLOG;
      return kFALSE;
    }
    fCount = fData[fPosition++];
    fPad = fData[fPosition++];
    fRow = fData[fPosition++];
    fSector = fData[fPosition++];
    fBunchLength = 0;
  }

  if (fBunchLength == 0) {
    if (fPosition >= fDataSize) {
      LOG(AliHLTLog::kError,"AliHLTDDLTPCRawStream::Next","Data")
      <<"Could not read bunch length"<<ENDLOG;
      return kFALSE;
    }
    fBunchLength = fData[fPosition++] - 2;
    fCount--;

    if (fPosition >= fDataSize) {
      LOG(AliHLTLog::kError,"AliHLTDDLTPCRawStream::Next","Data")
      <<"Could not read time bin"<<ENDLOG;
      return kFALSE;
    }
    fTime = fData[fPosition++] - fBunchLength;
    fCount--;
  }

  fTime++;
  if (fPosition >= fDataSize) {
    LOG(AliHLTLog::kError,"AliHLTDDLTPCRawStream::Next","Data")
      <<"Could not read sample amplitude"<<ENDLOG;
    return kFALSE;
  }

  fSignal = fData[fPosition++] + fgkOffset;
  fCount--;
  fBunchLength--;

  return kTRUE;
}

UShort_t AliHLTDDLTPCRawStream::Get10BitWord(UChar_t* buffer, Int_t position) const
{
  // return a word in a 10 bit array as an UShort_t
  Int_t iBit = position * 10;
  Int_t iByte = iBit / 8;
  Int_t shift = iBit % 8;

  // recalculate the byte numbers and the shift because
  // the raw data is written as integers where the high bits are filled first
  // -> little endian is assumed here !
  Int_t iByteHigh = 4 * (iByte / 4) + 3 - (iByte % 4);
  iByte++;
  Int_t iByteLow  = 4 * (iByte / 4) + 3 - (iByte % 4);
  shift = 6 - shift;
  return ((buffer[iByteHigh] * 256 + buffer[iByteLow]) >> shift) & 0x03FF;
}
