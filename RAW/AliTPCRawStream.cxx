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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
///
/// This class provides access to TPC digits in raw data.
///
/// It loops over all TPC digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliTPCRawStream.h"
#include "AliTPCHNode.h"
#include "AliRawReader.h"

ClassImp(AliTPCRawStream)


AliTPCHNode** AliTPCRawStream::fgRootNode = NULL;


AliTPCRawStream::AliTPCRawStream(AliRawReader* rawReader)
{
// create an object to read TPC raw digits

  fRawReader = rawReader;
  fRawReader->Select(0);
  fData = new UShort_t[fgkDataMax];
  fDataSize = fPosition = 0;
  fCount = fBunchLength = 0;

  if (!fgRootNode) {
    fgRootNode = new AliTPCHNode*[fgkNumTables];
    fCompression.CreateTreesFromFile(fgRootNode, fgkNumTables);
  }

  fSector = fPrevSector = fRow = fPrevRow = fPad = fPrevPad = fTime = fSignal = -1;
}

AliTPCRawStream::AliTPCRawStream(const AliTPCRawStream& stream) :
  TObject(stream)
{
  Fatal("AliTPCRawStream", "copy constructor not implemented");
}

AliTPCRawStream& AliTPCRawStream::operator = (const AliTPCRawStream& 
					      /* stream */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

AliTPCRawStream::~AliTPCRawStream()
{
// clean up

  delete[] fData;
}


Bool_t AliTPCRawStream::Next()
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
	UInt_t size = 0;
	fCompression.Decompress(fgRootNode, fgkNumTables, 
				(char*) data, fRawReader->GetDataSize(),
				fData, size);
	fDataSize = size;

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
      Error("Next", "could not read trailer");
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
      Error("Next", "could not read bunch length");
      return kFALSE;
    }
    fBunchLength = fData[fPosition++] - 2;
    fCount--;

    if (fPosition >= fDataSize) {
      Error("Next", "could not read time bin");
      return kFALSE;
    }
    fTime = fData[fPosition++] - fBunchLength;
    fCount--;
  }

  fTime++;
  if (fPosition >= fDataSize) {
    Error("Next", "could not read sample amplitude");
    return kFALSE;
  }
  fSignal = fData[fPosition++] + fgkOffset;
  fCount--;
  fBunchLength--;

  return kTRUE;
}


UShort_t AliTPCRawStream::Get10BitWord(UChar_t* buffer, Int_t position) const
{
// return a word in a 10 bit array as an UShort_t

  Int_t iBit = position * 10;
  Int_t iByte = iBit / 8;
  Int_t shift = iBit % 8;
//  return ((buffer[iByte+1] * 256 + buffer[iByte]) >> shift) & 0x03FF;

  // recalculate the byte numbers and the shift because
  // the raw data is written as integers where the high bits are filled first
  // -> little endian is assumed here !
  Int_t iByteHigh = 4 * (iByte / 4) + 3 - (iByte % 4);
  iByte++;
  Int_t iByteLow  = 4 * (iByte / 4) + 3 - (iByte % 4);
  shift = 6 - shift;
  return ((buffer[iByteHigh] * 256 + buffer[iByteLow]) >> shift) & 0x03FF;
}
