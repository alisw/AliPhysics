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
/// This class provides access to Altro digits in raw data.
///
/// It loops over all Altro digits in the raw data given by the AliRawReader.
/// The Next method goes to the next digit. If there are no digits left
/// it returns kFALSE.
/// Several getters provide information about the current digit.
///
///////////////////////////////////////////////////////////////////////////////

#include "AliAltroRawStream.h"
#include "AliRawReader.h"

ClassImp(AliAltroRawStream)


//_____________________________________________________________________________
AliAltroRawStream::AliAltroRawStream(AliRawReader* rawReader) :
  fSector(-1),
  fPrevSector(-1),
  fRow(-1),
  fPrevRow(-1),
  fPad(-1),
  fPrevPad(-1),
  fTime(-1),
  fSignal(-1),
  fRawReader(rawReader),
  fData(NULL),
  fPosition(0),
  fCount(0),
  fBunchLength(0)
{
// create an object to read Altro raw digits

}

//_____________________________________________________________________________
AliAltroRawStream::AliAltroRawStream(const AliAltroRawStream& stream) :
  TObject(stream),
  fSector(-1),
  fPrevSector(-1),
  fRow(-1),
  fPrevRow(-1),
  fPad(-1),
  fPrevPad(-1),
  fTime(-1),
  fSignal(-1),
  fRawReader(NULL),
  fData(NULL),
  fPosition(0),
  fCount(0),
  fBunchLength(0)
{
  Fatal("AliAltroRawStream", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliAltroRawStream& AliAltroRawStream::operator = (const AliAltroRawStream& 
					      /* stream */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliAltroRawStream::~AliAltroRawStream()
{
// clean up

}


//_____________________________________________________________________________
Bool_t AliAltroRawStream::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  fPrevSector = fSector;
  fPrevRow = fRow;
  fPrevPad = fPad;

  while (fCount == 0) {  // next trailer
    if (fPosition <= 0) {  // next payload
      do {
	if (!fRawReader->ReadNextData(fData)) return kFALSE;
      } while (fRawReader->GetDataSize() == 0);

      fPosition = (fRawReader->GetDataSize() * 8) / 10;
      while (Get10BitWord(fData, fPosition-1) == 0x2AA) fPosition--;
    }

    if (fPosition > 0) {
      // read the trailer
      if (fPosition <= 4) {
	Error("Next", "could not read trailer");
	return kFALSE;
      }
      fSector = Get10BitWord(fData, --fPosition);
      fRow    = Get10BitWord(fData, --fPosition);
      fPad    = Get10BitWord(fData, --fPosition);
      fCount  = Get10BitWord(fData, --fPosition);

      fPosition -= (4 - (fCount % 4)) % 4;  // skip fill words
      fBunchLength = 0;
    }
  }

  if (fBunchLength == 0) {
    if (fPosition <= 0) {
      Error("Next", "could not read bunch length");
      return kFALSE;
    }
    fBunchLength = Get10BitWord(fData, --fPosition) - 2;
    fTimeBunch = fBunchLength;
    fCount--;

    if (fPosition <= 0) {
      Error("Next", "could not read time bin");
      return kFALSE;
    }
    fTime = Get10BitWord(fData, --fPosition);
    fCount--;
  } else {
    fTime--;
  }

  if (fPosition <= 0) {
    Error("Next", "could not read sample amplitude");
    return kFALSE;
  }
  fSignal = Get10BitWord(fData, --fPosition);
  fCount--;
  fBunchLength--;

  return kTRUE;
}


//_____________________________________________________________________________
UShort_t AliAltroRawStream::Get10BitWord(UChar_t* buffer, Int_t position) const
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
