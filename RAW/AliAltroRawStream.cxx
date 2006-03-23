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
#include "AliLog.h"

ClassImp(AliAltroRawStream)


//_____________________________________________________________________________
AliAltroRawStream::AliAltroRawStream(AliRawReader* rawReader) :
  fSector(-1),
  fPrevSector(-1),
  fRow(-1),
  fPrevRow(-1),
  fPad(-1),
  fPrevPad(-1),
  fHWAddress(-1),
  fPrevHWAddress(-1),
  fTime(-1),
  fSignal(-1),
  fRawReader(rawReader),
  fData(NULL),
  fNoAltroMapping(kTRUE),
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
  fHWAddress(-1),
  fPrevHWAddress(-1),
  fTime(-1),
  fSignal(-1),
  fRawReader(NULL),
  fData(NULL),
  fNoAltroMapping(kTRUE),
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
void AliAltroRawStream::Reset()
{
// reset altro raw stream params

  fPosition = fCount = fBunchLength = 0;

  fSector = fPrevSector = fRow = fPrevRow = fPad = fPrevPad = fHWAddress = fPrevHWAddress = fTime = fSignal = -1;

  if (fRawReader) fRawReader->Reset();
}

//_____________________________________________________________________________
Bool_t AliAltroRawStream::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  fPrevSector = fSector;
  fPrevRow = fRow;
  fPrevPad = fPad;
  fPrevHWAddress = fHWAddress;

  while (fCount == 0) {  // next trailer
    if (fPosition <= 0) {  // next payload
      do {
	if (!fRawReader->ReadNextData(fData)) return kFALSE;
      } while (fRawReader->GetDataSize() == 0);

      fPosition = GetPosition();
    }

    if (!ReadTrailer())
      AliFatal("Incorrect trailer information !");

    fBunchLength = 0;
  }

  if (fBunchLength == 0) ReadBunch();
  else fTime--;

  ReadAmplitude();

  return kTRUE;
}

//_____________________________________________________________________________
UShort_t AliAltroRawStream::GetNextWord()
{
  // Read the next 10 bit word in backward direction
  // The input stream access is given by fData and fPosition

  fPosition--;

  Int_t iBit = fPosition * 10;
  Int_t iByte = iBit / 8;
  Int_t shift = iBit % 8;

  // recalculate the byte numbers and the shift because
  // the raw data is written as integers where the high bits are filled first
  // -> little endian is assumed here !
  Int_t iByteHigh = 4 * (iByte / 4) + 3 - (iByte % 4);
  iByte++;
  Int_t iByteLow  = 4 * (iByte / 4) + 3 - (iByte % 4);
  shift = 6 - shift;
  return ((fData[iByteHigh] * 256 + fData[iByteLow]) >> shift) & 0x03FF;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStream::ReadTrailer()
{
  //Read a trailer of 40 bits in the backward reading mode
  //In case of no mapping is provided, read a dummy trailer
  if (fNoAltroMapping) {
    AliError("No ALTRO mapping information is loaded! Reading a dummy trailer!");
    return ReadDummyTrailer();
  }

  //First reading filling words
  UShort_t temp;
  Int_t nFillWords = 0;
  while ((temp = GetNextWord()) == 0x2AA) nFillWords++;
  if (nFillWords == 0)
    AliFatal("Incorrect trailer found ! Expected 0x2AA not found !");

  //Then read the trailer
  if (fPosition <= 4)
    AliFatal(Form("Incorrect raw data size ! Expected at lest 4 words but found %d !",fPosition));

  fCount = (temp << 4) & 0x3FF;
  if ((temp >> 6) != 0xA)
    AliFatal(Form("Incorrect trailer found ! Expecting 0xA but found %x !",temp >> 6));

  temp = GetNextWord();
  fHWAddress = (temp & 0x3) << 10;
  if (((temp >> 2) & 0xF) != 0xA)
    AliFatal(Form("Incorrect trailer found ! Expecting second 0xA but found %x !",(temp >> 2) & 0xF));
  fCount |= ((temp & 0x3FF) >> 6);
  if (fCount == 0) return kFALSE;

  temp = GetNextWord();
  fHWAddress |= temp;

  fPosition -= (4 - (fCount % 4)) % 4;  // skip fill words

  ApplyAltroMapping();

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStream::ReadDummyTrailer()
{
  //Read a trailer of 40 bits in the backward reading mode
  //In case of no mapping is provided, read a dummy trailer
  UShort_t temp;
  while ((temp = GetNextWord()) == 0x2AA);

  fSector = temp;
  fRow = GetNextWord();
  fPad = GetNextWord();
  fCount = GetNextWord();
  if (fCount == 0) return kFALSE;
  fHWAddress = -1;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAltroRawStream::ReadBunch()
{
  // Read altro payload in 
  // backward direction
  if (fPosition <= 0)
    AliFatal("Could not read bunch length !");

  fBunchLength = GetNextWord() - 2;
  fTimeBunch = fBunchLength;
  fCount--;

  if (fPosition <= 0)
    AliFatal("Could not read time bin !");

  fTime = GetNextWord();
  fCount--;

  return;
}

//_____________________________________________________________________________
void AliAltroRawStream::ReadAmplitude()
{
  // Read next time bin amplitude
  if (fPosition <= 0)
    AliFatal("Could not read sample amplitude !");

  fSignal = GetNextWord();
  fCount--;
  fBunchLength--;

  return;
}

//_____________________________________________________________________________
Int_t AliAltroRawStream::GetPosition()
{
  // Sets the position in the
  // input stream
  Int_t position = (fRawReader->GetDataSize() * 8) / 10;
  if (position <= 4)
    AliFatal(Form("Incorrect raw data size ! Expected at lest 4 words but found %d !",position));

  return position;
}
