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

#include "AliAltroRawStreamOld.h"
#include "AliRawReader.h"
#include "AliLog.h"

ClassImp(AliAltroRawStreamOld)


//_____________________________________________________________________________
AliAltroRawStreamOld::AliAltroRawStreamOld(AliRawReader* rawReader) :
  fNoAltroMapping(kTRUE),
  fDDLNumber(-1),
  fPrevDDLNumber(-1),
  fHWAddress(-1),
  fPrevHWAddress(-1),
  fTime(-1),
  fPrevTime(-1),
  fSignal(-1),
  fTimeBunch(-1),
  fRawReader(rawReader),
  fData(NULL),
  fPosition(0),
  fCount(0),
  fBunchLength(0)
{
// create an object to read Altro raw digits
  fSegmentation[0] = fSegmentation[1] = fSegmentation[2] = -1;
}

//_____________________________________________________________________________
AliAltroRawStreamOld::AliAltroRawStreamOld(const AliAltroRawStreamOld& stream) :
  TObject(stream),
  fNoAltroMapping(kTRUE),
  fDDLNumber(-1),
  fPrevDDLNumber(-1),
  fHWAddress(-1),
  fPrevHWAddress(-1),
  fTime(-1),
  fPrevTime(-1),
  fSignal(-1),
  fTimeBunch(-1),
  fRawReader(NULL),
  fData(NULL),
  fPosition(0),
  fCount(0),
  fBunchLength(0)
{
  Fatal("AliAltroRawStreamOld", "copy constructor not implemented");
}

//_____________________________________________________________________________
AliAltroRawStreamOld& AliAltroRawStreamOld::operator = (const AliAltroRawStreamOld& 
					      /* stream */)
{
  Fatal("operator =", "assignment operator not implemented");
  return *this;
}

//_____________________________________________________________________________
AliAltroRawStreamOld::~AliAltroRawStreamOld()
{
// clean up

}

//_____________________________________________________________________________
void AliAltroRawStreamOld::Reset()
{
// reset altro raw stream params

  fPosition = fCount = fBunchLength = 0;

  fDDLNumber = fPrevDDLNumber = fHWAddress = fPrevHWAddress = fTime = fPrevTime = fSignal = fTimeBunch = -1;

  if (fRawReader) fRawReader->Reset();

  fSegmentation[0] = fSegmentation[1] = fSegmentation[2] = -1;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStreamOld::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  fPrevDDLNumber = fDDLNumber;
  fPrevHWAddress = fHWAddress;
  fPrevTime = fTime;

  while (fCount == 0) {  // next trailer
    if (fPosition <= 0) {  // next payload
      do {
	if (!fRawReader->ReadNextData(fData)) return kFALSE;
      } while (fRawReader->GetDataSize() == 0);

      fDDLNumber = fRawReader->GetDDLID();

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
void AliAltroRawStreamOld::SelectRawData(Int_t detId)
{
  // Select the raw data for specific
  // detector id
  AliDebug(1,Form("Selecting raw data for detector %d",detId));
  fRawReader->Select(detId);
}

//_____________________________________________________________________________
UShort_t AliAltroRawStreamOld::GetNextWord()
{
  // Read the next 10 bit word in backward direction
  // The input stream access is given by fData and fPosition

  fPosition--;

  Int_t iBit = fPosition * 10;
  Int_t iByte = iBit / 8;
  Int_t shift = iBit % 8;

  // the raw data is written as integers where the low bits are filled first
  // -> little endian is assumed here !
  Int_t iByteLow = iByte;
  iByte++;
  Int_t iByteHigh  = iByte;
  return ((fData[iByteHigh] * 256 + fData[iByteLow]) >> shift) & 0x03FF;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStreamOld::ReadTrailer()
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

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStreamOld::ReadDummyTrailer()
{
  //Read a trailer of 40 bits in the backward reading mode
  //In case of no mapping is provided, read a dummy trailer
  UShort_t temp;
  while ((temp = GetNextWord()) == 0x2AA);

  fSegmentation[0] = temp;
  fSegmentation[1] = GetNextWord();
  fSegmentation[2] = GetNextWord();
  fCount = GetNextWord();
  if (fCount == 0) return kFALSE;
  fHWAddress = -1;

  fPosition -= (4 - (fCount % 4)) % 4;  // skip fill words

  return kTRUE;
}

//_____________________________________________________________________________
void AliAltroRawStreamOld::ReadBunch()
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
void AliAltroRawStreamOld::ReadAmplitude()
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
Int_t AliAltroRawStreamOld::GetPosition()
{
  // Sets the position in the
  // input stream

  // Get the payload size from the RCU trailer
  // The trailer is actually one 32-bit word
  Int_t position = (fData[fRawReader->GetDataSize()-1]) << 24;
  position |= (fData[fRawReader->GetDataSize()-2]) << 16;
  position |= (fData[fRawReader->GetDataSize()-3]) << 8;
  position |= (fData[fRawReader->GetDataSize()-4]);
  // The size is specified in a number of 40bits
  // Therefore we need to transform it to number of bytes
  position *= 5;

  // Check the consistency of the header and trailer
  if ((fRawReader->GetDataSize() - 32) != position)
    AliFatal(Form("Inconsistent raw data size ! Expected %d bytes (from the header), found %d words (in the RCU trailer)!",
		  fRawReader->GetDataSize()-32,
		  position));

  // Skip the Common Data Header which contains
  // only 7 (!) words
  fData += 28;

  // Return the position in units of 10-bit words
  return position*8/10;
}
