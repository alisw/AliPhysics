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
  fNoAltroMapping(kTRUE),
  fIsOldRCUFormat(kFALSE),
  fIsShortDataHeader(kFALSE),
  fDDLNumber(-1),
  fPrevDDLNumber(-1),
  fRCUId(-1),
  fPrevRCUId(-1),
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
  fBunchLength(0),
  fRCUTrailerData(NULL),
  fRCUTrailerSize(0)
{
// create an object to read Altro raw digits
  fSegmentation[0] = fSegmentation[1] = fSegmentation[2] = -1;
}

//_____________________________________________________________________________
AliAltroRawStream::AliAltroRawStream(const AliAltroRawStream& stream) :
  TObject(stream),
  fNoAltroMapping(kTRUE),
  fIsOldRCUFormat(kFALSE),
  fIsShortDataHeader(kFALSE),
  fDDLNumber(-1),
  fPrevDDLNumber(-1),
  fRCUId(-1),
  fPrevRCUId(-1),
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
  fBunchLength(0),
  fRCUTrailerData(NULL),
  fRCUTrailerSize(0)
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

  fRCUTrailerData = NULL;
  fRCUTrailerSize = 0;

  fDDLNumber = fPrevDDLNumber = fRCUId = fPrevRCUId = fHWAddress = fPrevHWAddress = fTime = fPrevTime = fSignal = fTimeBunch = -1;

  if (fRawReader) fRawReader->Reset();

  fSegmentation[0] = fSegmentation[1] = fSegmentation[2] = -1;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStream::Next()
{
// read the next raw digit
// returns kFALSE if there is no digit left

  fPrevDDLNumber = fDDLNumber;
  fPrevRCUId = fRCUId;
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
void AliAltroRawStream::SelectRawData(Int_t detId)
{
  // Select the raw data for specific
  // detector id
  AliDebug(1,Form("Selecting raw data for detector %d",detId));
  fRawReader->Select(detId);
}

//_____________________________________________________________________________
void AliAltroRawStream::SelectRawData(const char *detName)
{
  // Select the raw data for specific
  // detector name
  AliDebug(1,Form("Selecting raw data for detector %s",detName));
  fRawReader->Select(detName);
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

  // the raw data is written as integers where the low bits are filled first
  // -> little endian is assumed here !
  Int_t iByteLow = iByte;
  iByte++;
  Int_t iByteHigh  = iByte;
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
  if (nFillWords == 0) {
    PrintDebug();
    AliFatal("Incorrect trailer found ! Expected 0x2AA not found !");
  }

  //Then read the trailer
  if (fPosition <= 4) {
    PrintDebug();
    AliFatal(Form("Incorrect raw data size ! Expected at lest 4 words but found %d !",fPosition));
  }

  fCount = (temp << 4) & 0x3FF;
  if ((temp >> 6) != 0xA) {
    PrintDebug();
    AliFatal(Form("Incorrect trailer found ! Expecting 0xA but found %x !",temp >> 6));
  }

  temp = GetNextWord();
  fHWAddress = (temp & 0x3) << 10;
  if (((temp >> 2) & 0xF) != 0xA) {
    PrintDebug();
    AliFatal(Form("Incorrect trailer found ! Expecting second 0xA but found %x !",(temp >> 2) & 0xF));
  }
  fCount |= ((temp & 0x3FF) >> 6);
  if (fCount == 0) return kFALSE;

  temp = GetNextWord();
  fHWAddress |= temp;

  fPosition -= (4 - (fCount % 4)) % 4;  // skip fill words

  return kTRUE;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStream::ReadDummyTrailer()
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
void AliAltroRawStream::ReadBunch()
{
  // Read altro payload in 
  // backward direction
  if (fPosition <= 0) {
    PrintDebug();
    AliFatal("Could not read bunch length !");
  }

  fBunchLength = GetNextWord() - 2;
  fTimeBunch = fBunchLength;
  fCount--;

  if (fPosition <= 0) {
    PrintDebug();
    AliFatal("Could not read time bin !");
  }

  fTime = GetNextWord();
  fCount--;

  return;
}

//_____________________________________________________________________________
void AliAltroRawStream::ReadAmplitude()
{
  // Read next time bin amplitude
  if (fPosition <= 0) {
    PrintDebug();
    AliFatal("Could not read sample amplitude !");
  }

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
  // Read the RCU trailer
  // This includes the trailer size,
  // RCU identifier and raw data payload.
  // The RCU trailer format is described
  // in details in the RCU manual.

  if (!fIsOldRCUFormat) {
    // First read 32-bit word with the
    // trailer size (22 bits) and RCU ID (the rest)
    Int_t index = fRawReader->GetDataSize();
    UInt_t word = Get32bitWord(index);
    fRCUId = (Int_t)(word >> 22);
    Int_t trailerSize = (word & 0x3FFFFF);

    // Now read the beginning of the trailer
    // where the payload size is written
    if (trailerSize < 2) {
      PrintDebug();
      AliFatal(Form("Invalid trailer size found (%d bytes) !",trailerSize*4));
    }
    fRCUTrailerSize = (trailerSize-2)*4;
    index -= fRCUTrailerSize;
    if (index < 4) {
      PrintDebug();
      AliFatal(Form("Invalid trailer size found (%d bytes) ! The size is bigger than the raw data size (%d bytes)!",
		    trailerSize*4,
		    fRawReader->GetDataSize()));
    }
    fRCUTrailerData = fData + index;
    Int_t position = Get32bitWord(index);
    // The size is specified in a number of 40bits
    // Therefore we need to transform it to number of bytes
    position *= 5;

    // Check the consistency of the header and trailer
    if ((fRawReader->GetDataSize() - trailerSize*4) != position) {
      PrintDebug();
      AliFatal(Form("Inconsistent raw data size ! Raw data size - %d bytes (from the header), RCU trailer - %d bytes, raw data paylod - %d bytes !",
		    fRawReader->GetDataSize(),
		    trailerSize*4,
		    position));
    }

    return position * 8 / 10;
  }
  else {
    // In case of the Old RCU trailer format
    // we have to read just the size of altro payload
    // in units of 40-bit words
    Int_t index = fRawReader->GetDataSize();
    Int_t position = Get32bitWord(index);

    fRCUId = -1;
    fRCUTrailerSize = 0;
    fRCUTrailerData = NULL;

    // The size is specified in a number of 40bits
    // Therefore we need to transform it to number of bytes
    position *= 5;

    if (!fIsShortDataHeader) {

      // Check the consistency of the header and trailer
      if ((fRawReader->GetDataSize() - 4) != position) {
	PrintDebug();
	AliFatal(Form("Inconsistent raw data size ! Expected %d bytes (from the header), found %d bytes (in the RCU trailer)!",
		      fRawReader->GetDataSize()-4,
		      position));
      }
    }
    else {
      // Check the consistency of the header and trailer
      // In this case the header is shorter by 4 bytes
      if (fRawReader->GetDataSize() != position) {
	PrintDebug();
	AliFatal(Form("Inconsistent raw data size ! Expected %d bytes (from the header), found %d bytes (in the RCU trailer)!",
		      fRawReader->GetDataSize(),
		      position));
      }

      // 7 32-bit words Common Data Header
      // therefore we have to shift back by 4 bytes
      // the pointer to the raw data payload
      fData -= 4;
    }
     
    // Return the position in units of 10-bit words
    return position*8/10;
  }
}

//_____________________________________________________________________________
UInt_t AliAltroRawStream::Get32bitWord(Int_t &index)
{
  // This method returns the 32 bit word at a given
  // position inside the raw data payload.
  // The 'index' points to the beginning of the next word.
  // The method is supposed to be endian (platform)
  // independent.
  if (!fData) {
    PrintDebug();
    AliFatal("Raw data paylod buffer is not yet initialized !");
  }

  if (index < 4) {
    PrintDebug();
    AliFatal(Form("Invalid raw data payload index (%d) !",index));
  }

  UInt_t word = 0;
  word  = fData[--index] << 24;
  word |= fData[--index] << 16;
  word |= fData[--index] << 8;
  word |= fData[--index];

  return word;
}

//_____________________________________________________________________________
Bool_t AliAltroRawStream::GetRCUTrailerData(UChar_t*& data) const
{
  // Return a pointer to the RCU trailer
  // data. Should be called always after
  // the RCU trailer was already processed
  // in the GetPosition() method
  if (!fRCUTrailerSize || !fRCUTrailerData) {
    AliError("No valid RCU trailer data is found !");
    data = NULL;
    return kFALSE;
  }

  data = fRCUTrailerData;

  return kTRUE;
}

//_____________________________________________________________________________
void AliAltroRawStream::PrintDebug() const
{
  // The method prints all the available
  // debug information.
  // Its is used in case of decoding errors.

  AliError("Start of debug printout\n--------------------");

  Dump();
  if (fRawReader) fRawReader->Dump();

  AliError("End of debug printout\n--------------------");
}

//_____________________________________________________________________________
Int_t AliAltroRawStream::GetBranch() const
{
  // The method provides the RCU branch index (0 or 1)
  // for the current hardware address.
  // In case the hardware address has not been yet
  // initialized, the method returns -1
  if (fHWAddress == -1) return -1;

  return ((fHWAddress >> 11) & 0x1);
}

//_____________________________________________________________________________
Int_t AliAltroRawStream::GetFEC() const
{
  // The method provides the front-end card index
  // for the current hardware address.
  // In case the hardware address has not been yet
  // initialized, the method returns -1
  if (fHWAddress == -1) return -1;

  return ((fHWAddress >> 7) & 0xF);
}

//_____________________________________________________________________________
Int_t AliAltroRawStream::GetAltro() const
{
  // The method provides the altro chip index
  // for the current hardware address.
  // In case the hardware address has not been yet
  // initialized, the method returns -1
  if (fHWAddress == -1) return -1;

  return ((fHWAddress >> 4) & 0x7);
}

//_____________________________________________________________________________
Int_t AliAltroRawStream::GetChannel() const
{
  // The method provides the channel index
  // for the current hardware address.
  // In case the hardware address has not been yet
  // initialized, the method returns -1
  if (fHWAddress == -1) return -1;

  return (fHWAddress & 0xF);
}
