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

//////////////////////////////////////////////////////////////////////
//                                                                          
// Buffer to read RAW ALTRO FMD format from a AliRawReader 
// 
// Currently, I had to overload the Next member function and introduce
// my own members fMyData, fMyPosition, fMyCount, and fMyBunchLength.
// The reason is, that we can use the fMyCount to determine the
// sampling rate used in the ALTRO.   However,
// AliAltroBuffer::fCount is private, so we can not access it. 
//
// If it wasn't I'd implement Next as 
//
//     Bool_t
//     AliFMDRawStreamer::Next()
//     {
//        if (AliAltroRawStreamer::Next()) {
//          if (fPrevPad != fPad) 
//             fSampleRate = (fCount - 2) / 128;
//          return kTRUE;
//        }
//        return kFALSE;
//     } 
//
//////////////////////////////////////////////////////////////////////
#ifndef ALIFMDRAWSTREAM_H
# include "AliFMDRawStream.h"
#endif
#ifndef ALIRAWREADER_H
# include "AliRawReader.h"
#endif
#ifndef __IOSTREAM__
# include <iostream>
#endif

//____________________________________________________________________
ClassImp(AliFMDRawStream);

//____________________________________________________________________
AliFMDRawStream::AliFMDRawStream(AliRawReader* reader) 
  : AliAltroRawStream(reader), 
    // fMyData(0), 
    // fMyPosition(0), 
    // fMyCount(0),
    //    fMyBunchLength(0), 
    fPrevTime(-1)
{}



//_____________________________________________________________________________
Bool_t 
AliFMDRawStream::Next()
{
  // read the next raw digit
  // returns kFALSE if there is no digit left
  fPrevTime = fTime;
  if (AliAltroRawStream::Next()) {
    if (fPrevPad != fPad) {
      fSampleRate = fTimeBunch / 128;
#if 1
      std::cout << "Setting Sample rate to (" << fTimeBunch << "/" 
		<< 128 << "=" << fSampleRate << std::endl;
#endif
    }
    return kTRUE;
  }
  return kFALSE;
}

#if 0
//_____________________________________________________________________________
Bool_t 
AliFMDRawStream::Next()
{
  // read the next raw digit
  // returns kFALSE if there is no digit left
  fPrevSector = fSector;
  fPrevRow    = fRow;
  fPrevPad    = fPad;
  fPrevTime   = fTime;
  
  while (fMyCount == 0) {  // next trailer
    if (fMyPosition <= 0) {  // next payload
      do {
        if (!fRawReader->ReadNextData(fMyData)) return kFALSE;
      } while (fRawReader->GetDataSize() == 0);

      fMyPosition = (fRawReader->GetDataSize() * 8) / 10;
      while (Get10BitWord(fMyData, fMyPosition-1) == 0x2AA) fMyPosition--;
    }

    if (fMyPosition > 0) {
      // read the trailer
      if (fMyPosition <= 4) {
        Error("Next", "could not read trailer");
        return kFALSE;
      }
      fSector   = Get10BitWord(fMyData, --fMyPosition);
      fRow      = Get10BitWord(fMyData, --fMyPosition);
      fPad      = Get10BitWord(fMyData, --fMyPosition);
      fMyCount  = Get10BitWord(fMyData, --fMyPosition);

      fMyPosition -= (4 - (fMyCount % 4)) % 4;  // skip fill words
      fMyBunchLength = 0;

      // Set the sample rate, based on the number of samples in the
      // channel. 
      fSampleRate = (fMyCount - 2) / 128;
#if 0
      std::cout << "Setting Sample rate to (" << fMyCount << " - 2)/" 
		<< 128 << "=" << fSampleRate << std::endl;
#endif
    }
  }

  if (fMyBunchLength == 0) {
    if (fMyPosition <= 0) {
      Error("Next", "could not read bunch length");
      return kFALSE;
    }
    fMyBunchLength = Get10BitWord(fMyData, --fMyPosition) - 2;
    fMyCount--;
    

    if (fMyPosition <= 0) {
      Error("Next", "could not read time bin");
      return kFALSE;
    }
    fTime = Get10BitWord(fMyData, --fMyPosition);
    fMyCount--;
  } else {
    fTime--;
  }

  if (fMyPosition <= 0) {
    Error("Next", "could not read sample amplitude");
    return kFALSE;
  }
  fSignal = Get10BitWord(fMyData, --fMyPosition);
  fMyCount--;
  fMyBunchLength--;

  return kTRUE;
}

//_____________________________________________________________________________
UShort_t 
AliFMDRawStream::Get10BitWord(UChar_t* buffer, Int_t position) const
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
#endif 
//_____________________________________________________________________________
//
// EOF
//
