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

//____________________________________________________________________
//                                                                          
// Buffer to read RAW ALTRO FMD format from a AliRawReader 
// 
// This class derives from AliAltroBuffer, but overloads the memer
// function Next to do some extra processing.  In particular, it tries
// to autodetect the sample rate.  If zero-suppression was used when
// writing the raw data, then the automatic discovery will not work,
// and the sample rate should be set explicitly. 
//
#include "AliFMDRawStream.h"		// ALIFMDRAWSTREAM_H
#include <AliRawReader.h>		// ALIRAWREADER_H
#include <AliLog.h>
#include <iomanip>
#include <iostream>

//____________________________________________________________________
ClassImp(AliFMDRawStream)
#if 0
  ; // This is here to keep Emacs for indenting the next line
#endif

//____________________________________________________________________
AliFMDRawStream::AliFMDRawStream(AliRawReader* reader, UShort_t sampleRate) 
  : AliAltroRawStream(reader), 
    fSampleRate(sampleRate),
    fPrevTime(-1), 
    fExplicitSampleRate(kFALSE), 
    fPos(0),
    fCur(0),
    fRead(0)
{
  if (fSampleRate > 0) fExplicitSampleRate = kTRUE;
}

//_____________________________________________________________________________
Int_t
AliFMDRawStream::ReadTrailer(UInt_t& addr, UInt_t& len)
{
  if (fPos <= 0) return 0;
  if (fPos <  4) {
    AliError("could not read trailer");
    return -1;
  }
  AliDebug(1, Form("Reading a trailer at %d", fPos));
  Int_t temp = Get10BitWord();
  if (temp != 0x2AA) {
    AliError(Form("Incorrect trailer! Expected 0x2AA but got %x!",temp));
    return -1;
  }
  temp = Get10BitWord();
  if ((temp >> 6) != 0xA) {
    AliError(Form("Incorrect trailer! Expected 0xA but got %x!",temp >> 6));
    return -1;
  }

  len  =  (temp << 4) & 0x3FF;
  temp =  Get10BitWord();
  len  |= (temp >> 6);
  if (((temp >> 2) & 0xF) != 0xA) {
    AliError(Form("Incorrect trailer! Expected 0xA but got %x!",temp >> 6));
    return -1;
  }
  addr =  (temp & 0x3) << 10;
  temp =  Get10BitWord();
  addr |= temp;

  return 4;
}

//_____________________________________________________________________________
Int_t 
AliFMDRawStream::ReadFillWords(UInt_t len)
{
  if (len % 4 == 0) return 0;
  Int_t nFill = (4 - (len % 4)) % 4;
  AliDebug(1, Form("Reading %d fill words", nFill));
  for (Int_t i = 0; i < nFill; i++) {
    UInt_t fill = Get10BitWord();
    if (fill != 0x2AA) {
      AliError(Form("Invalid fill! Expected 0x2AA, but got %X!", fill));
      return -1;
    }
  }
  return nFill;
}

//_____________________________________________________________________________
Int_t 
AliFMDRawStream::ReadBunch(UShort_t* data)
{
  AliDebug(1, "Reading a bunch");
  if (fPos <= 0) {
    AliError("could not read bunch length");
    return -1;
  }
  UShort_t len  = Get10BitWord();
  if (fPos <= 0) {
    AliError("could not read bunch length");
    return -1;
  }
  UShort_t time = Get10BitWord();
  
  AliDebug(1, Form("Bunch is %d long and ends at t=%d", len, time));
  for (UInt_t i = 2; i < len; i++) {
    Int_t amp = Get10BitWord();
    if (amp < 0) { 
      AliError(Form("Bad adc value (%X) !", amp));
      return -1;
    }
    data[time - (i-2)] = amp;
  }
  return len;
}

//_____________________________________________________________________________
Int_t 
AliFMDRawStream::ReadIntoBuffer()
{
  if (fPos > 0) return kTRUE;
  do {
    AliDebug(1, Form("Reading into the buffer"));
    if (!fRawReader->ReadNextData(fRead)) return -1;
  } while (fRawReader->GetDataSize() == 0);
  fPos = (fRawReader->GetDataSize() * 8) / 10;
  // Skip trailing `0x2AA's - is this needed?  Won't it break the
  // trailer? 
#if 0
  UShort_t skip;
  while ((skip = Get10BitWord()) != 0x2AA) 
    AliDebug(1,Form("Skipping one %x", skip));
#endif
  fPos++;
  return fPos;
}

//_____________________________________________________________________________
Bool_t 
AliFMDRawStream::ReadChannel(UInt_t& addr, UInt_t& len, UShort_t* data)
{
  Int_t ret = 0;
  AliDebug(1, "Reading a channel");
  if ((ret = ReadIntoBuffer())       < 0) return kFALSE;
  if ((ret = ReadTrailer(addr, len)) < 0) return kFALSE;
  if ((ret = ReadFillWords(len))     < 0) return kFALSE;
  Int_t toRead = len;
  while (toRead > 0) {
    if ((ret = ReadBunch(data)) < 0) return kFALSE;
    toRead -= ret;
  }
  len -= 2;
  return kTRUE;
}

//_____________________________________________________________________________
Bool_t 
AliFMDRawStream::DumpData()
{
  Int_t ret;
  if ((ret = ReadIntoBuffer())       < 0) return kFALSE;
  UShort_t data;
  Int_t i = 0;
  while ((data = Get10BitWord()) != 0xffff) {
    if (i % 4 == 0) {
      if (i != 0) std::cout << "\n";
      std::cout << std::setw(6) << i << ":";
    }
    std::cout << "  0x" << std::setfill('0') << std::setw(3) 
	      << std::hex << data << std::dec << std::setfill(' ')
	      << std::flush;
    i++;
  }
  return kTRUE;
}

//_____________________________________________________________________________
UShort_t 
AliFMDRawStream::Get10BitWord()
{
  // return a word in a 10 bit array as an UShort_t
  --fPos;
  if (fPos < 0) { 
    AliWarning("At high water mark");
    return 0xFFFF;
  }
  Int_t iBit  = fPos * 10;
  Int_t iByte = iBit / 8;
  Int_t shift = iBit % 8;
  // return ((buffer[iByte+1] * 256 + buffer[iByte]) >> shift) & 0x03FF;

  // recalculate the byte numbers and the shift because
  // the raw data is written as integers where the high bits are filled first
  // -> little endian is assumed here !
  Int_t iByteHigh = 4 * (iByte / 4) + 3 - (iByte % 4);
  iByte++;
  Int_t iByteLow  = 4 * (iByte / 4) + 3 - (iByte % 4);
  shift = 6 - shift;
  return ((fRead[iByteHigh] * 256 + fRead[iByteLow]) >> shift) & 0x03FF;
}

//_____________________________________________________________________________
Bool_t 
AliFMDRawStream::Next()
{
  // read the next raw digit
  // returns kFALSE if there is no digit left
  fPrevTime = fTime;
  if (AliAltroRawStream::Next()) {
    if (!fExplicitSampleRate && fPrevPad != fPad) 
      fSampleRate = fTimeBunch / 128;
    return kTRUE;
  }
  return kFALSE;
}

//_____________________________________________________________________________
//
// EOF
//
