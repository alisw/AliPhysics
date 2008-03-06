// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTAltroEncoder.cxx
    @author Matthias Richter
    @date   
    @brief  Encoder class for 10/40bit Altro Data format
*/

#include <cassert>
#include <cerrno>
#include "AliHLTAltroEncoder.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAltroEncoder)

AliHLTAltroEncoder::AliHLTAltroEncoder()
  :
  fpBuffer(NULL),
  fBufferSize(0),
  fPrevTimebin(AliHLTUInt16MAX),
  fBunchLength(0),
  fChannelStart(0),
  fChannel(AliHLTUInt16MAX),
  fChannels(),
  fOffset(0),
  f10bitWords(0),
  fOrder(kUnknownOrder)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTAltroEncoder::AliHLTAltroEncoder(AliHLTUInt8_t* pBuffer, int iSize)
  :
  fpBuffer(pBuffer),
  fBufferSize(iSize),
  fPrevTimebin(AliHLTUInt16MAX),
  fBunchLength(0),
  fChannelStart(0),
  fChannel(AliHLTUInt16MAX),
  fChannels(),
  fOffset(0),
  f10bitWords(0),
  fOrder(kUnknownOrder)
{
  // see header file for class documentation
}

AliHLTAltroEncoder::~AliHLTAltroEncoder()
{
  // see header file for class documentation
}

int AliHLTAltroEncoder::SetBuffer(AliHLTUInt8_t* pBuffer, int iSize)
{
  // see header file for class documentation
  fpBuffer=pBuffer;
  fBufferSize=iSize;

  return 0;
}

int AliHLTAltroEncoder::AddSignal(AliHLTUInt16_t signal, AliHLTUInt16_t timebin)
{
  // see header file for class documentation
  int iResult=0;
  if (fPrevTimebin!=AliHLTUInt16MAX) {
    assert(fPrevTimebin!=timebin);
    if (fOrder==kUnknownOrder) {
      if (fPrevTimebin+1==timebin) fOrder=kAscending;
      else if (fPrevTimebin==timebin+1) fOrder=kDescending;
    }
    if ((fOrder!=kAscending || fPrevTimebin+1!=timebin) &&
	(fOrder!=kDescending || fPrevTimebin!=timebin+1)) {
      // Finalize bunch and start new one
      iResult=SetBunch();
    }
  }

  if (iResult>=0 && (iResult=Add10BitValue(signal))>=0) {
    fBunchLength++;
  }
  assert(fOffset*4<=f10bitWords*5);
  fPrevTimebin=timebin;
  return iResult;
}

int AliHLTAltroEncoder::SetChannel(AliHLTUInt16_t hwaddress)
{
  // see header file for class documentation
  int iResult=0;
  int added10BitWords=0;
  if (!fpBuffer) return -ENODEV;
  if (fOffset+5>=fBufferSize) {
    HLTWarning("buffer too small too finalize channel: %d of %d byte(s) already used", fOffset, fBufferSize);
    return -ENOSPC;
  }

  if (iResult>=0 && 
      (iResult=SetBunch())>=0) {
    AliHLTUInt16_t length=f10bitWords-fChannelStart;
    if ((iResult=Pad40Bit())<0) return iResult;
    // 2 words for the SetBunch (end time and length) and the
    // padded words to fill 40bit word
    added10BitWords=iResult+2;
    assert((length+iResult)%4==0);
    //HLTInfo("%d %x", hwaddress, hwaddress);
    fpBuffer[fOffset++]=hwaddress&0xff;
    fpBuffer[fOffset++]=0xa0 | (hwaddress>>8)&0xf;
    fpBuffer[fOffset++]=length&0xff;
    fpBuffer[fOffset++]=0xa8 | (length>>8)&0x3;
    fpBuffer[fOffset++]=0xaa;
    f10bitWords+=4;
    fChannelStart=f10bitWords;
    fChannels.push_back(hwaddress);
    fPrevTimebin=AliHLTUInt16MAX;
  }
  if (iResult<0) return iResult;
  return added10BitWords;
}

int AliHLTAltroEncoder::AddChannelSignal(AliHLTUInt16_t signal, AliHLTUInt16_t timebin, AliHLTUInt16_t hwaddress)
{
  // see header file for class documentation
  int iResult=0;
  int added10BitWords=0;
  if (fChannel==AliHLTUInt16MAX) {
    fChannel=hwaddress;
  } else if (fChannel!=hwaddress) {
    iResult=SetChannel(fChannel);
    added10BitWords=iResult;
    fChannel=hwaddress;
  }

  if (iResult>=0) {
    if ((iResult=AddSignal(signal, timebin))>=0)
      added10BitWords++;
  }

  if (iResult<0) return iResult;
  return added10BitWords;
}

int AliHLTAltroEncoder::GetTotal40bitWords()
{
  // see header file for class documentation
  if (fChannelStart!=f10bitWords) {
    HLTWarning("unterminated channel found, check calling sequence");
  }
  assert(fChannelStart%4==0);

  return fChannelStart;
}

int AliHLTAltroEncoder::SetBunch()
{
  // see header file for class documentation
  int iResult=0;

  // return if the bunch has already been set
  if (fBunchLength==0) return 0;

  // fill time bin and bunch length
  if ((iResult=Add10BitValue(fPrevTimebin))>=0) {
    iResult=Add10BitValue(fBunchLength+2);
    fBunchLength=0;
    iResult=2;
  }
  return iResult;
}

int AliHLTAltroEncoder::Add10BitValue(AliHLTUInt16_t value)
{
  // see header file for class documentation
  int iResult=0;
  if (!fpBuffer) return -ENODEV;
  if (fOffset+2>=fBufferSize) {
    HLTWarning("buffer too small too add 10bit word: %d of %d byte(s) already used", fOffset, fBufferSize);
    return -ENOSPC;
  }

  int bit=(f10bitWords%4)*10;
  int shift=bit%8;
  unsigned short maskLow=~((0xff<<shift)>>8);
  //unsigned short maskHigh=~((0xff<<((bit+10)%8))>>8);
  fpBuffer[fOffset++]|=maskLow&(value<<shift);
  fpBuffer[fOffset]=(value&0x3ff)>>(8-shift);
  f10bitWords++;
  if (f10bitWords%4==0) fOffset++;

  return iResult;
}

int AliHLTAltroEncoder::Pad40Bit()
{
  // see header file for class documentation
  int iResult=0;
  int added10BitWords=0;
  while (iResult>=0 && f10bitWords%4!=0) {
    if ((iResult=Add10BitValue(0x2aa))>=0) {
      added10BitWords++;
    }
  }
  if (iResult<0) return iResult;
  return added10BitWords;
}
