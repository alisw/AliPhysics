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

/** @file   AliHLTAltroGenerator.cxx
    @author Matthias Richter
    @date   
    @brief  Simulation class of 10/40bit Altro Data.
*/

#include <cassert>
#include <cerrno>
#include "AliHLTAltroGenerator.h"
#include "TArrayS.h"
#include "TArrayC.h"
#include "TRandom.h"
#include "TDatime.h"
#include "AliRawDataHeader.h"
#include "AliHLTAltroEncoder.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAltroGenerator)

AliHLTAltroGenerator::AliHLTAltroGenerator(int maxChannels,
					   int maxBunches,
					   int maxBunchLength,
					   int maxTimebin,
					   int maxSignal)
  :
  fpData(NULL),
  fpSimData(NULL),
  fChannelPositions(),
  fNof10BitWords(0),
  fpCDH(NULL),
  fCDHSize(0),
  fpTrailer(NULL),
  fTrailerSize(0),
  fMaxChannels(maxChannels),
  fMaxBunches(maxBunches),
  fMaxBunchLength(maxBunchLength),
  fMaxTimebin(maxTimebin),
  fMaxSignal(maxSignal),
  fpRand(NULL),
  fDirection(kBackwards),
  fCurrentPosition(-1),
  fCurrentBunch(-1),
  fCurrentTimeOffset(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTAltroGenerator::~AliHLTAltroGenerator()
{
  // see header file for class documentation
  if (fpTrailer) delete[] fpTrailer;
  if (fpCDH) delete[] fpCDH;
  if (fpSimData) delete fpSimData;
  if (fpData) delete fpData;
}

int AliHLTAltroGenerator::Generate()
{
  // see header file for class documentation
  int iResult=0;

  if (!fpSimData) fpSimData=new TArrayS;
  if (!fpSimData) {
    return -ENOMEM;
  }

  Reset();

  int nofChannels=GetRandom(1, fMaxChannels);
  if (nofChannels==0) nofChannels=1;

  HLTDebug("number of channels: %d", nofChannels);
  int channelAddress=-1;
  int lastChannel=-1;
  int dataPos=0;
  int repetitions=0;
  for (int channel=0; channel<nofChannels; channel++) {
    channelAddress=GetRandom(0, fMaxChannels);
    //HLTDebug("channel %d: address %d, %d bunch(es)", channel, channelAddress, nofBunches);
    if (channelAddress==lastChannel) {
      channel--;
      if (repetitions++>5) break;
      continue;
    }
    int nofBunches=GetRandom(1, fMaxBunches);
    if (nofBunches==0) {
      channel--;
      if (repetitions++>5) break;
      continue;
    }
    repetitions=0;
    int totalBunches=0;
    int bunchEndTime=0;
    int lastBunchEndTime=0;

    HLTDebug("simulate channel %d: address %d", channel, channelAddress);

    // save beginning of this channel for better navigation
    AliChannelPosition position={channelAddress, dataPos, 0};

    // add channel address and bunch count at the beginning
    if (fpSimData->GetSize()<dataPos+2) fpSimData->Set(dataPos+2);
    (*fpSimData)[dataPos++]=channelAddress;
    dataPos++; // placeholder for number of bunches

    int bunch=0;

    // Matthias Oct 2008:
    // Data was simulated in the wrong order: bunches for the higher timebins
    // first. But real data is the other way round: bunches are written in the order
    // of ascending timebins. A new consistency check was added to AliAltroDecoder
    // (AliAltroBunch) in revision 29090. Now the channels are checked for overlapping
    // bunches, the time of a bunch must be smaller than time of the previous one
    // minus its length.
    // Data is now simulated in the right order in order to fullfil this check.
    for (bunch=0; bunch<nofBunches && bunchEndTime<fMaxTimebin-3; bunch++) {
      while ((bunchEndTime+=GetRandom(0, fMaxTimebin-bunchEndTime))-lastBunchEndTime<3) {/*empty body*/};
      int bunchLength=GetRandom(0, bunchEndTime-lastBunchEndTime<fMaxBunchLength?bunchEndTime-lastBunchEndTime:fMaxBunchLength);
      if (bunchLength==0) continue;
      totalBunches++;

      HLTDebug("       bunch %d, length %d, end time %d ", bunch, bunchLength, bunchEndTime);

      if (fpSimData->GetSize()<dataPos+bunchLength+4) fpSimData->Set(dataPos+bunchLength+4);
      // write bunch length and time at both ends
      (*fpSimData)[dataPos++]=bunchLength;
      int time=bunchEndTime-bunchLength+1;
      (*fpSimData)[dataPos++]=time;
      for (; time<=bunchEndTime; time++) {	
	int signal=GetRandom(0, fMaxSignal);
	(*fpSimData)[dataPos++]=signal;
      }
      (*fpSimData)[dataPos++]=bunchEndTime;
      (*fpSimData)[dataPos++]=bunchLength;
      fNof10BitWords+=bunchLength+2;
      lastBunchEndTime=bunchEndTime;
      bunchEndTime+=bunchLength;
    }
    if (totalBunches>0) {
      (*fpSimData)[position.fPosition+1]=totalBunches;
      if (fpSimData->GetSize()<dataPos+2) fpSimData->Set(dataPos+2);
      (*fpSimData)[dataPos++]=totalBunches;
      position.fEnd=dataPos;
      (*fpSimData)[dataPos++]=channelAddress;
      lastChannel=channelAddress;
      fNof10BitWords=(fNof10BitWords+7)/4; fNof10BitWords*=4; // align to 4 and add 4
      fChannelPositions.push_back(position);
      assert((*fpSimData)[position.fPosition]==(*fpSimData)[position.fEnd]);
      HLTDebug("       channel %d added: address %d, %d bunch(es)", channel, channelAddress, totalBunches);
    } else {
      dataPos-=2;
      HLTDebug("       channel %d skipped: address %d, %d bunch(es)", channel, channelAddress, totalBunches);
    }
  }

  assert(fNof10BitWords%4==0);
  if (iResult<0) {
    fpSimData->Set(0);
    return iResult;
  }
  fpSimData->Set(dataPos);
  return GetDataSize();
}

int AliHLTAltroGenerator::GetNof40BitAltroWords() const
{
  // see header file for class documentation
  assert(fNof10BitWords%4==0);
  return fNof10BitWords/4;
}

int AliHLTAltroGenerator::GetDataSize()
{
  // see header file for class documentation
  int iResult=0;
  iResult=(fNof10BitWords*5)/4;
  if (fpTrailer) {
    *(reinterpret_cast<AliHLTUInt32_t*>(fpTrailer))=GetNof40BitAltroWords();
    iResult+=fTrailerSize;
  }
  if (fpCDH) iResult+=fCDHSize;
  return iResult;
}

int AliHLTAltroGenerator::GetData(AliHLTUInt8_t* &pBuffer)
{
  // see header file for class documentation
  int iResult=GetDataSize();
  if (iResult>0) {
    if (!fpData) fpData=new TArrayC(iResult);
    if (fpData) {
      if (fpData->GetSize()<iResult) fpData->Set(iResult);
      if ((iResult=GetData(reinterpret_cast<AliHLTUInt8_t*>(fpData->GetArray()), fpData->GetSize()))>=0) {
	pBuffer=reinterpret_cast<AliHLTUInt8_t*>(fpData->GetArray());
      }
    } else {
      iResult=-ENOMEM;
    }
  }
  return iResult;
}

int AliHLTAltroGenerator::GetData(AliHLTUInt8_t* pBuffer, int size)
{
  // see header file for class documentation
  int iResult=0;
  int dataPos=0;

  if (size<GetDataSize()) return -ENOSPC;

  // copy Common Data Header
  if (fpCDH) {
    fpCDH->fSize=GetDataSize();
    memcpy(pBuffer+dataPos, fpCDH, fCDHSize);
    dataPos+=fCDHSize;
  }

  // encode simulated data
  if ((iResult=EncodeData(pBuffer+dataPos, size-dataPos))>=0) {
    dataPos+=iResult;
  }

  // copy trailer
  if (fpTrailer) {
    memcpy(pBuffer+dataPos, fpTrailer, fTrailerSize);
    AliHLTUInt32_t* pLast=reinterpret_cast<AliHLTUInt32_t*>(fpTrailer+fTrailerSize-sizeof(AliHLTUInt32_t));
    *pLast=GetNof40BitAltroWords();
    dataPos+=fTrailerSize;
  }

  if (iResult<0) return iResult;
  assert(fpCDH==NULL || (int)fpCDH->fSize==dataPos);
  return dataPos;
}

int AliHLTAltroGenerator::SetCDH(AliRawDataHeader* pCDH, int size)
{
  // see header file for class documentation
  int iResult=0;
  if (pCDH && size>0) {
    if (fpCDH) delete[] fpCDH;
    fpCDH=new AliRawDataHeader;
    if (fpCDH) {
      memcpy(fpCDH, pCDH, size);
      fCDHSize=size;
    } else {
      iResult=-ENOMEM;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTAltroGenerator::SetRCUTrailer(AliHLTUInt8_t* pTrailer, int size)
{
  // see header file for class documentation
  int iResult=0;
  if (pTrailer && size>=(int)sizeof(AliHLTUInt32_t)) {
    AliHLTUInt32_t* pLast=reinterpret_cast<AliHLTUInt32_t*>(pTrailer+size-sizeof(AliHLTUInt32_t));
    if (size!=sizeof(AliHLTUInt32_t)) {
      // if more than one trailer words, the last one is the trailer length (# 32bit words)
      if (*pLast!=size/sizeof(AliHLTUInt32_t)) {
	HLTError("invalid trailer: trailer length (last 32bit word) does not match trailer size (bytes)");
	return -EBADF;
      }
    }
    if (fpTrailer) delete[] fpTrailer;
    fpTrailer=new AliHLTUInt8_t[size];
    if (fpTrailer) {
      memcpy(fpTrailer, pTrailer, size);
      fTrailerSize=size;
    } else {
      iResult=-ENOMEM;
    }
  } else {
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTAltroGenerator::GetChannels(vector<AliHLTUInt16_t> list)
{
  // see header file for class documentation
  int iResult=0;
  list.clear();
  for (vector<AliChannelPosition>::iterator element=fChannelPositions.begin();
       element!=fChannelPositions.end();
       element++) {
    list.push_back(element->fChannel);
  }
  iResult=list.size();
  return iResult;
}

int AliHLTAltroGenerator::SetSorting(AliHLTUInt16_t */*array*/, int /*arraySize*/)
{
  // see header file for class documentation
  int iResult=0;
  HLTError("function not yet implemented");
  return iResult;
}

int AliHLTAltroGenerator::EncodeData(AliHLTUInt8_t* pBuffer, int size)
{
  // see header file for class documentation
  int iResult=0;
  if (!pBuffer) return -EINVAL;

  AliHLTAltroEncoder encoder;
  encoder.SetBuffer(pBuffer, size);

  Short_t channelAddress=-1;
  for (vector<AliChannelPosition>::iterator element=fChannelPositions.begin();
       element!=fChannelPositions.end() && iResult>=0;
       element++) {
    if (!fpSimData ||
	fpSimData->GetSize()<=element->fPosition ||
	fNof10BitWords==0) {
      iResult=-ENODATA;
      break;
    }
    channelAddress=element->fChannel;
    assert(fpSimData->At(element->fPosition)==channelAddress);
    if (fpSimData->At(element->fPosition)!=channelAddress) {
      iResult=-ENODATA;
      break;
    }
    int dataPos=element->fPosition+1;
    int nofBunches=fpSimData->At(dataPos++);
    int bunch=0;
    for (; bunch<nofBunches; bunch++) {
      int bunchLength=fpSimData->At(dataPos++);
      int startTime=fpSimData->At(dataPos++);
      int time=startTime;
      for (; time<startTime+bunchLength; time++) {
	iResult=encoder.AddSignal(fpSimData->At(dataPos++), time);
      }
      assert(time-1==fpSimData->At(dataPos));
      dataPos++; // DO NOT PUT INTO ASSERT
      assert(bunchLength==fpSimData->At(dataPos));
      dataPos++; // DO NOT PUT INTO ASSERT
    }

    if (iResult>=0 && channelAddress>=0) {
      assert(nofBunches==fpSimData->At(dataPos));
      dataPos++; // DO NOT PUT INTO ASSERT
      assert(channelAddress==fpSimData->At(dataPos));
      dataPos++; // DO NOT PUT INTO ASSERT
    }
    encoder.SetChannel(channelAddress);
  }

  if (iResult>=0) {
    iResult=(encoder.GetTotal40bitWords()*5)/4;
  }

  return iResult;
}

int AliHLTAltroGenerator::GetRandom(int min, int max)
{
  // see header file for class documentation
  if (max-min<2) return min;
  bool setTheSeed=fpRand!=NULL;
  if (fpRand==NULL) {
    fpRand=new TRandom;
  }
  if (fpRand==NULL) return min;
  if (!setTheSeed) {
    TDatime dt;
    fpRand->SetSeed(dt.Get());
  }
  return fpRand->Integer(max-min);
}

int AliHLTAltroGenerator::Reset()
{
  // see header file for class documentation
  fChannelPositions.clear();
  fNof10BitWords=0;
  Rewind();
  return 0;
}

int AliHLTAltroGenerator::Rewind()
{
  // see header file for class documentation
  fCurrentPosition=-1;
  fCurrentBunch=-1;
  fCurrentTimeOffset=-1;
  return 0;
}

bool AliHLTAltroGenerator::Next()
{
  // see header file for class documentation
  bool haveData=false;
  if (!fpSimData) return false;
  do {
    if (haveData=(fCurrentTimeOffset>=0 &&
		  fCurrentBunch>0 &&
		  ++fCurrentTimeOffset<GetBunchSize())) {
      break;
    }

    if (haveData=(NextBunch() && GetBunchSize()>0)) {
      fCurrentTimeOffset=0;
      break;
    }
  } while (NextChannel());
  return haveData;
}

AliHLTUInt16_t AliHLTAltroGenerator::GetSignal()
{
  // see header file for class documentation
  if (!fpSimData || fCurrentTimeOffset<0) return 0;
  assert(fCurrentPosition>=0 && fCurrentPosition<(int)fChannelPositions.size());
  assert(fCurrentBunch>=0);
  if (fDirection==kForwards) {
    return fpSimData->At(fCurrentBunch+2+fCurrentTimeOffset);
  } else if (fDirection==kBackwards){
    return fpSimData->At(fCurrentBunch-(2+fCurrentTimeOffset));
  }
  return 0;
}

bool AliHLTAltroGenerator::NextChannel()
{
  // see header file for class documentation
  bool haveData=false;
  if (fpSimData && fChannelPositions.size()==0) return false;
  fpSimData->GetArray();
  if (fCurrentPosition==-1) {
    if (fDirection==kForwards) fCurrentPosition=0;
    else fCurrentPosition=fChannelPositions.size()-1;
    haveData=true;
  } else {
    if (fDirection==kForwards && (haveData=(fCurrentPosition+1<(int)fChannelPositions.size()))) {
      fCurrentPosition++;
    } else if (fDirection==kBackwards && (haveData=(fCurrentPosition>0))) {
      fCurrentPosition--;
    }
  }
  
  fCurrentBunch=-1;
  fCurrentTimeOffset=-1;
  return haveData;
}

AliHLTUInt16_t AliHLTAltroGenerator::GetHwAddress()
{
  // see header file for class documentation
  if (fCurrentPosition>=0 && fCurrentPosition<(int)fChannelPositions.size())
    return fChannelPositions[fCurrentPosition].fChannel;
  return ~((AliHLTUInt16_t)0);
}

int AliHLTAltroGenerator::GetBunchCount()
{
  // see header file for class documentation
  if (fpSimData && fCurrentPosition>=0 && fCurrentPosition<(int)fChannelPositions.size()) {
    if (fDirection==kForwards)
      return fpSimData->At(fChannelPositions[fCurrentPosition].fPosition+1);
    else if (fDirection==kBackwards)
      return fpSimData->At(fChannelPositions[fCurrentPosition].fEnd-1);
  }
  return ~((AliHLTUInt16_t)0);
}

bool AliHLTAltroGenerator::NextBunch()
{
  // see header file for class documentation
  bool haveData=false;
  if (fpSimData && fCurrentPosition>=0 && fCurrentPosition<(int)fChannelPositions.size()) {
    if (fDirection==kBackwards) {
      if (fCurrentBunch<0) {
	// bunch count in channel end - 1
	if (haveData=(fpSimData->At(fChannelPositions[fCurrentPosition].fEnd-1))>0) {
	  // first bunch length at channel end - 2
	  fCurrentBunch=fChannelPositions[fCurrentPosition].fEnd-2;
	}
      } else if (fCurrentBunch>fChannelPositions[fCurrentPosition].fPosition+1) {
	fCurrentBunch-=fpSimData->At(fCurrentBunch)+4;
	haveData=fCurrentBunch>fChannelPositions[fCurrentPosition].fPosition+1;
      }
      // cross check
      if (haveData) {
	assert(fpSimData->At(fCurrentBunch)==fpSimData->At(fCurrentBunch-fpSimData->At(fCurrentBunch)-3));
	haveData=fpSimData->At(fCurrentBunch)==fpSimData->At(fCurrentBunch-fpSimData->At(fCurrentBunch)-3);
      }
    } else if (fDirection==kForwards) {
      if (fCurrentBunch<0) {
	// bunch count in channel start + 1
	if (haveData=(fpSimData->At(fChannelPositions[fCurrentPosition].fPosition+1))>0) {
	  // first bunch length at channel start + 2
	  fCurrentBunch=fChannelPositions[fCurrentPosition].fPosition+2;
	}
      } else if (fCurrentBunch<fChannelPositions[fCurrentPosition].fEnd-1) {
	fCurrentBunch+=fpSimData->At(fCurrentBunch)+4;
	haveData=fCurrentBunch<fChannelPositions[fCurrentPosition].fEnd-1;
      }
      // cross check
      if (haveData) {
	assert(fpSimData->At(fCurrentBunch)==fpSimData->At(fCurrentBunch+fpSimData->At(fCurrentBunch)+3));
	haveData=fpSimData->At(fCurrentBunch)==fpSimData->At(fCurrentBunch+fpSimData->At(fCurrentBunch)+3);
      }
    }
  }
  fCurrentTimeOffset=-1;
  return haveData;
}

AliHLTUInt16_t AliHLTAltroGenerator::GetBunchSize()
{
  // see header file for class documentation
  if (fCurrentBunch<0) return 0;
  if (fCurrentBunch<fChannelPositions[fCurrentPosition].fPosition+2) return 0;
  if (fCurrentBunch>fChannelPositions[fCurrentPosition].fEnd-2) return 0;
  return fpSimData->At(fCurrentBunch);
}

AliHLTUInt16_t AliHLTAltroGenerator::GetStartTime()
{
  // see header file for class documentation
  if (!fpSimData || GetBunchSize()==0) return 0;
  AliHLTUInt16_t startTime=0;
  if (fDirection==kForwards) {
    startTime=fpSimData->At(fCurrentBunch+1);
  } else if (fDirection==kBackwards) {
    startTime=fpSimData->At(fCurrentBunch-fpSimData->At(fCurrentBunch)-2);
  }
  if (fCurrentTimeOffset>=0) {
    if (fDirection==kForwards) {
      startTime+=fCurrentTimeOffset;
    } else if (fDirection==kBackwards) {
      startTime=GetEndTime();
    }
  }
  return startTime;
}

AliHLTUInt16_t AliHLTAltroGenerator::GetEndTime()
{
  // see header file for class documentation
  if (!fpSimData || GetBunchSize()==0) return 0;
  AliHLTUInt16_t endTime=0;
  if (fDirection==kForwards) {
    endTime=fpSimData->At(fCurrentBunch+fpSimData->At(fCurrentBunch)+2);
  } else if (fDirection==kBackwards) {
    endTime=fpSimData->At(fCurrentBunch-1);
  }
  if (fCurrentTimeOffset>=0) {
    assert(fCurrentTimeOffset<fpSimData->At(fCurrentBunch));
    if (fDirection==kForwards) {
      endTime=GetStartTime();
    } else if (fDirection==kBackwards) {
      endTime-=fCurrentTimeOffset;
    }
  }
  return endTime;
}

const Short_t* AliHLTAltroGenerator::GetSignals()
{
  // see header file for class documentation
  if (!fpSimData || GetBunchSize()==0) return NULL;
  if (fDirection==kForwards) {
    return fpSimData->GetArray()+fCurrentBunch+2;
  } else if (fDirection==kBackwards) {
    return fpSimData->GetArray()+(fCurrentBunch-fpSimData->At(fCurrentBunch)-1);
  }
  return NULL;
}

void AliHLTAltroGenerator::Print()
{
  // see header file for class documentation
  cout << *this << endl;
}

ostream &operator<<(ostream &stream, AliHLTAltroGenerator &generator)
{
  // see header file for class documentation
  int iResult=0;
  Short_t channelAddress=-1;
  for (vector<AliHLTAltroGenerator::AliChannelPosition>::iterator element=generator.fChannelPositions.begin();
       element!=generator.fChannelPositions.end() && iResult>=0;
       element++) {
    if (!generator.fpSimData ||
	generator.fpSimData->GetSize()<=element->fPosition ||
	generator.fNof10BitWords==0) {
      stream << "AliHLTAltroGenerator: no data available" << endl;;
      break;
    }
    channelAddress=element->fChannel;
    assert(generator.fpSimData->At(element->fPosition)==channelAddress);
    if (generator.fpSimData->At(element->fPosition)!=channelAddress) {
      stream << "AliHLTAltroGenerator: internal data mismatch" << endl;;
      iResult=-ENODATA;
      break;
    }
    int dataPos=element->fPosition+1;
    int nofBunches=generator.fpSimData->At(dataPos++);
    stream << "***************************************************************" << endl;
    stream << "channel address: " << channelAddress << "    " << nofBunches << " bunch(es)" << endl;
    int bunch=0;
    for (; bunch<nofBunches; bunch++) {
      int bunchLength=generator.fpSimData->At(dataPos++);
      int startTime=generator.fpSimData->At(dataPos++);
      int time=startTime;
      stream << "   length " << bunchLength << " start time " << startTime << ":     ";
      for (; time<startTime+bunchLength; time++) {
	stream << " " << generator.fpSimData->At(dataPos++);
      }
      assert(time-1==generator.fpSimData->At(dataPos));
      dataPos++; // DO NOT PUT INTO ASSERT
      assert(bunchLength==generator.fpSimData->At(dataPos));
      dataPos++; // DO NOT PUT INTO ASSERT
      stream << "      -> end time " << time-1 << endl;
    }

    if (iResult>=0 && channelAddress>=0) {
      assert(nofBunches==generator.fpSimData->At(dataPos));
      dataPos++; // DO NOT PUT INTO ASSERT
      assert(channelAddress==generator.fpSimData->At(dataPos));
      dataPos++; // DO NOT PUT INTO ASSERT
    }
  }

  return stream;
}
