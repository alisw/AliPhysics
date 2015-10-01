// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
 *                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTTPCDigitReader.cxx
    @author Timm Steinbeck, Jochen Thaeder, Matthias Richter, Kenneth Aamodt
    @date   
    @brief  An abstract reader class for TPC data.
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTPCDigitReader.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTStdIncludes.h"

using namespace std;

ClassImp(AliHLTTPCDigitReader)

AliHLTTPCDigitReader::AliHLTTPCDigitReader()
  :
  fFlags(0),
  fLckRow(-1),
  fLckPad(-1)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCDigitReader::~AliHLTTPCDigitReader()
{
  // see header file for class documentation
}

int AliHLTTPCDigitReader::InitBlock(void* ptr,unsigned long size,Int_t firstrow,Int_t lastrow, Int_t patch, Int_t slice)
{
  // see header file for class documentation
  if (patch<0 || patch>=AliHLTTPCGeometry::GetNumberOfPatches()) {
    HLTError("invalid readout partition number %d", patch);
    return -EINVAL;
  }
  if (firstrow!=AliHLTTPCGeometry::GetFirstRow(patch)) {
    HLTWarning("The firstrow parameter does not match the layout of the readout partition %d "
	       "(firstrow=%d). Parameter is ignored", patch, AliHLTTPCGeometry::GetFirstRow(patch));
  }
  if (lastrow!=AliHLTTPCGeometry::GetLastRow(patch)) {
    HLTWarning("The lastrow parameter does not match the layout of the readout partition %d "
	       "(lastrow=%d). Parameter is ignored", patch, AliHLTTPCGeometry::GetLastRow(patch));
  }
  return InitBlock(ptr, size, patch, slice);
}

void AliHLTTPCDigitReader::SetOldRCUFormat(Bool_t /*oldrcuformat*/)
{
  // default method of the base class
}

void AliHLTTPCDigitReader::SetUnsorted(Bool_t /*unsorted*/)
{
  // default method of the base class
  HLTWarning("common sorting functionality has not yet been implemented");
}

bool AliHLTTPCDigitReader::Next(int /*type*/)
{
  // see header file for class documentation
  if (!CheckFlag(kLocked)) return NextSignal();

  bool haveData=false;
  if (!CheckFlag(kChannelOverwrap))
    haveData=NextSignal();

  if (haveData && (fLckRow!=GetRow() || fLckPad!=GetPad())) {
    SetFlag(kChannelOverwrap);
    haveData=false;
  }

  return haveData;
}

bool AliHLTTPCDigitReader::NextChannel()
{
  // see header file for class documentation
  PrintWarningOnce(kWarnMissFastAccess,"\n"
		   "      !!! This digit reader does not implement the methods for       !!!\n"
		   "      !!! fast data access on channel/bunch basis. Data is discarded !!!");
  return false;
}

int AliHLTTPCDigitReader::NextBunch()
{
  // see header file for class documentation
  PrintWarningOnce(kWarnMissFastAccess,"\n"
		   "      !!! This digit reader does not implement the methods for       !!!\n"
		   "      !!! fast data access on channel/bunch basis. Data is discarded !!!");
  return false;
}

const UInt_t* AliHLTTPCDigitReader::GetSignals()
{
  // see header file for class documentation
  PrintWarningOnce(kWarnMissFastAccess,"\n"
		   "      !!! This digit reader does not implement the methods for       !!!\n"
		   "      !!! fast data access on channel/bunch basis. Data is discarded !!!");
  return 0;
}

const UShort_t* AliHLTTPCDigitReader::GetSignalsShort()
{
  // see header file for class documentation
  PrintWarningOnce(kWarnMissFastAccess,"\n"
		   "      !!! This digit reader does not implement the methods for       !!!\n"
		   "      !!! fast data access on channel/bunch basis. Data is discarded !!!");
  return 0;
}

void AliHLTTPCDigitReader::EnableCaching(bool bCache)
{
  // see header file for class documentation
  if (bCache) SetFlag(kChannelCaching);
  else ClearFlag(kChannelCaching);
}

int AliHLTTPCDigitReader::RewindChannel()
{
  // see header file for class documentation
  int iResult=0;
  
  return iResult;
}

unsigned int AliHLTTPCDigitReader::SetFlag(unsigned int flag)
{
  // see header file for class documentation
  return fFlags|=flag;
}
	
unsigned int AliHLTTPCDigitReader::ClearFlag(unsigned int flag)
{
  // see header file for class documentation
  return fFlags&=~flag;
}

// int operator[](int timebin)
// {
//   return -1;
// }

int AliHLTTPCDigitReader::RewindCurrentChannel()
{
  // see header file for class documentation
  SetFlag(kNoRewind);
  if (!CheckFlag(kChannelCaching)) return -ENODATA;
  return -ENOSYS;
}

int AliHLTTPCDigitReader::RewindToPrevChannel()
{
  // see header file for class documentation
  SetFlag(kNoRewind);
  if (!CheckFlag(kChannelCaching)) return -ENODATA;
  return -ENOSYS;
}

int AliHLTTPCDigitReader::GetBunchSize()
{
  // see header file for class documentation
  PrintWarningOnce(kWarnMissFastAccess,"\n"
		   "      !!! This digit reader does not implement the methods for       !!!\n"
		   "      !!! fast data access on channel/bunch basis. Data is discarded !!!");
  return 0;
}

int AliHLTTPCDigitReader::GetRowOffset() const
{
  // see header file for class documentation
  return 0;
}

AliHLTUInt32_t AliHLTTPCDigitReader::GetAltroBlockHWaddr() const
{
  // see header file for class documentation
  return 0;
}

AliHLTUInt32_t AliHLTTPCDigitReader::GetAltroBlockHWaddr(Int_t /*row*/, Int_t /*pad*/) const
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCDigitReader::GetRCUTrailerSize()
{
  // see header file for class documentation
  PrintWarningOnce(kWarnMissTrailerGetters,"\n"
		   "      !!! This digit reader does not implement the Getters       !!!\n"
		   "      !!! for RCU trailer. Ignoring call.                        !!!");
  return 0;
}

bool AliHLTTPCDigitReader::GetRCUTrailerData(UChar_t*& trData)
{
  // see header file for class documentation
  PrintWarningOnce(kWarnMissTrailerGetters,"\n"
		   "      !!! This digit reader does not implement the Getters       !!!\n"
		   "      !!! for RCU trailer. Ignoring call.                        !!!");
  if (trData) trData=NULL;
  return 0;
}


void AliHLTTPCDigitReader::PrintWarningOnce(int type, const char* message)
{
  // see header file for class documentation
  if (CheckFlag(type)) return;
  SetFlag(type);
  HLTWarning(message);
}

