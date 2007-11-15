// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTOUT.cxx
    @author Matthias Richter
    @date   
    @brief  The control class for HLTOUT data.                            */

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include <cerrno>
#include "AliHLTOUT.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUT)

AliHLTOUT::AliHLTOUT()
  :
  fSearchDataType(kAliHLTVoidDataType),
  fSearchSpecification(kAliHLTVoidDataSpec),
  fFlags(0),
  fBlockDescList(),
  fCurrent(fBlockDescList.begin()),
  fpBuffer(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

// definitions from ALICE internal note ALICE-INT-2002-010
const unsigned char AliHLTOUT::fgkCDHStatusWord=4;
const unsigned char AliHLTOUT::fgkCDHStatusFlagsOffset=12;

// definitions from ALICE internal note ALICE-INT-2006-XXX
const unsigned char AliHLTOUT::fgkCDHFlagsHLTDecision=6;
const unsigned char AliHLTOUT::fgkCDHFlagsHLTPayload=7;

AliHLTOUT::~AliHLTOUT()
{
  // see header file for class documentation
}

int AliHLTOUT::GetNofDataBlocks()
{
  // see header file for class documentation
  return fBlockDescList.size();
}

int AliHLTOUT::SelectFirstDataBlock(AliHLTComponentDataType dt, AliHLTUInt32_t spec)
{
  // see header file for class documentation
  if (CheckStatusFlag(kLocked)) return -EPERM;
  fCurrent=fBlockDescList.begin();
  fSearchDataType=dt;
  fSearchSpecification=spec;
  return FindAndSelectDataBlock();
}

int AliHLTOUT::SelectNextDataBlock()
{
  // see header file for class documentation
  if (CheckStatusFlag(kLocked)) return -EPERM;
  fCurrent++;
  return FindAndSelectDataBlock();
}

int AliHLTOUT::FindAndSelectDataBlock()
{
  // see header file for class documentation
  if (CheckStatusFlag(kLocked)) return -EPERM;
  int iResult=-ENOENT;
  while (fCurrent!=fBlockDescList.end() && iResult==-ENOENT) {
    if ((fSearchDataType==kAliHLTAnyDataType || (*fCurrent)==fSearchDataType) &&
	fSearchSpecification==kAliHLTVoidDataSpec || (*fCurrent)==fSearchSpecification) {
      iResult=0;
      // TODO: check the byte order on the current system and the byte order of the
      // data block, print warning when missmatch and user did not check
      //AliHLTOUTByteOrder_t blockBO=CheckByteOrder();
      /*
	if (blockBO!=fByteOrder) {
	SetStatusFlag(kByteOrderWarning);

	}
       */
      ClearStatusFlag(kByteOrderChecked);

      // TODO: check the alignment on the current system and the alignment of the
      // data block, print warning when missmatch and user did not check
      ClearStatusFlag(kAlignmentChecked);
    }
    fCurrent++;
  }
  return iResult;
}

int AliHLTOUT::GetDataBlockDescription(AliHLTComponentDataType& dt, AliHLTUInt32_t& spec)
{
  // see header file for class documentation
  int iResult=-ENOENT;
  if (fCurrent!=fBlockDescList.end()) {
    iResult=0;
    dt=(*fCurrent);
    spec=(*fCurrent);
  }
  return iResult;
}

int AliHLTOUT::GetDataBuffer(const AliHLTUInt8_t* &pBuffer, AliHLTUInt32_t& size)
{
  // see header file for class documentation
  int iResult=-ENOENT;
  pBuffer=NULL;
  size=0;
  if (fCurrent!=fBlockDescList.end()) {
    if ((iResult=GetDataBuffer((*fCurrent).GetIndex(), pBuffer, size))>=0) {
      fpBuffer=pBuffer;
    }
  }
  return iResult;  
}

int AliHLTOUT::ReleaseDataBuffer(const AliHLTUInt8_t* pBuffer)
{
  // see header file for class documentation
  int iResult=0;
  if (pBuffer==fpBuffer) {
    fpBuffer=NULL;
  } else {
    HLTWarning("buffer %p does not match the provided one %p", pBuffer, fpBuffer);
  }
  return iResult;  
}

int AliHLTOUT::AddBlockDescriptor(const AliHLTOUTBlockDescriptor desc)
{
  // see header file for class documentation
  if (!CheckStatusFlag(kCollecting)) return -EPERM;
  int iResult=0;
  fBlockDescList.push_back(desc);
  return iResult;  
}

AliHLTOUT::AliHLTOUTByteOrder_t AliHLTOUT::CheckByteOrder()
{
  if (fCurrent!=fBlockDescList.end()) {
    SetStatusFlag(kByteOrderChecked);
    AliHLTOUT::AliHLTOUTByteOrder_t order=CheckBlockByteOrder((*fCurrent).GetIndex());
    return order;
  }
  return kInvalidByteOrder;
}

int AliHLTOUT::CheckAlignment(AliHLTOUT::AliHLTOUTDataType_t type)
{
  if (fCurrent!=fBlockDescList.end()) {
    SetStatusFlag(kAlignmentChecked);
    int alignment=CheckBlockAlignment((*fCurrent).GetIndex(), type);
    return alignment;
  }
  return -ENOENT;
}
