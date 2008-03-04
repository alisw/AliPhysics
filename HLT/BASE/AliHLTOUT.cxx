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
#include <cassert>
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
  fpBuffer(NULL),
  fDataHandlers()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTOUT::~AliHLTOUT()
{
  // see header file for class documentation
}

int AliHLTOUT::Init()
{
  // see header file for class documentation
  int iResult=0;
  SetStatusFlag(kCollecting);
  if ((iResult=GenerateIndex())>=0) {
    if ((iResult=InitHandlers())>=0) {
    }
  }
  ClearStatusFlag(kCollecting);
  return iResult;
}

int AliHLTOUT::GetNofDataBlocks()
{
  // see header file for class documentation
  return fBlockDescList.size();
}

int AliHLTOUT::SelectFirstDataBlock(AliHLTComponentDataType dt, AliHLTUInt32_t spec,
				    AliHLTModuleAgent::AliHLTOUTHandlerType handlerType)
{
  // see header file for class documentation
  if (CheckStatusFlag(kLocked)) return -EPERM;
  fCurrent=fBlockDescList.begin();
  fSearchDataType=dt;
  fSearchSpecification=spec;
  //fSearchHandlerType=handlerType;
  return FindAndSelectDataBlock();
}

int AliHLTOUT::SelectNextDataBlock()
{
  // see header file for class documentation
  if (CheckStatusFlag(kLocked)) return -EPERM;
  if (fCurrent==fBlockDescList.end()) return -ENOENT;
  fCurrent++;
  return FindAndSelectDataBlock();
}

int AliHLTOUT::FindAndSelectDataBlock()
{
  // see header file for class documentation
  if (CheckStatusFlag(kLocked)) return -EPERM;
  int iResult=-ENOENT;
  while (fCurrent!=fBlockDescList.end() && iResult==-ENOENT) {
    if ((*fCurrent)==fSearchDataType &&
	fSearchSpecification==kAliHLTVoidDataSpec || (*fCurrent)==fSearchSpecification &&
	1/*fSearchHandlerType==AliHLTModuleAgent::kUnknownOutput*/) {
      iResult=fCurrent->GetIndex();
      // TODO: check the byte order on the current system and the byte order of the
      // data block, print warning when missmatch and user did not check
      //AliHLTOUTByteOrder blockBO=CheckByteOrder();
      CheckByteOrder();
      /*
	if (blockBO!=fByteOrder) {
	SetStatusFlag(kByteOrderWarning);

	}
       */
      ClearStatusFlag(kByteOrderChecked);

      // TODO: check the alignment on the current system and the alignment of the
      // data block, print warning when missmatch and user did not check
      ClearStatusFlag(kAlignmentChecked);

      break;
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

AliHLTUInt32_t AliHLTOUT::GetDataBlockIndex()
{
  // see header file for class documentation
  if (fCurrent==fBlockDescList.end()) return AliHLTOUTInvalidIndex;
  return fCurrent->GetIndex();
}

int AliHLTOUT::GetDataBuffer(const AliHLTUInt8_t* &pBuffer, AliHLTUInt32_t& size)
{
  // see header file for class documentation
  int iResult=-ENOENT;
  pBuffer=NULL;
  size=0;
  if (fCurrent!=fBlockDescList.end()) {
    if ((iResult=GetDataBuffer(fCurrent->GetIndex(), pBuffer, size))>=0) {
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

AliHLTOUTHandler* AliHLTOUT::GetHandler()
{
  // see header file for class documentation
  AliHLTOUTHandler* pHandler=NULL;
  pHandler=FindHandlerDesc(GetDataBlockIndex());
  return pHandler;
}

int AliHLTOUT::AddBlockDescriptor(const AliHLTOUTBlockDescriptor desc)
{
  // see header file for class documentation
  if (!CheckStatusFlag(kCollecting)) return -EPERM;
  int iResult=0;
  fBlockDescList.push_back(desc);
  return iResult;  
}

AliHLTOUT::AliHLTOUTByteOrder AliHLTOUT::CheckByteOrder()
{
  // see header file for class documentation
  if (fCurrent!=fBlockDescList.end()) {
    SetStatusFlag(kByteOrderChecked);
    AliHLTOUT::AliHLTOUTByteOrder order=CheckBlockByteOrder((*fCurrent).GetIndex());
    return order;
  }
  return kInvalidByteOrder;
}

int AliHLTOUT::CheckAlignment(AliHLTOUT::AliHLTOUTDataType type)
{
  // see header file for class documentation
  if (fCurrent!=fBlockDescList.end()) {
    SetStatusFlag(kAlignmentChecked);
    int alignment=CheckBlockAlignment((*fCurrent).GetIndex(), type);
    return alignment;
  }
  return -ENOENT;
}

int AliHLTOUT::InitHandlers()
{
  // see header file for class documentation
  int iResult=0;
  AliHLTOUTIndexList remnants;
  for (int havedata=SelectFirstDataBlock(kAliHLTAnyDataType, kAliHLTVoidDataSpec); havedata>=0; havedata=SelectNextDataBlock()) {
    remnants.push_back(GetDataBlockIndex());
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    if (GetDataBlockDescription(dt, spec)<0) break;
    for (AliHLTModuleAgent* pAgent=AliHLTModuleAgent::GetFirstAgent(); pAgent && iResult>=0; pAgent=AliHLTModuleAgent::GetNextAgent()) {
      AliHLTModuleAgent::AliHLTOUTHandlerDesc handlerDesc;
      if (pAgent->GetHandlerDescription(dt, spec, handlerDesc)>0) {
	AliHLTOUTHandlerListEntry entry(pAgent->GetOutputHandler(dt, spec), handlerDesc, pAgent, GetDataBlockIndex());
	InsertHandler(entry);
	remnants.pop_back();
	break;
      }
    }
  }
  if (remnants.size()>0) {
    HLTWarning("no handlers found for %d data blocks out of %d", remnants.size(), GetNofDataBlocks());
    AliHLTOUTBlockDescriptorVector::iterator block=fBlockDescList.begin();
    for (AliHLTOUTIndexList::iterator element=remnants.begin(); element!=remnants.end(); element++) {
      for (int trials=0; trials<2; trials++) {
	do {
	  // we start searching the index from the current position in the block list
	  if ((*block).GetIndex()==*element) break;
	} while ((++block)!=fBlockDescList.end());
	if (block==fBlockDescList.end()) {
	  // rewind and try again
	  block=fBlockDescList.begin();
	}
      }
      assert(block!=fBlockDescList.end());
      if (block!=fBlockDescList.end()) {
	HLTDebug("   %s", AliHLTComponent::DataType2Text((AliHLTComponentDataType)*block).c_str());
      }
    }
  }
  return iResult;
}

int AliHLTOUT::InsertHandler(const AliHLTOUTHandlerListEntry &entry)
{
  // see header file for class documentation
  int iResult=0;
  AliHLTOUTHandlerListEntryVector::iterator element=fDataHandlers.begin();
  while (element!=fDataHandlers.end()) {
    if (entry==(*element)) break;
    element++;
  }
  if (element==fDataHandlers.end()) {
    fDataHandlers.push_back(entry);
  } else {
    element->AddIndex(const_cast<AliHLTOUTHandlerListEntry&>(entry));
  }
  return iResult;
}

AliHLTOUT::AliHLTOUTHandlerListEntry AliHLTOUT::FindHandlerDesc(AliHLTUInt32_t blockIndex)
{
  // see header file for class documentation
  AliHLTOUTHandlerListEntryVector::iterator element=fDataHandlers.begin();
  while (element!=fDataHandlers.end()) {
    if (element->HasIndex(blockIndex)) {
      return *element;
    }
    element++;
  }
  return AliHLTOUT::AliHLTOUTHandlerListEntry::fgkVoidHandlerListEntry;
}

AliHLTOUT::AliHLTOUTHandlerListEntry::AliHLTOUTHandlerListEntry()
  :
  fpHandler(NULL),
  fpHandlerDesc(NULL),
  fpAgent(NULL),
  fBlocks()
{
  // see header file for class documentation
}

AliHLTOUT::AliHLTOUTHandlerListEntry::AliHLTOUTHandlerListEntry(AliHLTOUTHandler* pHandler, 
								AliHLTModuleAgent::AliHLTOUTHandlerDesc& handlerDesc,
								AliHLTModuleAgent* pAgent,
								AliHLTUInt32_t index)
  :
  fpHandler(pHandler),
  fpHandlerDesc(new AliHLTModuleAgent::AliHLTOUTHandlerDesc),
  fpAgent(pAgent),
  fBlocks()
{
  // see header file for class documentation
  *fpHandlerDesc=handlerDesc;
  fBlocks.push_back(index);
}

AliHLTOUT::AliHLTOUTHandlerListEntry::AliHLTOUTHandlerListEntry(const AliHLTOUTHandlerListEntry& src)
  :
  fpHandler(src.fpHandler),
  fpHandlerDesc(new AliHLTModuleAgent::AliHLTOUTHandlerDesc),
  fpAgent(src.fpAgent),
  fBlocks()
{
  // see header file for class documentation
  *fpHandlerDesc=*src.fpHandlerDesc;
  fBlocks.assign(src.fBlocks.begin(), src.fBlocks.end());
}

AliHLTOUT::AliHLTOUTHandlerListEntry::~AliHLTOUTHandlerListEntry()
{
  // see header file for class documentation
  if (fpHandlerDesc) delete fpHandlerDesc;
  fpHandlerDesc=NULL;
}

AliHLTOUT::AliHLTOUTHandlerListEntry& AliHLTOUT::AliHLTOUTHandlerListEntry::operator=(const AliHLTOUTHandlerListEntry& src)
{
  // see header file for class documentation
  fpHandler=src.fpHandler;
  *fpHandlerDesc=*src.fpHandlerDesc;
  fpAgent=src.fpAgent;
  fBlocks.assign(src.fBlocks.begin(), src.fBlocks.end());
  return *this;
}

AliHLTUInt32_t AliHLTOUT::AliHLTOUTHandlerListEntry::operator[](int i) const
{
  // see header file for class documentation
  return (int)fBlocks.size()>i?fBlocks[i]:AliHLTOUTInvalidIndex;
}

bool AliHLTOUT::AliHLTOUTHandlerListEntry::operator==(const AliHLTOUTHandlerListEntry& entry) const
{
  // see header file for class documentation
  if (entry.fpHandler!=fpHandler || fpHandler==NULL) return false;
  assert(entry.fpAgent==fpAgent);
  if (entry.fpAgent!=fpAgent) return false;
  return true;
}

void AliHLTOUT::AliHLTOUTHandlerListEntry::AddIndex(AliHLTOUT::AliHLTOUTHandlerListEntry &desc)
{
  // see header file for class documentation
  AliHLTOUTIndexList::iterator element;
  for (element=desc.fBlocks.begin(); element!=desc.fBlocks.end(); element++) {
    AddIndex(*element);
  }  
}

void AliHLTOUT::AliHLTOUTHandlerListEntry::AddIndex(AliHLTUInt32_t index)
{
  // see header file for class documentation
  fBlocks.push_back(index);
}

bool AliHLTOUT::AliHLTOUTHandlerListEntry::HasIndex(AliHLTUInt32_t index)
{
  // see header file for class documentation
  AliHLTOUTIndexList::iterator element;
  for (element=fBlocks.begin(); element!=fBlocks.end(); element++) {
    if (*element==index) return true;
  }
  return false;
}

const AliHLTOUT::AliHLTOUTHandlerListEntry AliHLTOUT::AliHLTOUTHandlerListEntry::fgkVoidHandlerListEntry;
