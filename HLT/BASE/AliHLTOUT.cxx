// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
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

/// @file   AliHLTOUT.cxx
/// @author Matthias Richter
/// @date   
/// @brief  The control class for HLTOUT data.
///

#include <cerrno>
#include <cassert>
#include "AliHLTOUT.h"
#include "AliHLTMessage.h"
#include "AliHLTMisc.h"
#include "TSystem.h"
#include "TClass.h"
#include "TROOT.h"
#include <sstream>
#include <iomanip>

using std::cout;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUT)

AliHLTOUT::AliHLTOUT()
  :
  fSearchDataType(kAliHLTVoidDataType),
  fSearchSpecification(kAliHLTVoidDataSpec),
  fSearchHandlerType(AliHLTModuleAgent::kUnknownOutput),
  fFlags(kSkipProcessed),
  fBlockDescList(),
  fCurrent(0),
  fpBuffer(NULL),
  fDataHandlers(),
  fbVerbose(false),
  fLog()
  , fpDataObject(NULL)
  , fpObjectBuffer(NULL)
  , fObjectBufferSize(0)
  , fCurrentEventId(kAliHLTVoidEventID)
{
  // constructor
  // 
  // The control class for HLTOUT data
  // see header file for class documentation
  // author Matthias Richter
}

AliHLTOUT::~AliHLTOUT()
{
  // destructor
  if (CheckStatusFlag(kIsSubCollection)) {
    fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "~AliHLTOUT" , __FILE__ , __LINE__ , "severe internal error: collection has not been released, potential crash due to invalid pointer");
  }

  if (fpDataObject) {
    fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "GetDataObject" , __FILE__ , __LINE__ , "data object has not been released, potential memory leak");
  }
  fpDataObject=NULL;
}
AliHLTOUT* AliHLTOUT::fgGlobalInstance=NULL;

int AliHLTOUT::Init()
{
  // see header file for class documentation
  int iResult=0;

  // ignore if already initialized
  if (fBlockDescList.size()>0) {
    return 0;
  }

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
  // get number of data blocks
  return fBlockDescList.size();
}

int AliHLTOUT::SelectFirstDataBlock(AliHLTComponentDataType dt, AliHLTUInt32_t spec,
				    AliHLTModuleAgent::AliHLTOUTHandlerType handlerType,
				    bool skipProcessed)
{
  // select the first data block according to data type, specification and
  // handler type
  fCurrent=0;
  fSearchDataType=dt;
  fSearchSpecification=spec;
  fSearchHandlerType=handlerType;
  if (skipProcessed) SetStatusFlag(kSkipProcessed);
  else ClearStatusFlag(kSkipProcessed);
  return FindAndSelectDataBlock();
}

int AliHLTOUT::SelectNextDataBlock()
{
  // select next data block according to selection criteria specified
  // for SelectFirstDataBlock
  if (fCurrent>=fBlockDescList.size()) return -ENOENT;
  fCurrent++;
  return FindAndSelectDataBlock();
}

int AliHLTOUT::FindAndSelectDataBlock()
{
  // Select data block according to data type and specification, internal function
  // invoked by SelectFirstDataBlock/SelectNextDataBlock
  if (CheckStatusFlag(kLocked)) return -EPERM;
  int iResult=-ENOENT;
  while (fCurrent<fBlockDescList.size() && iResult==-ENOENT) {
    if (fBlockDescList[fCurrent]==fSearchDataType &&
	(fSearchSpecification==kAliHLTVoidDataSpec || fBlockDescList[fCurrent]==fSearchSpecification) &&
	(fSearchHandlerType==AliHLTModuleAgent::kUnknownOutput || FindHandlerDesc(fCurrent)==fSearchHandlerType) &&
	(!CheckStatusFlag(kBlockSelection) || fBlockDescList[fCurrent].IsSelected()) &&
	(!CheckStatusFlag(kSkipProcessed) || !fBlockDescList[fCurrent].IsProcessed())) {
      iResult=fBlockDescList[fCurrent].GetIndex();
      // TODO: check the byte order on the current system and the byte order of the
      // data block, print warning when mismatch and user did not check
      //AliHLTOUTByteOrder blockBO=CheckByteOrder();
      CheckByteOrder();
      /*
	if (blockBO!=fByteOrder) {
	SetStatusFlag(kByteOrderWarning);

	}
       */
      ClearStatusFlag(kByteOrderChecked);

      // TODO: check the alignment on the current system and the alignment of the
      // data block, print warning when mismatch and user did not check
      ClearStatusFlag(kAlignmentChecked);

      break;
    }
    fCurrent++;
  }
  return iResult;
}

int AliHLTOUT::GetDataBlockDescription(AliHLTComponentDataType& dt, AliHLTUInt32_t& spec)
{
  // fill data type and specification
  int iResult=-ENOENT;
  if (fCurrent<fBlockDescList.size()) {
    iResult=0;
    dt=fBlockDescList[fCurrent];
    spec=fBlockDescList[fCurrent];
  }
  return iResult;
}

const AliHLTOUT::AliHLTOUTHandlerListEntry& AliHLTOUT::GetDataBlockHandlerDesc()
{
  // Get handler descriptor of the selected data block.
  return FindHandlerDesc(fCurrent);
}

AliHLTModuleAgent::AliHLTOUTHandlerType AliHLTOUT::GetDataBlockHandlerType()
{
  // Get handler type of the selected data block.
  AliHLTModuleAgent::AliHLTOUTHandlerDesc desc=FindHandlerDesc(fCurrent);
  AliHLTModuleAgent::AliHLTOUTHandlerType type=desc;
  return type;
}

AliHLTUInt32_t AliHLTOUT::GetDataBlockIndex()
{
  // Get the index of the current data block.
  if (fCurrent>=fBlockDescList.size()) return AliHLTOUTInvalidIndex;
  return fBlockDescList[fCurrent].GetIndex();
}

int AliHLTOUT::GetDataBuffer(AliHLTComponentBlockData& desc)
{
  // fill block data descriptor and select the current data buffer
  // buffer has to be released using ReleaseDataBuffer
  int iResult=-ENOENT;
  if (fCurrent<fBlockDescList.size()) {
    AliHLTComponent::FillBlockData(desc);
    if ((iResult=fBlockDescList[fCurrent].GetDataBuffer(fpBuffer, desc.fSize))>=0) {
      desc.fPtr=const_cast<void*>(reinterpret_cast<const void*>(fpBuffer));
      desc.fDataType=fBlockDescList[fCurrent];
      desc.fSpecification=fBlockDescList[fCurrent];
    }
  }
  return iResult;
}

int AliHLTOUT::GetDataBuffer(const AliHLTUInt8_t* &pBuffer, AliHLTUInt32_t& size)
{
  // select and return the current data buffer
  // buffer has to be released using ReleaseDataBuffer
  int iResult=-ENOENT;
  pBuffer=NULL;
  size=0;
  if (fCurrent<fBlockDescList.size()) {
    if ((iResult=fBlockDescList[fCurrent].GetDataBuffer(fpBuffer, size))>=0) {
      pBuffer=fpBuffer;
    }
  }
  return iResult;  
}

int AliHLTOUT::ReleaseDataBuffer(const AliHLTUInt8_t* pBuffer)
{
  // release data buffer, previously returned by GetDataBuffer
  int iResult=0;
  if (pBuffer==fpBuffer) {
    fpBuffer=NULL;
  } else {
    fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "ReleaseDataBuffer" , __FILE__ , __LINE__ , "buffer %p does not match the provided one %p", pBuffer, fpBuffer);
  }
  return iResult;  
}

AliHLTModuleAgent* AliHLTOUT::GetAgent()
{
  // get module agent of the selected data block
  AliHLTModuleAgent* pAgent=NULL;
  pAgent=FindHandlerDesc(fCurrent);
  return pAgent;
}

AliHLTOUTHandler* AliHLTOUT::GetHandler()
{
  // get HLTOUT handler of the selected data block
  AliHLTOUTHandler* pHandler=NULL;
  pHandler=FindHandlerDesc(fCurrent);
  return pHandler;
}

int AliHLTOUT::WriteESD(const AliHLTUInt8_t* /*pBuffer*/, AliHLTUInt32_t /*size*/, AliHLTComponentDataType /*dt*/, AliESDEvent* /*tgtesd*/) const
{
  // default function, child must overload
  fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "WriteESD" , __FILE__ , __LINE__ , "method not implemented in base class");
  return -ENOSYS;
}

int AliHLTOUT::AddBlockDescriptor(const AliHLTOUTBlockDescriptor desc)
{
  // add new block descriptor
  if (!CheckStatusFlag(kCollecting)) return -EPERM;
  int iResult=0;
  fBlockDescList.push_back(desc);
  return iResult;  
}

AliHLTOUT::AliHLTOUTByteOrder AliHLTOUT::CheckByteOrder()
{
  // check the byte order of the current data block
  // NOTE: this functionality has not been tested and development was hardly finished
  if (fCurrent<fBlockDescList.size()) {
    SetStatusFlag(kByteOrderChecked);
    AliHLTOUT::AliHLTOUTByteOrder order=CheckBlockByteOrder(fBlockDescList[fCurrent].GetIndex());
    return order;
  }
  return kInvalidByteOrder;
}

int AliHLTOUT::CheckAlignment(AliHLTOUT::AliHLTOUTDataType type)
{
  // check alignment of the current data block
  // NOTE: this functionality has not been tested and development was hardly finished
  if (fCurrent<fBlockDescList.size()) {
    SetStatusFlag(kAlignmentChecked);
    int alignment=CheckBlockAlignment(fBlockDescList[fCurrent].GetIndex(), type);
    return alignment;
  }
  return -ENOENT;
}

int AliHLTOUT::InitHandlers()
{
  // init handlers for all registered blocks
  int iResult=0;
  AliHLTOUTIndexList remnants;
  int iCount=0;
  
  // -- Sergey Gorbunov  
  // Clear list of handler pointers, created in previous event, because some of them  may not exist anymore
  // --
  fDataHandlers.clear(); 

  for (int havedata=SelectFirstDataBlock(kAliHLTAnyDataType, kAliHLTVoidDataSpec); havedata>=0; havedata=SelectNextDataBlock()) {
    iCount++;
    remnants.push_back(GetDataBlockIndex());
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    if (GetDataBlockDescription(dt, spec)<0) break;
    bool bHaveHandler=false;
    for (AliHLTModuleAgent* pAgent=AliHLTModuleAgent::GetFirstAgent(); pAgent && iResult>=0; pAgent=AliHLTModuleAgent::GetNextAgent()) {
      AliHLTModuleAgent::AliHLTOUTHandlerDesc handlerDesc;
      if (pAgent->GetHandlerDescription(dt, spec, handlerDesc)>0) {
	AliHLTOUTHandlerListEntry entry(pAgent->GetOutputHandler(dt, spec), handlerDesc, pAgent, GetDataBlockIndex());
	InsertHandler(fDataHandlers, entry);
	remnants.pop_back();
	bHaveHandler=true;
	if (fbVerbose) {
	  stringstream sout;
	  sout << "adding handler for block " << AliHLTComponent::DataType2Text(dt).c_str()
	       << " 0x" << setfill('0') << setw(8) << hex << spec;
	  cout << sout.str() << endl;
	}
	break;
      }
    }
    if (!bHaveHandler && (dt==kAliHLTDataTypeESDObject || dt==kAliHLTDataTypeESDTree)) {
      // ESDs are handled by the framework
      remnants.pop_back();
    }
  }

  // warning if some of the data blocks are not selected by the kAliHLTAnyDataType
  // criterion
  if (GetNofDataBlocks()>iCount) {
    fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "InitHandlers" , __FILE__ , __LINE__ , "incomplete data type in %d out of %d data block(s)", GetNofDataBlocks()-iCount, GetNofDataBlocks());
  }

  // warning if handler not found
  if (remnants.size()>0) {
    fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "InitHandlers" , __FILE__ , __LINE__ , "no handlers found for %d data blocks out of %d", remnants.size(), iCount);
    AliHLTOUTBlockDescriptorVector::iterator block=fBlockDescList.begin();
    for (AliHLTOUTIndexList::iterator element=remnants.begin();
	 element!=remnants.end() && block!=fBlockDescList.end();
	 element++) {
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
    }
  }

  if (fbVerbose) Print();

  return iResult;
}

int AliHLTOUT::InsertHandler(AliHLTOUTHandlerListEntryVector& list, const AliHLTOUTHandlerListEntry &entry)
{
  // insert handler into list, called from child implementations
  int iResult=0;
  AliHLTOUTHandlerListEntryVector::iterator element=list.begin();
  for (; element!=list.end();
	 element++) {
    if (entry==(*element)) break;
  }
  if (element==list.end()) {
    list.push_back(entry);
  } else {
    element->AddIndex(const_cast<AliHLTOUTHandlerListEntry&>(entry));
  }
  return iResult;
}

int AliHLTOUT::FillHandlerList(AliHLTOUTHandlerListEntryVector& list, AliHLTModuleAgent::AliHLTOUTHandlerType handlerType)
{
  // fill a list according to specified handler type
  int iResult=0;
  for (iResult=SelectFirstDataBlock(kAliHLTAnyDataType, kAliHLTVoidDataSpec, handlerType);
       iResult>=0;
       iResult=SelectNextDataBlock()) {
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    GetDataBlockDescription(dt, spec);
    AliHLTOUTHandler* pHandler=GetHandler();
    if (!pHandler) {
      fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "FillHandlerList" , __FILE__ , __LINE__ , 
			 "missing HLTOUT handler for block of type kChain: agent %s, data type %s, specification %#x, ... skipping data block",
			 GetAgent()?GetAgent()->GetModuleId():"invalid",
			 AliHLTComponent::DataType2Text(dt).c_str(), spec);
    } else {
      InsertHandler(list, GetDataBlockHandlerDesc());
    }
  }
  // TODO: the return value of SelectFirst/NextDataBlock must be
  // changed in order to avoid this check
  if (iResult==-ENOENT) iResult=0;

  return iResult;
}

int AliHLTOUT::RemoveEmptyDuplicateHandlers(AliHLTOUTHandlerListEntryVector& list)
{
  // remove empty handlers from list
  int iResult=0;
  AliHLTOUTHandlerListEntryVector::iterator element=list.begin();
  while (element!=list.end()) {
    if (element->IsEmpty()) {
      AliHLTOUTHandler* pHandler=*element;
      AliHLTModuleAgent* pAgent=*element;
      AliHLTModuleAgent::AliHLTOUTHandlerDesc desc=*element;
      if (FindHandler(list, desc)>=0) {
	element=list.erase(element);
	if (pAgent) {
	  // -- Sergey Gorbunov
	  // Do not delete the handler here, only remove it from the input list. 
	  // When the pointer is deleted form AliHLTSystem::fpChainHandlers list, it is still kept in AliHLTOut::fDataHandlers list and may be used later
	  // -- 
	  // pAgent->DeleteOutputHandler(pHandler); 
	}
	// we are already at the next element
	continue;
      }
    }
    element++;
  }
  return iResult;
}

int AliHLTOUT::FindHandler(AliHLTOUTHandlerListEntryVector& list, const AliHLTModuleAgent::AliHLTOUTHandlerDesc desc)
{
  // find handler according to descriptor
  for (int i=0; i<(int)list.size(); i++) {
    if (list[i]==desc) return i;
  }
  return -ENOENT;
}

int AliHLTOUT::InvalidateBlocks(AliHLTOUTHandlerListEntryVector& list)
{
  // invalidate all handlers in a list
  for (AliHLTOUTHandlerListEntryVector::iterator element=list.begin();
	 element!=list.end();
	 element++) {
    element->InvalidateBlocks();
  }
  return 0;
}

const AliHLTOUT::AliHLTOUTHandlerListEntry& AliHLTOUT::FindHandlerDesc(AliHLTUInt32_t blockIndex)
{
  // get handler description
  if (blockIndex<fBlockDescList.size()) {
    return fBlockDescList[blockIndex].GetHandlerDesc();
  }
  return const_cast<AliHLTOUT::AliHLTOUTHandlerListEntry&>(AliHLTOUT::AliHLTOUTHandlerListEntry::VoidHandlerListEntry());
}

AliHLTOUT::AliHLTOUTHandlerListEntry::AliHLTOUTHandlerListEntry()
  :
  fpHandler(NULL),
  fpHandlerDesc(NULL),
  fpAgent(NULL),
  fBlocks()
{
  // default constructor
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
  // constructor
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
  // copy constructor
  *fpHandlerDesc=*src.fpHandlerDesc;
  fBlocks.assign(src.fBlocks.begin(), src.fBlocks.end());
}

AliHLTOUT::AliHLTOUTHandlerListEntry::~AliHLTOUTHandlerListEntry()
{
  // destructor
  if (fpHandlerDesc) delete fpHandlerDesc;
  fpHandlerDesc=NULL;
}

AliHLTOUT::AliHLTOUTHandlerListEntry& AliHLTOUT::AliHLTOUTHandlerListEntry::operator=(const AliHLTOUTHandlerListEntry& src)
{
  // assignment operator
  if (this==&src) return *this;  
  fpHandler=src.fpHandler;
  if (src.fpHandlerDesc)
    *fpHandlerDesc=*src.fpHandlerDesc;
  fpAgent=src.fpAgent;
  fBlocks.assign(src.fBlocks.begin(), src.fBlocks.end());
  return *this;
}

AliHLTUInt32_t AliHLTOUT::AliHLTOUTHandlerListEntry::operator[](int i) const
{
  // access operator
  return (int)fBlocks.size()>i?fBlocks[i]:AliHLTOUTInvalidIndex;
}

bool AliHLTOUT::AliHLTOUTHandlerListEntry::operator==(const AliHLTOUTHandlerListEntry& entry) const
{
  // comparison operator
  if (entry.fpHandler!=fpHandler || fpHandler==NULL) return false;
  assert(entry.fpAgent==fpAgent);
  if (entry.fpAgent!=fpAgent) return false;
  return true;
}

bool AliHLTOUT::AliHLTOUTHandlerListEntry::operator==(const AliHLTModuleAgent::AliHLTOUTHandlerType handlerType) const
{
  // comparison operator
  if (!fpHandlerDesc) return false;
  return *fpHandlerDesc==handlerType;
}

bool AliHLTOUT::AliHLTOUTHandlerListEntry::operator==(const AliHLTModuleAgent::AliHLTOUTHandlerDesc desc) const
{
  // comparison operator
  if (!fpHandlerDesc) return false;
  return *fpHandlerDesc==desc;
}

void AliHLTOUT::AliHLTOUTHandlerListEntry::AddIndex(const AliHLTOUT::AliHLTOUTHandlerListEntry &desc)
{
  // add block index, a handler can serve multiple blocks
  AliHLTOUTIndexList::const_iterator element;
  for (element=desc.fBlocks.begin(); element!=desc.fBlocks.end(); element++) {
    AddIndex(*element);
  }  
}

void AliHLTOUT::AliHLTOUTHandlerListEntry::AddIndex(AliHLTUInt32_t index)
{
  // add block index, a handler can serve multiple blocks
  fBlocks.push_back(index);
}

bool AliHLTOUT::AliHLTOUTHandlerListEntry::HasIndex(AliHLTUInt32_t index) const
{
  // check if handler serves the specified block
  AliHLTOUTIndexList::iterator element;
  for (unsigned int i=0; i<fBlocks.size(); i++) {
    if (fBlocks[i]==index) return true;
  }
  return false;
}

const AliHLTOUT::AliHLTOUTHandlerListEntry AliHLTOUT::AliHLTOUTHandlerListEntry::fgkVoidHandlerListEntry;

AliHLTUInt64_t AliHLTOUT::ByteSwap64(AliHLTUInt64_t src)
{
  // swap a 64 bit number
  return ((src & 0xFFULL) << 56) | 
    ((src & 0xFF00ULL) << 40) | 
    ((src & 0xFF0000ULL) << 24) | 
    ((src & 0xFF000000ULL) << 8) | 
    ((src & 0xFF00000000ULL) >> 8) | 
    ((src & 0xFF0000000000ULL) >> 24) | 
    ((src & 0xFF000000000000ULL) >>  40) | 
    ((src & 0xFF00000000000000ULL) >> 56);
}

AliHLTUInt32_t AliHLTOUT::ByteSwap32(AliHLTUInt32_t src)
{
  // swap a 32 bit number
  return ((src & 0xFFULL) << 24) | 
    ((src & 0xFF00ULL) << 8) | 
    ((src & 0xFF0000ULL) >> 8) | 
    ((src & 0xFF000000ULL) >> 24);
}

AliHLTOUT* AliHLTOUT::New(AliRawReader* pRawReader)
{
  // transparently create HLTOUT implementation for AliRawReader
  AliHLTOUT* instance=AliHLTMisc::LoadInstance((AliHLTOUT*)0, "AliHLTOUTRawReader", "libHLTrec.so");
  if (instance) {
    instance->SetParam(pRawReader);
  }
  return instance;
}

AliHLTOUT* AliHLTOUT::New(TTree* pDigitTree, int event)
{
  // transparently create HLTOUT implementation for digit tree
  AliHLTOUT* instance=AliHLTMisc::LoadInstance((AliHLTOUT*)0, "AliHLTOUTDigitReader", "libHLTrec.so");
  if (instance) {
    instance->SetParam(pDigitTree, event);
  }
  return instance;
}

AliHLTOUT* AliHLTOUT::New(const char* filename, int event)
{
  // transparently create HLTOUT implementation for raw file
  AliHLTOUT* instance=AliHLTMisc::LoadInstance((AliHLTOUT*)0, "AliHLTOUTDigitReader", "libHLTrec.so");
  if (instance) {
    instance->SetParam(filename, event);
  }
  return instance;
}

void AliHLTOUT::Delete(AliHLTOUT* pInstance)
{
  // delete the HLTOUT instance
  // check if the library is still there in order to have the
  // destructor available

  if (!pInstance) return;
  if (pInstance==fgGlobalInstance) return;

  TClass* pCl1=TClass::GetClass("AliHLTOUTRawReader");
  TClass* pCl2=TClass::GetClass("AliHLTOUTDigitReader");
  if (!pCl1 && !pCl2) {
    AliHLTLogging log;
    log.Logging(kHLTLogError, "AliHLTOUT::Delete", "HLTOUT handling", "potential memory leak: libHLTrec library not available, skipping destruction %p", pInstance);    
    return;
  }

  delete pInstance;  
}

void AliHLTOUT::SetParam(AliRawReader* /*pRawReader*/)
{
  // see header file for class documentation
  // default implementation, we should never get here
  // this function can only be called from the class itsself and
  // is intended to be used with the New functions. If we get into
  // the default implementation there is a class mismatch.
  assert(0);
  fLog.LoggingVarargs(kHLTLogFatal, "AliHLTOUT", "SetParam" , __FILE__ , __LINE__ , "severe internal error: class mismatch");
}

void AliHLTOUT::SetParam(TTree* /*pDigitTree*/, int /*event*/)
{
  // see header file for class documentation
  // default implementation, we should never get here
  // this function can only be called from the class itsself and
  // is intended to be used with the New functions. If we get into
  // the default implementation there is a class mismatch.
  assert(0);
  fLog.LoggingVarargs(kHLTLogFatal, "AliHLTOUT", "SetParam" , __FILE__ , __LINE__ , "severe internal error: class mismatch");
}

void AliHLTOUT::SetParam(const char* /*filename*/, int /*event*/)
{
  // see header file for class documentation
  // default implementation, we should never get here
  // this function can only be called from the class itsself and
  // is intended to be used with the New functions. If we get into
  // the default implementation there is a class mismatch.
  assert(0);
  fLog.LoggingVarargs(kHLTLogFatal, "AliHLTOUT", "SetParam" , __FILE__ , __LINE__ , "severe internal error: class mismatch");
}

int AliHLTOUT::SelectDataBlock()
{
  // mark the current data block for processing
  int iResult=0;
  if (fCurrent>=fBlockDescList.size()) return 0;
  fBlockDescList[fCurrent].Select(true);
  EnableBlockSelection();
  return iResult;
}

int AliHLTOUT::SelectDataBlocks(const AliHLTOUTHandlerListEntry* pHandlerEntry)
{
  // mark all data blocks served by specified handler for processing
  int iResult=0;
  if (!pHandlerEntry) return 0;

  AliHLTModuleAgent* pAgent=*pHandlerEntry;
  AliHLTLogging log;
  log.Logging(kHLTLogDebug, "AliHLTOUT::SelectDataBlocks", "HLTOUT handling", "selecting blocks for handler %s", pAgent->GetModuleId());
  AliHLTOUTBlockDescriptorVector::iterator element;
  for (AliHLTOUTBlockDescriptorVector::iterator block=fBlockDescList.begin();
       block!=fBlockDescList.end();
       block++) {
    if (block->GetHandlerDesc()==*pHandlerEntry && pHandlerEntry->HasIndex(block->GetIndex())) {
      block->Select(true);
      log.Logging(kHLTLogDebug, "AliHLTOUT::SelectDataBlocks", "HLTOUT handling", "   select block %s", AliHLTComponent::DataType2Text(*block).c_str());
    } else {
      log.Logging(kHLTLogDebug, "AliHLTOUT::SelectDataBlocks", "HLTOUT handling", "   skip block %s", AliHLTComponent::DataType2Text(*block).c_str());
      block->Select(false);
    }
  }
  EnableBlockSelection();

  // Matthias 2009-07-03 bugfix: the fCurrent position was not reset at that
  // place. Also I think the data type and specification must be set in order
  // to make SelectFirst/NextDataBlock working on the selected collection
  // of data blocks
  AliHLTModuleAgent::AliHLTOUTHandlerDesc pHandlerDesc=*pHandlerEntry; 
  fSearchDataType=pHandlerDesc;
  fSearchSpecification=kAliHLTVoidDataSpec;
  fSearchHandlerType=pHandlerDesc;
  fCurrent=0;
  
  return iResult;
}

int AliHLTOUT::EnableBlockSelection()
{
  // enable block selection, in this mode only the blocks marked for
  // processing can be accessed
  SetStatusFlag(kBlockSelection);
  return 0;
}

int AliHLTOUT::DisableBlockSelection()
{
  // disable block selection
  ClearStatusFlag(kBlockSelection);
  return 0;
}

int AliHLTOUT::ResetBlockSelection()
{
  // reset the 'selected' flag for all blocks
  for (AliHLTOUTBlockDescriptorVector::iterator block=fBlockDescList.begin();
       block!=fBlockDescList.end();
       block++) {
    block->Select(false);
  }
  return 0;
}

int AliHLTOUT::MarkDataBlockProcessed()
{
  // mark the current data block as 'processed'
  int iResult=0;
  if (fCurrent>=fBlockDescList.size()) return 0;
  fBlockDescList[fCurrent].MarkProcessed();
  return iResult;
}

int AliHLTOUT::MarkDataBlocksProcessed(const AliHLTOUTHandlerListEntry* pHandlerDesc)
{
  // mark all data blocks served by handler as processed
  int iResult=0;
  if (!pHandlerDesc) return 0;

  AliHLTOUTBlockDescriptorVector::iterator element;
  for (AliHLTOUTBlockDescriptorVector::iterator block=fBlockDescList.begin();
       block!=fBlockDescList.end();
       block++) {
    if (block->GetHandlerDesc()==*pHandlerDesc && pHandlerDesc->HasIndex(block->GetIndex()))
      block->MarkProcessed();
  }
  
  return iResult;
}

int AliHLTOUT::AddSubCollection(AliHLTOUT* pCollection)
{
  // add a sub-collection to the HLTOUT instance
  // all blocks of the sub-collection are accessed transparently through the master instance
  int iResult=0;
  if (!pCollection) return 0;

  SetStatusFlag(kCollecting);  
  int index=-1;
  for (index=pCollection->SelectFirstDataBlock();
       index>=0;
       index=pCollection->SelectNextDataBlock()) {
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    pCollection->GetDataBlockDescription(dt, spec);  
    AliHLTOUTBlockDescriptor desc(dt, spec, index, pCollection);
    AddBlockDescriptor(desc);
    iResult++;
  }
  if (iResult>0) {
    if (CheckStatusFlag(kIsSubCollection)) {
      fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "AddSubCollection" , __FILE__ , __LINE__ , "HLTOUT object %p has already been added as sub-collection", pCollection);
    } else {
      pCollection->SetStatusFlag(kIsSubCollection);
    }
  }
  ClearStatusFlag(kCollecting);  

  return iResult;
}

int AliHLTOUT::ReleaseSubCollection(AliHLTOUT* pCollection)
{
  // release a sub-collection
  int iResult=0;
  if (!pCollection) return 0;

  AliHLTOUTBlockDescriptorVector::iterator block=fBlockDescList.begin();
  while (block!=fBlockDescList.end()) {
    if ((*block)==pCollection) {
      block=fBlockDescList.erase(block);
      continue;
    }
    block++;
  }
  pCollection->ClearStatusFlag(kIsSubCollection);

  return iResult;
}

int AliHLTOUT::Reset()
{
  // reset HLTOUT instance
  // clears all blocks and handler descriptions
  int iResult=0;
  AliHLTOUTPVector subCollections;
  AliHLTOUTBlockDescriptorVector::iterator block=fBlockDescList.begin();
  while (block!=fBlockDescList.end()) {
    if (!((*block)==this)) {
      AliHLTOUTPVector::iterator collection=subCollections.begin();
      for (; collection!=subCollections.end(); collection++)
	if((*block)==*collection) break;
      if (collection==subCollections.end())
	subCollections.push_back(block->GetCollection());
    }
    block=fBlockDescList.erase(block);
  }

  for (AliHLTOUTPVector::iterator collection=subCollections.begin(); 
       collection!=subCollections.end(); collection++) {
    (*collection)->Reset();
    (*collection)->ClearStatusFlag(kIsSubCollection);
  }

  ResetInput();
  fCurrentEventId=kAliHLTVoidEventID;

  return iResult;
}

int AliHLTOUT::ResetInput()
{
  // default implementation, nothing to do
  return 0;
}

const AliHLTOUT::AliHLTOUTHandlerListEntry& AliHLTOUT::AliHLTOUTBlockDescriptor::GetHandlerDesc()
{
  // see header file for class documentation
  if (fpCollection) {
    AliHLTOUTHandlerListEntryVector::iterator element=fpCollection->fDataHandlers.begin();
    while (element!=fpCollection->fDataHandlers.end()) {
      if (element->HasIndex(GetIndex())) {
	return *element;
      }
      element++;
    }
  }
  return const_cast<AliHLTOUT::AliHLTOUTHandlerListEntry&>(AliHLTOUT::AliHLTOUTHandlerListEntry::VoidHandlerListEntry());
}

TObject* AliHLTOUT::GetDataObject()
{
  // check if the current block encodes a ROOT object and expand it
  if (fpDataObject) {
    fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "GetDataObject" , __FILE__ , __LINE__ , "data object has not been released, potential memory leak");
    ReleaseDataBuffer(fpObjectBuffer);
  }
  fpObjectBuffer=NULL;
  fObjectBufferSize=0;
  fpDataObject=NULL;

  if (GetDataBuffer(fpObjectBuffer, fObjectBufferSize)>=0) {
    fpDataObject=AliHLTMessage::Extract(fpObjectBuffer, fObjectBufferSize);
  } else {
    fLog.LoggingVarargs(kHLTLogError, "AliHLTOUT", "GetDataObject" , __FILE__ , __LINE__ , "can not fetch data buffer");    
  }

  return fpDataObject;
}

int AliHLTOUT::ReleaseDataObject(TObject* pObject)
{
  // release a ROOT object previously expanded from the currentr data block
  if (!pObject) return -EINVAL;
  if (pObject!=fpDataObject) {
    fLog.LoggingVarargs(kHLTLogError, "AliHLTOUT", "GetDataObject" , __FILE__ , __LINE__ , "attempt to release wrong data object %p, expected %p", pObject, fpDataObject);
    return -EINVAL;
  }

  delete fpDataObject;
  fpDataObject=NULL;
  ReleaseDataBuffer(fpObjectBuffer);
  fpObjectBuffer=NULL;
  fObjectBufferSize=0;

  return 0;
}

void AliHLTOUT::SetEventId(AliHLTUInt64_t id)
{
  // set event id
  if (fCurrentEventId!=kAliHLTVoidEventID && fCurrentEventId!=id) {
    fLog.LoggingVarargs(kHLTLogWarning, "AliHLTOUT", "SetEventId" , __FILE__ , __LINE__ , "event id was already set to 0x%llx, setting now to 0x%llx", fCurrentEventId, id);
  }
  fCurrentEventId=id;
}

void AliHLTOUT::Print(const char* option) const
{
  // print info
  {
    for (AliHLTOUTBlockDescriptorVector::const_iterator  i=fBlockDescList.begin();
	 i!=fBlockDescList.end(); i++)
      i->Print(option);
  }
  {
    for (AliHLTOUTHandlerListEntryVector::const_iterator  i=fDataHandlers.begin();
	 i!=fDataHandlers.end(); i++)
      i->Print(option);
  }
}

void AliHLTOUT::AliHLTOUTBlockDescriptor::Print(const char* /*option*/) const
{
  // print info
  stringstream sout;
  sout << "AliHLTOUTBlockDescriptor index 0x" << setfill('0') << setw(8) << hex << right << fIndex 
       << ":   " << AliHLTComponent::DataType2Text(fDataType).c_str()
       << "  0x" << setfill('0') << setw(8) << hex << fSpecification
       << "  processed " << dec << fProcessed;
  cout << sout.str() << endl;
}

void AliHLTOUT::AliHLTOUTHandlerListEntry::Print(const char* /*option*/) const
{
  // print info
  stringstream sout;
  AliHLTModuleAgent::AliHLTOUTHandlerType type=AliHLTModuleAgent::kUnknownOutput;
  AliHLTComponentDataType dt=kAliHLTVoidDataType;
  if (this->fpHandlerDesc) {
    type=*(this->fpHandlerDesc);
    dt=*(this->fpHandlerDesc);
  }
  const char* stype="";
  switch(type) {
  case AliHLTModuleAgent::kEsd: stype="ESD"; break;
  case AliHLTModuleAgent::kRawReader: stype="RawReader"; break;
  case AliHLTModuleAgent::kRawStream: stype="RawStream"; break;
  case AliHLTModuleAgent::kChain: stype="Chain"; break;
  case AliHLTModuleAgent::kProprietary: stype="Proprietary"; break;
  default: stype="unknown";
  }
  sout << "HLTOUT handler: "
       << "  " << type << " (" << stype << ")"
       << "  " << AliHLTComponent::DataType2Text(dt).c_str();
  cout << sout.str() << endl;
}
