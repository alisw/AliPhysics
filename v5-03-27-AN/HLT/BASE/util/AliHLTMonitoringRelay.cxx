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

/// @file   AliHLTMonitoringRelay.cxx
/// @author Matthias Richter
/// @date   2009-11-11
/// @brief  Relay components for monitoring objects.
///

#include <cstdlib>
#include <cassert>
#include "AliHLTMonitoringRelay.h"
#include "AliHLTMessage.h"
#include "TArrayC.h"
#include "TObject.h"
#include "TDatime.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTMonitoringRelay)

AliHLTMonitoringRelay::AliHLTMonitoringRelay()
  : AliHLTProcessor()
  , fItems()
  , fOutputSize()
  , fFlags(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTMonitoringRelay::~AliHLTMonitoringRelay()
{
  // see header file for class documentation
}

void AliHLTMonitoringRelay::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTMonitoringRelay::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTAnyDataType;
}

void AliHLTMonitoringRelay::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase=fOutputSize;
  inputMultiplier=1.0;
}

int AliHLTMonitoringRelay::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;
  fOutputSize=0;

  iResult=ConfigureFromArgumentString(argc, argv);

  return iResult;
}

int AliHLTMonitoringRelay::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -verbose
  if (argument.CompareTo("-verbose")==0) {
    // 
    return 1;
  } else if (argument.CompareTo("-check-object")==0) { // check the objects in the blocks
    SetFlag(kCheckObject);
    return 1;
  }

  return 0;
}

int AliHLTMonitoringRelay::DoDeinit()
{
  // see header file for class documentation
  int iResult=0;
  AliHLTMonitoringItemPList::iterator element=fItems.begin();
  while (element!=fItems.end()) {
    AliHLTMonitoringItem* pItem=*element;
    element=fItems.erase(element);
    if (pItem) {
      delete pItem;
    }
  }
  return iResult;
}

int AliHLTMonitoringRelay::DoEvent(const AliHLTComponentEventData& /*evtData*/,
				   AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  int iResult=0;
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock();
       pBlock!=NULL; 
       pBlock=GetNextInputBlock()) {
    // ignore private blocks
    if (pBlock->fDataType==(kAliHLTAnyDataType|kAliHLTDataOriginPrivate)) continue;
    TObject* pObject=NULL;
    if (CheckFlag(kCheckObject)) pObject=AliHLTMessage::Extract(pBlock->fPtr, pBlock->fSize);
    
    AliHLTMonitoringItem* pItem=FindItem(pBlock, pObject);
    if (pItem) {
      HLTInfo("found block %s 0x%0lx %s %s", DataType2Text(pItem->GetDataType()).c_str(), pItem->GetSpecification(), pObject?pObject->GetName():"", pObject?pObject->GetName():"");
      if (pItem->GetSize()<pBlock->fSize) {
	// update with the new maximum
	assert(fOutputSize>=pItem->GetSize());
	fOutputSize-=pItem->GetSize();
	fOutputSize+=pBlock->fSize;
      }
      pItem->SetData(pBlock->fPtr, pBlock->fSize);
      HLTInfo("setting item size %d, total size %d", pItem->GetSize(), fOutputSize);
    } else {
      pItem=new AliHLTMonitoringItem(pBlock, pObject);
      fItems.push_back(pItem);
      fOutputSize+=pBlock->fSize;
      HLTInfo("new item size %d (%d),  %s 0x%0lx %s %s", pItem->GetSize(), fOutputSize, DataType2Text(pItem->GetDataType()).c_str(), pItem->GetSpecification(), pObject?pObject->GetName():"", pObject?pObject->GetName():"");
    }
    if (pObject) delete pObject;
  }

  int nofObjects=0;
  for (AliHLTMonitoringItemPList::iterator element=fItems.begin();
       element!=fItems.end(); element++) {
    AliHLTMonitoringItem* pItem=*element;
    if (pItem) {
      HLTInfo("push back item size %d (%d),  %s 0x%0lx", pItem->GetSize(), fOutputSize, DataType2Text(pItem->GetDataType()).c_str(), pItem->GetSpecification());
      PushBack(pItem->GetBuffer(), pItem->GetSize(), pItem->GetDataType(), pItem->GetSpecification());
      if (!pItem->GetObjectName().IsNull()) nofObjects++;
    }
  }
  // info output once every 5 seconds
  const TDatime time;
  static UInt_t lastTime=0;
  if (time.Get()-lastTime>5) {
      lastTime=time.Get();
      HLTBenchmark("accumulated %d items containing %d TObjects", fItems.size(), nofObjects);
  }

  return iResult;
}

AliHLTMonitoringRelay::AliHLTMonitoringItem* AliHLTMonitoringRelay::FindItem(const AliHLTComponentBlockData* pBlock, const TObject* pObject) const
{
  // find an item by data type, specification, name and title
  for (unsigned i=0; i<fItems.size(); i++) {
    AliHLTMonitoringItem* pItem=fItems[i];
    if (pItem &&
	(*pItem)==(*pBlock) &&
	(pObject==NULL || (*pItem)==(*pObject))) {
      return pItem;
    }
  }
  return NULL;  
}

AliHLTMonitoringRelay::AliHLTMonitoringItem::AliHLTMonitoringItem()
  : fDt(kAliHLTVoidDataType)
  , fSpecification(kAliHLTVoidDataSpec)
  , fName()
  , fTitle()
  , fData(new TArrayC)
  , fDataSize(0)
{
  // standard constructor
}

AliHLTMonitoringRelay::AliHLTMonitoringItem::AliHLTMonitoringItem(const AliHLTComponentBlockData* pBlock, const TObject* pObject)
  : fDt(kAliHLTVoidDataType)
  , fSpecification(kAliHLTVoidDataSpec)
  , fName()
  , fTitle()
  , fData(new TArrayC)
  , fDataSize(0)
{
  // constructor
  if (pBlock) {
    fDt=pBlock->fDataType;
    fSpecification=pBlock->fSpecification;
    if (fData) {
      fData->Set(pBlock->fSize, reinterpret_cast<const Char_t*>(pBlock->fPtr));
      fDataSize=pBlock->fSize;
    }
  }
  if (pObject) {
    fName=pObject->GetName();
    fTitle=pObject->GetTitle();
  }
}

AliHLTMonitoringRelay::AliHLTMonitoringItem::~AliHLTMonitoringItem()
{
  // desstructor
  if (fData) delete fData;
}

int AliHLTMonitoringRelay::AliHLTMonitoringItem::SetData(void* pBuffer, int size)
{
  // copy the data buffer
  if (!fData) {
    fData=new TArrayC(size, reinterpret_cast<const Char_t*>(pBuffer));
  }
  if (!fData) {
    return -ENOMEM;
  }

  if (fData->GetSize()<size) {
    fData->Set(size, reinterpret_cast<const Char_t*>(pBuffer));
  } else {
    memcpy(fData->GetArray(), pBuffer, size);
  }

  fDataSize=size;
  return 0;
}


void* AliHLTMonitoringRelay::AliHLTMonitoringItem::GetBuffer() const
{
  // get buffer pointer of the current data
  return fData!=NULL?fData->GetArray():NULL;
}

unsigned AliHLTMonitoringRelay::AliHLTMonitoringItem::GetSize() const
{
  // get size of the current data
  return fDataSize;
}

const AliHLTComponentDataType& AliHLTMonitoringRelay::AliHLTMonitoringItem::GetDataType() const
{
  // get data type
  return fDt;
}

AliHLTUInt32_t AliHLTMonitoringRelay::AliHLTMonitoringItem::GetSpecification() const
{
  // get specification
  return fSpecification;
}

bool AliHLTMonitoringRelay::AliHLTMonitoringItem::operator==(const AliHLTComponentBlockData& bd) const
{
  // equal to data type and specification
  if (bd.fDataType!=fDt) return false;
  if (bd.fSpecification!=fSpecification) return false;
  return true;
}

bool AliHLTMonitoringRelay::AliHLTMonitoringItem::operator!=(const AliHLTComponentBlockData& bd) const
{
  // not equal to data type and specification
  return not operator==(bd);
}

bool AliHLTMonitoringRelay::AliHLTMonitoringItem::operator==(const TObject& object) const
{
  // equal to name and title
  if (fName.CompareTo(object.GetName())!=0) return false;
  if (!fTitle.IsNull() && fTitle.CompareTo(object.GetTitle())!=0) return false;
  return true;
}

bool AliHLTMonitoringRelay::AliHLTMonitoringItem::operator!=(const TObject& object) const
{
  // not equal to name and title
  return not operator==(object);
}
