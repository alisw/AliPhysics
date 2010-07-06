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

/// @file   AliHLTOUTHandlerEsdBranch.cxx
/// @author Matthias Richter
/// @date   01.07.2010
/// @brief  HLTOUT handler of type kEsd to merge objects into the hltEsd.

#include "AliHLTOUTHandlerEsdBranch.h"
#include "AliHLTOUT.h"
#include "AliHLTMessage.h"
#include "AliHLTErrorGuard.h"
#include "AliESDEvent.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TArrayC.h"
#include <cassert>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTHandlerEsdBranch)

AliHLTOUTHandlerEsdBranch::AliHLTOUTHandlerEsdBranch(const char* branchname)
  : AliHLTOUTHandler() 
  , fBranch(branchname)
  , fESD(NULL)
  , fpData(NULL)
  , fSize(0)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTOUTHandlerEsdBranch::~AliHLTOUTHandlerEsdBranch()
{
  // see header file for class documentation
}

int AliHLTOUTHandlerEsdBranch::ProcessData(AliHLTOUT* pData)
{
  // see header file for class documentation
  if (!pData) return -EINVAL;
  int iResult=0;

  if (CheckStatus(kHandlerError)) {
    HLTWarning("kEsd handler for ESD branch '%s' in error state, skipping processing of associated HLTOUT blocks", fBranch.Data());
    return -EPERM;
  }

  if (!fESD) {
    // create the ESD container, but without std content
    fESD = new AliESDEvent;
  }
  if (!fpData) fpData=new TArrayC;
  if (fESD && fpData) {
    fESD->Reset();
    iResult=ExtractAndAddObjects(pData);
  }

  AliHLTMessage* pMsg=AliHLTMessage::Stream(fESD);
  if (pMsg) {
    if (!pMsg->CompBuffer()) {
      fSize=pMsg->Length();
      fpData->Set(fSize, pMsg->Buffer());
    } else {
      fSize=pMsg->CompLength();
      fpData->Set(fSize, pMsg->CompBuffer());
    }
    delete pMsg;
    pMsg=NULL;
  } else {
    HLTError("streaming of object failed");
  }

  return iResult;
}

int AliHLTOUTHandlerEsdBranch::ExtractAndAddObjects(AliHLTOUT* pData)
{
  // Default method
  // Extract streamed object from the HLTOUT and add to ESD
  // The default method works only for single blocks in the HLTOUT,
  // A specific child class is required if multiple blocks should be
  // treated.

  int iResult=0;
  iResult=pData->SelectFirstDataBlock(); 
  if (iResult<0) return iResult;

  TObject* pObject=pData->GetDataObject();
  if (pObject) {
    TString bname=fBranch;
    if (bname.IsNull()) {
      bname=pObject->GetName();
      if (bname.CompareTo(pObject->ClassName())==0) {
	ALIHLTERRORGUARD(5, "no branch name specified for unnamed object %s, not added to ESD", bname.Data());
	bname="";
      }
    }
    if (!bname.IsNull()) {
      iResult=Add(pObject, bname.Data());
    } else {
      iResult=-EBADF;
    }
    pData->ReleaseDataObject(pObject);
    pObject=NULL;
  } else {
    AliHLTComponentDataType dt=kAliHLTVoidDataType;
    AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
    pData->GetDataBlockDescription(dt, spec);
    HLTError("Can not get TObject from HLTOUT buffer for block %s 0x%x", AliHLTComponent::DataType2Text(dt).c_str(), spec);
    iResult=-ENODATA;
  }

  if (pData->SelectNextDataBlock()>=0) {
    ALIHLTERRORGUARD(5, "the default function can only handle one single data block/object");
  }

  return iResult;
}

int AliHLTOUTHandlerEsdBranch::Add(TObject* object, const char* branchname)
{
  // add object to the specified branch in the internal ESD
  if (fESD && object) {
    TObject* pESDObject=fESD->FindListObject(branchname);
    if (pESDObject) {
      // copy the content to the already existing object
      object->Copy(*pESDObject);
    } else {
      // add a new object
      fESD->AddObject(object->Clone());
    }
  }
  return 0;
}

int AliHLTOUTHandlerEsdBranch::GetProcessedData(const AliHLTUInt8_t* &pData)
{
  // see header file for class documentation
  if (!fpData) {
    pData=NULL;
    return 0;
  }

  pData=reinterpret_cast<AliHLTUInt8_t*>(fpData->GetArray());
  return fSize;
}

int AliHLTOUTHandlerEsdBranch::ReleaseProcessedData(const AliHLTUInt8_t* pData, int size)
{
  // see header file for class documentation
  int iResult=0;
  if (!fpData || size != fSize ||
      const_cast<AliHLTUInt8_t*>(pData) != reinterpret_cast<AliHLTUInt8_t*>(fpData->GetArray())) {
    HLTError("attempt to release to wrong data buffer %p size %d, expected %p size %d", pData, size, fpData?fpData->GetArray():NULL, fSize);
  }
  fSize=0;
  return iResult;
}
