// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no         *
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

//  @file   AliHLTGlobalPreprocessor.cxx
//  @author Matthias Richter
//  @date   2010-08-20
//  @brief  HLT Preprocessor plugin for global HLT
// 

#include "AliHLTGlobalPreprocessor.h"
#include "AliHLTRootSchemaEvolutionComponent.h"
#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliGRPObject.h"
#include "AliCDBMetaData.h"
#include "AliCDBId.h"
#include "AliLog.h"
#include <cassert>
#include <cerrno>

#include <TObjArray.h>
#include <TStreamerInfo.h>
//#include <AliDCSValue.h>
//#include <TTimeStamp.h>


ClassImp(AliHLTGlobalPreprocessor)

AliHLTGlobalPreprocessor::AliHLTGlobalPreprocessor()
  : AliHLTModulePreprocessor()
{
  // constructor
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

const char* AliHLTGlobalPreprocessor::fgkStreamerInfoAlias="StreamerInfo";
const char* AliHLTGlobalPreprocessor::fgkStreamerInfoName="StreamerInfo";
const char* AliHLTGlobalPreprocessor::fgkStreamerInfoType="Calib";

AliHLTGlobalPreprocessor::~AliHLTGlobalPreprocessor()
{
  // destructor
}


void AliHLTGlobalPreprocessor::Initialize(Int_t /*run*/, UInt_t /*startTime*/, 
					  UInt_t /*endTime*/)
{
  // initializes AliHLTGlobalPreprocessor

}


UInt_t AliHLTGlobalPreprocessor::Process(TMap* dcsAliasMap)
{
  // processes the DCS value map
  
  if (!dcsAliasMap) return -EINVAL;
  if (dcsAliasMap->GetEntries() == 0 ) return 0;
  
  TObject* streamerinfo=GetFromMap(dcsAliasMap, fgkStreamerInfoAlias);
  if (streamerinfo) ProcessStreamerInfo(streamerinfo);

  return 0;
}

Int_t AliHLTGlobalPreprocessor::GetModuleNumber() {
  Int_t modulenumber = 0;
  //modulenumber = AliHLTModulePreprocessor::DetectorBitMask("GRP");
  return modulenumber;
};

int AliHLTGlobalPreprocessor::ProcessStreamerInfo(TObject* object)
{
  /// process the StreamerInfo object
  int iResult=0;
  if (!object) return -EINVAL;

  TObjArray* streamerinfos=dynamic_cast<TObjArray*>(object);
  if (!streamerinfos) {
    AliError(Form("StreamerInfo object has wrong class type %s, expecting TObjArray", object->ClassName()));
    return -EINVAL;
  }
  if (streamerinfos->GetEntriesFast()==0) return 0;

  bool bStore=false;
  AliCDBEntry* entry = GetFromOCDB(fgkStreamerInfoType, fgkStreamerInfoName);
  TObjArray* clone=NULL;
  if (entry && entry->GetObject()) {
    TObject* cloneObj=entry->GetObject()->Clone();
    if (cloneObj) clone=dynamic_cast<TObjArray*>(cloneObj);
    bStore=AliHLTRootSchemaEvolutionComponent::MergeStreamerInfo(clone, streamerinfos)>0;
  } else {
    TObject* cloneObj=streamerinfos->Clone();
    if (cloneObj) clone=dynamic_cast<TObjArray*>(cloneObj);
    bStore=true;
  }

  if (clone) {
    AliCDBMetaData* metaData=entry?entry->GetMetaData():NULL;
    AliCDBMetaData* newMetaData=NULL;
    if (!metaData) {
      newMetaData=new AliCDBMetaData;
      if (newMetaData) {
	metaData=newMetaData;
	metaData->SetBeamPeriod(0);
	metaData->SetResponsible("ALICE HLT Matthias.Richter@cern.ch");
	metaData->SetComment("Streamer info for HLTOUT payload");
	//metaData->SetAliRootVersion(ALIROOT_SVN_BRANCH);
      } else {
	iResult=-ENOMEM;
      }
    }
    Store(fgkStreamerInfoType, fgkStreamerInfoName, clone, metaData, GetRun(), kFALSE);
    delete clone;
    if (newMetaData) delete newMetaData;
    newMetaData=NULL;
    metaData=NULL;
    // TODO
    // - what to do with variable 'entry', to be deleted?

  } else {
    AliError("failed to clone streamer info object array");
    return -ENOMEM;
  }

  return iResult;
}
