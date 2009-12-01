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

/** @file   AliHLTRootSchemaEvolutionComponent.cxx
    @author Matthias Richter
    @date   2009-10-18
    @brief  Handler component for ROOT schema evolution of streamed objects
*/

#include "AliHLTRootSchemaEvolutionComponent.h"
#include "AliHLTMessage.h"
#include "TObjArray.h"
#include "TStreamerInfo.h"
#include "TList.h"
#include "TFile.h"

#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRootSchemaEvolutionComponent)

AliHLTRootSchemaEvolutionComponent::AliHLTRootSchemaEvolutionComponent()
  : AliHLTProcessor()
  , fFlags(0)
  , fpStreamerInfos(NULL)
  , fFXSPrescaler(0)
  , fFileName()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}
const char* AliHLTRootSchemaEvolutionComponent::fgkConfigurationObject=NULL;

AliHLTRootSchemaEvolutionComponent::~AliHLTRootSchemaEvolutionComponent()
{
  // see header file for class documentation
  if (fpStreamerInfos) {
    fpStreamerInfos->Clear();
    delete fpStreamerInfos;
  }
  fpStreamerInfos=NULL;
}

void AliHLTRootSchemaEvolutionComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTRootSchemaEvolutionComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeStreamerInfo;
}

void AliHLTRootSchemaEvolutionComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation

  // this is nothing more than an assumption, in fact it's very difficult to predict how
  // much output the component produces
  constBase=100*1024;
  inputMultiplier=3;
}

int AliHLTRootSchemaEvolutionComponent::DoInit(int argc, const char** argv)
{
  // see header file for class documentation

  int iResult=0;

  // default configuration from CDB
  if (iResult>=0 && fgkConfigurationObject!=NULL) iResult=ConfigureFromCDBTObjString(fgkConfigurationObject);

  // custom configuration from command line arguments
  if (iResult>=0 && argc>0) iResult=ConfigureFromArgumentString(argc, argv);

  if (iResult>=0) {
    fpStreamerInfos=new TObjArray();
    if (!fpStreamerInfos) iResult=-ENOMEM;
  }
  return 0;
}

int AliHLTRootSchemaEvolutionComponent::DoDeinit()
{
  // see header file for class documentation
  if (fFileName.IsNull()==0) {
    WriteToFile(fFileName, fpStreamerInfos);
    fFileName.Clear();
  }

  if (fpStreamerInfos) {
    fpStreamerInfos->Clear();
    delete fpStreamerInfos;
  }
  fpStreamerInfos=NULL;

  return 0;
}

int AliHLTRootSchemaEvolutionComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/,
						 AliHLTComponentTriggerData& /*trigData*/ )
{
  // see header file for class documentation
  int iResult=0;
  AliHLTUInt32_t eventType=gkAliEventTypeUnknown;
  if (!IsDataEvent(&eventType) && 
      eventType==gkAliEventTypeStartOfRun) {
    return 0;
  }

  AliHLTMessage msg(kMESS_OBJECT);
  msg.EnableSchemaEvolution();
  for (const TObject* pObj=GetFirstInputObject();
       pObj && iResult>=0;
       pObj=GetNextInputObject()) {
    msg.WriteObject(pObj);
    iResult=UpdateStreamerInfos(msg.GetStreamerInfos(), fpStreamerInfos);
  }

  if (iResult>=0) {
    if ((TestBits(kHLTOUTatFirstEvent) && GetEventCount()==0) ||
	(TestBits(kHLTOUTatAllEvents)) || 
	(TestBits(kHLTOUTatEOR) && eventType==gkAliEventTypeEndOfRun)) {
      PushBack(fpStreamerInfos, kAliHLTDataTypeStreamerInfo);
    }
  }

  if (TestBits(kFXS) && eventType==gkAliEventTypeEndOfRun) {
    // push to FXS, needs to be implemented in the base class
  }

  if (fFileName.IsNull()==0 && eventType==gkAliEventTypeEndOfRun) {
    WriteToFile(fFileName, fpStreamerInfos);
    fFileName.Clear();
  }

  return iResult;
}

int AliHLTRootSchemaEvolutionComponent::UpdateStreamerInfos(const TList* list, TObjArray* infos) const
{
  // see header file for class documentation
  int iResult=0;
  if (!list || !infos) {
    return -EINVAL;
  }

  TObject* element=NULL;
  TIter next((TList*)list);
  while ((element = next())) {
    TStreamerInfo* pInfo=dynamic_cast<TStreamerInfo*>(element);
    if (!pInfo) continue;
    TString name=pInfo->GetName();
    int i=0;
    for (; i<infos->GetEntriesFast(); i++) {
      if (name.CompareTo(infos->At(i)->GetName())==0 &&
	  pInfo->GetClassVersion() == ((TStreamerInfo*)infos->At(i))->GetClassVersion()) {
	// have it already
	break;
      }
    }

    // Add streamer info
    infos->Add(pInfo);
  }

  return iResult;
}

int AliHLTRootSchemaEvolutionComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // see header file for class documentation
  int iResult=0;
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -hltout=[all,first,eor,off]
  if (argument.Contains("-hltout")) {
    argument.ReplaceAll("-hltout", "");
    argument.ReplaceAll("=", "");
    if (argument.IsNull() || argument.CompareTo("all")) {
      SetBits(kHLTOUTatAllEvents);
    } else if (argument.CompareTo("first")) {
      SetBits(kHLTOUTatFirstEvent);
    } else if (argument.CompareTo("eor")) {
      SetBits(kHLTOUTatEOR);
    } else if (argument.CompareTo("off")) {
      ClearBits(kHLTOUTatAllEvents | kHLTOUTatFirstEvent | kHLTOUTatEOR);
    } else {
      HLTError("invalid parameter for argument -hltout= : %s", argument.Data());
      return -EINVAL;
    }
    return 1;
  }

  // -fxs=[n,off]
  if (argument.Contains("-fxs")) {
    argument.ReplaceAll("-fxs", "");
    argument.ReplaceAll("=", "");
    if (argument.IsNull()) {
      SetBits(kFXS);
    } else if (argument.CompareTo("off")) {
      ClearBits(kFXS);
    } else if (argument.IsDigit()) {
      fFXSPrescaler=argument.Atoi();
    } else {
      HLTError("invalid parameter for argument -fxs= : %s", argument.Data());
      return -EINVAL;
    }
    return 1;
  }

  // -file=<filename>
  if (argument.Contains("-file=")) {
    argument.ReplaceAll("-file=", "");
    if (!argument.IsNull()) {
      fFileName=argument;
    } else {
      HLTError("argument -file= expects file name");
      return -EINVAL;
    }
    return 1;
  }
  
  return iResult;
}

int AliHLTRootSchemaEvolutionComponent::WriteToFile(const char* filename, const TObjArray* infos) const
{
  // write aray of streamer infos to file
  if (!filename || !infos) return -EINVAL;

  TFile out(filename, "RECREATE");
  if (out.IsZombie()) {
    HLTError("failed to open file %s", filename);
    return -EBADF;
  }

  const char* entrypath="HLT/Calib/StreamerInfo";
  int version = -1;
  AliCDBStorage* store = AliCDBManager::Instance()->GetDefaultStorage();
  if (store) {
    version = store->GetLatestVersion(entrypath, GetRunNo());
  }
  version++;

  AliCDBPath cdbPath(entrypath);
  AliCDBId cdbId(cdbPath, GetRunNo(), AliCDBRunRange::Infinity(), version, 0);
  AliCDBMetaData* cdbMetaData=new AliCDBMetaData;
  cdbMetaData->SetResponsible("ALICE HLT");
  cdbMetaData->SetComment("Streamer info for HLTOUT payload");
  AliCDBEntry* entry=new AliCDBEntry(infos->Clone(), cdbId, cdbMetaData, kTRUE);

  out.cd();
  entry->Write();
  // this is a small memory leak
  // seg fault in ROOT object handling if the two objects are deleted
  // investigate later
  //delete entry;
  //delete cdbMetaData;
  out.Close();

  return 0;
}
