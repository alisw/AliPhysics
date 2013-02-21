// $Id$

//**************************************************************************
//* This file is property of and copyright by the                          * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTRootSchemaEvolutionComponent.cxx
/// @author Matthias Richter
/// @date   2009-10-18
/// @brief  Handler component for ROOT schema evolution of streamed objects
///

#include "AliHLTRootSchemaEvolutionComponent.h"
#include "AliHLTMessage.h"
#include "AliHLTReadoutList.h"
#include "AliHLTMisc.h"
#include "TObjArray.h"
#include "TStreamerInfo.h"
#include "TList.h"
#include "TFile.h"
#include "TStopwatch.h"
#include "TTimeStamp.h"
#include "TDatime.h"

#include "AliCDBStorage.h"
#include "AliCDBManager.h"
#include "AliCDBPath.h"
#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCDBEntry.h"

#include <numeric>
using namespace std;

namespace
{
  // Helper class for std::accumulate algorithm.
  class AliTimeSum {
  public:
    typedef int first_argument_type;
    typedef AliHLTRootSchemaEvolutionComponent::AliHLTDataBlockItem second_argument_type;
    typedef bool result_type;
    int operator() (int a, AliHLTRootSchemaEvolutionComponent::AliHLTDataBlockItem b) {
      return a+b.GetTotalTime();
    }
  };
} // end of namespace

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTRootSchemaEvolutionComponent)

AliHLTRootSchemaEvolutionComponent::AliHLTRootSchemaEvolutionComponent()
  : AliHLTCalibrationProcessor()
  , fList()
  , fPropertyFlags(kFXS)
  , fpStreamerInfos(NULL)
  , fpEventTimer(NULL)
  , fpCycleTimer(NULL)
  , fMaxEventTime(500)
  , fFXSPrescaler(0)
  , fFileName()
{
  // Collects streamer info for all input objects and produces the corresponding
  // calibration object for reconstruction of HLT. The component runs with a
  // configurable rate constraint and skips the processing of known data blocks
  // for the sake of performance. New data blocks are always processed and added
  // to the list.
  //
  // Component ID: \b ROOTSchemaEvolutionComponent                        <br>
  // Library: \b libAliHLTUtil.so						<br>
  // Input Data Types: ::kAliHLTAnyDataType				<br>
  // Output Data Types: none						<br>
}

// FIXME: read below when defining an OCDB object here
const char* AliHLTRootSchemaEvolutionComponent::fgkConfigurationObject=NULL;
const AliHLTUInt32_t AliHLTRootSchemaEvolutionComponent::fgkTimeScale=1000000;

AliHLTRootSchemaEvolutionComponent::~AliHLTRootSchemaEvolutionComponent()
{
  // destructor
  if (fpStreamerInfos) {
    fpStreamerInfos->Clear();
    delete fpStreamerInfos;
  }
  fpStreamerInfos=NULL;
}

void AliHLTRootSchemaEvolutionComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // overloaded from AliHLTComponent
  list.push_back(kAliHLTAnyDataType);
}

AliHLTComponentDataType AliHLTRootSchemaEvolutionComponent::GetOutputDataType()
{
  // overloaded from AliHLTComponent
  return kAliHLTDataTypeStreamerInfo;
}

void AliHLTRootSchemaEvolutionComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // overloaded from AliHLTComponent

  // this is nothing more than an assumption, in fact it's very difficult to predict how
  // much output the component produces
  constBase=100*1024;
  inputMultiplier=3;
}

int AliHLTRootSchemaEvolutionComponent::InitCalibration()
{
  // overloaded from AliHLTCalibrationProcessor: initialization

  int iResult=0;

  // default configuration from CDB
  // FIXME: has to be called from AliHLTCalibrationProcessor::DoInit in order to set
  // the default parameters from OCDB before the custom argument scan
  // not valid at the moment because fgkConfigurationObject==NULL
  if (iResult>=0 && fgkConfigurationObject!=NULL) iResult=ConfigureFromCDBTObjString(fgkConfigurationObject);

  if (iResult>=0) {
    fpStreamerInfos=new TObjArray();
    if (!fpStreamerInfos) iResult=-ENOMEM;

    fpEventTimer=new TStopwatch;
    if (fpEventTimer) {
      fpEventTimer->Reset();
    }
    fpCycleTimer=new TStopwatch;
    if (fpCycleTimer) {
      fpCycleTimer->Reset();
    }
  }

  return 0;
}

int AliHLTRootSchemaEvolutionComponent::DeinitCalibration()
{
  // overloaded from AliHLTCalibrationProcessor: termination and cleanup
  if (fFileName.IsNull()==0) {
    WriteToFile(fFileName, fpStreamerInfos);
    fFileName.Clear();
  }

  if (fpStreamerInfos) {
    fpStreamerInfos->Clear();
    delete fpStreamerInfos;
  }
  fpStreamerInfos=NULL;

  if (fpEventTimer) {
    delete fpEventTimer;
    fpEventTimer=NULL;
  }
  if (fpCycleTimer) {
    delete fpCycleTimer;
    fpCycleTimer=NULL;
  }
  return 0;
}

int AliHLTRootSchemaEvolutionComponent::ProcessCalibration( const AliHLTComponentEventData& /*evtData*/,
							    AliHLTComponentTriggerData& /*trigData*/ )
{
  // overloaded from AliHLTCalibrationProcessor: event processing
  int iResult=0;
  AliHLTUInt32_t eventType=gkAliEventTypeUnknown;
  if (!IsDataEvent(&eventType) && 
      eventType==gkAliEventTypeStartOfRun) {
    return 0;
  }

  AliHLTUInt32_t listtime=accumulate(fList.begin(), fList.end(), int(0), AliTimeSum());
  AliHLTUInt32_t averageEventTime=0;
  AliHLTUInt32_t averageCycleTime=0;

  AliHLTUInt32_t proctime=0;
  if (fpEventTimer) {
    averageEventTime=AliHLTUInt32_t(fpEventTimer->RealTime()*fgkTimeScale)/(GetEventCount()+1);
    proctime=AliHLTUInt32_t(fpEventTimer->RealTime()*fgkTimeScale);
    fpEventTimer->Start(kFALSE);
  }
  if (fpCycleTimer) {
    fpCycleTimer->Stop();
    averageCycleTime=AliHLTUInt32_t(fpCycleTimer->RealTime()*fgkTimeScale)/(GetEventCount()+1);
  }

  // scale down the event processing according to the required rate
  // and average processing time.
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock();
       pBlock && iResult>=0;
       pBlock=GetNextInputBlock()) {
    bool processBlock=true;
    AliHLTDataBlockItem* item=FindItem(pBlock->fDataType, pBlock->fSpecification);
    if (item) {
      // TODO: do a selection of blocks on basis of the time spent in its processing
      // for now only the global processing time is checked
      // process if the average event time is smaller then the cycle time, i.e.
      // the time is spent outside the component
      // apply a factor 4 margin
      processBlock=4*averageEventTime<fMaxEventTime || 2*averageEventTime<averageCycleTime;
    } else {
      // always process new incoming blocks
      processBlock=true;
      fList.push_back(AliHLTDataBlockItem(pBlock->fDataType, pBlock->fSpecification));
      item=&fList[fList.size()-1];
    }
    if (processBlock) {
      TObject* pObj=item->Extract(pBlock);
      if (pObj) {
	AliHLTMessage msg(kMESS_OBJECT);
	msg.EnableSchemaEvolution();
	if ((iResult=item->Stream(pObj, msg))>=0) {
	  iResult=UpdateStreamerInfos(msg.GetStreamerInfos(), fpStreamerInfos);
	} else {
	  HLTError("failed to stream object %s of type %s", pObj->GetName(), pObj->ClassName());
	}
	delete pObj;
	pObj=NULL;
      }
    }
  }

  if (fpEventTimer) {
    fpEventTimer->Stop();
    proctime=AliHLTUInt32_t(fpEventTimer->RealTime()*fgkTimeScale)-proctime;
    averageEventTime=AliHLTUInt32_t(fpEventTimer->RealTime()*fgkTimeScale)/(GetEventCount()+1);

    // info output once every 2 seconds
    static UInt_t lastTime=0;
    TDatime time;
    if (time.Get()-lastTime>2) {
      lastTime=time.Get();
      HLTInfo("event time %d, average time %d, list time %d, cycle time %d", proctime, averageEventTime, listtime, averageCycleTime);
    }
  }
  if (fpCycleTimer) {
    fpCycleTimer->Start(kFALSE);
  }

  if (iResult>=0) {
    if ((TestBits(kHLTOUTatFirstEvent) && GetEventCount()==0) ||
	(TestBits(kHLTOUTatAllEvents))) {
      PushBack(fpStreamerInfos, kAliHLTDataTypeStreamerInfo);
    }
  }

  if (TestBits(kFXS) && fFXSPrescaler>0 && (GetEventCount()%fFXSPrescaler)==0) {
    // push to FXS
    AliHLTReadoutList rdList(AliHLTReadoutList::kHLT);
    PushToFXS((TObject*)fpStreamerInfos, "HLT", "StreamerInfo", &rdList );
  }

  return iResult;
}

int AliHLTRootSchemaEvolutionComponent::ShipDataToFXS( const AliHLTComponentEventData& /*evtData*/,
						       AliHLTComponentTriggerData& /*trigData*/)
{
  // overloaded from AliHLTCalibrationProcessor: ship data
  if (TestBits(kFXS)) {
    // push to FXS
    AliHLTReadoutList rdList(AliHLTReadoutList::kHLT);
    PushToFXS((TObject*)fpStreamerInfos, "HLT", "StreamerInfo", &rdList );
  }

  if (fFileName.IsNull()==0) {
    WriteToFile(fFileName, fpStreamerInfos);
    fFileName.Clear();
  }

    if (TestBits(kHLTOUTatEOR)) {
      PushBack(fpStreamerInfos, kAliHLTDataTypeStreamerInfo);
    }

    for (unsigned i=0; i<fList.size(); i++) {
      if (CheckFilter(kHLTLogDebug)) fList[i].Print("short");
      else if (fList[i].IsObject()) {
	HLTInfo("AliHLTDataBlockItem %s %08x\n"
		"   average extraction time: %d usec\n"
		"   average streaming time: %d usec"
		, AliHLTComponent::DataType2Text(fList[i]).c_str()
		, fList[i].GetSpecification()
		, fList[i].GetExtractionTime()
		, fList[i].GetStreamingTime());
      }
    }

  return 0;
}

int AliHLTRootSchemaEvolutionComponent::UpdateStreamerInfos(const TList* list, TObjArray* infos) const
{
  // update streamer infos
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
    if (pInfo->GetClassVersion()==0) continue; // skip classes which are not for storage
    for (; i<infos->GetEntriesFast(); i++) {
      if (name.CompareTo(infos->At(i)->GetName())==0 &&
	  pInfo->GetClassVersion() == ((TStreamerInfo*)infos->At(i))->GetClassVersion()) {
	// have it already
	break;
      }
    }

    // Add streamer info if not yet there
    if (i>=infos->GetEntriesFast()) {
      infos->Add(pInfo);
    }
  }

  return iResult;
}

int AliHLTRootSchemaEvolutionComponent::ScanConfigurationArgument(int argc, const char** argv)
{
  // overloaded from AliHLTComponent
  int iResult=0;
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -hltout=[all,first,eor,off]
  if (argument.Contains("-hltout")) {
    argument.ReplaceAll("-hltout", "");
    argument.ReplaceAll("=", "");
    if (argument.IsNull() || argument.CompareTo("all")==0) {
      SetBits(kHLTOUTatAllEvents|kHLTOUTatEOR);
    } else if (argument.CompareTo("first")==0) {
      SetBits(kHLTOUTatFirstEvent);
    } else if (argument.CompareTo("eor")==0) {
      SetBits(kHLTOUTatEOR);
    } else if (argument.CompareTo("off")==0) {
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
    SetBits(kFXS);
    if (argument.IsNull()) {
    } else if (argument.CompareTo("off")==0) {
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

  if (argument.Contains("-rate=")) {
    argument.ReplaceAll("-rate=", "");
    AliHLTUInt32_t rate=argument.Atoi();
    if (rate>0 && rate<fgkTimeScale) {
      fMaxEventTime=fgkTimeScale/rate;
    } else {
      HLTError("argument -file= expects number [Hz]");
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
  AliCDBEntry* existingEntry=AliHLTMisc::Instance().LoadOCDBEntry(entrypath);
  if (existingEntry) {
    version=existingEntry->GetId().GetVersion();
  }
  version++;

  TObjArray* clone=NULL;

  if (existingEntry && existingEntry->GetObject()) {
    TObject* cloneObj=existingEntry->GetObject()->Clone();
    if (cloneObj) clone=dynamic_cast<TObjArray*>(cloneObj);
    if (AliHLTMisc::Instance().MergeStreamerInfo(clone, infos)==0) {
      // no change, store with identical version
      version=existingEntry->GetId().GetVersion();
    }
  } else {
    TObject* cloneObj=infos->Clone();
    if (cloneObj) clone=dynamic_cast<TObjArray*>(cloneObj);
  }
  if (!clone) {
    HLTError("failed to clone streamer info object array");
    return -ENOMEM;
  }

  AliCDBPath cdbPath(entrypath);
  AliCDBId cdbId(cdbPath, AliCDBManager::Instance()->GetRun(), AliCDBRunRange::Infinity(), version, 0);
  AliCDBMetaData* cdbMetaData=new AliCDBMetaData;
  cdbMetaData->SetResponsible("ALICE HLT Matthias.Richter@cern.ch");
  cdbMetaData->SetComment("Streamer info for HLTOUT payload");
  AliCDBEntry* entry=new AliCDBEntry(clone, cdbId, cdbMetaData, kTRUE);

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

AliHLTRootSchemaEvolutionComponent::AliHLTDataBlockItem*
AliHLTRootSchemaEvolutionComponent::FindItem(AliHLTComponentDataType dt,
					     AliHLTUInt32_t spec)
{
  /// find item in the list
  // vector<AliHLTDataBlockItem>::iterator element=std::find(fList.begin(), fList.end(), AliHLTDataBlockItem(dt,spec));
  // if (element!=fList.end()) return &(*element);
  for (unsigned i=0; i<fList.size(); i++) {
    if (fList[i]==dt && fList[i]==spec) return &fList[i];
  }
  return NULL;
}

AliHLTRootSchemaEvolutionComponent::AliHLTDataBlockItem::AliHLTDataBlockItem(AliHLTComponentDataType dt,
									     AliHLTUInt32_t spec)
  : fDt(dt)
  , fSpecification(spec)
  , fIsObject(false)
  , fNofExtractions(0)
  , fExtractionTimeUsec(0)
  , fLastExtraction(0)
  , fNofStreamings(0)
  , fStreamingTimeUsec(0)
  , fLastStreaming(0)
{
  // helper class to keep track of input data blocks
  // in the AliHLTRootSchemaEvolutionComponent
  //
}

AliHLTRootSchemaEvolutionComponent::AliHLTDataBlockItem::~AliHLTDataBlockItem()
{
  // destructor
}

TObject* AliHLTRootSchemaEvolutionComponent::AliHLTDataBlockItem::Extract(const AliHLTComponentBlockData* bd)
{
  /// extract data block to root object, and update performance parameters
  /// object needs to be deleted externally
  if (!bd || !bd->fPtr || bd->fSize<8) return NULL;

  AliHLTUInt32_t firstWord=*((AliHLTUInt32_t*)bd->fPtr);
  if (!(fIsObject=(firstWord==bd->fSize-sizeof(AliHLTUInt32_t)))) return NULL;

  TStopwatch sw;
  sw.Start();
  AliHLTMessage msg(bd->fPtr, bd->fSize);
  TClass* objclass=msg.GetClass();
  if (!(fIsObject=(objclass!=NULL))) return NULL;
  TObject* pObj=msg.ReadObject(objclass);
  if (!(fIsObject=(pObj!=NULL))) return NULL;
  sw.Stop();
  AliHLTUInt32_t usec=AliHLTUInt32_t(sw.RealTime()*fgkTimeScale);
  fNofExtractions++;
  fExtractionTimeUsec+=usec;
  TTimeStamp ts;
  fLastExtraction=(ts.GetSec()%1000)*fgkTimeScale + ts.GetNanoSec()/1000;
  return pObj;
}

int AliHLTRootSchemaEvolutionComponent::AliHLTDataBlockItem::Stream(const TObject* obj, AliHLTMessage& msg)
{
  /// stream object and update performance parameters
  if (!obj) return -EINVAL;
  TStopwatch sw;
  sw.Start();
  msg.WriteObject(obj);

  AliHLTUInt32_t usec=AliHLTUInt32_t(sw.RealTime()*fgkTimeScale);
  fNofStreamings++;
  fStreamingTimeUsec+=usec;
  TTimeStamp ts;
  fLastStreaming=(ts.GetSec()%1000)*fgkTimeScale + ts.GetNanoSec()/1000;
  return 0;
}

void AliHLTRootSchemaEvolutionComponent::AliHLTDataBlockItem::Print(const char* option) const
{
  /// print status
  if (fIsObject || !(strcmp(option, "short")==0))
    cout << "AliHLTDataBlockItem: " << AliHLTComponent::DataType2Text(fDt).c_str() << " " << hex << fSpecification << dec << endl;
  if (fIsObject) {
    if (fNofExtractions>0) cout << "   average extraction time: " << fExtractionTimeUsec/fNofExtractions << " usec" << endl;
    else cout << "   never extracted" << endl;
    if (fNofStreamings>0) cout << "   average streaming time: " << fStreamingTimeUsec/fNofStreamings << " usec" << endl;
    else cout << "   never streamed" << endl;
  }
}
