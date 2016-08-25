//-*- Mode: C++ -*-
// $Id: AliHLTAnalysisManagerComponent.cxx $
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: David Rohr, Jens Wiechula, C. Zampolli, M. Krzewicki  *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file    AliHLTAnalysisManagerComponent.cxx
    @author  David Rohr, Jens Wiechula, C. Zampolli, I. Vorobyev, M.Krzewicki
    @brief   Runs AliAnalysisTasks in a HLT component
*/

#include "TMap.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TStopwatch.h"
#include "TObjString.h"
#include "TH1F.h"
#include "TList.h"

#include <map>
#include <string>
#include "TString.h"
#include "TArrayI.h"

#include "AliESDEvent.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTAnalysisManagerComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "AliHLTAnalysisManager.h"
#include "AliHLTVEventInputHandler.h"
#include "AliAnalysisDataContainer.h"
#include "TTree.h"
#include "TChain.h"
#include "AliFlatESDEvent.h"
#include "AliFlatESDFriend.h"
#include "AliVEvent.h"
#include "AliVfriendEvent.h"
#include "AliSysInfo.h"
#include "AliHLTObjArray.h"

#include "TPRegexp.h"

#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliGRPObject.h"

#include "TROOT.h"
#include "TTreeStream.h"

#include "TGeoManager.h"
#include "AliCDBPath.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliHLTLogging.h"
#include "AliHLTCTPData.h"
#include "TClass.h"
#include "TDataMember.h"

using namespace std;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTAnalysisManagerComponent)

/*
 * ---------------------------------------------------------------------------------
 *                            Constructor / Destructor
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
AliHLTAnalysisManagerComponent::AliHLTAnalysisManagerComponent() :
  AliHLTProcessor(),
  fUID(0),
  fQuickEndRun(false),
  fAnalysisInitialized(false),
  fAnalysisManager(NULL),
  fInputHandler(NULL),
  fAddTaskMacro(""),
  fWriteAnalysisToFile(kFALSE),
  fInitializeGeometry(kTRUE),
  fEnableDebug(kFALSE),
  fResetAfterPush(kTRUE),
  fPushEventModulo(0),
  fNEvents(0),
  fMinTracks(0),
  fNumEvents(0),
  fQueueDepth(0),
  fAsyncProcess(0),
  fForceKillAsyncProcess(-1),
  fPushRequestOngoing(kFALSE),
  fAsyncProcessor(),
  fAnalysisOutputContainer(NULL)
{
  // an example component which implements the ALICE HLT processor
  // interface and does some analysis on the input raw data
  //
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  //
  // NOTE: all helper classes should be instantiated in DoInit()
}

// #################################################################################
AliHLTAnalysisManagerComponent::~AliHLTAnalysisManagerComponent() {
  // see header file for class documentation
}

/*
 * ---------------------------------------------------------------------------------
 * Public functions to implement AliHLTComponent's interface.
 * These functions are required for the registration process
 * ---------------------------------------------------------------------------------
 */

// #################################################################################
const Char_t* AliHLTAnalysisManagerComponent::GetComponentID() {
  // see header file for class documentation
  return "HLTAnalysisManagerComponent";
}

// #################################################################################
void AliHLTAnalysisManagerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginAny);
  list.push_back(kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD);
  list.push_back(kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO);
  list.push_back(kAliHLTDataTypeESDfriendObject|kAliHLTDataOriginAny);
  list.push_back(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
  list.push_back(kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut);
  list.push_back(kAliHLTDataTypeEOR|kAliHLTDataOriginAny);
}

// #################################################################################
AliHLTComponentDataType AliHLTAnalysisManagerComponent::GetOutputDataType() {
  // see header file for class documentation
  return kAliHLTDataTypeTObject|kAliHLTDataOriginHLT;
}

// #################################################################################
void AliHLTAnalysisManagerComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
  // see header file for class documentation
  constBase = 600000; //fixed size output objects
  inputMultiplier = 0; //do not scale with the input size
}

// #################################################################################
void AliHLTAnalysisManagerComponent::GetOCDBObjectDescription( TMap* const targetMap) {
  // see header file for class documentation

  if (!targetMap) return;
  /*  targetMap->Add(new TObjString("HLT/ConfigGlobal/MultiplicityCorrelations"),
		 new TObjString("configuration object"));
  targetMap->Add(new TObjString("HLT/ConfigGlobal/MultiplicityCorrelationsCentrality"),
		 new TObjString("centrality configuration object"));
  */
  return;
}

// #################################################################################
AliHLTComponent* AliHLTAnalysisManagerComponent::Spawn() {
  // see header file for class documentation
  return new AliHLTAnalysisManagerComponent;
}

/*
 * ---------------------------------------------------------------------------------
 * Protected functions to implement AliHLTComponent's interface.
 * These functions provide initialization as well as the actual processing
 * capabilities of the component. 
 * ---------------------------------------------------------------------------------
 */

void* AliHLTAnalysisManagerComponent::AnalysisManagerInit(void*)
{
  if (!gGeoManager && fInitializeGeometry)
  {
    AliCDBPath path("GRP","Geometry","Data");
    if(path.GetPath())
    {
      //      HLTInfo("configure from entry %s", path.GetPath());
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
      if (pEntry) 
	{
	  gGeoManager = (TGeoManager*) pEntry->GetObject();
	  
	  if(!gGeoManager)
	    {
	       HLTFatal("can not get gGeoManager from OCDB");
	    }
	}
      else
	{
	    HLTFatal("can not fetch object \"%s\" from OCDB", path.GetPath().Data());
	}
    }
  }

  fAnalysisManager = new AliHLTAnalysisManager();
  fInputHandler    = new AliHLTVEventInputHandler("HLTinputHandler","HLT input handler");
  fAnalysisManager->SetInputEventHandler(fInputHandler);

  if (fAddTaskMacro.Length()>0) 
  {
    HLTInfo("Executing the macro: %s\n",fAddTaskMacro.Data());
    TString addTaskMacro(fAddTaskMacro);
    if (fQueueDepth!=0) { addTaskMacro+="\\;"; }
    gROOT->Macro(addTaskMacro);
  }

  if (fAnalysisManager->InitAnalysis() == kFALSE)
  {
    return((void*) -1);
  }

  ////this disables streaming of the fProducer and fConsumers data members of
  ////AliAnalysisDataContainer.
  ////It is needed to avoid memory leaks downstream.
  ////Proper fix in the container itself may be too dangerous
  //TClass::GetClass("AliAnalysisDataContainer")->
  //  GetDataMember("fProducer")->
  //  SetBit(BIT(2),0);
  //TClass::GetClass("AliAnalysisDataContainer")->
  //  GetDataMember("fConsumers")->
  //  SetBit(BIT(2),0);

  return(NULL);
}

// #################################################################################
Int_t AliHLTAnalysisManagerComponent::DoInit( Int_t /*argc*/, const Char_t** /*argv*/ ) {
  // see header file for class documentation
  HLTInfo("----> AliHLTAnalysisManagerComponent::DoInit");

  //process arguments
  ProcessOptionString(GetComponentArgs());

  if (!fEnableDebug)
  {
    //shut up AliSysInfo globally, prevents syswatch.log files from being created
    //dont open any debug files
    AliSysInfo::SetDisabled(kTRUE);
    TTreeSRedirector::SetDisabled(kTRUE);
  }

  HLTInfo("AliHLTAnalysisManagerComponent::DoInit (with QueueDepth %d)", fQueueDepth);
  if (fAsyncProcessor.Initialize(fQueueDepth, fAsyncProcess > 0, fAsyncProcess)) return(1);

  void* initRetVal;
  if (fAsyncProcessor.InitializeAsyncMemberTask(this, &AliHLTAnalysisManagerComponent::AnalysisManagerInit, NULL, &initRetVal) == 0)
  {
    if (initRetVal == 0) fAnalysisInitialized = kTRUE;
  }
  if (!fAnalysisInitialized)
  {
    HLTWarning("Error initializing analysis task");
  }

  //Init the CTP data
  if (SetupCTPData() == -ENOMEM) 
  {
    HLTError("could not SetupCTPData(); ENOMEM");
    return -ENOMEM;
  }

  return 0;
}

void* AliHLTAnalysisManagerComponent::AnalysisManagerExit(void*)
{
  if (fWriteAnalysisToFile && fAnalysisManager) 
    fAnalysisManager->WriteAnalysisToFile();
  delete fAnalysisManager;
  return(NULL);
}

// #################################################################################
Int_t AliHLTAnalysisManagerComponent::DoDeinit() {
  // see header file for class documentation

  if (fForceKillAsyncProcess != -1) fAsyncProcessor.ForceChildExit(fForceKillAsyncProcess);
  if (fAsyncProcessor.GetNumberOfAsyncTasksInQueue())
  {
    HLTError("Cannot deinitialize AsyncProcessor - Still tasks in queue");
    //We wait for all tasks in the queue, fetch the results, and drop them.
    //This might result in a memory leak but will at least shut down the thread properly.
    fAsyncProcessor.WaitForTasks(0);
    while (fAsyncProcessor.IsQueuedTaskCompleted()) fAsyncProcessor.RetrieveQueuedTaskResult();
  }

  if (fForceKillAsyncProcess == -1)
  {
    fAsyncProcessor.InitializeAsyncMemberTask(this, &AliHLTAnalysisManagerComponent::AnalysisManagerExit, NULL);
  }

  fAsyncProcessor.Deinitialize();
  return 0;
}

//_________________________________________________________________________________________________
void AliHLTAnalysisManagerComponent::CleanEventData(AnalysisManagerQueueData* eventData)
{
  if (fAsyncProcess)
  {
    fAsyncProcessor.FreeBuffer(eventData);
  }
  else
  {
    delete eventData->fEvent;
    delete eventData->fFriend;
    delete eventData;
  }
}

//_________________________________________________________________________________________________
TCollection* AliHLTAnalysisManagerComponent::GetOutputs()
{
  //if not yet populated, populate the output container with pointers
  //to output data contained in TCollections instead of AliAnalysisDataContainers
  //which don't support proper streaming.
  //otherwise just return the pointer.
  if (!fAnalysisOutputContainer) {
    fAnalysisOutputContainer = new AliHLTObjArray(1);
    TIter i(fAnalysisManager->GetOutputs());
    while (AliAnalysisDataContainer* adc = (AliAnalysisDataContainer*)i.Next())
    {
      AliHLTObjArray* ersatzCont = new AliHLTObjArray(1);
      ersatzCont->SetName(adc->GetName());
      ersatzCont->AddAt(adc->GetData(),0);
      fAnalysisOutputContainer->Add(ersatzCont);
    }
  }
  return fAnalysisOutputContainer;
}

//_________________________________________________________________________________________________
void AliHLTAnalysisManagerComponent::ResetOutputData()
{
  fAnalysisManager->ResetOutputData();
  delete fAnalysisOutputContainer; //this will delete the ersatz containers but not the data
  fAnalysisOutputContainer = NULL;
}

//_________________________________________________________________________________________________
void* AliHLTAnalysisManagerComponent::AnalysisManagerDoEvent(void* tmpEventData)
{
  AnalysisManagerQueueData* eventData = (AnalysisManagerQueueData*) tmpEventData;
  
  void* retVal = NULL;

  if (!*((volatile bool*) &fQuickEndRun))
  {
    TStopwatch stopwatch;
    if (fQueueDepth == 0) stopwatch.Start();
  
    fInputHandler->InitTaskInputData(eventData->fEvent, eventData->fFriend, fAnalysisManager->GetTasks());
    fAnalysisManager->ExecAnalysis();
    fInputHandler->FinishEvent();
  
    //pushes once every n seconds if
    //configured with -pushback-period=n
    //fAnalysisManager->GetOutputs() is an AliHLTObjArray of AliHLTObjArray objects
    //containing user output data
    fNEvents++;
    if (fPushEventModulo == 0 || fNEvents % fPushEventModulo == 0)
    {
      retVal = GetOutputs();
    }
  
    if (fQueueDepth == 0)
    {
      stopwatch.Stop();
      AliSysInfo::AddStamp("analysisTiming",eventData->fEvent->GetNumberOfTracks(),stopwatch.RealTime()*1000,stopwatch.CpuTime()*1000);
    }
  }

  bool requestPush = eventData->fRequestPush;

  //In AsyncMode, we copy the content to a temporary buffer and reset the Analysis Manager directly
  if (retVal && fQueueDepth)
  {
    //If we are an async process, we cannot access the pushback-period of the parent process, so we use this flag
    if (!(fAsyncProcess ? requestPush : CheckPushbackPeriod()))
    {
      retVal = NULL; 
    }
    else
    {
      retVal = fAsyncProcessor.SerializeIntoBuffer((TObject*) retVal, this);
      if (fResetAfterPush) {ResetOutputData();}
    }
  }

  if (fQueueDepth)
  {
    CleanEventData(eventData);
  }

  //In synchronous mode, we can just return the object
  return(retVal);
}

// #################################################################################
Int_t AliHLTAnalysisManagerComponent::DoEvent(const AliHLTComponentEventData& evtData,
    AliHLTComponentTriggerData& /*trigData*/) {
  // see header file for class documentation
  
  if (!fAnalysisInitialized) return(0);

  TStopwatch stopwatch;
  stopwatch.Start();
  
  HLTInfo("AliHLTAnalysisManagerComponent::DoEvent");
  Int_t iResult=0;

  if( fUID == 0 ){
    TTimeStamp t;
    fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
  }

  AnalysisManagerQueueData* eventData = NULL;
  if (IsDataEvent())
  {
    // -- Get ESD object
    // -------------------
    AliVEvent* vEvent=NULL;
    AliVfriendEvent* vFriend=NULL;
    if (fAsyncProcess)
    {
      eventData = (AnalysisManagerQueueData*) fAsyncProcessor.AllocateBuffer();
    }
    else
    {
      eventData = new AnalysisManagerQueueData;
    }
    if (eventData == NULL)
    {
      HLTError("Memory Allocation Error");
      return(-ENOSPC);
    }
    iResult = ReadInput(eventData);
    if (iResult != 0 || !eventData->fEvent)
    {
      CleanEventData(eventData);
      return(iResult);
    }

    if (eventData->fEvent) {HLTInfo("----> event %p has %d tracks: \n", eventData->fEvent, eventData->fEvent->GetNumberOfTracks());}
    if (eventData->fFriend) {HLTInfo("----> friend %p has %d tracks: \n", eventData->fFriend, eventData->fFriend->GetNumberOfTracks());}
    const AliHLTCTPData* ctpData = CTPData();
    std::string activeTriggers;
    activeTriggers.reserve(1000);
    ctpData->GetFiredTriggerClasses(activeTriggers);
    HLTInfo("active triggers: %s", activeTriggers.c_str());

    if (eventData->fEvent->GetNumberOfTracks() >= fMinTracks)
    {
      if (fMinTracks) HLTInfo("Event has %d tracks, running AnalysisManager", eventData->fEvent->GetNumberOfTracks());
      eventData->fRequestPush = CheckPushbackPeriod() && !fPushRequestOngoing;
      if (fAsyncProcessor.QueueAsyncMemberTask(this, &AliHLTAnalysisManagerComponent::AnalysisManagerDoEvent, eventData))
      {
        CleanEventData(eventData); //Could not queue this task, clean up
      }
      else
      {
        fNumEvents++;
        if (eventData->fRequestPush) fPushRequestOngoing = kTRUE;
      }
    }
    else
    {
      CleanEventData(eventData); //Buffers have been allocated, clean up
    }
  }
  
  if (!IsDataEvent() && GetFirstInputBlock(kAliHLTDataTypeSOR | kAliHLTDataOriginAny))
  {
    fQuickEndRun = false;
  }
  if (!IsDataEvent() && GetFirstInputBlock(kAliHLTDataTypeEOR | kAliHLTDataOriginAny))
  {
    fQuickEndRun = true;
    if (fForceKillAsyncProcess != -1) fAsyncProcessor.ForceChildExit(fForceKillAsyncProcess);
    fAsyncProcessor.WaitForTasks(0);
  }

  while (fAsyncProcessor.IsQueuedTaskCompleted())
  {
    void* retVal = fAsyncProcessor.RetrieveQueuedTaskResult();
    if (retVal)
    {
      if (fQueueDepth)
      {
        AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* ret = (AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer*) retVal;
	int pushResult = PushBack(ret->fPtr, ret->fSize, kAliHLTDataTypeTObject|kAliHLTDataOriginHLT,fUID);
	fPushRequestOngoing = kFALSE;
	if (pushResult)
	{
	  HLTInfo("HLT Analysis Manager pushing output: %p (%d bytes, %d events)", retVal, pushResult, fNumEvents);
	  fNumEvents = 0;
	}

	fAsyncProcessor.FreeBuffer(ret);
      }
      else
      {
        TObject* retObj = (TObject*) retVal;
        int pushResult = PushBack(retObj, kAliHLTDataTypeTObject|kAliHLTDataOriginHLT,fUID);
        if (pushResult > 0)
        {
           if (fResetAfterPush) ResetOutputData();
           HLTInfo("HLT Analysis Manager pushing output: %p (%d bytes, %d events)", retVal, pushResult, fNumEvents);
           fNumEvents = 0;
        }
      }
    }
  }

  if (fQueueDepth == 0 && eventData)
  {
    CleanEventData(eventData);
  }

  return 0;
}

// #################################################################################
Int_t AliHLTAnalysisManagerComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
  // see header file for class documentation

  Int_t iResult=0;
  TString cdbPath;
  if (cdbEntry) {
    cdbPath=cdbEntry;
  } else {
    cdbPath="HLT/ConfigGlobal/";
    cdbPath+=GetComponentID();
  }

  AliInfoClass(Form("reconfigure '%s' from entry %s%s", chainId, cdbPath.Data(), cdbEntry?"":" (default)"));
  iResult=ConfigureFromCDBTObjString(cdbPath);

  return iResult;
}

// #################################################################################
Int_t AliHLTAnalysisManagerComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
  // see header file for class documentation
  ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
  return 0;
}

// #################################################################################
Int_t AliHLTAnalysisManagerComponent::ReadInput(AnalysisManagerQueueData* eventData)
{
  //read input, return code 0 for normal events, 1 for EOR
  
  AliVEvent* vEvent = NULL;
  AliVfriendEvent* vFriend = NULL;

  Int_t iResult=0;
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {
    vEvent = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
    if( !vEvent ){ 
      HLTWarning("Wrong ESDEvent object received");
      iResult = -1;
      continue;
    }
    if (fAsyncProcess)
    {
      HLTFatal("Invalid configuration: HLTAnalysisManager must process FlatESDs not ESDs in async-process mode!");
    }
    if (RemoveInputObjectFromCleanupList(vEvent) == NULL)
    {
      HLTError("Error taking ownership for esdEvent, cannot queue async calibration task.");
    }
    vEvent->GetStdContent();
  }
  for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDfriendObject); iter != NULL; iter = GetNextInputObject() ) {
    vFriend = dynamic_cast<AliESDfriend*>(const_cast<TObject*>( iter ) );
    if( !vFriend ){ 
      HLTWarning("Wrong ESDFriend object received");
      iResult = -1;
      continue;
    }
    if (fAsyncProcess)
    {
      HLTFatal("Invalid configuration: HLTAnalysisManager must process FlatESDs not ESDs in async-process mode!");
    }
    if (RemoveInputObjectFromCleanupList(vFriend) == NULL)
    {
      HLTError("Error taking ownership for esdEvent-friends, cannot queue async calibration task.");
    }
  }

    size_t flatEsdSize;
    if (!vEvent)
    {
	const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESD|kAliHLTDataOriginOut);
	if (pBlock)
	{
	    AliFlatESDEvent* tmpFlatEvent=reinterpret_cast<AliFlatESDEvent*>( pBlock->fPtr );
	    if (tmpFlatEvent)
	    {
		AliFlatESDEvent* tmpCopy;
		if (fAsyncProcess)
		{
			if (pBlock->fSize + sizeof(AnalysisManagerQueueData) > fAsyncProcessor.BufferSize())
			{
				HLTError("Insufficient buffer size for Flat-ESD");
				tmpCopy = NULL;
			}
			else
			{
				tmpCopy = (AliFlatESDEvent*) (((char*) eventData) + sizeof(AnalysisManagerQueueData));
			}
		}
		else
		{
			tmpCopy = (AliFlatESDEvent*) new Byte_t[pBlock->fSize];
		}

		if (tmpCopy == NULL)
		{
			HLTError("Memory allocation error");
			tmpFlatEvent = NULL;
		}
		else
		{
			memcpy((void*) tmpCopy, (void*) tmpFlatEvent, pBlock->fSize);
			tmpFlatEvent = tmpCopy;
			tmpFlatEvent->Reinitialize();

			if (tmpFlatEvent->GetSize() != pBlock->fSize)
			{
				tmpFlatEvent = NULL;
				iResult=-1;
    				HLTWarning("data mismatch in block %s (0x%08x): size %d -> ignoring flatESD information",
				DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize);
			}
			else
			{
				vEvent = tmpFlatEvent;
				flatEsdSize = pBlock->fSize;
				if (flatEsdSize % 32) flatEsdSize += 32 - flatEsdSize % 32;
			}
		}
	    }
	}
    }

    if( vEvent )
    {
	const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeFlatESDFriend|kAliHLTDataOriginOut);
	if (pBlock)
	{
	    AliFlatESDFriend* tmpFlatFriend = reinterpret_cast<AliFlatESDFriend*>( pBlock->fPtr );
	    if (tmpFlatFriend)
	    {
		AliFlatESDFriend* tmpCopy;
		if (fAsyncProcess)
		{
			if (pBlock->fSize + flatEsdSize + sizeof(AnalysisManagerQueueData) > fAsyncProcessor.BufferSize())
			{
				HLTError("Insufficient buffer size for Flat-ESD Friend");
				tmpCopy = NULL;
			}
			else
			{
				tmpCopy = (AliFlatESDFriend*) (((char*) vEvent) + flatEsdSize);
			}
		}
		else
		{
			tmpCopy  = (AliFlatESDFriend*) new Byte_t[pBlock->fSize];
		}

		if (tmpCopy == NULL)
		{
			HLTError("Memory allocation error");
			tmpFlatFriend = NULL;
		}
		else
		{
			memcpy((void*) tmpCopy, (void*) tmpFlatFriend, pBlock->fSize);
			tmpFlatFriend = tmpCopy;
			tmpFlatFriend->Reinitialize();

			if (tmpFlatFriend->GetSize() != pBlock->fSize)
			{
				tmpFlatFriend = NULL;
				iResult=-1;
				HLTWarning("data mismatch in block %s (0x%08x): size %d -> ignoring flatESDfriend information", 
				DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, pBlock->fSize);
			}
			else
			{
				vFriend = tmpFlatFriend;
			}
		}
	    }
	}
    }
    
    eventData->fEvent = vEvent;
    eventData->fFriend = vFriend;
    return iResult;
}

// #################################################################################
int AliHLTAnalysisManagerComponent::ProcessOption(TString option, TString value)
{
  //process option
  //to be implemented by the user
  if (option.EqualTo("WriteAnalysisToFile")) 
  {
    fWriteAnalysisToFile=(value.Contains("0"))?kFALSE:kTRUE;
    HLTInfo("fWriteAnalysisToFile=%i\n",fWriteAnalysisToFile?1:0);
  }
  else if (option.EqualTo("AddTaskMacro"))
  {
    fAddTaskMacro=value;
    HLTInfo("fAddTaskMacro=%s\n",fAddTaskMacro.Data());
  }
  else if (option.EqualTo("PushEventModulo"))
  {
    fPushEventModulo=atoi(value.Data());
    HLTInfo("fPushEventModulo=%d\n",fPushEventModulo);
  }
  else if (option.EqualTo("MinTracks"))
  {
    fMinTracks=atoi(value.Data());
    HLTInfo("fMinTracks=%d\n",fMinTracks);
  }
  else if (option.EqualTo("QueueDepth"))
  {
    fQueueDepth=atoi(value.Data());
    HLTInfo("fQueueDepth=%d\n",fQueueDepth);
  }
  else if (option.EqualTo("AsyncProcess"))
  {
    fAsyncProcess=atoi(value.Data());
    if (fAsyncProcess == 1) fAsyncProcess = 100000000;
    if (fAsyncProcess) HLTInfo("AsyncProcess (%d bytes buffer)\n", fAsyncProcess);
  }
  else if (option.EqualTo("ForceKillAsyncProcess"))
  {
    fForceKillAsyncProcess=atoi(value.Data());
    if (fAsyncProcess) HLTInfo("ForceKillAsyncProcess set to %d\n", fForceKillAsyncProcess);
  }
  else if (option.EqualTo("ResetAfterPush"))
  {
    fResetAfterPush=(value.Contains("0")?kFALSE:kTRUE);
    HLTInfo("fResetAfterPush=%i\n",fResetAfterPush?1:0);
  }
  else if (option.EqualTo("NoFullQueueWarning"))
  {
    fAsyncProcessor.SetFullQueueWarning(0);
  }
  else if (option.EqualTo("InitializeGeometry"))
  {
    fInitializeGeometry =  value.EqualTo("1");
  }
  else if (option.EqualTo("EnableDebug"))
  {
    fEnableDebug=value.Contains("1");
    HLTInfo("fEnableDebug=%s",fEnableDebug?"1":"0");
  }
  return 1; 
}

// #################################################################################
int AliHLTAnalysisManagerComponent::ProcessOptionString(TString arguments)
{
  //process passed options
  HLTInfo("Argument string: %s\n", arguments.Data());
  stringMap* options = TokenizeOptionString(arguments);
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    HLTInfo("  %s : %s\n", i->first.data(), i->second.data());
    ProcessOption(i->first,i->second);
  }
  delete options; //tidy up

  return 1; 
}

// #################################################################################
AliHLTAnalysisManagerComponent::stringMap* AliHLTAnalysisManagerComponent::TokenizeOptionString(const TString str)
{
  //options have the form:
  // -option value
  // -option=value
  // -option
  // --option value
  // --option=value
  // --option
  // option=value
  // option value
  // (value can also be a string like 'some string')
  //
  // options can be separated by ' ' or ',' arbitrarily combined, e.g:
  //"-option option1=value1 --option2 value2, -option4=\'some string\'"
  
  //optionRE by construction contains a pure option name as 3rd submatch (without --,-, =)
  //valueRE does NOT match options
  TPRegexp optionRE("(?:(-{1,2})|((?='?[^,=]+=?)))"
                    "((?(2)(?:(?(?=')'(?:[^'\\\\]++|\\.)*+'|[^, =]+))(?==?))"
                    "(?(1)[^, =]+(?=[= ,$])))");
  TPRegexp valueRE("(?(?!(-{1,2}|[^, =]+=))"
                   "(?(?=')'(?:[^'\\\\]++|\\.)*+'"
                   "|[^, =]+))");

  stringMap* options = new stringMap;

  TArrayI pos;
  const TString mods="";
  Int_t start = 0;
  while (1) {
    Int_t prevStart=start;
    TString optionStr="";
    TString valueStr="";
    
    //check if we have a new option in this field
    Int_t nOption=optionRE.Match(str,mods,start,10,&pos);
    if (nOption>0)
    {
      optionStr = str(pos[6],pos[7]-pos[6]);
      optionStr=optionStr.Strip(TString::kBoth,'\'');
      start=pos[1]; //update the current character to the end of match
    }

    //check if the next field is a value
    Int_t nValue=valueRE.Match(str,mods,start,10,&pos);
    if (nValue>0)
    {
      valueStr = str(pos[0],pos[1]-pos[0]);
      valueStr=valueStr.Strip(TString::kBoth,'\'');
      start=pos[1]; //update the current character to the end of match
    }
    
    //skip empty entries
    if (nOption>0 || nValue>0)
    {
      (*options)[optionStr.Data()] = valueStr.Data();
    }
    
    if (start>=str.Length()-1 || start==prevStart ) break;
  }

/*
  for (stringMap::iterator i=options->begin(); i!=options->end(); ++i)
  {
    HLTInfo("%s : %s", i->first.data(), i->second.data());
  }
*/
  return options;
}
