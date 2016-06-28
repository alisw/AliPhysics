
//***************************************************************************
//* This file is property of and copyright by the ALICE HLT Project         *
//* ALICE Experiment at CERN, All rights reserved.                          *
//*                                                                         *
//* Primary Authors: Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de *
//*                  for The ALICE HLT Project.                             *
//*                                                                         *
//* Permission to use, copy, modify and distribute this software and its    *
//* documentation strictly for non-commercial purposes is hereby granted    *
//* without fee, provided that the above copyright notice appears in all    *
//* copies and that both the copyright notice and this permission notice    *
//* appear in the supporting documentation. The authors make no claims      *
//* about the suitability of this software for any purpose. It is           *
//* provided "as is" without express or implied warranty.                   *
//***************************************************************************

/** @file   AliHLTTPCOfflinePreprocessorWrapperComponent.cxx
    @author Sergey Gorbunov
    @date   
    @brief 
*/

#include "AliHLTTPCOfflinePreprocessorWrapperComponent.h"
#include "AliHLTTPCClusterTransformation.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTTPCFastTransformObject.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "TObjArray.h"
#include "AliTPCPreprocessorOffline.h"
#include "AliTPCcalibTime.h"
#include "AliHLTMessage.h"
#include "TDirectory.h"
#include "TH1.h"

#include "TMath.h"
#include "TObjString.h" 
#include "TStopwatch.h"
#include <cstdlib>
#include <cerrno>
#include <sys/time.h>

class AliTPCcalibTime;

using namespace std;

ClassImp(AliHLTTPCOfflinePreprocessorWrapperComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCOfflinePreprocessorWrapperComponent::AliHLTTPCOfflinePreprocessorWrapperComponent() :
  fAsyncProcessor()
  {}

AliHLTTPCOfflinePreprocessorWrapperComponent::~AliHLTTPCOfflinePreprocessorWrapperComponent()
{ 
  // destructor
}

const char* AliHLTTPCOfflinePreprocessorWrapperComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCOfflinePreprocessorWrapper";
}

void AliHLTTPCOfflinePreprocessorWrapperComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
  // see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeTObject | kAliHLTDataOriginHLT );
}

AliHLTComponentDataType AliHLTTPCOfflinePreprocessorWrapperComponent::GetOutputDataType() { 
  // see header file for class documentation

  return kAliHLTDataTypeTObject;
}

void AliHLTTPCOfflinePreprocessorWrapperComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
  // see header file for class documentation
  constBase = 1000000;
  inputMultiplier = 1.0;
}

AliHLTComponent* AliHLTTPCOfflinePreprocessorWrapperComponent::Spawn() { 
  // see header file for class documentation

  return new AliHLTTPCOfflinePreprocessorWrapperComponent();
}

AliCDBEntry* AliHLTTPCOfflinePreprocessorWrapperComponent::RunPreprocessor(TObjArray* dataContainer)
{
	if (fAsyncProcess)
	{
		AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* buffer = (AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer*) dataContainer;
		TObject* tmpObject = AliHLTMessage::Extract(buffer->fPtr, buffer->fSize);
		if (tmpObject == NULL)
		{
			HLTError("Error deserializing object");
			return(NULL);
		}
		fAsyncProcessor.FreeBuffer(buffer);
		dataContainer = GetDataContainer((TObject*) tmpObject);
	}
	AliTPCcalibTime* timeObj = dynamic_cast<AliTPCcalibTime*>(dataContainer->At(0));
	AliCDBEntry* retVal = NULL;
	if (timeObj == NULL)
	{
		HLTError("Error obtaining AliTPCcalibTime Object from received calibration object, cannot reinitialize transformation map");
	}
	else
	{
		HLTImportant("TPC Offline Preprocessor preparing drift time cdb object");
		AliCDBManager* cdbManager = AliCDBManager::Instance();
		
		static AliTPCPreprocessorOffline* preprocessor = new AliTPCPreprocessorOffline;
		preprocessor->SetTimeGainRange(0.5,5.0);
		preprocessor->SetMaxVDriftCorr(0.2);
		preprocessor->CalibTimeVdrift(timeObj, GetRunNo(), GetRunNo());
		preprocessor->CreateDriftCDBentryObject(GetRunNo(), GetRunNo());
		retVal = preprocessor->GetDriftCDBentry();
		preprocessor->TakeOwnershipDriftCDBEntry();
	}
	delete dataContainer;
	//TODO: Memory Leaks: We can NOT delete fNewCalibObject. The CDBEntry created by the Preprocessor contains a TObjArray that links to some data in fNewCalibObject.
	//We have no control over where in AliRoot someone queries that CDBEntry and we do not know when it is no longer needed!
	
	if (fAsyncProcess)
	{
		void* serRetVal = fAsyncProcessor.SerializeIntoBuffer(retVal, this);
		delete retVal;
		return((AliCDBEntry*) serRetVal);
	}
	
	return(retVal);
}

void* AliHLTTPCOfflinePreprocessorWrapperComponent::AsyncRunPreprocessor(void* obj)
{
	return((void*) RunPreprocessor((TObjArray*) obj));
}

int AliHLTTPCOfflinePreprocessorWrapperComponent::DoInit( int argc, const char** argv ) 
{ 
  // see header file for class documentation
  
  int iResult=0;
  //!! iResult = ConfigureFromCDBTObjString(fgkOCDBEntryClusterTransformation);

  //don't keep track of root objects
  if (getenv("HLT_ONLINE_MODE") && strcmp(getenv("HLT_ONLINE_MODE"), "on") == 0)
  {
    TDirectory::AddDirectory(kFALSE);
    TH1::AddDirectory(kFALSE);
  }

  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);
  
  if (fAsyncProcessor.Initialize(fAsyncProcessorQueueDepth, fAsyncProcess, fAsyncProcess ? 100000000 : 0)) return(1);

  return iResult;
} // end DoInit()

int AliHLTTPCOfflinePreprocessorWrapperComponent::DoDeinit() { 
  // see header file for class documentation
  if (fAsyncProcessor.GetNumberOfAsyncTasksInQueue())
  {
    //We just dump all remaining tasks, since we do not need the transformation maps anymore
	fAsyncProcessor.WaitForTasks(0);
	while (fAsyncProcessor.IsQueuedTaskCompleted()) fAsyncProcessor.RetrieveQueuedTaskResult();
  }
  fAsyncProcessor.Deinitialize();
  return 0;
}

int AliHLTTPCOfflinePreprocessorWrapperComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) { 
  // see header file for class documentation
  return 0;//!! ConfigureFromCDBTObjString(fgkOCDBEntryClusterTransformation);
}

int AliHLTTPCOfflinePreprocessorWrapperComponent::ScanConfigurationArgument(int argc, const char** argv){

  // see header file for class documentation

  if (argc<=0) return 0;
  int iRet = 0;
  for( int i=0; i<argc; i++ ){
    TString argument=argv[i];  
    if (argument.CompareTo("-QueueDepth")==0){
	  if (++i >= argc)
	  {
	    HLTError("Value missing for -QueueDepth parameter");
		return(-EINVAL);
	  }
	  fAsyncProcessorQueueDepth = atoi(argv[i]);
      HLTInfo("Queue Depth set to %d.", fAsyncProcessorQueueDepth);
      iRet+=2;
          }
    else if (argument.CompareTo("-AsyncProcess")==0){
	  fAsyncProcess = 1;
      iRet++;
    } else {
      iRet = -EINVAL;
      HLTError("Unknown argument %s",argv[i]);     
    }
  } 
  return iRet;
}

TObjArray* AliHLTTPCOfflinePreprocessorWrapperComponent::GetDataContainer(TObject* obj)
{
	TObjArray* tmpContainer = NULL;
	TObjArray* tmpArray = (TObjArray*) dynamic_cast<const TObjArray*>(obj);
	tmpContainer = tmpArray ? NULL : (TObjArray*) dynamic_cast<const TObjArray*>(obj);
	if (tmpContainer)
	{
		if (!fAsyncProcess) RemoveInputObjectFromCleanupList(tmpContainer);
	}
	else if (tmpArray)
	{
		for (int i = 0;i <= tmpArray->GetLast();i++)
		{
			tmpContainer = (TObjArray*) dynamic_cast<const TObjArray*>((*tmpArray)[i]);
			if (tmpContainer != NULL && strcmp(tmpContainer->GetName(), "calibTime") == 0)
			{
				if (!fAsyncProcess) RemoveInputObjectFromCleanupList(tmpArray);
				tmpArray->Remove(tmpContainer);
				break;
			}
			else
			{
				tmpContainer = NULL;
			}
		}
	}
	
	if (!tmpContainer)
	{
		HLTImportant("TPC Offline Preprocessor component received object that is no TObjArray!");
	}
	return(tmpContainer);
}

Int_t AliHLTTPCOfflinePreprocessorWrapperComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/) {
	// see header file for class documentation

	// -- Only use data event
	if (IsDataEvent())
	{
		TObjArray* tmpContainer = NULL;
		if (fAsyncProcess)
		{
			const AliHLTComponentBlockData* blk = GetFirstInputBlock(kAliHLTDataTypeTObject);
			if (blk)
			{
				AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* myBuffer = fAsyncProcessor.AllocateBuffer(blk->fSize);
				if (myBuffer == NULL)
				{
					HLTError("No shared buffer available (Buf Size %d, Requested %d)", (int) fAsyncProcessor.BufferSize(), blk->fSize);
				}
				else
				{
					memcpy(myBuffer->fPtr, blk->fPtr, blk->fSize);
					myBuffer->fSize = blk->fSize;
					if (fAsyncProcessor.QueueAsyncMemberTask(this, &AliHLTTPCOfflinePreprocessorWrapperComponent::AsyncRunPreprocessor, myBuffer))
					{
						fAsyncProcessor.FreeBuffer(myBuffer);
					}
				}
			}
		}
		else
		{
			for (const TObject *iter = GetFirstInputObject(kAliHLTDataTypeTObject); iter != NULL; iter = GetNextInputObject())
			{
				if ((tmpContainer = GetDataContainer((TObject*) iter))) break;
			}
			if (tmpContainer)
			{
				if (fAsyncProcessor.QueueAsyncMemberTask(this, &AliHLTTPCOfflinePreprocessorWrapperComponent::AsyncRunPreprocessor, tmpContainer))
				{
					delete tmpContainer;
				}
			}
		}
	}
	
	if (!IsDataEvent() && GetFirstInputBlock(kAliHLTDataTypeEOR | kAliHLTDataOriginAny))
	{
		fAsyncProcessor.WaitForTasks(0);
	}

	while (fAsyncProcessor.IsQueuedTaskCompleted())
	{
		AliCDBEntry* retVal;
		int pushResult;
		if (fAsyncProcess)
		{
			AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* buf = (AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer*) fAsyncProcessor.RetrieveQueuedTaskResult();
			pushResult = PushBack(buf->fPtr, buf->fSize, kAliHLTDataTypeTObject | kAliHLTDataOriginHLT);
			fAsyncProcessor.FreeBuffer(buf);
		}
		else
		{
			retVal = (AliCDBEntry*) fAsyncProcessor.RetrieveQueuedTaskResult();
			pushResult = PushBack(dynamic_cast<TObject*>(retVal), kAliHLTDataTypeTObject | kAliHLTDataOriginHLT);
			delete retVal;
		}
		if (pushResult) HLTImportant("Offline Preprocessor pushed CDB entry (%d bytes)", pushResult);
	}

	return 0;
}

void AliHLTTPCOfflinePreprocessorWrapperComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description needed for the particular component
  if (!targetMap) return;
  
  // OCDB entries for component arguments

  //!! targetMap->Add(new TObjString(fgkOCDBEntryClusterTransformation), new TObjString("component argument for the charge threshold"));
  
  // OCDB entries to be fetched by the TAXI (access via the AliTPCcalibDB class)
  targetMap->Add(new TObjString("TPC/Calib/Parameters"),    new TObjString("unknown content"));
  targetMap->Add(new TObjString("TPC/Calib/TimeDrift"),     new TObjString("drift velocity calibration"));
  targetMap->Add(new TObjString("TPC/Calib/TimeGain"),      new TObjString("time gain  calibration"));
  targetMap->Add(new TObjString("TPC/Calib/Temperature"),   new TObjString("temperature map"));
  targetMap->Add(new TObjString("TPC/Calib/PadGainFactor"), new TObjString("gain factor pad by pad"));
  targetMap->Add(new TObjString("TPC/Calib/ClusterParam"),  new TObjString("cluster parameters"));
  targetMap->Add(new TObjString("TPC/Calib/Correction"),  new TObjString("coreection"));
  targetMap->Add(new TObjString("TPC/Calib/RecoParam"),  new TObjString("reconstruction parameters"));
 
  // OCDB entries needed to be fetched by the Pendolino
  targetMap->Add(new TObjString("TPC/Calib/AltroConfig"), new TObjString("contains the altro config, e.g. info about the L0 trigger timing"));
  targetMap->Add(new TObjString("GRP/CTP/CTPtiming"),     new TObjString("content used in the cluster coordinate transformation in relation to the L0 trigger timing"));

  // OCDB entries necessary for replaying data on the HLT cluster
  targetMap->Add(new TObjString("GRP/GRP/Data"), new TObjString("contains magnetic field info"));  
 
  // OCDB entries needed to suppress fatals/errors/warnings during reconstruction
  targetMap->Add(new TObjString("TPC/Calib/Distortion"),  new TObjString("distortion map"));
  targetMap->Add(new TObjString("TPC/Calib/GainFactorDedx"), new TObjString("gain factor dedx"));
  targetMap->Add(new TObjString("TPC/Calib/PadTime0"),    new TObjString("time0 offset pad by pad"));
  targetMap->Add(new TObjString("TPC/Calib/PadNoise"),    new TObjString("pad noise values"));
  targetMap->Add(new TObjString("TPC/Calib/Pedestals"),   new TObjString("pedestal info"));
  targetMap->Add(new TObjString("TPC/Calib/Pulser"),      new TObjString("pulser info"));
  targetMap->Add(new TObjString("TPC/Calib/CE"),          new TObjString("CE laser calibration result"));
  targetMap->Add(new TObjString("TPC/Calib/Raw"),         new TObjString("unknown content"));
  targetMap->Add(new TObjString("TPC/Calib/QA"),          new TObjString("not important"));
  targetMap->Add(new TObjString("TPC/Calib/Mapping"),     new TObjString("unknown content"));
  targetMap->Add(new TObjString("TPC/Calib/Goofie"),      new TObjString("Goofie values, not used at the moment (05.03.2010)"));
  targetMap->Add(new TObjString("TPC/Calib/HighVoltage"), new TObjString("high voltage values, not used"));
  targetMap->Add(new TObjString("TPC/Calib/Ref"),         new TObjString("unknown content"));
}
