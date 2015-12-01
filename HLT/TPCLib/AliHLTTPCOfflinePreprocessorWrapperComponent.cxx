
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
#include "AliAnalysisDataContainer.h"
#include "AliTPCPreprocessorOffline.h"
#include "AliTPCcalibTime.h"

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

AliCDBEntry* AliHLTTPCOfflinePreprocessorWrapperComponent::RunPreprocessor(AliAnalysisDataContainer* dataContainer)
{

	AliTPCcalibTime* timeObj = dynamic_cast<AliTPCcalibTime*>(dataContainer->GetData());
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
		preprocessor->CalibTimeVdrift(timeObj, GetRunNo(), GetRunNo());
		preprocessor->CreateDriftCDBentryObject(GetRunNo(), GetRunNo());
		retVal = preprocessor->GetDriftCDBentry();
		preprocessor->TakeOwnershipDriftCDBEntry();
	}
	delete dataContainer;
	//TODO: Memory Leaks: We can NOT delete fNewCalibObject. The CDBEntry created by the Preprocessor contains a TObjArray that links to some data in fNewCalibObject.
	//We have no control over where in AliRoot someone queries that CDBEntry and we do not know when it is no longer needed!
	
	return(retVal);
}

void* AliHLTTPCOfflinePreprocessorWrapperComponent::AsyncRunPreprocessor(void* obj)
{
	return((void*) RunPreprocessor((AliAnalysisDataContainer*) obj));
}

int AliHLTTPCOfflinePreprocessorWrapperComponent::DoInit( int argc, const char** argv ) 
{ 
  // see header file for class documentation
  
  int iResult=0;
  //!! iResult = ConfigureFromCDBTObjString(fgkOCDBEntryClusterTransformation);

  if (iResult>=0 && argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);
  
  if (fAsyncProcessor.Initialize(fAsyncProcessorQueueDepth)) return(1);

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
    } else {
      iRet = -EINVAL;
      HLTInfo("Unknown argument %s",argv[i]);     
    }
  } 
  return iRet;
}


Int_t AliHLTTPCOfflinePreprocessorWrapperComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/) {
	// see header file for class documentation

	// -- Only use data event
	if (IsDataEvent())
	{
		AliAnalysisDataContainer* tmpContainer = NULL;
		for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeTObject); iter != NULL; iter = GetNextInputObject() )
		{
			TObjArray* tmpArray = (TObjArray*) dynamic_cast<const TObjArray*>(iter);
			tmpContainer = tmpArray ? NULL : (AliAnalysisDataContainer*) dynamic_cast<const AliAnalysisDataContainer*>(iter);
			if (tmpContainer)
			{
				RemoveInputObjectFromCleanupList(tmpContainer);
			}
			else if (tmpArray)
			{
				for (int i = 0;i <= tmpArray->GetLast();i++)
				{
					tmpContainer = (AliAnalysisDataContainer*) dynamic_cast<const AliAnalysisDataContainer*>((*tmpArray)[i]);
					if (tmpContainer != NULL && strcmp(tmpContainer->GetName(), "calibTime") == 0)
					{
						RemoveInputObjectFromCleanupList(tmpArray);
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
				HLTImportant("TPC Offline Preprocessor component received object that is no AliAnalysisDataContainer!");
			}
		}
		
		if (tmpContainer)
		{
			fAsyncProcessor.QueueAsyncMemberTask(this, &AliHLTTPCOfflinePreprocessorWrapperComponent::AsyncRunPreprocessor, tmpContainer);
		}
	}
	
	if (!IsDataEvent() && GetFirstInputBlock(kAliHLTDataTypeEOR | kAliHLTDataOriginAny))
	{
		fAsyncProcessor.WaitForTasks(0);
	}
	
	//If a new transform map is available from an async creation task, we ship the newest one.
	while (fAsyncProcessor.IsQueuedTaskCompleted())
	{
		AliCDBEntry* retVal = (AliCDBEntry*) fAsyncProcessor.RetrieveQueuedTaskResult();
		int pushResult = PushBack(dynamic_cast<TObject*>(retVal), kAliHLTDataTypeTObject | kAliHLTDataOriginHLT);
		if (pushResult) HLTImportant("Offline Preprocessor pushed CDB entry (%d bytes)", pushResult);
		delete retVal;
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
