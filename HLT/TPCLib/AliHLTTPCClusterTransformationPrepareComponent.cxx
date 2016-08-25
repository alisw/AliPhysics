
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

/** @file   AliHLTTPCClusterTransformationPrepareComponent.cxx
    @author Sergey Gorbunov
    @date   
    @brief 
*/

#include "AliHLTTPCClusterTransformationPrepareComponent.h"
#include "AliHLTTPCClusterTransformation.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCRawCluster.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTTPCFastTransformObject.h"

#include "AliCDBManager.h"
#include "AliCDBEntry.h"
#include "AliTPCcalibDB.h"
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

ClassImp(AliHLTTPCClusterTransformationPrepareComponent) //ROOT macro for the implementation of ROOT specific class methods

const char* AliHLTTPCClusterTransformationPrepareComponent::fgkOCDBEntryClusterTransformation="HLT/ConfigTPC/TPCClusterTransformation";

AliHLTTPCClusterTransformation AliHLTTPCClusterTransformationPrepareComponent::fgTransform;
Bool_t AliHLTTPCClusterTransformationPrepareComponent::fgTimeInitialisedFromEvent = 0;

AliHLTTPCClusterTransformationPrepareComponent::AliHLTTPCClusterTransformationPrepareComponent() :
  fMinInitSec(-1),
  fMaxInitSec(-1),
  fNoInitialObject(false),
  fTmpFastTransformObject(NULL),
  fAsyncProcessor(),
  fAsyncProcessorQueueDepth(0),
  fNewCalibObject(NULL)
  {}

AliHLTTPCClusterTransformationPrepareComponent::~AliHLTTPCClusterTransformationPrepareComponent()
{ 
  // destructor
}

const char* AliHLTTPCClusterTransformationPrepareComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCClusterTransformationPrepare";
}

void AliHLTTPCClusterTransformationPrepareComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
  // see header file for class documentation

  list.clear(); 
  list.push_back( kAliHLTDataTypeTObject | kAliHLTDataOriginHLT );
}

AliHLTComponentDataType AliHLTTPCClusterTransformationPrepareComponent::GetOutputDataType() { 
  // see header file for class documentation

  return AliHLTTPCDefinitions::fgkTPCFastTransformDataObjectDataType;
}

void AliHLTTPCClusterTransformationPrepareComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
  // see header file for class documentation
  constBase = 18000000;
  inputMultiplier = 0.0;
}

AliHLTComponent* AliHLTTPCClusterTransformationPrepareComponent::Spawn() { 
  // see header file for class documentation

  return new AliHLTTPCClusterTransformationPrepareComponent();
}

AliHLTTPCFastTransformObject* AliHLTTPCClusterTransformationPrepareComponent::GenerateFastTransformObject()
{
	AliTPCcalibDB *calib=AliTPCcalibDB::Instance();
	if (!calib)
	{
		HLTError("AliTPCcalibDB does not exist");
		return(NULL);
	}
	
	if (fNewCalibObject)
	{
		AliCDBManager* cdbManager = AliCDBManager::Instance();
		cdbManager->PromptCacheEntry("TPC/Calib/TimeDrift", fNewCalibObject);

		//TODO: Memory Leaks: We can NOT delete fNewCalibObject. The CDBEntry created by the Preprocessor contains a TObjArray that links to some data in fNewCalibObject.
		//We have no control over where in AliRoot someone queries that CDBEntry and we do not know when it is no longer needed!
		fNewCalibObject = NULL;
	}
	else
	{
		HLTImportant("No updated calibration data available, creating transformation map with last calibration data available.");
	}

	calib->SetRun(GetRunNo());
	calib->UpdateRunInformations(GetRunNo());
	calib->Update();

	TStopwatch timer;
	timer.Start();
	//Workaround for offline simulation. In offline simulation there is one shared static instance of AliHLTTPCFastTransform, so we have to reset the sector borders every time.
	if (fMinInitSec != -1 && fMaxInitSec != -1)
	{
		fgTransform.SetInitSec(fMinInitSec, fMaxInitSec);
	}
	int err = fgTransform.Init( GetBz(), GetTimeStamp() );

	timer.Stop();
	HLTImportant("Initialization time: %f / %f", timer.CpuTime(), timer.RealTime());
	   
	AliHLTTPCFastTransformObject* obj = new AliHLTTPCFastTransformObject;
	
	fgTransform.GetFastTransformNonConst().WriteToObject(*obj);
	fgTransform.DeInit();
	
	return(obj);
}

void* AliHLTTPCClusterTransformationPrepareComponent::AsyncGenerateFastTransformObject(void*)
{
	return((void*) GenerateFastTransformObject());
}

int AliHLTTPCClusterTransformationPrepareComponent::DoInit( int argc, const char** argv ) 
{ 
  // see header file for class documentation
  
  int iResult=0;

  if (argc>0)
    iResult=ConfigureFromArgumentString(argc, argv);
  
  if (fMinInitSec != -1 || fMaxInitSec != -1)
  {
    if (fMinInitSec != -1 && fMaxInitSec != -1)
	{
	  fgTransform.SetInitSec(fMinInitSec, fMaxInitSec);
	}
	else
	{
	  HLTError("Both MinSector and MaxSector must be provided!");
	  return(1);
	}
  }

  if (fNoInitialObject)
  {
	HLTInfo("Skipping creation of initial transform object");
  }
  else
  {
	fTmpFastTransformObject = GenerateFastTransformObject();
  }
  
  if (fAsyncProcessor.Initialize(fAsyncProcessorQueueDepth)) return(1);

  return iResult;
} // end DoInit()

int AliHLTTPCClusterTransformationPrepareComponent::DoDeinit() { 
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

int AliHLTTPCClusterTransformationPrepareComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) { 
  // see header file for class documentation
  return 0;//!! ConfigureFromCDBTObjString(fgkOCDBEntryClusterTransformation);
}

int AliHLTTPCClusterTransformationPrepareComponent::ScanConfigurationArgument(int argc, const char** argv){

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
    else if (argument.CompareTo("-MinSector")==0){
	  if (++i >= argc)
	  {
	    HLTError("Value missing for -MinSector parameter");
		return(-EINVAL);
	  }
	  fMinInitSec = atoi(argv[i]);
      HLTInfo("Min Sector set to %d.", fAsyncProcessorQueueDepth);
      iRet+=2;
	}
    else if (argument.CompareTo("-MaxSector")==0){
	  if (++i >= argc)
	  {
	    HLTError("Value missing for -MaxSector parameter");
		return(-EINVAL);
	  }
	  fMaxInitSec = atoi(argv[i]) + 1; //If we want to process sector n, the for loop with < comparison has to go to n+1
      HLTInfo("Max Sector set to %d.", fAsyncProcessorQueueDepth);
      iRet+=2;
	}
    else if (argument.CompareTo("-NoInitialObject")==0){
      fNoInitialObject = true;
      HLTInfo("Not creating initial transform object");
      iRet+=1;
    } else {
      iRet = -EINVAL;
      HLTError("Unknown argument %s",argv[i]);     
    }
  } 
  return iRet;
}


Int_t AliHLTTPCClusterTransformationPrepareComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/) {
	// see header file for class documentation

	AliHLTTPCFastTransformObject* transformMap = NULL;
	if (fTmpFastTransformObject)
	{
		HLTImportant("Shipping initial transformation map");
		//If we have prepared a first transform map in DoInit, we ship this as soon as possible
		transformMap = fTmpFastTransformObject;
		fTmpFastTransformObject = NULL;
	}

	// -- Only use data event
	if (IsDataEvent())
	{
		if (fNewCalibObject == NULL)
		{
			for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeTObject); iter != NULL; iter = GetNextInputObject() )
			{
				AliCDBEntry* tmpEntry = dynamic_cast<AliCDBEntry*>((TObject*) iter);
				
				if (tmpEntry)
				{
					RemoveInputObjectFromCleanupList(tmpEntry);
					fNewCalibObject = tmpEntry;
					break;
				}
				else
				{
					HLTImportant("Transformation Prepare component received object that is no CDBEntry!");
				}
			}
			
			if (fNewCalibObject)
			{
				fAsyncProcessor.QueueAsyncMemberTask(this, &AliHLTTPCClusterTransformationPrepareComponent::AsyncGenerateFastTransformObject, NULL);
			}
		}
		
		//If a new transform map is available from an async creation task, we ship the newest one.
		while (fAsyncProcessor.IsQueuedTaskCompleted())
		{
			if (transformMap) delete transformMap;
			transformMap = (AliHLTTPCFastTransformObject*) fAsyncProcessor.RetrieveQueuedTaskResult();
		}
	}
	
	//If there is a new map, ship it
	if (transformMap)
	{
		PushBack(dynamic_cast<TObject*>(transformMap), AliHLTTPCDefinitions::fgkTPCFastTransformDataObjectDataType | kAliHLTDataOriginTPC);
		delete transformMap;
	}

	return 0;
}

void AliHLTTPCClusterTransformationPrepareComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
  // Get a list of OCDB object description needed for the particular component
  if (!targetMap) return;
  
  // OCDB entries for component arguments

  //!! targetMap->Add(new TObjString(fgkOCDBEntryClusterTransformation), new TObjString("component argument for the charge threshold"));
  
  // OCDB entries to be fetched by the TAXI (access via the AliTPCcalibDB class)
  targetMap->Add(new TObjString("TPC/Calib/Parameters"),    new TObjString("unknown content"));
  targetMap->Add(new TObjString("TPC/Calib/TimeDrift"),     new TObjString("drift velocity calibration"));
  targetMap->Add(new TObjString("TPC/Calib/TimeGain"),     new TObjString("time gain  calibration"));
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
