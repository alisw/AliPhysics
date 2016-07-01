/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliHLTAsyncCalibrationComponent.cxx
@author  David Rohr (drohr@cern.ch)
*/

#include "TMap.h"
#include "TSystem.h"
#include "TTimeStamp.h"
#include "TObjString.h"
#include "TString.h"
#include "TH1F.h"
#include "TList.h"
#include "AliESDtrackCuts.h"
#include "AliESDEvent.h"
#include "AliHLTErrorGuard.h"
#include "AliHLTDataTypes.h"
#include "AliHLTAsyncCalibrationComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "TTree.h"

using namespace std;

ClassImp(AliHLTAsyncCalibrationComponent)

AliHLTAsyncCalibrationComponent::AliHLTAsyncCalibrationComponent() :
AliHLTProcessor(),
fAsyncProcessor(),
fAsyncProcessorQueueDepth(0),
fNEvents(0)
{
}

void AliHLTAsyncCalibrationComponent::GetOCDBObjectDescription(TMap*) {}

AliHLTAsyncCalibrationComponent::~AliHLTAsyncCalibrationComponent() {
}

const Char_t* AliHLTAsyncCalibrationComponent::GetComponentID() { 
	return "AsyncCalibration";
}

void AliHLTAsyncCalibrationComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
	list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginAny);
	list.push_back(kAliHLTDataTypeEOR|kAliHLTDataOriginAny);
}

AliHLTComponentDataType AliHLTAsyncCalibrationComponent::GetOutputDataType() {
	return kAliHLTDataTypeTObject;
}

void AliHLTAsyncCalibrationComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
	constBase = 10000000;
	inputMultiplier = 0.;
}

AliHLTComponent* AliHLTAsyncCalibrationComponent::Spawn() {
	// see header file for class documentation
	return new AliHLTAsyncCalibrationComponent;
}

/*
* ---------------------------------------------------------------------------------
* Protected functions to implement AliHLTComponent's interface.
* These functions provide initialization as well as the actual processing
* capabilities of the component. 
* ---------------------------------------------------------------------------------
*/

int AliHLTAsyncCalibrationComponent::Configure(const char* cdbEntry, const char* chainId, const char *commandLine)
{
	if (commandLine && commandLine[0] != '\0')
	{
		HLTInfo("received configuration string from HLT framework: \"%s\"", commandLine);
		int retVal = ReadConfigurationString( commandLine );
		if (retVal) return(retVal);
	}

	return(0);
} 

// #################################################################################
Int_t AliHLTAsyncCalibrationComponent::DoInit(Int_t argc, const Char_t** argv)
{
	// see header file for class documentation

	// Configure the component
	TString arguments = "";
	for ( int i = 0; i < argc; i++ ) {
		if ( !arguments.IsNull() ) arguments += " ";
		arguments += argv[i];
	}

	int retVal = Configure(NULL, NULL, arguments.Data()); 
	if (retVal) return(retVal);

	HLTImportant("AliHLTAsyncCalibrationComponent::DoInit (with QueueDepth %d)", fAsyncProcessorQueueDepth);
	if (fAsyncProcessor.Initialize(fAsyncProcessorQueueDepth)) return(1);

	//We do not need an initializer for the async part yet
	/*
	int* initRetVal;
	retVal = fAsyncProcessor.InitializeAsyncMemberTask(this, &AliHLTAsyncCalibrationComponent::MemberInitializer, NULL, (void**) &initRetVal);
	if (retVal) return(1);
	if (initRetVal == NULL) return(1);
	retVal = *initRetVal;
	delete initRetVal;
	*/

	return retVal;
}


// #################################################################################
Int_t AliHLTAsyncCalibrationComponent::DoDeinit() {
	// see header file for class documentation
	if (fAsyncProcessor.GetNumberOfAsyncTasksInQueue())
	{
		HLTError("Cannot deinitialize AsyncProcessor - Still tasks in queue");
		//We wait for all tasks in the queue, fetch the results, and drop them.
		//This might result in a memory leak but will at least shut down the thread properly.
		fAsyncProcessor.WaitForTasks(0);
		while (fAsyncProcessor.IsQueuedTaskCompleted()) fAsyncProcessor.RetrieveQueuedTaskResult();
	}

	fAsyncProcessor.Deinitialize();
	return 0;
}

// #################################################################################
Int_t AliHLTAsyncCalibrationComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/) {
	// see header file for class documentation

	// -- Only use data event
	if (IsDataEvent())
	{
		HLTImportant("AliHLTAsyncCalibrationComponent::DoEvent");
		// -- Get ESD object
		// -------------------
		AliESDEvent *esdEvent = NULL;
		for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() )
		{
			esdEvent = dynamic_cast<AliESDEvent*>(const_cast<TObject*>(iter));
			if( !esdEvent ){ 
				HLTWarning("Wrong ESDEvent object received");
				continue;
			}
			esdEvent->GetStdContent();
			break;
		}
		if (!esdEvent) return(-1);

		if (esdEvent->GetNumberOfTracks() > 10)
		{
			HLTImportant("Event contains %d tracks, running calibration", esdEvent->GetNumberOfTracks());
			if (RemoveInputObjectFromCleanupList(esdEvent) == NULL)
			{
				HLTError("Error taking ownership for esdEvent, cannot queue async calibration task.");
			}
			else
			{
				if (fNEvents++ % 3 == 0)
				{
					fAsyncProcessor.QueueAsyncMemberTask(this, &AliHLTAsyncCalibrationComponent::CalibrationTask, esdEvent);
				}
				else
				{
					delete esdEvent;
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
		TObject* retObj = (TObject*) fAsyncProcessor.RetrieveQueuedTaskResult();
		if (PushBack(retObj, kAliHLTDataTypeTObject ))
		{
			HLTError("Error pushing out calibration object");
		}
		delete retObj;
	}

	return 0;
}

void* AliHLTAsyncCalibrationComponent::CalibrationTask(void* data) {
	AliESDEvent* esdEvent = dynamic_cast<AliESDEvent*>((TObject*) data);
	if (fAsyncProcessorQueueDepth) usleep(rand() % 1000 * 500);
	
	TH1F* retVal = new TH1F("calib", "calib", 100, 0, 100);
	retVal->Fill(esdEvent->GetNumberOfTracks());
	return(retVal);
}

// #################################################################################
Int_t AliHLTAsyncCalibrationComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
	// see header file for class documentation

	return(0);
}

// #################################################################################
Int_t AliHLTAsyncCalibrationComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
	// see header file for class documentation
	ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
	return 0;
}

int AliHLTAsyncCalibrationComponent::ReadConfigurationString(const char* arguments)
{
	// Set configuration parameters for the component from the string

	int iResult = 0;
	if (!arguments) return iResult;

	TString allArgs = arguments;
	TString argument;
	int bMissingParam = 0;

	TObjArray* pTokens = allArgs.Tokenize(" ");

	int nArgs =  pTokens ? pTokens->GetEntries() : 0;

	for (int i = 0; i < nArgs; i++)
	{
		argument = ( ( TObjString* )pTokens->At( i ) )->GetString();
		if ( argument.IsNull() ) continue;

		/*if (argument.CompareTo( "-GlobalTracking" ) == 0) {
		fGlobalTracking = 1;
		HLTImportant( "Global Tracking Activated" );
		continue;
		}*/

		if (argument.CompareTo( "-QueueDepth" ) == 0)
		{
			if ((bMissingParam = (++i >= pTokens->GetEntries()))) break;
			fAsyncProcessorQueueDepth = ((TObjString*) pTokens->At(i))->GetString().Atoi();
			HLTInfo("AsyncProcessor Queue Depth set to: %d", fAsyncProcessorQueueDepth);
			continue;
		}

		HLTError("Unknown option \"%s\"", argument.Data());
		iResult = -EINVAL;
	}
	delete pTokens;

	if ( bMissingParam )
	{
		HLTError("Specifier missed for parameter \"%s\"", argument.Data());
		iResult = -EINVAL;
	}

	return iResult;
}
