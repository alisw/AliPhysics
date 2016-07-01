/* This file is property of and copyright by the ALICE HLT Project        * 
* ALICE Experiment at CERN, All rights reserved.                         *
* See cxx source for full Copyright notice                               */

/** @file    AliHLTAsyncTestComponent.cxx
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
#include "AliHLTAsyncTestComponent.h"
#include "AliHLTITSClusterDataFormat.h"
#include "TTree.h"

using namespace std;

ClassImp(AliHLTAsyncTestComponent)

AliHLTAsyncTestComponent::AliHLTAsyncTestComponent() :
AliHLTProcessor(),
fUID(0),
fAsyncProcessor(),
fAsyncProcessorQueueDepth(0),
fAsyncTaskData(0)
{
}

void AliHLTAsyncTestComponent::GetOCDBObjectDescription(TMap*) {}

AliHLTAsyncTestComponent::~AliHLTAsyncTestComponent() {
}

const Char_t* AliHLTAsyncTestComponent::GetComponentID() { 
	return "AsyncTest";
}

void AliHLTAsyncTestComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
	list.push_back(kAliHLTDataTypeESDObject|kAliHLTDataOriginAny);
	list.push_back(kAliHLTDataTypeEOR|kAliHLTDataOriginAny);
}

AliHLTComponentDataType AliHLTAsyncTestComponent::GetOutputDataType() {
	return kAliHLTDataTypeTObject|kAliHLTDataOriginHLT;
}

void AliHLTAsyncTestComponent::GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier ) {
	constBase = 1;
	inputMultiplier = 0.;
}

AliHLTComponent* AliHLTAsyncTestComponent::Spawn() {
	// see header file for class documentation
	return new AliHLTAsyncTestComponent;
}

/*
* ---------------------------------------------------------------------------------
* Protected functions to implement AliHLTComponent's interface.
* These functions provide initialization as well as the actual processing
* capabilities of the component. 
* ---------------------------------------------------------------------------------
*/

int AliHLTAsyncTestComponent::Configure(const char* cdbEntry, const char* chainId, const char *commandLine)
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
Int_t AliHLTAsyncTestComponent::DoInit(Int_t argc, const Char_t** argv)
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

	HLTImportant("AliHLTAsyncTestComponent::DoInit (with QueueDepth %d)", fAsyncProcessorQueueDepth);
	if (fAsyncProcessor.Initialize(fAsyncProcessorQueueDepth)) return(1);

	int* initRetVal;
	retVal = fAsyncProcessor.InitializeAsyncMemberTask(this, &AliHLTAsyncTestComponent::MemberInitializer, NULL, (void**) &initRetVal);
	if (retVal) return(1);
	if (initRetVal == NULL) return(1);
	retVal = *initRetVal;
	delete initRetVal;

	return retVal;
}


// #################################################################################
Int_t AliHLTAsyncTestComponent::DoDeinit() {
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
Int_t AliHLTAsyncTestComponent::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& /*trigData*/) {
	// see header file for class documentation

	HLTImportant("AliHLTAsyncTestComponent::DoEvent");

	// -- Only use data event
	if (IsDataEvent())
	{
		// -- Get ESD object
		// -------------------
		AliESDEvent *esdEvent = NULL;
		for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {
			esdEvent = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
			if( !esdEvent ){ 
				HLTWarning("Wrong ESDEvent object received");
				continue;
			}
			esdEvent->GetStdContent();
			break;
		}
		if (!esdEvent) return(-1);

		static int nEvents=0;
		nEvents++;

		if (nEvents % 2 == 0)
		{
			TString* processObj = new TString("Data Object for Async Task ");
			(*processObj) += nEvents;
			if (nEvents % 3)
			{
				HLTImportant("Starting ASYNC member task");
				fAsyncProcessor.QueueAsyncMemberTask(this, &AliHLTAsyncTestComponent::MemberTestTask, processObj);
			}
			else
			{
				fAsyncProcessor.QueueAsyncTask(&StaticTestTask, processObj);
			}

			if (true)
			{
				TString* output = new TString("Test Output");
				PushBack( dynamic_cast<TObject*>(output), kAliHLTDataTypeTObject | kAliHLTDataOriginHLT,fUID);
				delete output;
			}
		}
	}

	if (!IsDataEvent() && GetFirstInputBlock(kAliHLTDataTypeEOR | kAliHLTDataOriginAny))
	{
		fAsyncProcessor.WaitForTasks(0);
	}

	while (fAsyncProcessor.IsQueuedTaskCompleted())
	{
		TString* retObj = (TString*) fAsyncProcessor.RetrieveQueuedTaskResult();
		HLTImportant("Returned object from async task: %s", retObj->Data());
		delete retObj;
	}

	return 0;
}

void* AliHLTAsyncTestComponent::MemberTestTask(void* data) {
	TString* dataObj = (TString*) data;
	printf("Running ASYNC Member Test Task (%s)\n", dataObj->Data());
	delete dataObj;
	usleep(rand() % 1000 * 500);
	fAsyncProcessor.LockMutex();
	fAsyncTaskData++;	//Secure some pseudo critical section with a mutex
	fAsyncProcessor.UnlockMutex();
	TString* retVal = new TString("Return Object (Member Task");
	*retVal += fAsyncTaskData;
	*retVal += ")";
	return(retVal);
}

//Example of async static function
void* AliHLTAsyncTestComponent::StaticTestTask(void* data) {
	TString* dataObj = (TString*) data;
	printf("Running ASYNC Static Test Task (%s)\n", dataObj->Data());
	delete dataObj;
	usleep(rand() % 1000 * 500);
	TString* retVal = new TString("Return Object (Static Task)");
	return(retVal);
}

//Example function for the initialization inside the async task
void* AliHLTAsyncTestComponent::MemberInitializer(void*) {
	HLTImportant("Running ASYNC Initializer via member function interface");
	fAsyncTaskData = 10; //Initialize different to 0 to check initialization
	//We must communicate via a pointer, so we return a zero integer to say ok
	return(new int(0));
}

// #################################################################################
Int_t AliHLTAsyncTestComponent::Reconfigure(const Char_t* cdbEntry, const Char_t* chainId) {
	// see header file for class documentation

	return(0);
}

// #################################################################################
Int_t AliHLTAsyncTestComponent::ReadPreprocessorValues(const Char_t* /*modules*/) {
	// see header file for class documentation
	ALIHLTERRORGUARD(5, "ReadPreProcessorValues not implemented for this component");
	return 0;
}

int AliHLTAsyncTestComponent::ReadConfigurationString(const char* arguments)
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
