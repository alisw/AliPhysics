
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

/** @file   AliHLTRootObjectMergerComponent.cxx
    @author Sergey Gorbunov
    @date   
    @brief 
*/

#include "AliHLTRootObjectMergerComponent.h"
#include "transform/AliHLTTPCFastTransformObject.h"
#include "AliHLTTPCDefinitions.h"
#include "TList.h"
#include "TClass.h"

using namespace std;

ClassImp(AliHLTRootObjectMergerComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTRootObjectMergerComponent::AliHLTRootObjectMergerComponent() :
fCumulative(0), fTotalInputs(0), fObj(NULL), fQueueDepth(0), fAsyncProcessor()
{}

AliHLTRootObjectMergerComponent::~AliHLTRootObjectMergerComponent()
{ 
  // destructor
}

const char* AliHLTRootObjectMergerComponent::GetComponentID() { 
// see header file for class documentation

  return "RootObjectMerger";
}

AliHLTComponentDataType AliHLTRootObjectMergerComponent::GetOutputDataType() { 
  // see header file for class documentation
  return kAliHLTAllDataTypes;
}

void AliHLTRootObjectMergerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
  // see header file for class documentation

  list.clear(); 
  list.push_back(kAliHLTAllDataTypes|kAliHLTDataOriginAny);
}

void AliHLTRootObjectMergerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
  // see header file for class documentation
  constBase = 0;
  inputMultiplier = 1.0;
}

AliHLTComponent* AliHLTRootObjectMergerComponent::Spawn() { 
  return new AliHLTRootObjectMergerComponent();
}
	
int AliHLTRootObjectMergerComponent::DoInit( int argc, const char** argv ) 
{
  ConfigureFromArgumentString(argc, argv);

  HLTImportant("AliHLTRootObjectMergerComponent::DoInit (with QueueDepth %d)", fQueueDepth);
  if (fAsyncProcessor.Initialize(fQueueDepth)) return(1);

  return 0;
}

int AliHLTRootObjectMergerComponent::DoDeinit() {

	if (fAsyncProcessor.GetNumberOfAsyncTasksInQueue())
	{
		HLTError("Cannot deinitialize AsyncProcessor - Still tasks in queue");
		//We wait for all tasks in the queue, fetch the results, and drop them.
		//This might result in a memory leak but will at least shut down the thread properly.
		fAsyncProcessor.WaitForTasks(0);
		while (fAsyncProcessor.IsQueuedTaskCompleted()) fAsyncProcessor.RetrieveQueuedTaskResult();
	}

  fAsyncProcessor.Deinitialize();

  if (fObj) delete fObj;
  fObj = NULL;
  return 0;
}

int AliHLTRootObjectMergerComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) { 
  return 0;
}

int AliHLTRootObjectMergerComponent::ScanConfigurationArgument(int argc, const char** argv){

  // see header file for class documentation

  if (argc<=0) return 0;
  int iRet = 0;
  for( int i=0; i<argc; i++ ){
    TString argument=argv[i];  
    if (argument.CompareTo("-cumulative")==0){
	  fCumulative = 1;
      HLTInfo("Cumulative object merging activated");
      iRet++;
    }
 	else if (argument.CompareTo( "-QueueDepth" ) == 0)
	{
		if (++i >= argc)
		{
			HLTError("QueueDepth value missing");
			continue;
		}
		TString val = argv[i];  
		fQueueDepth = val.Atoi();
		HLTInfo("AliHLTRootObjectMergerComponent Queue Depth set to: %d", fQueueDepth);
		iRet+=2;
		if (fQueueDepth) HLTFatal("AliHLTRootObjectMergerComponent cannot run with QueueDepth != 0 yet, must be synchronized properly!");
	}
    else
	{
      iRet = -EINVAL;
      HLTInfo("Unknown argument %s",argv[i]);     
    }
  } 
  return iRet;
}

void* AliHLTRootObjectMergerComponent::MergeObjects(void* tmpObjects)
{
	MergeObjectStruct* objects = (MergeObjectStruct*) tmpObjects;
	TList* mergeList = objects->fList;
	TObject* returnObj = objects->fObject;
	delete objects;
	
	Int_t error = 0;
	TString listHargs;
	listHargs.Form("((TCollection*)0x%lx)", (ULong_t) mergeList);
	returnObj->Execute("Merge", listHargs.Data(), &error);
	if (error)
	{
		HLTError("Error running merge!");
		return(NULL);
	}
	
	mergeList->Delete();

	return(returnObj);
}

int AliHLTRootObjectMergerComponent::DoEvent(const AliHLTComponentEventData& evtData, 
					          const AliHLTComponentBlockData* blocks, 
					          AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					          AliHLTUInt32_t& size, 
					          vector<AliHLTComponentBlockData>& outputBlocks )
{
	TObject* returnObj = fObj;
	int nInputs = 0;
	TList* mergeList = new TList;
	bool writeOutput = false; //For cumulative mode we should make sure we do not send the same merged object again and again

	for (const TObject *obj = GetFirstInputObject(kAliHLTAllDataTypes); obj != NULL; obj = GetNextInputObject())
	{
		writeOutput = true;
		TObject* nonConstObj;
		if ((nonConstObj = RemoveInputObjectFromCleanupList(obj)) == NULL)
		{
			HLTError("Error taking ownership of object.");
			return(-1);
		}

		if (returnObj == NULL)
		{
			returnObj = nonConstObj;
			if (fCumulative)
			{
				fObj = returnObj;
			}
			nInputs = 1;
		}
		else
		{
			mergeList->Add(nonConstObj);
		}
	}
	
	if (mergeList->GetEntries())
	{
		if (!returnObj->IsA()->GetMethodWithPrototype("Merge", "TCollection*"))
		{
			HLTError("Object does not implement a merge function!");
			return(-1);
		}
		
		MergeObjectStruct* tmp = new MergeObjectStruct;
		tmp->fObject = returnObj;
		tmp->fList = mergeList;
		nInputs += mergeList->GetEntries();
		
		fAsyncProcessor.QueueAsyncMemberTask(this, &AliHLTRootObjectMergerComponent::MergeObjects, tmp);
	}
	
	if (!IsDataEvent() && GetFirstInputBlock(kAliHLTDataTypeEOR | kAliHLTDataOriginAny))
	{
		fAsyncProcessor.WaitForTasks(0);
	}

	while (fAsyncProcessor.IsQueuedTaskCompleted())
	{
		TObject* returnObj = (TObject*) fAsyncProcessor.RetrieveQueuedTaskResult();
		if (writeOutput && returnObj)
		{
			if (fCumulative)
			{
				fTotalInputs += nInputs;
				HLTImportant("Cluster tranformation map merged cumulatively: %d new inputs, %d total inputs", nInputs, fTotalInputs);
			}
			else
			{
				HLTImportant("Cluster tranformation map merged from %d inputs", nInputs);
			}
			PushBack(dynamic_cast<TObject*>(returnObj), GetDataType());
			char tmpType[100];
			GetDataType().PrintDataType(tmpType, 100);
			HLTImportant("Merger Component pushing data type %s (Class name %s)", tmpType, returnObj->ClassName());
		}
		if (!fCumulative) delete returnObj;
	}
	
	return 0;
} // end DoEvent()

void AliHLTRootObjectMergerComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
}
