
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
#include "AliHLTMessage.h"
#include "TList.h"
#include "TClass.h"
#include "AliAnalysisDataContainer.h"

using namespace std;

ClassImp(AliHLTRootObjectMergerComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTRootObjectMergerComponent::AliHLTRootObjectMergerComponent() :
fCumulative(0), fTotalInputs(0), fObj(NULL), fQueueDepth(0), fAsyncProcess(0), fDataType(kAliHLTAllDataTypes|kAliHLTDataOriginAny), fDataTypeSet(false), fAsyncProcessor()
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
  constBase = 1000000;
  inputMultiplier = 1.5;
}

AliHLTComponent* AliHLTRootObjectMergerComponent::Spawn() { 
  return new AliHLTRootObjectMergerComponent();
}
	
int AliHLTRootObjectMergerComponent::DoInit( int argc, const char** argv ) 
{
  ConfigureFromArgumentString(argc, argv);
  
  HLTInfo("AliHLTRootObjectMergerComponent::DoInit (with QueueDepth %d)", fQueueDepth);
  if (fAsyncProcessor.Initialize(fQueueDepth, fAsyncProcess, fAsyncProcess ? 100000000 : 0)) return(1);

  return 0;
}

void* AliHLTRootObjectMergerComponent::cleanup(void*)
{
  if (fObj) delete fObj;
  fObj = NULL;
  return(NULL);
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

  fAsyncProcessor.InitializeAsyncMemberTask(this, &AliHLTRootObjectMergerComponent::cleanup, NULL);
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
			return(-1);
		}
		TString val = argv[i];  
		fQueueDepth = val.Atoi();
		HLTInfo("AliHLTRootObjectMergerComponent Queue Depth set to: %d", fQueueDepth);
		iRet+=2;
	}
	else if (argument.CompareTo( "-DataType" ) == 0)
	{
		if (i + 2 > argc)
		{
			HLTError("DataType missing");
			return(-1);
		}
		
		fDataType = AliHLTComponentDataTypeInitializerWithPadding(argv[i + 1], argv[i + 2]);
		fDataTypeSet = true;
		char tmpType[32];
		fDataType.PrintDataType(tmpType, 32);
		HLTInfo("AliHLTRootObjectMergerComponent Data Type set to: %s", tmpType);
		i += 2;
		iRet += 3;
	}
 	else if (argument.CompareTo( "-AsyncProcess" ) == 0)
	{
		fAsyncProcess = 1;
		iRet++;
	}
    else
	{
      iRet = -EINVAL;
      HLTError("Unknown argument %s",argv[i]);     
    }
  } 
  return iRet;
}

void* AliHLTRootObjectMergerComponent::MergeObjects(void* input)
{
	MergeObjectStruct* objects;
	if (fAsyncProcess)
	{
		AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer* buf = (AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer*) input;
		objects = new MergeObjectStruct;
		objects->fObject = fObj;
		objects->fList = new TList;
		for (int i = 0;i < buf->fNumberOfEntries;i++)
		{
			AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* entry = fAsyncProcessor.GetEntry(buf, i);
			BuildMergeList(objects->fObject, objects->fList, AliHLTMessage::Extract(entry->fPtr, entry->fSize));
		}
	}
	else
	{
		objects = (MergeObjectStruct*) input;
	}
	
	TList* mergeList = objects->fList;
	TObject* returnObj = objects->fObject;
	
	if (returnObj && mergeList->GetEntries())
	{
		Int_t error = 0;
		TString listHargs;
		listHargs.Form("((TCollection*)0x%lx)", (ULong_t) mergeList);
		returnObj->Execute("Merge", listHargs.Data(), &error);
		if (error)
		{
			HLTError("Error running merge!");
			return(NULL);
		}
	}
	
	AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* ptr;
	if (returnObj && (fAsyncProcess ? fAsyncProcessor.PushRequested() : CheckPushbackPeriod()))
	{
		ptr = fAsyncProcessor.SerializeIntoBuffer(returnObj, this);
	}
	else
	{
		ptr = NULL;
	}
	
	ClearBuffers(input);
	if (fAsyncProcess) ClearBuffers(objects, true);
	return(ptr);
}

void AliHLTRootObjectMergerComponent::ClearBuffers(void* buffer, bool isMergeObjectStruct)
{
	if (fAsyncProcess && !isMergeObjectStruct)
	{
		fAsyncProcessor.FreeBuffer((AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer*) buffer);
	}
	else
	{
		MergeObjectStruct* tmp = (MergeObjectStruct*) buffer;
		
		TIter next(tmp->fList);
		while (TObject *obj = next())
		{
			TObjArray* objArray = (TObjArray*) dynamic_cast<const TObjArray*>(obj);
			if (objArray != NULL)
			{
				for (int i = 0;i <= objArray->GetLast();i++)
				{
					AliAnalysisDataContainer* tmpContainer = (AliAnalysisDataContainer*) dynamic_cast<const AliAnalysisDataContainer*>((*objArray)[i]);
					if (tmpContainer != NULL)
					{
						tmpContainer->SetDataOwned(kTRUE);
						tmpContainer->DeleteData();
					}
				}
				objArray->Delete();
			}
		}
		
		tmp->fList->Delete();
		delete tmp->fList;
		if (!fCumulative) delete tmp->fObject;
		delete tmp;
	}
}

int AliHLTRootObjectMergerComponent::BuildMergeList(TObject*& returnObj, TList*& mergeList, TObject* obj)
{
	if (obj == NULL)
	{
		HLTError("NULL Ptr to object received");
		return(-1);
	}

	if (!obj->IsA()->GetMethodWithPrototype("Merge", "TCollection*"))
	{
		HLTError("Object does not implement a merge function!");
		return(-1);
	}

	if (returnObj == NULL)
	{
		returnObj = obj;
		if (fCumulative)
		{
			fObj = returnObj;
		}
	}
	else
	{
		mergeList->Add(obj);
	}

	return(0);
}

int AliHLTRootObjectMergerComponent::DoEvent(const AliHLTComponentEventData& evtData, 
					          const AliHLTComponentBlockData* blocks, 
					          AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					          AliHLTUInt32_t& size, 
					          vector<AliHLTComponentBlockData>& outputBlocks )
{
	if (evtData.fBlockCnt == 0) return(0);
	
	TList* mergeList = new TList;
	size_t inputSize = 0;
	int nInputs = 0;
	
	void* objectForAsyncProcessor;
	if (fAsyncProcess)
	{
		AliHLTAsyncProcessor::AliHLTAsyncProcessorMultiBuffer* buf = NULL;
		for (const AliHLTComponentBlockData* blk = GetFirstInputBlock(fDataType); blk != NULL; blk = GetNextInputBlock())
		{
			if (blk->GetDataType() == (kAliHLTAnyDataType | kAliHLTDataOriginPrivate)) continue;
			if (buf == NULL) buf = fAsyncProcessor.AllocateMultiBuffer();
			if (buf == NULL)
			{
				HLTError("Error obtaining multi-part buffer");
				return(1);
			}
			nInputs++;
			inputSize += blk->fSize;
			if (fAsyncProcessor.AddBuffer(buf, blk->fSize, blk->fPtr) == NULL)
			{
				HLTError("Error adding input to asynchronous buffer");
			}
		}
		objectForAsyncProcessor = buf;
	}
	else
	{
		MergeObjectStruct* tmp = NULL;
		for (const TObject *obj = GetFirstInputObject(fDataType); obj != NULL; obj = GetNextInputObject())
		{
			if (tmp == NULL)
			{
				tmp = new MergeObjectStruct;
				tmp->fObject = fObj;
				tmp->fList = new TList;
			}

			nInputs++;
			inputSize += blocks[GetCurrentInputBlockIndex()].fSize;
			TObject* nonConstObj;
			if ((nonConstObj = RemoveInputObjectFromCleanupList(obj)) == NULL)
			{
				HLTError("Error taking ownership of object.");
				return(-1);
			}
			BuildMergeList(tmp->fObject, tmp->fList, nonConstObj);
		}
		objectForAsyncProcessor = tmp;
	}
	
	if (nInputs)
	{
		if (!fDataTypeSet)
		{
			fDataType = GetDataType();
			char tmpType[32];
			fDataType.PrintDataType(tmpType, 32);
			HLTInfo("Automatic data type setting set to %s", tmpType);
			fDataTypeSet = true;
		}
		
		if (fCumulative)
		{
			fTotalInputs += nInputs;
			HLTInfo("Root objects merging cumulatively: %d new inputs (%d bytes), %d total inputs", nInputs, inputSize, fTotalInputs);
		}
		else
		{
			HLTInfo("Root objects merging from %d inputs", nInputs);
		}
		if (fAsyncProcessor.QueueAsyncMemberTask(this, &AliHLTRootObjectMergerComponent::MergeObjects, objectForAsyncProcessor))
		{
			ClearBuffers(objectForAsyncProcessor);
		}
	}
	
	if (!IsDataEvent() && GetFirstInputBlock(kAliHLTDataTypeEOR | kAliHLTDataOriginAny))
	{
		fAsyncProcessor.WaitForTasks(0);
	}

	while (fAsyncProcessor.IsQueuedTaskCompleted())
	{
		AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer* returnObj = (AliHLTAsyncProcessor::AliHLTAsyncProcessorBuffer*) fAsyncProcessor.RetrieveQueuedTaskResult();
		if (returnObj)
		{
			int pushresult = PushBack(returnObj->fPtr, returnObj->fSize, fDataType);
			fAsyncProcessor.FreeBuffer(returnObj);
			if (pushresult > 0)
			{
				char tmpType[32];
				fDataType.PrintDataType(tmpType, 32);
				HLTInfo("Merger Component pushing (%s, %d bytes)", tmpType, pushresult);
			}
		}
		else if (fAsyncProcess && CheckPushbackPeriod())
		{
			fAsyncProcessor.RequestPush();
		}
	}
	
	return 0;
}

void AliHLTRootObjectMergerComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
}
