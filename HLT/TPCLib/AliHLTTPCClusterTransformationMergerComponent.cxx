
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

/** @file   AliHLTTPCClusterTransformationMergerComponent.cxx
    @author Sergey Gorbunov
    @date   
    @brief 
*/

#include "AliHLTTPCClusterTransformationMergerComponent.h"
#include "transform/AliHLTTPCFastTransformObject.h"
#include "AliHLTTPCDefinitions.h"
#include "TList.h"
#include "TClass.h"

using namespace std;

ClassImp(AliHLTTPCClusterTransformationMergerComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCClusterTransformationMergerComponent::AliHLTTPCClusterTransformationMergerComponent() :
fCumulative(0), fTotalInputs(0), fObj(NULL)
{}

AliHLTTPCClusterTransformationMergerComponent::~AliHLTTPCClusterTransformationMergerComponent()
{ 
  // destructor
}

const char* AliHLTTPCClusterTransformationMergerComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCClusterTransformationMerger";
}

AliHLTComponentDataType AliHLTTPCClusterTransformationMergerComponent::GetOutputDataType() { 
  // see header file for class documentation
  return kAliHLTAllDataTypes;
}

void AliHLTTPCClusterTransformationMergerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
  // see header file for class documentation

  list.clear(); 
  list.push_back(kAliHLTAllDataTypes|kAliHLTDataOriginAny);
}

void AliHLTTPCClusterTransformationMergerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
  // see header file for class documentation
  constBase = 0;
  inputMultiplier = 1.0;
}

AliHLTComponent* AliHLTTPCClusterTransformationMergerComponent::Spawn() { 
  return new AliHLTTPCClusterTransformationMergerComponent();
}
	
int AliHLTTPCClusterTransformationMergerComponent::DoInit( int argc, const char** argv ) 
{
  ConfigureFromArgumentString(argc, argv);
  return 0;
}

int AliHLTTPCClusterTransformationMergerComponent::DoDeinit() {
  if (fObj) delete fObj;
  fObj = NULL;
  return 0;
}

int AliHLTTPCClusterTransformationMergerComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/) { 
  return 0;
}

int AliHLTTPCClusterTransformationMergerComponent::ScanConfigurationArgument(int argc, const char** argv){

  // see header file for class documentation

  if (argc<=0) return 0;
  int iRet = 0;
  for( int i=0; i<argc; i++ ){
    TString argument=argv[i];  
    if (argument.CompareTo("-cumulative")==0){
	  fCumulative = 1;
      HLTInfo("Cumulative object merging activated");
      iRet++;
    } else {
      iRet = -EINVAL;
      HLTInfo("Unknown argument %s",argv[i]);     
    }
  } 
  return iRet;
}


int AliHLTTPCClusterTransformationMergerComponent::DoEvent(const AliHLTComponentEventData& evtData, 
					          const AliHLTComponentBlockData* blocks, 
					          AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, 
					          AliHLTUInt32_t& size, 
					          vector<AliHLTComponentBlockData>& outputBlocks )
{
	TObject* returnObj = fObj;
	int nInputs = 0;
	TList mergeList;
	bool writeOutput = false; //For cumulative mode we should make sure we do not send the same merged object again and again

	for ( const TObject *obj = GetFirstInputObject(kAliHLTAllDataTypes); obj != NULL; obj = GetNextInputObject() )
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
			mergeList.Add(nonConstObj);
		}
	}
	
	if (mergeList.GetEntries())
	{
		if (!returnObj->IsA()->GetMethodWithPrototype("Merge", "TCollection*"))
		{
			HLTError("Object does not implement a merge function!");
			return(-1);
		}
		Int_t error = 0;
		TString listHargs;
		listHargs.Form("((TCollection*)0x%lx)", (ULong_t) &mergeList);
		returnObj->Execute("Merge", listHargs.Data(), &error);
		if (error)
		{
			HLTError("Error running merge!");
			return(-1);
		}
		
		nInputs += mergeList.GetEntries();
	}
	
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
		if (!fCumulative) delete returnObj;
	}
	mergeList.Delete();
	return 0;
} // end DoEvent()

void AliHLTTPCClusterTransformationMergerComponent::GetOCDBObjectDescription( TMap* const targetMap)
{
}
