// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
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

/** @file   AliHLTTPCOfflineClustererComponent.cxx
    @author Matthias Richter
    @date   
    @brief  Wrapper component to the TPC offline cluster finder
*/

#include "AliHLTTPCOfflineClustererComponent.h"
#include "AliHLTTPCDefinitions.h"
#include "AliRawReaderMemory.h"
#include "AliTPCParam.h"
#include "AliTPCParamSR.h"
#include "AliTPCclustererMI.h"
#include "AliDAQ.h"
#include "TString.h"
#include "TObjArray.h"
#include "TObjString.h"
#include "TTree.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCOfflineClustererComponent)

AliHLTTPCOfflineClustererComponent::AliHLTTPCOfflineClustererComponent()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCOfflineClustererComponent::~AliHLTTPCOfflineClustererComponent()
{
  // see header file for class documentation
}

const char* AliHLTTPCOfflineClustererComponent::GetComponentID()
{
  // see header file for class documentation
  return "TPCOfflineClusterer";
}

void AliHLTTPCOfflineClustererComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  // see header file for class documentation
  list.push_back(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC);
}

AliHLTComponentDataType AliHLTTPCOfflineClustererComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeAliTreeR|kAliHLTDataOriginTPC;
}

void AliHLTTPCOfflineClustererComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // see header file for class documentation
  constBase = 0;inputMultiplier = 1;
}

AliHLTComponent* AliHLTTPCOfflineClustererComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCOfflineClustererComponent;
}

int AliHLTTPCOfflineClustererComponent::DoInit( int argc, const char** argv )
{
  // see header file for class documentation
  int iResult=0;

  TString argument="";
  TString configuration=""; 
  int bMissingParam=0;
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (argument.IsNull()) continue;

  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }

  if (iResult>=0 && !configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  } else {
    iResult=Reconfigure(NULL, NULL);
  }

  return iResult;
}

int AliHLTTPCOfflineClustererComponent::DoDeinit()
{
  // see header file for class documentation
  return 0;
}

int AliHLTTPCOfflineClustererComponent::DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  // see header file for class documentation
  int iResult=0;
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginTPC);
       pBlock!=NULL; 
       pBlock=GetNextInputBlock()) {
    int slice=AliHLTTPCDefinitions::GetMinSliceNr(pBlock->fSpecification);
    int patch=AliHLTTPCDefinitions::GetMinPatchNr(pBlock->fSpecification);
    if (slice!=AliHLTTPCDefinitions::GetMaxSliceNr(pBlock->fSpecification) ||
	patch!=AliHLTTPCDefinitions::GetMaxPatchNr(pBlock->fSpecification)) {
      HLTError("ambiguous readout partition (specification 0x%08x), skipping input block", pBlock->fSpecification);
      break;
    }
    if (slice<0 || slice>35 || patch<0 || patch>5) {
      HLTError("invalid readout partition %d/%d (specification 0x%08x, skipping input block", slice, patch,  pBlock->fSpecification);
      break;
    }

    TTree* pTreeR=new TTree("TreeR", "Reconstructed Points Container");
    AliRawReaderMemory* pRawReader=new AliRawReaderMemory;

    // TODO: choose the right parameter
    AliTPCParam* pParam = new AliTPCParamSR;
    AliTPCclustererMI* pClusterer = new AliTPCclustererMI(pParam);

    if (pTreeR && pRawReader && pClusterer) {
      // setup raw reader and cluster finder
      pRawReader->SetMemory( reinterpret_cast<UChar_t*>( pBlock->fPtr ), pBlock->fSize );
      int ddlId=AliDAQ::DdlIDOffset("TPC");
      if (patch<2) {
	ddlId+=2*slice+patch;
      } else {
	ddlId+=72;
	ddlId+=4*slice+patch;	  
      }
      pRawReader->SetEquipmentID(ddlId);

      // TODO: find out from the data
      //pClusterer->SetOldRCUFormat(kTRUE);

      // run the cluster finder
      pClusterer->SetOutput(pTreeR);
      pClusterer->Digits2Clusters(pRawReader);

      // insert tree into output stream
      PushBack(pTreeR, kAliHLTDataTypeAliTreeR|kAliHLTDataOriginTPC, pBlock->fSpecification);

      if (pClusterer) delete pClusterer; pClusterer=NULL;
      if (pParam)     delete pParam;	 pParam=NULL;
      if (pRawReader) delete pRawReader; pRawReader=NULL;
      if (pTreeR)     delete pTreeR;	 pTreeR=NULL;
    } else {
      iResult=-ENOMEM;
    }
  }

  return iResult;
}

int AliHLTTPCOfflineClustererComponent::Configure(const char* arguments)
{
  // see header file for class documentation
  int iResult=0;
  if (!arguments) return iResult;

  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;

      if (argument.CompareTo("-something")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;

      } else {
	HLTError("unknown argument %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}

int AliHLTTPCOfflineClustererComponent::Reconfigure(const char* /*cdbEntry*/, const char* /*chainId*/)
{
  // see header file for class documentation
  int iResult=0;
  // CDB stuff needs to be implemented
  return iResult;
}
