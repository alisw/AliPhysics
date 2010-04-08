//-*- Mode: C++ -*-
 /**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Svein Lindal, Oeystein Djuvsland                      * 
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** 
 * @file   AliHLTCaloHistoComponent.cxx
 * @author Svein Lindal
 * @date   
 * @brief  A physics histogram producer component for Calo HLT
*/

#if __GNUC__>= 3
using namespace std;
#endif


#include "AliHLTCaloHistoComponent.h"
#include "AliHLTCaloHistoCellEnergy.h"
#include "AliHLTCaloHistoClusterEnergy.h"
#include "AliHLTCaloHistoInvMass.h"
#include "AliHLTCaloHistoMatchedTracks.h"
#include "AliESDEvent.h"
#include "TRefArray.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"
#include "TObjArray.h"

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


AliHLTCaloHistoComponent gAliHLTCaloHistoComponent;

ClassImp(AliHLTCaloHistoComponent);

AliHLTCaloHistoComponent::AliHLTCaloHistoComponent() :
  AliHLTProcessor(),
  fClusterReader(NULL),
  fEmcalClustersArray(NULL),
  fPhosClustersArray(NULL),
  fPhosHistogramArray(NULL),
  fEmcalHistogramArray(NULL),
  fDoEmcal(kFALSE),
  fDoPhos(kFALSE)
{
  //see header file for documentation
}

AliHLTCaloHistoComponent::~AliHLTCaloHistoComponent() 
{
  //see header file for documentation
  //Deinit();
}

Int_t AliHLTCaloHistoComponent::DoInit(int argc, const char** argv ) {
  //see header file for documentation
  
  
  fEmcalHistogramArray = new TObjArray();
  fEmcalHistogramArray->SetOwner(kTRUE);
  fPhosHistogramArray = new TObjArray();
  fPhosHistogramArray->SetOwner(kTRUE);
  
  bool doPhos = true;
  bool doEmcal = true;

  
  for(int i = 0; i < argc; i++) {
    
    if(!strcmp("-phos", argv[i])) {
      fDoPhos = kTRUE;
      doPhos = true;
      doEmcal = false;
    }
    if(!strcmp("-emcal", argv[i])) {
      fDoEmcal = kTRUE;
      doEmcal = true;
      doPhos = false;
    }
    if(!strcmp("-both", argv[i])) {
      fDoEmcal = kTRUE;
      fDoPhos = kTRUE;
      doEmcal = true;
      doPhos = true;
    }
    
    if(!strcmp("-clusterenergy", argv[i])){
      if(doEmcal){
	AliHLTCaloHistoClusterEnergy * histo = new AliHLTCaloHistoClusterEnergy("EMCAL");
	fEmcalHistogramArray->AddLast(dynamic_cast<TObject*>(histo));
	HLTImportant("Adding EMCAL cluster energy histogram");
      }	
      if(doPhos){
	AliHLTCaloHistoClusterEnergy * histo = new AliHLTCaloHistoClusterEnergy("PHOS");
	fPhosHistogramArray->AddLast(dynamic_cast<TObject*>(histo));
	HLTImportant("Adding PHOS cluster energy histogram");
      }
    } 

    if(!strcmp("-invariantmass", argv[i])){
      if(doEmcal){
	AliHLTCaloHistoInvMass * histo = new AliHLTCaloHistoInvMass("EMCAL");
	fEmcalHistogramArray->AddLast(dynamic_cast<TObject*>(histo));
	HLTImportant("Adding EMCAL invariant mass histogram");
      }	
      if(doPhos){
	AliHLTCaloHistoInvMass * histo = new AliHLTCaloHistoInvMass("PHOS");
	fPhosHistogramArray->AddLast(dynamic_cast<TObject*>(histo));
	HLTImportant("Adding PHOS invariant mass histogram");
      }
    } 
   
    if(!strcmp("-matchedtracks", argv[i])) {
      if(doEmcal){
	AliHLTCaloHistoMatchedTracks * histo = new AliHLTCaloHistoMatchedTracks("EMCAL");
	fEmcalHistogramArray->AddLast(dynamic_cast<TObject*>(histo));
	HLTImportant("Adding EMCAL track-matching histograms");
      }	
      if(doPhos){
	AliHLTCaloHistoMatchedTracks * histo = new AliHLTCaloHistoMatchedTracks("PHOS");
	fPhosHistogramArray->AddLast(dynamic_cast<TObject*>(histo));
	HLTImportant("Adding PHOS track-matching histograms");
      }
    }
  }
    
  if(fDoPhos)
    fPhosClustersArray = new TRefArray();
  if(fDoEmcal)
    fEmcalClustersArray = new TRefArray();
  
  fClusterReader = new AliHLTCaloClusterReader();
  
  return 0;
}


Int_t AliHLTCaloHistoComponent::DoDeinit()
{ 
  //see header file for documentation

  if(fEmcalClustersArray)
    delete fEmcalClustersArray;
  fEmcalClustersArray = NULL;

  if(fPhosClustersArray)
    delete fPhosClustersArray;
  fPhosClustersArray = NULL;

  //Deleting these should also destroy histograms!!?
  if(fEmcalHistogramArray)
    delete fEmcalHistogramArray;
  fEmcalHistogramArray = NULL;
 
 if(fPhosHistogramArray)
    delete fPhosHistogramArray;
  fPhosHistogramArray = NULL;
 
  return 0;
}

const char* AliHLTCaloHistoComponent::GetComponentID()
{
  //see header file for documentation
  return "CaloPhysicsHistos";
}


void
AliHLTCaloHistoComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{ 
  //see header file for documentation
  list.clear();
  list.push_back( kAliHLTDataTypeESDObject | kAliHLTDataOriginOut );
  list.push_back( kAliHLTDataTypeCaloCluster | kAliHLTDataOriginEMCAL );
  list.push_back( kAliHLTDataTypeCaloCluster | kAliHLTDataOriginPHOS );

  //   list.push_back(AliHLTPHOSDefinitions::fgkClusterDataType);
  //   list.push_back(AliHLTPHOSDefinitions::fgkESDCaloClusterDataType);
  //   list.push_back(AliHLTPHOSDefinitions::fgkESDCaloCellsDataType);

}

AliHLTComponentDataType AliHLTCaloHistoComponent::GetOutputDataType()
{
  //see header file for documentation
  return kAliHLTDataTypeHistogram  | kAliHLTDataOriginAny ;
}


void AliHLTCaloHistoComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //see header file for documentation
  constBase = 30;
  inputMultiplier = 5;
}

AliHLTComponent* AliHLTCaloHistoComponent::Spawn() {
  //see header file for documentation
  return new AliHLTCaloHistoComponent();
}

Int_t AliHLTCaloHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/) {

  //see header file for documentation
  Int_t iResult = 0;


  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  
  if (fDoEmcal) {
    //    HLTInfo("Processing EMCAL blocks");
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock( kAliHLTDataTypeCaloCluster | kAliHLTDataOriginEMCAL ); pBlock!=NULL; pBlock=GetNextInputBlock()) {
      ProcessBlocks(pBlock, fEmcalHistogramArray);
    }
  }

  if (fDoPhos) {
    //HLTInfo("Processing PHOS blocks");
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock( kAliHLTDataTypeCaloCluster | kAliHLTDataOriginPHOS ); pBlock!=NULL; pBlock=GetNextInputBlock()) {
      ProcessBlocks(pBlock, fPhosHistogramArray);
    }
  }

  //Push histos
  for(int ih = 0; ih < fPhosHistogramArray->GetEntriesFast(); ih++) {
    //HLTInfo("Pushing PHOS histograms");
    PushBack(static_cast<AliHLTCaloHistoProducer*>(fPhosHistogramArray->At(ih))->GetHistograms(), kAliHLTDataTypeHistogram | kAliHLTDataOriginPHOS );
  }
  for(int ih = 0; ih < fEmcalHistogramArray->GetEntriesFast(); ih++) {
    //HLTInfo("Pushing EMCAL histograms");
    PushBack(static_cast<AliHLTCaloHistoProducer*>(fEmcalHistogramArray->At(ih))->GetHistograms(), kAliHLTDataTypeHistogram | kAliHLTDataOriginEMCAL );
  }

  return iResult;

}


Int_t AliHLTCaloHistoComponent::ProcessBlocks(const AliHLTComponentBlockData * pBlock, TObjArray * histoArray) {

  Int_t iResult = 0;

  AliHLTCaloClusterDataStruct * clusterStruct;
  vector<AliHLTCaloClusterDataStruct*> clustersVector;
  
  AliHLTCaloClusterHeaderStruct *clusterHeader = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(pBlock->fPtr);
  fClusterReader->SetMemory(clusterHeader);
     
  if ( (clusterHeader->fNClusters) < 0) {
    HLTError("Event has negative number of clusters: %d! Very bad for vector resizing", (Int_t) (clusterHeader->fNClusters));
    return -1;
  }
  
  clustersVector.resize((int) (clusterHeader->fNClusters)); 
  Int_t nClusters = 0;
  
  while( (clusterStruct = fClusterReader->NextCluster()) != 0) {
    clustersVector[nClusters++] = clusterStruct;  
  }
  
  nClusters = clusterHeader->fNClusters;
  
  for(int ih = 0; ih < histoArray->GetEntriesFast(); ih++) {
    iResult = static_cast<AliHLTCaloHistoProducer*>(histoArray->At(ih))->FillHistograms(nClusters, clustersVector);
  }
  
  return iResult;
}
