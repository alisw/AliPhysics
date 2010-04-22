//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Svein Lindal <svein.lindal@gmail.com>                 *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTGlobalTrackMatcherComponent.cxx
    @author Svein Lindal
    @brief  Component to match TPC tracks to Calo Clusters
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTProcessor.h"
#include "AliHLTGlobalTrackMatcherComponent.h"
#include "AliHLTGlobalTrackMatcher.h"
#include "TObjArray.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTCaloClusterDataStruct.h"
#include "AliHLTCaloClusterReader.h"


/** ROOT macro for the implementation of ROOT specific class methods */
AliHLTGlobalTrackMatcherComponent gAliHLTGlobalTrackMatcherComponent;

ClassImp(AliHLTGlobalTrackMatcherComponent);

AliHLTGlobalTrackMatcherComponent::AliHLTGlobalTrackMatcherComponent() :
  fTrackMatcher(NULL),
  fNEvents(0),
  fBz(-9999999),
  fClusterReader(NULL),
  fTrackArray(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTGlobalTrackMatcherComponent::~AliHLTGlobalTrackMatcherComponent()
{
  // see header file for class documentation
  if(fTrackMatcher)
    delete fTrackMatcher;
  fTrackMatcher = NULL;

  if(fClusterReader)
    delete fClusterReader;
  fClusterReader = NULL;

}


// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process
const char* AliHLTGlobalTrackMatcherComponent::GetComponentID()
{
  // see header file for class documentation
  return "TrackMatcher";
}

void AliHLTGlobalTrackMatcherComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( kAliHLTDataTypeTrack );
  list.push_back( kAliHLTDataTypeCaloCluster );
}

AliHLTComponentDataType AliHLTGlobalTrackMatcherComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeCaloCluster | kAliHLTDataOriginAny;
}

void AliHLTGlobalTrackMatcherComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 80000;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTGlobalTrackMatcherComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTGlobalTrackMatcherComponent;
}

Int_t AliHLTGlobalTrackMatcherComponent::DoInit( int argc, const char** argv ) {
 
  //BALLE TODO, use command line values to initialise matching vaules
  Int_t iResult = argc;

  
  if(!(argc > 0)){
    HLTError("Component must be configured with detector as first parameter, either \"PHOS\" or \"EMCAL\"");
    return -ENODEV;
  } else {
    TString allArgs = "";
    
    for ( int i = 0; i < argc; i++ ) {
      if ( !allArgs.IsNull() ) allArgs += " ";
      allArgs += argv[i];
    }
    
    HLTInfo("Configuring with string \"%s\"", allArgs.Data());
  }

  if(!strcmp("PHOS", argv[0])) {
    fTrackMatcher = new AliHLTGlobalTrackMatcher();
    fTrackMatcher->SetMaxX( 355. + 70. );
    fTrackMatcher->SetMaxZ( 64. + 70. );
    fTrackMatcher->SetYSign(kFALSE);
    fTrackMatcher->SetRadius(460.);
    fTrackMatcher->SetMatchDistance(40.*40);
  } else if(!strcmp("EMCAL", argv[0])) {
    fTrackMatcher = new AliHLTGlobalTrackMatcher();
    fTrackMatcher->SetMaxX( 500. + 70. );
    fTrackMatcher->SetMaxZ( 500. + 70. );
    fTrackMatcher->SetYSign(kTRUE);
    fTrackMatcher->SetRadius(470.);
    fTrackMatcher->SetMatchDistance(40.*40);
  } else {
    HLTError("Component must be configured with detector as first parameter, either \"PHOS\" or \"EMCAL\"");
    return -ENODEV;
  }

  if(!fTrackArray){
    fTrackArray = new TObjArray();
    fTrackArray->SetOwner(kFALSE);
  }
  
  if(!fClusterReader)
    fClusterReader = new AliHLTCaloClusterReader();
  
  fBz = GetBz();
  if(fBz == -999999) {
    HLTError("Magnetic field not properly set, current value: %d", fBz);
    return -ENODEV;
  }

  return 0; 
}

  
Int_t AliHLTGlobalTrackMatcherComponent::DoDeinit() {
  // see header file for class documentation
  Int_t iResult = 1;
  
  if(fTrackMatcher)
    delete fTrackMatcher;
  fTrackMatcher = NULL;
  
  if(fClusterReader)
    delete fClusterReader;
  fClusterReader = NULL;
  

  return iResult;
}


Int_t AliHLTGlobalTrackMatcherComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)  {
  
  //See header file for documentation

  Int_t iResult = 0;
  
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  fTrackArray->Clear();
  vector<AliHLTGlobalBarrelTrack> tracks;
   
  for (const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC); pBlock!=NULL; pBlock=GetNextInputBlock()) {
    
    if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>=0) {
      for(UInt_t it = 0; it < tracks.size(); it++) {
        	fTrackArray->AddLast(dynamic_cast<TObject*>(&(tracks.at(it))));
      }
    } else {
      HLTError("Converting tracks to vector failed");
    }
    
    ///Push the TPC blocks on, without changes
    PushBack(pBlock->fPtr, pBlock->fSize, pBlock->fDataType, pBlock->fSpecification);
  }

  //Get the clusters
  AliHLTCaloClusterDataStruct * caloClusterStruct;
  vector<AliHLTCaloClusterDataStruct*> clusterVec;

  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeCaloCluster | kAliHLTDataOriginAny); pBlock!=NULL; pBlock=GetNextInputBlock()) {
    AliHLTCaloClusterHeaderStruct *caloClusterHeader = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(pBlock->fPtr);
    fClusterReader->SetMemory(caloClusterHeader);

    
    if ( (caloClusterHeader->fNClusters) < 0) {
      HLTError("Event has negative number of clusters: %d! Very bad for vector resizing", (Int_t) (caloClusterHeader->fNClusters));
      continue;
    } else {    
      clusterVec.reserve(clusterVec.size()+(int) (caloClusterHeader->fNClusters)); 
      while( (caloClusterStruct = fClusterReader->NextCluster()) != 0) {
	clusterVec.push_back(caloClusterStruct);  
      }
    }
  }
          
  iResult = fTrackMatcher->Match(fTrackArray, clusterVec, fBz);
  
  for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeCaloCluster | kAliHLTDataOriginAny); pBlock!=NULL; pBlock=GetNextInputBlock()) {
    PushBack(pBlock->fPtr, pBlock->fSize, pBlock->fDataType, pBlock->fSpecification);
  }

  fTrackArray->Clear();
   
  return iResult;
   
}

// int AliHLTGlobalTrackMatcherComponent::Configure(const char* arguments)
// {
//   Int_t iResult = 1;
//   return iResult;
// }

// int AliHLTGlobalTrackMatcherComponent::Reconfigure(const char* cdbEntry, const char* chainId)
// {
//   Int_t iResult = 1;
//   return iResult;
// }
