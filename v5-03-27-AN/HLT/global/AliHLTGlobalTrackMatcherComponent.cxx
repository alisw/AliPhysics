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
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "TGeoManager.h"
#include "TRefArray.h"
#include "TString.h"
#include "TMap.h"

/** ROOT macro for the implementation of ROOT specific class methods */
AliHLTGlobalTrackMatcherComponent gAliHLTGlobalTrackMatcherComponent;

ClassImp(AliHLTGlobalTrackMatcherComponent);

AliHLTGlobalTrackMatcherComponent::AliHLTGlobalTrackMatcherComponent() :
  fOCDBEntry("HLT/ConfigHLT/GlobalTrackMatcher"), //TODO
  fMethod(1), //Method 1(PbPb) 2(pp)
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
  inputMultiplier = 1;
}

AliHLTComponent* AliHLTGlobalTrackMatcherComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTGlobalTrackMatcherComponent;
}

int AliHLTGlobalTrackMatcherComponent::DoInit( int argc, const char** argv ) 
{
  Int_t iResult=ConfigureFromCDBTObjString(fOCDBEntry); //MARCEL
    // configure from the command line parameters if specified
  if (iResult>=0 && argc>0) {
    iResult=ConfigureFromArgumentString(argc, argv);
    HLTImportant("Extrapolation Method from argument string:  %d", fMethod);   
  } else if ( iResult >=0 ) {
    HLTImportant("Extrapolation Method from OCDB database entry:  %d", fMethod);   
  } 
  
  
  //BALLE TODO, use command line values to initialise matching vaules
 // init
//   Int_t iResult = argc;
//   iResult = argc;
  
  if(argc > 0){
    HLTWarning("Ignoring all configuration args, starting with: argv %s", argv[0]);
  }

  if(!fClusterReader)
    fClusterReader = new AliHLTCaloClusterReader();

  fBz = GetBz();
  if(fBz == -999999) {
    HLTError("Magnetic field not properly set, current value: %d", fBz);
  }

  if(!fTrackMatcher)
    fTrackMatcher = new AliHLTGlobalTrackMatcher();

  fNEvents = 0;

  if(!fTrackArray){
    fTrackArray = new TObjArray();
    fTrackArray->SetOwner(kFALSE);
  }

  //*** GeoManager ***
   AliCDBPath path("GRP","Geometry","Data");
   AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path);
   if (pEntry) {
      if(!gGeoManager) {
      gGeoManager = (TGeoManager*)pEntry->GetObject();
      }
   }
   else {
      HLTError("can not fetch object \"%s\" from CDB",path.GetPath().Data());
   }
   // ****


  return iResult; 
}
  
int AliHLTGlobalTrackMatcherComponent::DoDeinit() 
{
  // see header file for class documentation
  Int_t iResult = 1;
  
  if(fTrackMatcher)
    delete fTrackMatcher;
  fTrackMatcher = NULL;
  
  if(fClusterReader)
    delete fClusterReader;
  fClusterReader = NULL;
  

  fNEvents = 0;

  return iResult;
}

int AliHLTGlobalTrackMatcherComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/) 
{
  
  //See header file for documentation
  Int_t iResult = 0;
  
  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;
  
  if(!IsDataEvent()){//marcel test
    return 0;//marcel test
  }//marcel test

  fNEvents++;

  //Loop over TPC blocks
  //BALLE TODO check that the tracks in the TObjArray are fine over several blocks
  
   fTrackArray->Clear();
   vector<AliHLTGlobalBarrelTrack> tracks;
   
   for (const AliHLTComponentBlockData* pBlock = GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC); pBlock!=NULL; pBlock=GetNextInputBlock()) {

     if ((iResult=AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(pBlock->fPtr), pBlock->fSize, tracks))>=0) {
       for(UInt_t it = 0; it < tracks.size(); it++) {
	 AliHLTGlobalBarrelTrack track = tracks.at(it);
 	fTrackArray->AddLast(dynamic_cast<TObject*>(&(tracks.at(it))));
	       }
     } else {
	      HLTWarning("Converting tracks to vector failed");
     }
    
     //     //Push the TPC block on, without any changes
     //PushBack(pBlock->fPtr, pBlock->fSize, pBlock->fDataType, pBlock->fSpecification);

   }
    
   AliHLTCaloClusterDataStruct * caloClusterStruct;
   //Get the PHOS Clusters
   vector<AliHLTCaloClusterDataStruct*> phosClustersVector;
   
   for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeCaloCluster | kAliHLTDataOriginPHOS); pBlock!=NULL; pBlock=GetNextInputBlock()) {
     AliHLTCaloClusterHeaderStruct *caloClusterHeader = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(pBlock->fPtr);
     fClusterReader->SetMemory(caloClusterHeader);
     if ( (caloClusterHeader->fNClusters) < 0) {
       HLTWarning("Event has negative number of clusters: %d! Very bad for vector resizing", (Int_t) (caloClusterHeader->fNClusters));
       continue;
     } else {    
       phosClustersVector.reserve( (int) (caloClusterHeader->fNClusters) + phosClustersVector.size() ); 
       while( (caloClusterStruct = fClusterReader->NextCluster()) != 0) {
	 phosClustersVector.push_back(caloClusterStruct);  
       }
     }
   }

    //Get the EMCAL Clusters
   vector<AliHLTCaloClusterDataStruct*> emcalClustersVector;
   for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeCaloCluster | kAliHLTDataOriginEMCAL); pBlock!=NULL; pBlock=GetNextInputBlock()) {
     AliHLTCaloClusterHeaderStruct *caloClusterHeader = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(pBlock->fPtr);
     fClusterReader->SetMemory(caloClusterHeader);
//           HLTInfo("\n EMCAL: estou aqui");//marcel
     if ( (caloClusterHeader->fNClusters) < 0) {
       HLTWarning("Event has negative number of clusters: %d! Very bad for vector resizing", (Int_t) (caloClusterHeader->fNClusters));
       continue;
     } else {    
       emcalClustersVector.reserve( (int) (caloClusterHeader->fNClusters) + emcalClustersVector.size() ); 
       while( (caloClusterStruct = fClusterReader->NextCluster()) != 0) {
	 emcalClustersVector.push_back(caloClusterStruct);  
       }
     }
   }

      iResult = fTrackMatcher->Match(fTrackArray, phosClustersVector, emcalClustersVector, fBz,fMethod); //With Method String


   //Push the blocks on
   for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeCaloCluster | kAliHLTDataOriginAny); pBlock!=NULL; pBlock=GetNextInputBlock()) {
     PushBack(pBlock->fPtr, pBlock->fSize, pBlock->fDataType, pBlock->fSpecification);
    }
   
   fTrackArray->Clear();

   return iResult;
}

int AliHLTGlobalTrackMatcherComponent::ScanConfigurationArgument(int argc, const char** argv) {
  // see header file for class documentation
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -maxpt
  if (argument.CompareTo("-method")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fMethod=argument.Atoi(); // 
    return 2;
  }    

// unknown argument
  return -EINVAL;
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

int AliHLTGlobalTrackMatcherComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/) {
  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) entry=fOCDBEntry;

  return ConfigureFromCDBTObjString(entry);
}

void AliHLTGlobalTrackMatcherComponent::GetOCDBObjectDescription( TMap* const targetMap) {
  
  // Get a list of OCDB object description.
  if (!targetMap) return;
  targetMap->Add(new TObjString(fOCDBEntry),
		 new TObjString(Form("Track-Matcher Method OCDB object") ) 
		 );
}


