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


// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


ClassImp(AliHLTCaloHistoComponent);

AliHLTCaloHistoComponent::AliHLTCaloHistoComponent() :
  AliHLTProcessor(),
  fEmcalClustersArray(NULL),
  fPhosClustersArray(NULL),
  fDoPhos(kFALSE), 
  fDoEmcal(kFALSE),
  fDoCellEnergy(kFALSE),
  fDoClusterEnergy(kFALSE),
  fDoInvariantMass(kFALSE),
  fDoMatchedTracks(kFALSE),
  fPhosCellEnergyHistProducer(NULL),
  fEmcalCellEnergyHistProducer(NULL),
  fPhosClusterEnergyHistProducer(NULL),
  fEmcalClusterEnergyHistProducer(NULL),
  fPhosInvariantMassHistProducer(NULL),
  fEmcalInvariantMassHistProducer(NULL),
  fPhosMatchedTracksHistProducer(NULL),
  fEmcalMatchedTracksHistProducer(NULL),
  fClusterReader(NULL)
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
  

  
  for(int i = 0; i < argc; i++) {
    if(!strcmp("-cellenergy", argv[i])) fDoCellEnergy = true;
    if(!strcmp("-clusterenergy", argv[i])) fDoClusterEnergy = true;
    if(!strcmp("-invariantmass", argv[i])) fDoInvariantMass= true;
    if(!strcmp("-matchedtracks", argv[i])) fDoMatchedTracks = true;
    if(!strcmp("-phos", argv[i])) fDoPhos = true;
    if(!strcmp("-emcal", argv[i])) fDoEmcal = true;
  }



  
  //PHOS
  if(fDoPhos){

    fPhosClustersArray = new TRefArray();
    
    if(fDoInvariantMass) 
      fPhosInvariantMassHistProducer = new AliHLTCaloHistoInvMass("PHOS");
    
    
    if(fDoMatchedTracks)
      fPhosMatchedTracksHistProducer = new AliHLTCaloHistoMatchedTracks("PHOS");
  
  
    if(fDoClusterEnergy)
      fPhosClusterEnergyHistProducer = new AliHLTCaloHistoClusterEnergy("PHOS");
  
  
    if(fDoCellEnergy) 
      fPhosCellEnergyHistProducer = new AliHLTCaloHistoCellEnergy("PHOS");

  }
  

  //EMCAL
  if(fDoEmcal) {

    fEmcalClustersArray = new TRefArray();
  
    if(fDoInvariantMass)
      fEmcalInvariantMassHistProducer = new AliHLTCaloHistoInvMass("EMCAL");

    
    if(fDoMatchedTracks)
      fEmcalMatchedTracksHistProducer = new AliHLTCaloHistoMatchedTracks("EMCAL");

  
    if(fDoClusterEnergy)
      fEmcalClusterEnergyHistProducer = new AliHLTCaloHistoClusterEnergy("EMCAL");
    

    if(fDoCellEnergy) 
      fEmcalCellEnergyHistProducer = new AliHLTCaloHistoCellEnergy("EMCAL");
    
  }
  
  fClusterReader = new AliHLTCaloClusterReader();

  return 0;
}


Int_t AliHLTCaloHistoComponent::DoDeinit()
{ 
  //see header file for documentation


  //Clusters Arrays
  if(fEmcalClustersArray)
    delete fEmcalClustersArray;
  fEmcalClustersArray = NULL;

  if(fPhosClustersArray)
    delete fPhosClustersArray;
  fPhosClustersArray = NULL;


  //CellEnergy
  if(fPhosCellEnergyHistProducer)
    delete fPhosCellEnergyHistProducer;
  fPhosCellEnergyHistProducer = NULL;


  if(fEmcalCellEnergyHistProducer)
    delete fEmcalCellEnergyHistProducer;
  fEmcalCellEnergyHistProducer = NULL;


  //ClusterEnergy
  if(fPhosClusterEnergyHistProducer)
    delete fPhosClusterEnergyHistProducer;
  fPhosClusterEnergyHistProducer = NULL;

  if(fEmcalClusterEnergyHistProducer)
    delete fEmcalClusterEnergyHistProducer;
  fEmcalClusterEnergyHistProducer = NULL;


  //Invariant mass histogram producers
  if(fPhosInvariantMassHistProducer)
    delete  fPhosInvariantMassHistProducer;
  fPhosInvariantMassHistProducer = NULL;

  if(fEmcalInvariantMassHistProducer)
    delete  fEmcalInvariantMassHistProducer;
  fEmcalInvariantMassHistProducer = NULL;


  //Matched track histogram producers
  if(fEmcalMatchedTracksHistProducer)
    delete fEmcalMatchedTracksHistProducer;
  fEmcalMatchedTracksHistProducer = NULL;

  if(fPhosMatchedTracksHistProducer)
    delete fPhosMatchedTracksHistProducer;
  fPhosMatchedTracksHistProducer = NULL;


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



  //Get the clusters from struct input
  AliHLTCaloClusterDataStruct * clusterStruct;
  vector<AliHLTCaloClusterDataStruct*> clustersVector;
  
  if (fDoEmcal) {
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock( kAliHLTDataTypeCaloCluster | kAliHLTDataOriginEMCAL ); pBlock!=NULL; pBlock=GetNextInputBlock()) {

      //Check for origin and continue if not do this detector

      AliHLTCaloClusterHeaderStruct *clusterHeader = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(pBlock->fPtr);
      fClusterReader->SetMemory(clusterHeader);
     
      if ( (clusterHeader->fNClusters) < 0) {
	HLTWarning("Event has negative number of clusters: %d! Very bad for vector resizing", (Int_t) (clusterHeader->fNClusters));
      } else {    
	HLTInfo("Event has positive number of clusters: %d", (Int_t ) (clusterHeader->fNClusters));

	//BALLE, TODO, make it able to do EMCAL as well!!!
	clustersVector.resize((int) (clusterHeader->fNClusters)); 
	Int_t nClusters = 0;
	while( (clusterStruct = fClusterReader->NextCluster()) != 0) {
	  clustersVector[nClusters++] = clusterStruct;  
	}
      
	nClusters = clusterHeader->fNClusters;

	//	if(fDoMatchedTracks)
	  //	  fEmcalMatchedTracksHistProducer->FillHistograms(nClusters, clustersVector);
       
	if(fDoInvariantMass)
	  fEmcalInvariantMassHistProducer->FillHistograms(nClusters, clustersVector);
       
	//	if(fDoClusterEnergy)
	  //fEmcalClusterEnergyHistProducer->FillHistograms(nClusters, clustersVector);
       
	//if(fDoCellEnergy)
	  //fEmcalCellEnergyHistProducer->FillHistograms(nClusters, clustersVector);
       
      }
 
      PushBack(pBlock->fPtr, pBlock->fSize, pBlock->fDataType, pBlock->fSpecification);
      //PushBack(pBlock->fPtr, pBlock->fSize, kAliHLTDataTypeCaloCluster | kAliHLTDataOriginAny );
 
    }
     
    if(iResult <0) {
      HLTWarning("Error in track matcher");
    }
  }


  clustersVector.clear();
  if (fDoPhos) {
    for (const AliHLTComponentBlockData* pBlock=GetFirstInputBlock( kAliHLTDataTypeCaloCluster | kAliHLTDataOriginPHOS ); pBlock!=NULL; pBlock=GetNextInputBlock()) {

      //Check for origin and continue if not do this detector

      AliHLTCaloClusterHeaderStruct *clusterHeader = reinterpret_cast<AliHLTCaloClusterHeaderStruct*>(pBlock->fPtr);
      fClusterReader->SetMemory(clusterHeader);
     
      if ( (clusterHeader->fNClusters) < 0) {
	HLTWarning("Event has negative number of clusters: %d! Very bad for vector resizing", (Int_t) (clusterHeader->fNClusters));
      } else {    
	HLTInfo("Event has positive number of clusters: %d", (Int_t ) (clusterHeader->fNClusters));

	//BALLE, TODO, make it able to do PHOS as well!!!
	clustersVector.resize((int) (clusterHeader->fNClusters)); 
	Int_t nClusters = 0;
	while( (clusterStruct = fClusterReader->NextCluster()) != 0) {
	  clustersVector[nClusters++] = clusterStruct;  
	}
      
	nClusters = clusterHeader->fNClusters;

	//	if(fDoMatchedTracks)
	//fPhosMatchedTracksHistProducer->FillHistograms(nClusters, clustersVector);
       
	if(fDoInvariantMass)
	  fPhosInvariantMassHistProducer->FillHistograms(nClusters, clustersVector);
       
// 	if(fDoClusterEnergy)
// 	  fPhosClusterEnergyHistProducer->FillHistograms(nClusters, clustersVector);
       
// 	if(fDoCellEnergy)
// 	  fPhosCellEnergyHistProducer->FillHistograms(nClusters, clustersVector);
       
      }
      PushBack(pBlock->fPtr, pBlock->fSize, pBlock->fDataType, pBlock->fSpecification);
    }
     
    if(iResult <0) {
      HLTWarning("Error in track matcher");
    }
    
  }




  
//   if( fUID == 0 ){
//     TTimeStamp t;
//     fUID = ( gSystem->GetPid() + t.GetNanoSec())*10 + evtData.fEventID;
//   }

   //TODO Set up for ESD input
//    AliESDCaloCluster * esdCluster;
//    vector<AliESDCaloCluster *> phosEsdClustersVector;
//    vector<AliESDCaloCluster *> EmcalEsdClustersVector;
//    for ( const TObject *iter = GetFirstInputObject(kAliHLTDataTypeESDObject); iter != NULL; iter = GetNextInputObject() ) {
     
//      AliESDEvent *event = dynamic_cast<AliESDEvent*>(const_cast<TObject*>( iter ) );
//      event->GetStdContent();
     

   
     

    
//     //EMCAL
//     if(fDoEmcal){
//       FillClustersToVector(nec, fEmcalClustersArray, emcalEsdClustersVector);
   
//       Int_t nec = event->GetEMCALClusters(fEmcalClustersArray);

      
      
//       if(fDoMatchedTracks)
// 	fEmcalMatchedTracksHistProducer->FillHistograms(nec, fEmcalClustersArray);
      
//       if(fDoInvariantMass)
// 	fEmcalInvariantMassHistProducer->FillHistograms(nec, fEmcalClustersArray);
      
//       if(fDoClusterEnergy)
// 	fEmcalClusterEnergyHistProducer->FillHistograms(nec, fEmcalClustersArray);
      
//       if(fDoCellEnergy)
// 	fEmcalCellEnergyHistProducer->FillHistograms(nec, fEmcalClustersArray);;
      
//     }
    
    
//     //PHOS
//     if(fDoPhos){
      
//       Int_t npc = event->GetPHOSClusters(fPhosClustersArray);
      
//       if(fDoMatchedTracks)
// 	fPhosMatchedTracksHistProducer->FillHistograms(npc, fPhosClustersArray);
      
//       if(fDoInvariantMass)
// 	fPhosInvariantMassHistProducer->FillHistograms(npc, fPhosClustersArray);
      
//       if(fDoClusterEnergy)
// 	fPhosClusterEnergyHistProducer->FillHistograms(npc, fPhosClustersArray);
      
//       if(fDoCellEnergy)
// 	fPhosCellEnergyHistProducer->FillHistograms(npc, fPhosClustersArray);
      
//     }
    
//   }
  
  
  //Push histos
  
  //PHOS
  if(fDoPhos){
    
    if(fDoInvariantMass)
      PushBack(fPhosInvariantMassHistProducer->GetHistograms(), kAliHLTDataTypeHistogram|kAliHLTDataOriginPHOS);
    
//     if(fDoMatchedTracks)
//       PushBack(fPhosMatchedTracksHistProducer->GetHistograms(), kAliHLTDataTypeHistogram, kAliHLTDataOriginPHOS);
    
//     if(fDoClusterEnergy) 
//       PushBack(fPhosClusterEnergyHistProducer->GetHistograms(), kAliHLTDataTypeHistogram, kAliHLTDataOriginPHOS);
    
//     if(fDoCellEnergy)
//       PushBack(fPhosCellEnergyHistProducer->GetHistograms(), kAliHLTDataTypeHistogram, kAliHLTDataOriginPHOS);
    
  }

  //EMCAL
  if(fDoEmcal) {

    if(fDoInvariantMass) 
      PushBack(fEmcalInvariantMassHistProducer->GetHistograms(), kAliHLTDataTypeHistogram);
        
//     if(fDoMatchedTracks) 
//       PushBack(fEmcalMatchedTracksHistProducer->GetHistograms(), kAliHLTDataOriginEMCALHLTDataTypeHistogram, kAliHLTDataOriginEMCAL);
  
//     if(fDoClusterEnergy) 
//       PushBack(fEmcalClusterEnergyHistProducer->GetHistograms(), kAliHLTDataOriginEMCALHLTDataTypeHistogram, kAliHLTDataOriginEMCAL);
    
//     if(fDoCellEnergy)
//       PushBack(fEmcalCellEnergyHistProducer->GetHistograms(), kAliHLTDataOriginEMCALHLTDataTypeHistogram, kAliHLTDataOriginEMCAL);
    
  }

  return 0;

}




