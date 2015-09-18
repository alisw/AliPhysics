// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
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

/** @file   AliHLTTPCCalibSeedMakerComponent.cxx
    @author Kalliopi Kanaki
    @date   2009-07-08
    @brief  
*/

#include <map>

#include "AliHLTTPCCalibSeedMakerComponent.h"
#include "AliHLTTPCGeometry.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCOfflineCluster.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTTrackMCLabel.h"

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliTPCcalibDB.h"
#include "AliTPCParam.h"

#include "AliRieman.h"

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliCDBStorage.h"

#include <cstdlib>
#include <cerrno>

#include "TObjArray.h"
#include "TClonesArray.h"
#include "TObject.h"
#include "TFile.h"
#include "TH2F.h"

#include <sys/time.h>

using namespace std;

ClassImp(AliHLTTPCCalibSeedMakerComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCCalibSeedMakerComponent::AliHLTTPCCalibSeedMakerComponent()
    :    
    fTPCGeomParam(0)
   ,fSeedArray(0x0)
   ,fdEdx(0x0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt  
  for( int i=0; i<fkNPartition; i++ ){
       fPartitionClusters[i]  = 0;    
       fNPartitionClusters[i] = 0;    
  }
}

AliHLTTPCCalibSeedMakerComponent::~AliHLTTPCCalibSeedMakerComponent() { 
// see header file for class documentation  
  
  for( int i=0; i<fkNPartition; i++ ){
       delete[] fPartitionClusters[i];
  }
}

const char* AliHLTTPCCalibSeedMakerComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCCalibSeedMaker";
}

void AliHLTTPCCalibSeedMakerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
// see header file for class documentation

  list.clear(); 
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType   );
  list.push_back( kAliHLTDataTypeTrack|kAliHLTDataOriginTPC   );
  list.push_back( kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC );
}

AliHLTComponentDataType AliHLTTPCCalibSeedMakerComponent::GetOutputDataType(){ 
// see header file for class documentation

  return kAliHLTMultipleDataType;
  //return kAliHLTDataTypeTObjArray;
}

int AliHLTTPCCalibSeedMakerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList){ 
// see header file for class documentation

  tgtList.clear();
  tgtList.push_back( kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC );
  tgtList.push_back( kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC );
  return tgtList.size();
}

void AliHLTTPCCalibSeedMakerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
// see header file for class documentation

  constBase=2000;
  inputMultiplier=2.0; // to be estimated
}

AliHLTComponent* AliHLTTPCCalibSeedMakerComponent::Spawn() { 
// see header file for class documentation

  return new AliHLTTPCCalibSeedMakerComponent();
}
	
int AliHLTTPCCalibSeedMakerComponent::DoInit( int /*argc*/, const char** /*argv*/ ) { 
// see header file for class documentation
 
  fTPCGeomParam = AliTPCcalibDB::Instance()->GetParameters();
  if(!fTPCGeomParam) HLTError("TPC Parameters are not loaded.");
  
  fSeedArray = new TObjArray(10000);  
  //fSeedArray->SetOwner(kTRUE);
  
  fdEdx = new TH2F("fdEdx","energy loss vs. momentum", 400, -200, 200, 300, 0, 300);
  
  return 0;

} // end DoInit()

int AliHLTTPCCalibSeedMakerComponent::DoDeinit() { 
// see header file for class documentation  
  
  if(fTPCGeomParam) delete fTPCGeomParam; fTPCGeomParam = NULL;	  
  if(fSeedArray)    delete fSeedArray;    fSeedArray	= NULL;	
  if(fdEdx)         delete fdEdx;         fdEdx         = NULL;

  return 0;
}

int AliHLTTPCCalibSeedMakerComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/){
// see header file for class documentation
 
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR)) return 0;

  int nInputClusters = 0;
  fSeedArray->Clear();
  
  for(Int_t i=0; i<fkNPartition; i++){
      delete[] fPartitionClusters[i];    
      fPartitionClusters[i]  = 0;
      fNPartitionClusters[i] = 0;    
  }


  // ---------- Access to clusters --------------------//
  
  for(const AliHLTComponentBlockData *iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock()){
      
      if(iter->fDataType != AliHLTTPCDefinitions::fgkClustersDataType) continue;
    
      Int_t slice     = AliHLTTPCDefinitions::GetMinSliceNr(iter->fSpecification);
      Int_t partition = AliHLTTPCDefinitions::GetMinPatchNr(iter->fSpecification);
    
      Int_t slicepartition = slice*6+partition;
      
      if(slicepartition > fkNPartition){
         HLTWarning("Wrong header of TPC cluster data, slice %d, partition %d", slice, partition );
         continue;
      }
      
      AliHLTTPCClusterData *inPtrSP = ( AliHLTTPCClusterData* )( iter->fPtr );
      nInputClusters += inPtrSP->fSpacePointCnt;

      delete[] fPartitionClusters[slicepartition];
      fPartitionClusters[slicepartition]  = new AliTPCclusterMI[inPtrSP->fSpacePointCnt];
      fNPartitionClusters[slicepartition] = inPtrSP->fSpacePointCnt;
    
      // create  offline clusters out of the HLT clusters
      // todo: check which cluster information is really needed for the dEdx

    for ( unsigned int i = 0; i < inPtrSP->fSpacePointCnt; i++ ) {
      AliHLTTPCSpacePointData *chlt = &( inPtrSP->fSpacePoints[i] );
      AliTPCclusterMI *c = fPartitionClusters[slicepartition]+i;
      c->SetX(chlt->fX);
      c->SetY(chlt->fY);
      c->SetZ(chlt->fZ);
      c->SetSigmaY2(chlt->fSigmaY2);
      c->SetSigmaYZ( 0 );
      c->SetSigmaZ2(chlt->fSigmaZ2);
      c->SetQ( chlt->fCharge );
      c->SetMax( chlt->fQMax );
      Int_t sector, row;
      Float_t padtime[3] = {0,chlt->fY,chlt->fZ};
      AliHLTTPCGeometry::Slice2Sector(slice,chlt->fPadRow, sector, row);
      AliHLTTPCGeometry::Local2Raw( padtime, sector, row);
      c->SetDetector( sector );
      c->SetRow( row );
      c->SetPad( (Int_t) padtime[1] );
      c->SetTimeBin( (Int_t) padtime[2] );
    }      
  } // end of loop over blocks of clusters

 
 
  
  // ------------ loop over the MC labels -----------------//
  
  std::map<int,int> mcLabels;
  for(const AliHLTComponentBlockData* pBlock=GetFirstInputBlock(kAliHLTDataTypeTrackMC|kAliHLTDataOriginTPC); pBlock!=NULL; pBlock=GetNextInputBlock()){
   
      AliHLTTrackMCData *dataPtr = reinterpret_cast<AliHLTTrackMCData*>( pBlock->fPtr );
    
      if(sizeof(AliHLTTrackMCData)+dataPtr->fCount*sizeof(AliHLTTrackMCLabel)==pBlock->fSize){     
        for(UInt_t il=0; il<dataPtr->fCount; il++){  
      	    AliHLTTrackMCLabel &lab = dataPtr->fLabels[il];
      	    mcLabels[lab.fTrackID]  = lab.fMCLabel;
      	    HLTDebug("MC labels, track ID:  %d, %d\n", lab.fMCLabel, lab.fTrackID);
        }
      } 
      else {
        HLTWarning("data mismatch in block %s (0x%08x): count %d, size %d -> ignoring track MC information", DataType2Text(pBlock->fDataType).c_str(), pBlock->fSpecification, 
		   dataPtr->fCount, pBlock->fSize);
      }
  } // end of loop over MC label blocks
  
  
  
 
 
  //------------------ loop over track data blocks --------------------//

  int nTracks = 0;
  for(const AliHLTComponentBlockData *pBlock = GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC); pBlock != NULL; pBlock = GetNextInputBlock()){
      
      AliHLTTracksData *dataPtr = (AliHLTTracksData*) pBlock->fPtr;
      //int nTracks = dataPtr->fCount;
      nTracks = dataPtr->fCount;
    
      AliHLTExternalTrackParam *currTrack = dataPtr->fTracklets;     
      
      for(Int_t itr=0; itr<nTracks && ( (AliHLTUInt8_t *)currTrack < ((AliHLTUInt8_t *) pBlock->fPtr)+pBlock->fSize); itr++){
    
         // create an offline track
         AliHLTGlobalBarrelTrack gb(*currTrack);
         AliTPCseed tTPC;
         tTPC.Set( gb.GetX(), gb.GetAlpha(), gb.GetParameter(), gb.GetCovariance() );
	            
         Int_t mcLabel = -1;
         if( mcLabels.find(gb.TrackID())!=mcLabels.end() ) mcLabel = mcLabels[gb.TrackID()];      
         tTPC.SetLabel(mcLabel);
            
         // set the clusters 
	     
         for(UInt_t ic=0; ic<currTrack->fNPoints; ic++){	
            
             tTPC.SetNumberOfClusters(currTrack->fNPoints);
          
	     UInt_t id      = currTrack->fPointIDs[ic];	     
	     int iSlice = AliHLTTPCSpacePointData::GetSlice(id);
	     int iPartition = AliHLTTPCSpacePointData::GetPatch(id);
	     int iCluster = AliHLTTPCSpacePointData::GetNumber(id);
	
	     if(iSlice<0 || iSlice>36 || iPartition<0 || iPartition>5){
	         HLTError("Corrupted TPC cluster Id: slice %d, partition %d, cluster %d", iSlice, iPartition, iCluster);
	         continue;
	     }
	
	     AliTPCclusterMI *patchClusters = fPartitionClusters[iSlice*6 + iPartition];
	     if(!patchClusters){
	        HLTError("Clusters are missed for slice %d, partition %d", iSlice, iPartition );
	        continue;
	     }
	
	     if(iCluster >= fNPartitionClusters[iSlice*6 + iPartition]){
	        HLTError("TPC slice %d, partition %d: ClusterID==%d >= N Cluaters==%d ", iSlice, iPartition,iCluster, fNPartitionClusters[iSlice*6 + iPartition] );
	        continue;
	     }
	
	     AliTPCclusterMI *c = &(patchClusters[iCluster]);	  	
	     int sec = c->GetDetector();
	     int row = c->GetRow();
	     if(sec >= 36) row = row + AliHLTTPCGeometry::GetNRowLow();
	
	     tTPC.SetClusterPointer(row, c);	
	
	     AliTPCTrackerPoint &point = *( tTPC.GetTrackPoint( row ) );
	     //tTPC.Propagate( TMath::DegToRad()*(sec%18*20.+10.), c->GetX(), fSolenoidBz );
	     Double_t angle2 = tTPC.GetSnp()*tTPC.GetSnp();
	     angle2 = (angle2<1) ?TMath::Sqrt(angle2/(1-angle2)) :10.; 
	     point.SetAngleY( angle2 );
	     point.SetAngleZ( tTPC.GetTgl() );
         } // end of associated cluster loop

      // Cook dEdx
           
      AliTPCseed *seed = &(tTPC);      
      fSeedArray->AddAt( seed, TMath::Abs(seed->GetLabel()) );
      fdEdx->Fill( seed->P()*seed->Charge(), seed->CookdEdx(0.02, 0.6) );
          
      unsigned int step = sizeof( AliHLTExternalTrackParam ) + currTrack->fNPoints * sizeof( unsigned int );
      currTrack = ( AliHLTExternalTrackParam* )( (( Byte_t * )currTrack) + step );  

      }// end of vector track loop           
  } // end of loop over blocks of merged tracks 
  
  HLTDebug("Number of reconstructed tracks %d, number of produced seeds %d\n", nTracks, fSeedArray->GetEntries());  
  
  PushBack((TObject*)fSeedArray, kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC, 0x0);
  PushBack((TObject*)fdEdx, kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC, 0x0);
         
  return 0;
} // end DoEvent()
