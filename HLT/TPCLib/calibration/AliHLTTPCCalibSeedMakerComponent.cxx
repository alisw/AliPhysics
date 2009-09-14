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

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCCalibSeedMakerComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCOfflineCluster.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTGlobalBarrelTrack.h"

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
#include "TObject.h"
#include <sys/time.h>

ClassImp(AliHLTTPCCalibSeedMakerComponent) //ROOT macro for the implementation of ROOT specific class methods

AliHLTTPCCalibSeedMakerComponent::AliHLTTPCCalibSeedMakerComponent()
    :    
    fSpecification(0),
    fMinSlice(35),
    fMaxSlice(0),
    fMinPartition(5),
    fMaxPartition(0),
    fTPCGeomParam(0),
    fSeedArray(0)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCCalibSeedMakerComponent::~AliHLTTPCCalibSeedMakerComponent() { 
// see header file for class documentation

}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCCalibSeedMakerComponent::GetComponentID() { 
// see header file for class documentation

  return "TPCCalibSeedMaker";
}

void AliHLTTPCCalibSeedMakerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list) { 
// see header file for class documentation

  list.clear(); 
  list.push_back(AliHLTTPCDefinitions::fgkClustersDataType);
  //list.push_back(AliHLTTPCDefinitions::fgkTrackSegmentsDataType);
  list.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);

}

AliHLTComponentDataType AliHLTTPCCalibSeedMakerComponent::GetOutputDataType(){ 
// see header file for class documentation

  //return kAliHLTMultipleDataType;
  return kAliHLTDataTypeTObjArray;
}

int AliHLTTPCCalibSeedMakerComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList){ 
// see header file for class documentation

  tgtList.clear();
  //tgtList.push_back(AliHLTTPCDefinitions::fgkAliHLTDataTypeTPCSeed);
  tgtList.push_back(kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC);
  return tgtList.size();
}

void AliHLTTPCCalibSeedMakerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) { 
// see header file for class documentation

  constBase=0;
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
  
  fSeedArray = new TObjArray;
  fSeedArray->SetOwner(kTRUE);
  return 0;

} // end DoInit()

int AliHLTTPCCalibSeedMakerComponent::DoDeinit() { 
// see header file for class documentation  
  
  if(fTPCGeomParam) delete fTPCGeomParam; fTPCGeomParam = NULL;          
  if(fSeedArray)   {fSeedArray->Clear();  delete fSeedArray; }   fSeedArray    = NULL;          
  return 0;
}

int AliHLTTPCCalibSeedMakerComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/){
// see header file for class documentation
 
  const AliHLTComponentBlockData *iter = NULL;      
  if(GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR)) return 0;

  Int_t totalSpacePoints = 0;
  Int_t usedSpacePoints  = 0;   
  
 
  // ---------- Access to clusters --------------------//
  
  for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock()){
  
      if(iter->fDataType!=AliHLTTPCDefinitions::fgkClustersDataType) continue;

      AliHLTUInt8_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter);
      AliHLTUInt8_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);
      
      const AliHLTTPCClusterData *clusterData = (const AliHLTTPCClusterData*)iter->fPtr;
      Int_t nSpacepoint = (Int_t)clusterData->fSpacePointCnt;	
      
      totalSpacePoints += nSpacepoint;
      
      AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*)clusterData->fSpacePoints; 
      
      if(fClustersArray[minSlice][minPartition] != NULL){
         //delete(fClustersArray[minSlice][minPartition]);
         fClustersArray[minSlice][minPartition] = NULL;
      }
      
      // fill the array with AliHLTTPCSpacePointData pointers
      // it will be used in the track loop to access information
      // for the used clusters only
      fClustersArray[minSlice][minPartition] = clusters;
      fNSpacePoints[minSlice][minPartition]  = nSpacepoint;
      
      if(nSpacepoint==0) fClustersArray[minSlice][minPartition] = NULL;
      
  } // end of loop over blocks of clusters

  HLTDebug("Total space points: %d", totalSpacePoints);
   
  
  //------------------ Access to track data blocks --------------------//
  
  TObjArray *offClusterArray = new TObjArray;
  //offClusterArray->SetOwner(kTRUE);
  offClusterArray->Clear();
  
  fSeedArray->Clear();
  
  for(iter = GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC); iter != NULL; iter = GetNextInputBlock()){ 
  
      if(iter->fDataType != (kAliHLTDataTypeTrack|kAliHLTDataOriginTPC)) continue; 
     
      AliHLTUInt8_t slice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter);
      AliHLTUInt8_t partition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);
      
      if( slice     < fMinSlice )     fMinSlice     = slice;
      if( slice     > fMaxSlice )     fMaxSlice     = slice;
      if( partition < fMinPartition ) fMinPartition = partition;
      if( partition > fMaxPartition ) fMaxPartition = partition;
      
      
      // create a vector of AliHLTGlobalBarrelTrack tracks from AliHLTTracksData
      vector<AliHLTGlobalBarrelTrack> tracks;
      AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(iter->fPtr), iter->fSize, tracks);      
    
      // loop over the elements(tracks) of the AliHLTGlobalBarrelTrack vector

      for(vector<AliHLTGlobalBarrelTrack>::iterator element=tracks.begin();  element!=tracks.end(); element++){
          
	  AliRieman rieman(element->GetNumberOfPoints());
	  rieman.Reset();    
          AliTPCseed *seed = 0x0;
          Double_t param[5]; for(Int_t i=0; i<5;  i++) param[i] = 0.;
          Double_t cov[15];  for(Int_t i=0; i<15; i++) cov[i]   = 0.;
          Double_t xmin  = 1000.;
          Double_t alpha = 0.;

    	  // calculate padrow radius (x - radius)
    	  Double_t xrow[160]; for(Int_t i=0; i<160; i++) xrow[i] = 0.;
    	  Int_t nrowlow = fTPCGeomParam->GetNRowLow();
    	  Int_t nrowup  = fTPCGeomParam->GetNRowUp();
	
    	  for (Int_t i=0;i<nrowlow;i++) xrow[i]         = fTPCGeomParam->GetPadRowRadiiLow(i);
    	  for (Int_t i=0;i<nrowup;i++)  xrow[i+nrowlow] = fTPCGeomParam->GetPadRowRadiiUp(i);
	  
    	  // sector rotation angles - only one angle is needed
    	  Double_t angle = fTPCGeomParam->GetInnerAngle();
    	    
	  const UInt_t *hitnum = element->GetPoints(); // store the clusters on each track in an array and loop over them
          for(UInt_t i=0; i<element->GetNumberOfPoints(); i++){
              
	      // the id of the cluster contains information about the slice and partition it belongs too
	      // as well as its index (pos)          
              UInt_t idTrack   = hitnum[i];
              Int_t sliceTrack = (idTrack>>25) & 0x7f;
              Int_t patchTrack = (idTrack>>22) & 0x7;
              UInt_t pos       = idTrack&0x3fffff;
	      
	      // use the fClustersArray that was filled in the cluster loop
	      if( !fClustersArray[sliceTrack][patchTrack]    ) continue;
    	      if(  fNSpacePoints[sliceTrack][patchTrack]<pos ) HLTError("Space point array out of boundaries!");
	      
	      // get the sector(detector) information
	      Int_t sector, row = -99;
	      AliHLTTPCTransform::Slice2Sector(sliceTrack, (fClustersArray[sliceTrack][patchTrack])[pos].fPadRow, sector, row);
	               
	      // next line recalculates rows in the sector(ROC) system
 	      //if(patchTrack>1) (fClustersArray[sliceTrack][patchTrack])[pos].fPadRow -= 63;//(Int_t)AliHLTTPCTransform::GetFirstRow(2);
	      
	      HLTDebug("slice %d, partition :%d, sector row: %d", sliceTrack, patchTrack, (fClustersArray[sliceTrack][patchTrack])[pos].fPadRow);
	       
	      // convert the HTL clusters to AliTPCclusterMI    	      
    	      AliHLTTPCOfflineCluster pConv;
    	      AliTPCclusterMI *offClus = pConv.ConvertHLTToOffline((fClustersArray[sliceTrack][patchTrack])[pos]);	      
    	      offClus->SetDetector(sector);
	      usedSpacePoints++;
	      
	      offClusterArray->Add(offClus);
	      	            
	      rieman.AddPoint( offClus->GetX(),offClus->GetY(),offClus->GetZ(),TMath::Sqrt(offClus->GetSigmaY2()),TMath::Sqrt(offClus->GetSigmaZ2()) );    
              
	      xmin = TMath::Min(xmin,xrow[offClus->GetRow()]); // min pad-row radius
	      	    
              alpha = 0.5*angle+angle*(sector%18); //sector rotation angle
              
              // HLTInfo("detector: %d, row: %d, xrow[row]: %f", sector, offClus->GetRow(), xrow[offClus->GetRow()]);  
	     
          } // end of cluster loop

          // creation of AliTPCseed by means of a Riemann fit
          rieman.Update();
          rieman.GetExternalParameters(xmin,param,cov);  
	  seed = new AliTPCseed(xmin,alpha,param,cov,0);
	 
	  // set up of the cluster pointers inside the seed
 	  for(Int_t j=0; j<offClusterArray->GetEntries(); j++){ 
              AliTPCclusterMI *cl = (AliTPCclusterMI*)offClusterArray->At(j);
              if(cl) seed->SetClusterPointer(cl->GetRow(),cl);
          }

	  offClusterArray->Clear();	
          HLTDebug("External track parameters: seed: 0x%08x, xmin: %f, alpha: %f, param[0]: %f, cov[0]: %f", seed, xmin, alpha, param[0], cov[0]);
	  fSeedArray->Add(seed);	 

      }// end of vector track loop           
  } // end of loop over blocks of merged tracks  
  
  HLTDebug("Used space points: %d", usedSpacePoints);
  HLTDebug("Number of entries in fSeedArray: %d", fSeedArray->GetEntries());
  
  if(offClusterArray) delete offClusterArray;
  offClusterArray = NULL;
 
  fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification( fMinSlice, fMaxSlice, fMinPartition, fMaxPartition );
  //PushBack((TObject*)offClusterArray, kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC, fSpecification);
  PushBack((TObject*)fSeedArray,       kAliHLTDataTypeTObjArray|kAliHLTDataOriginTPC, fSpecification);
 
  return 0;
} // end DoEvent()
