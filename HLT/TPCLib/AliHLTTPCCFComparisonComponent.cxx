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

/** @file   AliHLTTPCCFComparisonComponent.cxx
    @author Kalliopi Kanaki
    @date   
    @brief  A comparison component for FCF and SCF properties
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCCFComparisonComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTGlobalBarrelTrack.h"
#include "AliHLTExternalTrackParam.h"
#include "AliHLTDataTypes.h"

#include <TFile.h>
#include <TString.h>
#include "TNtuple.h"
#include "TH1F.h"
#include "TObjString.h"
#include "TObjArray.h"

#include "AliHLTTPCTrack.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCCFComparisonComponent)

AliHLTTPCCFComparisonComponent::AliHLTTPCCFComparisonComponent()
    :
    fEvtMod(20)
  , fBufferSize(5000)
  , fMultiplicity(NULL)
  , fClusters(NULL)
  , fTracks(NULL)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  for(int i=0; i<36; i++){
     for(int j=0; j<6; j++){
         fSCFClustersArray[i][j] = NULL;
         fSCFNSpacePoints[i][j]  = 0;     
         fFCFClustersArray[i][j] = NULL;
         fFCFNSpacePoints[i][j]  = 0;     
     }  
  }
}

const char* AliHLTTPCCFComparisonComponent::fgkOCDBEntry="HLT/ConfigTPC/TPCCFComparison";

AliHLTTPCCFComparisonComponent::~AliHLTTPCCFComparisonComponent(){
// see header file for class documentation
  for(int i=0; i<36; i++){
     for(int j=0; j<6; j++){
         delete[] fSCFClustersArray[i][j];
         fSCFClustersArray[i][j] = NULL;         
         delete[] fSCFClustersArray[i][j];
         fSCFClustersArray[i][j] = NULL;         
     }  
  }
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCCFComparisonComponent::GetComponentID(){
// see header file for class documentation
  
  return "TPCCFComparison";
}

void AliHLTTPCCFComparisonComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list){
// see header file for class documentation
  
  list.clear();
  list.push_back(AliHLTTPCDefinitions::fgkClustersDataType|kAliHLTDataOriginTPC);
  list.push_back(AliHLTTPCDefinitions::fgkAlterClustersDataType|kAliHLTDataOriginTPC);
  list.push_back(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC);  
}

AliHLTComponentDataType AliHLTTPCCFComparisonComponent::GetOutputDataType(){
// see header file for class documentation
  return kAliHLTMultipleDataType;
}

int AliHLTTPCCFComparisonComponent::GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList){
// see header file for class documentation
  tgtList.clear();
  tgtList.push_back(kAliHLTDataTypeTNtuple|kAliHLTDataOriginTPC);
  tgtList.push_back(kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC);
  return tgtList.size();
}

void AliHLTTPCCFComparisonComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ){
// see header file for class documentation
  
  constBase = 5000;
  inputMultiplier = 1;// XXX TODO: Find more realistic value
}

AliHLTComponent* AliHLTTPCCFComparisonComponent::Spawn(){
// see header file for class documentation
  return new AliHLTTPCCFComparisonComponent;
}

int AliHLTTPCCFComparisonComponent::DoInit( int argc, const char** argv ){
// see header file for class documentation
 
  fClusters = new TNtuple("fClusters", "fClusters", "charge:qmax:residualY:residualZ"); 
  fTracks   = new TNtuple("fTracks",  "fTracks",  "pt:eta:psi:nclusters"); 
 
  fClusters->SetCircular(fBufferSize);
  fTracks->SetCircular(fBufferSize);
 
  fMultiplicity = new TH1F("fMultiplicity","Track multiplicity per event", 1000, 0, 1000);
  
  
  // first configure the default
  int iResult=0;
  if (iResult>=0) iResult=ConfigureFromCDBTObjString(fgkOCDBEntry);

  // configure from the command line parameters if specified
  if (iResult>=0 && argc>0)  iResult=ConfigureFromArgumentString(argc, argv);
  
  return iResult;
}
  
int AliHLTTPCCFComparisonComponent::DoDeinit(){
// see header file for class documentation
  
  delete fClusters;
  delete fTracks;
  delete fMultiplicity;
  
  return 0;
}

int AliHLTTPCCFComparisonComponent::Reconfigure(const char* cdbEntry, const char* /*chainId*/){
// see header file for class documentation

  // configure from the specified antry or the default one
  const char* entry=cdbEntry;
  if (!entry || entry[0]==0) {
     entry=fgkOCDBEntry;
  }
  return ConfigureFromCDBTObjString(entry);
}


int AliHLTTPCCFComparisonComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/){
// see header file for class documentation

  if(GetFirstInputBlock(kAliHLTDataTypeSOR) || GetFirstInputBlock(kAliHLTDataTypeEOR)) return 0;  

  const AliHLTComponentBlockData *iter = NULL;
  
  for(int i=0; i<36; i++){
     for(int j=0; j<6; j++){
         fSCFClustersArray[i][j] = NULL;
         fSCFNSpacePoints[i][j]  = 0;     
         fFCFClustersArray[i][j] = NULL;
         fFCFNSpacePoints[i][j]  = 0;     
     }  
  }
 
 
  //----------------- loop over SCF output blocks ---------------------//
  
  //Int_t totalSpacePoints = 0;
  
  for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock()){
            
      if(iter->fDataType!=AliHLTTPCDefinitions::fgkClustersDataType) continue;

      AliHLTUInt8_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter);
      AliHLTUInt8_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);
      //HLTDebug("Input Data - TPC cluster - slice/partition: %d/%d.", minSlice, minPartition);

      const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*)iter->fPtr;
      Int_t nSpacepoint = (Int_t)clusterData->fSpacePointCnt;    
      //totalSpacePoints += nSpacepoint;
      HLTDebug("TPCCFComparison component found %d SCF spacepoints in slice %d partition %d", nSpacepoint, minSlice, minPartition);
      
      AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*)clusterData->fSpacePoints;
      
      if(fSCFClustersArray[minSlice][minPartition] != NULL){
         delete(fSCFClustersArray[minSlice][minPartition]); 
         fSCFClustersArray[minSlice][minPartition] = NULL;
      }      

      fSCFClustersArray[minSlice][minPartition] = clusters;
      fSCFNSpacePoints[minSlice][minPartition]  = nSpacepoint;
      
      HLTInfo("SCF cluster loop charge %d, slice %d partition %d\n", (UInt_t)clusters->fCharge, minSlice, minPartition);   
      
      // Here are all the cluster data you can use:
              
      //Float_t fX;	   // X coordinate in local coordinates
      //Float_t fY;	   // Y coordinate in local coordinates
      //Float_t fZ;	   // Z coordinate in local coordinates
      //UInt_t fID;	   // contains slice patch and number
      //UChar_t fPadRow;  // Pad row number
      //Float_t fSigmaY2; // error (former width) of the clusters
      //Float_t fSigmaZ2; // error (former width) of the clusters
      //UInt_t fCharge;   // total charge of cluster
      //UInt_t fQMax;     // QMax of cluster
  
      
      if(nSpacepoint==0) fSCFClustersArray[minSlice][minPartition] = NULL;

  } // end of loop over cluster data blocks
  
  //HLTDebug("TPCCFComparison found %d spacepoints",totalSpacePoints);
  
  
  //----------------- loop over FCF output blocks ---------------------//
  
  
  for(iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkAlterClustersDataType); iter != NULL; iter = GetNextInputBlock()){
            
      if(iter->fDataType!=AliHLTTPCDefinitions::fgkAlterClustersDataType) continue;

      AliHLTUInt8_t minSlice     = AliHLTTPCDefinitions::GetMinSliceNr(*iter);
      AliHLTUInt8_t minPartition = AliHLTTPCDefinitions::GetMinPatchNr(*iter);
      //HLTDebug("Input Data - TPC cluster - slice/partition: %d/%d.", minSlice, minPartition);

      const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*)iter->fPtr;
      Int_t nSpacepoint = (Int_t)clusterData->fSpacePointCnt;    
      //totalSpacePoints += nSpacepoint;
      HLTDebug("TPCCFComparison component found %d FCF spacepoints in slice %d partition %d", nSpacepoint, minSlice, minPartition);
      
      AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*)clusterData->fSpacePoints;
      
      if(fFCFClustersArray[minSlice][minPartition] != NULL){
         delete(fFCFClustersArray[minSlice][minPartition]); 
         fFCFClustersArray[minSlice][minPartition] = NULL;
      }      

      fFCFClustersArray[minSlice][minPartition] = clusters;
      fFCFNSpacePoints[minSlice][minPartition]  = nSpacepoint;
           
      HLTInfo("FCF cluster loop charge %d, slice %d partition %d\n", (UInt_t)clusters->fCharge, minSlice, minPartition);  
      
      if(nSpacepoint==0) fFCFClustersArray[minSlice][minPartition] = NULL;

  } // end of loop over cluster data blocks
  
  
  
  //----------------- loop over merged tracks -------------------//

//   Int_t totalTracks = 0;
//   
//   for(iter = GetFirstInputBlock(kAliHLTDataTypeTrack|kAliHLTDataOriginTPC); iter != NULL; iter = GetNextInputBlock()){    
//       if(iter->fDataType != (kAliHLTDataTypeTrack|kAliHLTDataOriginTPC)) continue; 
//       ReadTracks(iter,totalTracks);
//   }
//   
//   HLTDebug("TrackHisto found %d tracks", totalTracks);  
// 
//   fMultiplicity->Fill(totalTracks);
// 
//   PushHisto();
  
  return 0;
} // end of DoEvent
 
// void AliHLTTPCCFComparisonComponent::ReadTracks(const AliHLTComponentBlockData* iter,Int_t &tt){
// // see header file for class documentation
// 
//   //AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr(*iter);
//   
//   Int_t usedSpacePoints = 0;
//   
//   vector<AliHLTGlobalBarrelTrack> tracksVector;
//   AliHLTGlobalBarrelTrack::ConvertTrackDataArray(reinterpret_cast<const AliHLTTracksData*>(iter->fPtr), iter->fSize, tracksVector);
//   
//   tt = tracksVector.size();
//   
//   for(vector<AliHLTGlobalBarrelTrack>::iterator element=tracksVector.begin();  element!=tracksVector.end(); element++){
//        
//   
// 
//            
//        UInt_t nHits = element->GetNumberOfPoints();
//        fTracks->Fill( 1./element->OneOverPt(), element->GetSnp(), element->GetTgl(), nHits );  
//        //fdNdEta->Fill(element->GetSnp());
//  
//        Double_t totCharge = 0;
//        const UInt_t *hitnum = element->GetPoints();
//        for(UInt_t i=0; i<element->GetNumberOfPoints(); i++){
//            
//       	   UInt_t idTrack   = hitnum[i];
//            Int_t sliceTrack = (idTrack>>25) & 0x7f;
//            Int_t patchTrack = (idTrack>>22) & 0x7;
//            UInt_t pos	    = idTrack&0x3fffff;
// 	   
// 	   //printf("KKKK pos :%d\n", pos);
// 	   cout << "KKKK pos  " << pos << endl;
// 
//            if( !fClustersArray[sliceTrack][patchTrack] ) continue;	       
//       	   
//       	   if(sliceTrack<0 || sliceTrack>36 || patchTrack<0 || patchTrack>5 ){
//       	      HLTError("Corrupted TPC cluster Id: slice %d, patch %d, cluster %d", sliceTrack, patchTrack, idTrack);
//       	      continue;
//       	   }
// 
//            if(fNSpacePoints[sliceTrack][patchTrack]<=pos ){
//       	      HLTError("Space point array out of boundaries!");
//       	      continue;
//            }
//       	   
//       	   totCharge += (fClustersArray[sliceTrack][patchTrack])[pos].fCharge; 
//       	   
//       	   //Float_t xyz[3]; xyz[0] = xyz[1] = xyz[2] = 0.;
//         
//            //xyz[0] = (fClustersArray[sliceTrack][patchTrack])[pos].fX;
//            //xyz[1] = (fClustersArray[sliceTrack][patchTrack])[pos].fY;
//            //xyz[2] = (fClustersArray[sliceTrack][patchTrack])[pos].fZ;
//         
//            //AliHLTTPCTransform::Local2Global(xyz,slice); 
//       	   
//       	   //Double_t p[2]   = { xyz[1], xyz[2] };
//            //Double_t cov[3] = { (fClustersArray[sliceTrack][patchTrack])[pos].fSigmaY2, 0., (fClustersArray[sliceTrack][patchTrack])[pos].fSigmaZ2};  
//       	   //Double_t *res = element->GetResiduals(p,cov,kFALSE); 
//       	   
//       	   //HLTInfo("resy: %f, resz: %f", res[0], res[1]);
//       	   
//       	   //if(!res)  res[0] = res[1] = -1000.;
//       	   //else      fClusters->Fill( (fClustersArray[sliceTrack][patchTrack])[pos].fCharge, (fClustersArray[sliceTrack][patchTrack])[pos].fQMax, res[0], res[1]);
//       	  
//       	   fClusters->Fill( (fClustersArray[sliceTrack][patchTrack])[pos].fCharge, (fClustersArray[sliceTrack][patchTrack])[pos].fQMax, -1000., -1000.);
//             
//       	   usedSpacePoints++;	   
//        }	
// 
//   }
// }

// void AliHLTTPCCFComparisonComponent::PushHisto(){
// // see header file for class documentation
//     
//     PushBack( (TObject*)fTracks,       kAliHLTDataTypeTNtuple  |kAliHLTDataOriginTPC, 0x0);   
//     PushBack( (TObject*)fClusters,     kAliHLTDataTypeTNtuple  |kAliHLTDataOriginTPC, 0x0);
//     PushBack( (TObject*)fMultiplicity, kAliHLTDataTypeHistogram|kAliHLTDataOriginTPC, 0x0);
// }

int AliHLTTPCCFComparisonComponent::ScanConfigurationArgument(int argc, const char** argv){
// see header file for class documentation
 
  if (argc<=0) return 0;
  int i=0;
  TString argument=argv[i];

  // -event-modulo
  if (argument.CompareTo("-event-modulo")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fEvtMod=argument.Atof();
    return 2;
  }    

  // -buffer-size
  if (argument.CompareTo("-buffer-size")==0) {
    if (++i>=argc) return -EPROTO;
    argument=argv[i];
    fBufferSize=argument.Atof();
    return 2;
  }    
   
  return -EINVAL;
}

void AliHLTTPCCFComparisonComponent::GetOCDBObjectDescription( TMap* const targetMap){
// Get a list of OCDB object description needed for the particular component
  if (!targetMap) return;
  targetMap->Add(new TObjString("HLT/ConfigTPC/TPCCFComparison"), new TObjString("component arguments for setting the size of the filled ntuples and the event modulo for the mean multiplicity distribution"));
}
