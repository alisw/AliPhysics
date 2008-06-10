// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Gaute Ovrebekk <ovrebekk@ift.uib.no>                  *
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

/** @file   AliHLTTPCTrackHistoComponent.cxx
    @author Gaute Ovrebekk, Matthias Richter
    @date   
    @brief  The TPC conformal mapping tracker component.
*/


#if __GNUC__>= 3
using namespace std;
#endif

#include "AliHLTTPCTrackHistoComponent.h"
#include "AliHLTTPCTransform.h"
#include "AliHLTTPCClusterDataFormat.h"
#include "AliHLTTPCTrackletDataFormat.h"
#include "AliHLTTPCMemHandler.h"
#include "AliHLTTPCDefinitions.h"
#include <TFile.h>
#include <TString.h>
#include "TH1F.h"
#include "TObjString.h"
#include "TObjArray.h"

//#include "AliHLTTPC.h"
//#include <stdlib.h>
//#include <cerrno>

// this is a global object used for automatic component registration, do not use this
AliHLTTPCTrackHistoComponent gAliHLTTPCTrackHistoComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCTrackHistoComponent)

AliHLTTPCTrackHistoComponent::AliHLTTPCTrackHistoComponent()
:
fHistoNClustersOnTracks(NULL),                                  
  fHistoAllClusters(NULL),                                                                         
  fHistoUsedClusters(NULL),                                                                        
  fHistoPT(NULL),                                                                                  
  fHistoResidual(NULL),                                                                            
  fHistoTgl(NULL),                                                                                 
  fPlotAll(kFALSE),                                                
  fPlotNClustersOnTracks(kFALSE),                                                
  fPlotChargeClusters(kFALSE),                                                                                     
  fPlotChargeUsedClusters(kFALSE),                                                                                 
  fPlotPT(kFALSE),                                                                                                 
  fPlotResidual(kFALSE),                                                                                           
  fPlotTgl(kFALSE),
  fClusters(),
  fTracks()
{
  
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTTPCTrackHistoComponent::~AliHLTTPCTrackHistoComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTPCTrackHistoComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "TPCTrackHisto";
}

void AliHLTTPCTrackHistoComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back( AliHLTTPCDefinitions::fgkClustersDataType );
  list.push_back( AliHLTTPCDefinitions::fgkTrackSegmentsDataType );
  list.push_back( AliHLTTPCDefinitions::fgkTracksDataType );
}

AliHLTComponentDataType AliHLTTPCTrackHistoComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram;

}

void AliHLTTPCTrackHistoComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  // XXX TODO: Find more realistic values.
  constBase = 0;
  inputMultiplier = 1;
}

AliHLTComponent* AliHLTTPCTrackHistoComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTPCTrackHistoComponent;
}

int AliHLTTPCTrackHistoComponent::DoInit( int argc, const char** argv )
{
  fHistoNClustersOnTracks = new TH1F("fHistoNClustersOnTracks","Number of Clusters on Tracks",160,0,160);                                   
  fHistoAllClusters = new TH1F("fHistoAllClusters","Charge of All Clusters",4000,0,4000);                                                                    
  fHistoUsedClusters = new TH1F("fHistoUsedClusters","Charge of Clusters used on Tracks",4000,0,4000);                                                              
  fHistoPT = new TH1F("fHistoPT","pT of Tracks",100,0,10);                                                                    
  fHistoResidual = new TH1F("fHistoResidual","Resuduals",360,0,360);    //change. Testing                                                                       
  fHistoTgl = new TH1F("fHistoTgl","Tgl of Tracks",900,0,90);  

  fPlotAll=kFALSE;                                                
  fPlotNClustersOnTracks=kFALSE;                                                
  fPlotChargeClusters=kFALSE;                                                                                     
  fPlotChargeUsedClusters=kFALSE;                                                                                 
  fPlotPT=kFALSE;                                                                                                 
  fPlotResidual=kFALSE;                                                                                           
  fPlotTgl=kFALSE;                 

  int iResult=0;
  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }
  
  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  }  
  return iResult; 
}
  
int AliHLTTPCTrackHistoComponent::DoDeinit()
{
  // see header file for class documentation
  
  delete fHistoNClustersOnTracks;                                 
  delete fHistoAllClusters;                                                                         
  delete fHistoUsedClusters;                                                                        
  delete fHistoPT;                                                                                
  delete fHistoResidual;                                                                            
  delete fHistoTgl;        
  
  return 0;
}

int AliHLTTPCTrackHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  const AliHLTComponentBlockData* iter = NULL;

  if ( GetFirstInputBlock( kAliHLTDataTypeSOR ) || GetFirstInputBlock( kAliHLTDataTypeEOR ) )
    return 0;

  Int_t TotalTrack = 0;

  //Reading Merged Tracks
  for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTracksDataType); iter != NULL; iter = GetNextInputBlock() ) {
    if(iter->fDataType!=AliHLTTPCDefinitions::fgkTracksDataType){continue;}
    ReadTracks(iter,TotalTrack);  
  }

  //Reading Tracks form slice
  for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkTrackSegmentsDataType); iter != NULL; iter = GetNextInputBlock() ) {
    if(iter->fDataType!=AliHLTTPCDefinitions::fgkTrackSegmentsDataType){continue;}
    ReadTracks(iter,TotalTrack);
  }

  int TotalSpacePoint = 0;
  int nClustersUsed=0;
 
  for ( iter = GetFirstInputBlock(AliHLTTPCDefinitions::fgkClustersDataType); iter != NULL; iter = GetNextInputBlock() ) {
    
    if(iter->fDataType!=AliHLTTPCDefinitions::fgkClustersDataType){continue;}

    //AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
    //AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );

    //HLTDebug ( "Input Data - TPC cluster - Slice/Patch: %d/%d.", slice, patch );
    const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*) iter->fPtr;
    Int_t nSpacepoint = (Int_t) clusterData->fSpacePointCnt;    
    TotalSpacePoint += nSpacepoint;
    //HLTInfo("TrackHisto found %d Spacepoints in slice %d patch %d", nSpacepoint, slice, patch);
    AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*) clusterData->fSpacePoints;

    for(int i=0;i<nSpacepoint;i++){
      UInt_t idCluster = clusters[i].fID;
      Int_t sliceCl = (idCluster>>25) & 0x7f;
      Int_t patchCl = (idCluster>>22) & 0x7;
      UInt_t pos = idCluster&0x3fffff;
      if(fPlotChargeClusters || fPlotAll){fHistoAllClusters->Fill(clusters[i].fCharge);}
      for(UInt_t id=0;id<fTrackClusterID[sliceCl][patchCl].size();id++){
	if(fTrackClusterID[sliceCl][patchCl][id]==pos){
	  clusters[i].fUsed=kTRUE;
	  nClustersUsed++;
	  if(fPlotChargeUsedClusters || fPlotAll){fHistoUsedClusters->Fill(clusters[i].fCharge);}
	}
      } 
      fClusters.push_back(clusters[i]);
    }
  } 

  HLTInfo("TrackHisto found %d Spacepoints",TotalSpacePoint);
  HLTInfo("TrackHisto found %d Tracks",TotalTrack);
  
  PushHisto();

  fClusters.clear();
  fTracks.clear();

  return 0;
}
 
 int AliHLTTPCTrackHistoComponent::Configure(const char* arguments)
 {
   
   int iResult=0;
   if (!arguments) return iResult;
   
   TString allArgs=arguments;
   TString argument;
   
   TObjArray* pTokens=allArgs.Tokenize(" ");

   if (pTokens) {
     for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
       argument=((TObjString*)pTokens->At(i))->GetString();
       if (argument.IsNull()) continue;

       if (argument.CompareTo("-plot-All")==0) {
	HLTInfo("Ploting All Histograms for Tracks");
	fPlotAll = kTRUE;
	fPlotNClustersOnTracks=kTRUE;                                                
	fPlotChargeClusters=kTRUE;                                                                                     
	fPlotChargeUsedClusters=kTRUE;                                                                                 
	fPlotPT=kTRUE;                                                                                                 
	fPlotResidual=kTRUE;                                                                                           
	fPlotTgl=kTRUE;            
	continue;
       }
       else if (argument.CompareTo("-plot-nClusters")==0) {
	 HLTInfo("Ploting Number of clusters Used on Tracks");
	 fPlotNClustersOnTracks = kTRUE;
	 continue;
       }
       else if (argument.CompareTo("-plot-ChargeClusters")==0) {
	 HLTInfo("Ploting Charge of All Clusters");
	 fPlotChargeClusters = kTRUE;
	 continue;
       }
       else if (argument.CompareTo("-plot-ChargeUsedClusters")==0) {
	 HLTInfo("Ploting Charge of Clusters Used on Tracks");
	 fPlotChargeUsedClusters = kTRUE; 
	 continue;
       }
       else if (argument.CompareTo("-plot-pT")==0) {
	 HLTInfo("Ploting pT of Tracks");
	 fPlotPT=kTRUE;
	 continue;
       }
       else if (argument.CompareTo("-plot-Residuals")==0) {
	 HLTInfo("Ploting Residuals");
	 fPlotResidual=kTRUE; 
	 continue;
       }
       else if (argument.CompareTo("-plot-Tgl")==0) {
	 HLTInfo("Ploting Tgl of Tracks");
	 fPlotTgl=kTRUE;  
	 continue;
       }
       else {
	 HLTError("unknown argument %s", argument.Data());
	 iResult=-EINVAL;
	 break;
       }
     }
     delete pTokens;
   }
   
   return iResult;
 }
 
void AliHLTTPCTrackHistoComponent::ReadTracks(const AliHLTComponentBlockData* iter,Int_t &tt){

  //AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
  //AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
  
  //HLTDebug ( "Input Data - TPC cluster - Slice/Patch: %d/%d.", slice, patch );
  const AliHLTTPCTrackletData* trackData = (const AliHLTTPCTrackletData*) iter->fPtr;
  AliHLTUInt32_t nTracks = trackData->fTrackletCnt;
  tt += nTracks;
  //HLTInfo("TrackHisto found %d Tracks in slice %d patch %d", nTracks, slice, patch);
  AliHLTTPCTrackSegmentData *tracks = (AliHLTTPCTrackSegmentData*) trackData->fTracklets;
  
  for(AliHLTUInt32_t i=0;i<nTracks;i++){
    fTracks.push_back(tracks[i]);
    UInt_t nHits = tracks->fNPoints;
    if(fPlotNClustersOnTracks || fPlotAll){fHistoNClustersOnTracks->Fill(nHits);}
    if(fPlotPT || fPlotAll){fHistoPT->Fill(tracks[i].fPt);}
    if(fPlotTgl || fPlotAll){fHistoTgl->Fill(tracks[i].fTgl);}
    const UInt_t *hitnum = tracks->fPointIDs;
    for(UInt_t h=0; h<nHits; h++){
      UInt_t idTrack = hitnum[h];
      Int_t sliceTrack = (idTrack>>25) & 0x7f;
      Int_t patchTrack = (idTrack>>22) & 0x7;
      UInt_t pos = idTrack&0x3fffff;
      fTrackClusterID[sliceTrack][patchTrack].push_back(pos);
    }
    UChar_t *tmpP = (UChar_t*)tracks;
    tmpP += sizeof(AliHLTTPCTrackSegmentData)+tracks->fNPoints*sizeof(UInt_t);
    tracks = (AliHLTTPCTrackSegmentData*)tmpP;
  } 
}

void AliHLTTPCTrackHistoComponent::PushHisto(){

  if(fPlotNClustersOnTracks || fPlotAll){
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5);
    PushBack( (TObject*) fHistoNClustersOnTracks,kAliHLTDataTypeHistogram, fSpecification);   
  }
  if(fPlotChargeClusters || fPlotAll){
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5);
    PushBack( (TObject*) fHistoAllClusters,kAliHLTDataTypeHistogram, fSpecification);   
  }
  if(fPlotChargeUsedClusters || fPlotAll){
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5);
    PushBack( (TObject*) fHistoUsedClusters,kAliHLTDataTypeHistogram, fSpecification);   
  }
  if(fPlotPT || fPlotAll){
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5);
    PushBack( (TObject*) fHistoPT,kAliHLTDataTypeHistogram, fSpecification);   
  }
  if(fPlotResidual || fPlotAll){
    for(unsigned int i=0;i<fTracks.size();i++){
      fHistoResidual->Fill(fTracks[i].fPsi); //Not rigth. Change here and x in Histo. Just for test.
    }
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5);
    PushBack( (TObject*) fHistoResidual,kAliHLTDataTypeHistogram, fSpecification);   
  }
  if(fPlotTgl || fPlotAll){
    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5);
    PushBack( (TObject*) fHistoTgl,kAliHLTDataTypeHistogram, fSpecification);   
  }
  
  fHistoNClustersOnTracks->Reset();                                 
  fHistoAllClusters->Reset();                                                                         
  fHistoUsedClusters->Reset();                                                                        
  fHistoPT->Reset();                                                                                
  fHistoResidual->Reset();                                                                            
  fHistoTgl->Reset();        
}
