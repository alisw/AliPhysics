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
#include "TNtuple.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "AliHLTTPCTrackArray.h"
#include "AliHLTTPCTrack.h"

//#include "AliHLTTPC.h"
//#include <stdlib.h>
//#include <cerrno>

AliHLTTPCTrackHistoComponent gAliHLTTPCTrackHistoComponent;

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCTrackHistoComponent)

AliHLTTPCTrackHistoComponent::AliHLTTPCTrackHistoComponent()
:
fClusters(NULL),
  fTracks(NULL),
  fTracksArray(NULL)
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
  return kAliHLTDataTypeTNtuple;

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
 
  fClusters = new TNtuple("fCluster", "fCluster", "charge:qmax:residualY:residualZ:used:event"); 
  fTracks = new TNtuple("fTracks", "fTracks", "pt:eta:psi:nclusters:event"); 
  fTracksArray=new AliHLTTPCTrackArray();

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
  
  delete fClusters;
  delete fTracks;
  delete fTracksArray;

  return 0;
}

int AliHLTTPCTrackHistoComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/)
{
  if(!fTracksArray){fTracksArray=new AliHLTTPCTrackArray();}

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

    AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
    AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );

    //HLTDebug ( "Input Data - TPC cluster - Slice/Patch: %d/%d.", slice, patch );
    const AliHLTTPCClusterData* clusterData = (const AliHLTTPCClusterData*) iter->fPtr;
    Int_t nSpacepoint = (Int_t) clusterData->fSpacePointCnt;    
    TotalSpacePoint += nSpacepoint;
    //HLTInfo("TrackHisto found %d Spacepoints in slice %d patch %d", nSpacepoint, slice, patch);
    AliHLTTPCSpacePointData *clusters = (AliHLTTPCSpacePointData*) clusterData->fSpacePoints;
    
    if (fClustersArray[slice][patch]!=NULL) {
      delete(fClustersArray[slice][patch]);
      fClustersArray[slice][patch]=NULL;
    }
    Int_t arraysize=nSpacepoint*sizeof(AliHLTTPCSpacePointData);
    fClustersArray[slice][patch] = (AliHLTTPCSpacePointData*)new Byte_t[arraysize];
    if (fClustersArray[slice][patch]) {
      memcpy(fClustersArray[slice][patch], clusters, arraysize);
      fNcl[slice][patch]=nSpacepoint;
    } else {
      fNcl[slice][patch]=nSpacepoint;
      HLTError ( "Memory allocation failed!" );
    }

    for(int i=0;i<nSpacepoint;i++){
      UInt_t idCluster = clusters[i].fID;
      Int_t sliceCl = (idCluster>>25) & 0x7f;
      Int_t patchCl = (idCluster>>22) & 0x7;
      UInt_t pos = idCluster&0x3fffff;
      Int_t used = 0;
      Float_t resy = 0, resz = 0;
      for(UInt_t id=0;id<fTrackClusterID[sliceCl][patchCl].size();id++){
	if(fTrackClusterID[sliceCl][patchCl][id]==pos){
	  clusters[i].fUsed=kTRUE;
	  nClustersUsed++;
	  used=1;
	  FillResidual(pos,sliceCl,patchCl,resy,resz);
	}
      }
      if(used==1){
	fClusters->Fill(clusters[i].fCharge,clusters[i].fQMax,resy,resz,used,GetEventId()); 
      }
      else{
	fClusters->Fill(clusters[i].fCharge,clusters[i].fQMax,-100,-100,used,GetEventId()); 
      }
    }
  } 
  
  HLTInfo("TrackHisto found %d Spacepoints",TotalSpacePoint);
  HLTInfo("TrackHisto found %d Tracks",TotalTrack);
  
  PushHisto();

  delete fTracksArray;
  fTracksArray=NULL;

  for(UInt_t i=0;i<36;i++){
    for(UInt_t j=0;j<6;j++){ 
      fTrackClusterID[i][j].clear();
    }
  }

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

       HLTError("unknown argument %s", argument.Data());
       iResult=-EINVAL;
       break;
     }
     delete pTokens;
   }
   
   return iResult;
 }
 
void AliHLTTPCTrackHistoComponent::ReadTracks(const AliHLTComponentBlockData* iter,Int_t &tt){

  AliHLTUInt8_t slice = AliHLTTPCDefinitions::GetMinSliceNr( *iter );
  //AliHLTUInt8_t patch = AliHLTTPCDefinitions::GetMinPatchNr( *iter );
  
  //HLTDebug ( "Input Data - TPC cluster - Slice/Patch: %d/%d.", slice, patch );
  AliHLTTPCTrackletData* trackData = (AliHLTTPCTrackletData*) iter->fPtr;
  AliHLTUInt32_t nTracks = trackData->fTrackletCnt;
  fTracksArray->FillTracksChecked(trackData->fTracklets,trackData->fTrackletCnt,iter->fSize,slice,true);
  tt += nTracks;
  //HLTInfo("TrackHisto found %d Tracks in slice %d patch %d", nTracks, slice, patch);
  AliHLTTPCTrackSegmentData *tracks = (AliHLTTPCTrackSegmentData*) trackData->fTracklets;
  
  for(AliHLTUInt32_t i=0;i<nTracks;i++){
    UInt_t nHits = tracks->fNPoints;
    
    fTracks->Fill(tracks->fPt,tracks->fPsi,tracks->fTgl,nHits,GetEventId()); 

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

    AliHLTUInt32_t fSpecification = AliHLTTPCDefinitions::EncodeDataSpecification(0,35,0,5);
    PushBack( (TObject*) fTracks,kAliHLTDataTypeTNtuple, fSpecification);   
    PushBack( (TObject*) fClusters,kAliHLTDataTypeTNtuple, fSpecification);   
  
}
void AliHLTTPCTrackHistoComponent::FillResidual( UInt_t pos,AliHLTUInt8_t slice,AliHLTUInt8_t patch,Float_t& resy,Float_t& resz){

  AliHLTTPCSpacePointData *cl =  &fClustersArray[slice][patch][pos];
  if(!cl){return;}

  AliHLTTPCTrack *gtrack = NULL;

  for(int i;i<fTracksArray->GetNTracks();i++){
    AliHLTTPCTrack *tt = fTracksArray->GetCheckedTrack(i); 
    UInt_t *hitnum =tt->GetHitNumbers();
    Int_t nHits = tt->GetNHits();
    for(Int_t h=0; h<nHits; h++){
      UInt_t id=hitnum[h];
      Int_t Tslice = (id>>25) & 0x7f;
      Int_t Tpatch = (id>>22) & 0x7;
      UInt_t Tpos = id&0x3fffff; 
      if(Tslice==slice && Tpatch==patch && Tpos==pos) {
	gtrack = tt; 
	break;
      }
    }
  }
  
  if(!gtrack){return;}

  Int_t tslice = gtrack->GetSector();
  Double_t radius = gtrack->GetRadius();      // radius
  Double_t kappa = gtrack->GetKappa();        // curvature = 1/R , signed
  Double_t lambda = atan( gtrack->GetTgl() ); // dipAngle lambda

  // ------------------------------------
  // ++ Get first/last point of the track
  
  Double_t xyzL[3];      // lastpoint of track
  Double_t xyzF[3];      // firstpoint of track
  
  xyzF[0] = gtrack->GetFirstPointX();
  xyzF[1] = gtrack->GetFirstPointY();
  xyzF[2] = gtrack->GetFirstPointZ();
  
  xyzL[0] = gtrack->GetLastPointX();
  xyzL[1] = gtrack->GetLastPointY();
  xyzL[2] = gtrack->GetLastPointZ();

  // --------------------------
  // ++ Calculate length of the track
  
  Double_t s = 0.;       // length of the track
  if (  AliHLTTPCTransform::GetBFieldValue() == 0. || kappa == 0 ) 
    s = sqrt ( (xyzL[0] - xyzF[0])*(xyzL[0] - xyzF[0]) + (xyzL[1] - xyzF[1])*(xyzL[1] - xyzF[1]) ); 
  else {
    // Calculate the length of the track. If it is to flat in in s,z plane use sxy, otherwise use sz
    if (fabs(lambda) > 0.05){
      // length of track calculated out of z
      s = fabs( (xyzL[2] - xyzF[2]) / sin(lambda) ); // length of track calculated out of z
    }
    else {
      Double_t d = (xyzL[0] - xyzF[0])*(xyzL[0] - xyzF[0]) + (xyzL[1] - xyzF[1])*(xyzL[1] - xyzF[1]); 
      // length of track calculated out of xy
      s = fabs ( acos( 0.5 * (2 - (d / (radius*radius)))) / ( kappa * cos(lambda) ) ); 		
    }
  }
  
  gtrack->Rotate(tslice,kTRUE);
  
  Double_t padrows = 0;                   
    
  Float_t xyzC[3];       // cluster tmp
  Float_t xyzTtmp[3];    // track tmp

  xyzC[0] = cl->fX;
  xyzC[1] = cl->fY;
  xyzC[2] = cl->fZ;
 
  Int_t padrow = AliHLTTPCTransform::GetPadRow(cl->fX);

  xyzTtmp[0] = gtrack->GetFirstPointX();

  if(gtrack->GetCrossingPoint(padrow,xyzTtmp)) {
    // ----------------------
       // ++ Calculate Residuals
       
       Float_t deltaY = ( xyzC[1] - xyzTtmp[1] );
       Float_t deltaZ = ( xyzC[2] - xyzTtmp[2] );
       
       resy = deltaY;
       resz = deltaZ;
  }
  else{
    resy = -1000;
    resz = -1000;
  }
}
