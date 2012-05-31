/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Svein Lindal <slindal@fys.uio.no   >                  *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/// @file   AliHLTEvePhos.cxx
/// @author Svein Lindal <slindal@fys.uio.no>
/// @brief  Muon processor for the HLT EVE display

#include "AliHLTEveMuon.h"
#include "AliHLTMUONDataBlockReader.h"
#include "AliHLTHOMERBlockDesc.h"
#include "TCanvas.h"
#include "TEveStraightLineSet.h"
#include "TEvePointSet.h"
#include "AliEveHLTEventManager.h"
#include "TEveManager.h"


#include "TEveVSDStructs.h"
#include "TGeoGlobalMagField.h"
#include "AliMUONTrack.h"
#include "AliMUONTrackParam.h"
#include "AliMUONESDInterface.h"
#include "AliMUONVTrackReconstructor.h"
#include "AliEveMUONTrack.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONUtils.h"
#include "AliMUONVCluster.h"
#include "AliMUONVClusterStore.h"
#include "AliMUONConstants.h"
#include "TEveTrackPropagator.h"

using namespace std;

class AliHLTMUONUtils;
class AliEveMuonTrack;

ClassImp(AliHLTEveMuon);

AliHLTEveMuon::AliHLTEveMuon() : 
  AliHLTEveBase("Muon"),
  fFullTrackList(NULL),
  fTracks(NULL),
  fClusters(NULL)
{
  // Constructor.
  SetMaxHistograms(6);
}

AliHLTEveMuon::~AliHLTEveMuon()
{
  //Destructor
  if (fFullTrackList)
    delete fFullTrackList;
  fFullTrackList = NULL;
  
  if (fTracks)
    delete fTracks;
  fTracks = NULL;

  if(fClusters)
    delete fClusters;
  fClusters = NULL;
}


void AliHLTEveMuon::ProcessBlock(AliHLTHOMERBlockDesc * block) {
  //See header file for documentation
  if ( (block->GetDataType().CompareTo("RECHITS") == 0) || (block->GetDataType().CompareTo("TRIGRECS") == 0) ) {
    if(!fClusters) {
      fClusters = CreateClusters();
      AddElement(fClusters);
    }
    ProcessClusters( block, fClusters );
    
  }else if(block->GetDataType().CompareTo("MANTRACK") == 0){
    
    if ( !fTracks ) {
      fTracks = CreateTrackSet(); 
      AddElement(fTracks);
    }
    
    ProcessTracks( block, fTracks );

  }else if(block->GetDataType().CompareTo("TRACKS") == 0){
    
    if ( !fFullTrackList ) {
      fFullTrackList = CreateFullTrackList(); 
      AddElement(fFullTrackList);
    }
    
    ProcessFullTracks( block,  fFullTrackList );

  } else if(block->GetDataType().CompareTo("ROOTHIST") == 0) {
    ProcessHistogram(block);
  }
 
}

TEvePointSet * AliHLTEveMuon::CreateClusters() {
  //See header file for documentation
  TEvePointSet * ps = new TEvePointSet("MUON RecHits");
  ps->SetMainColor(kBlue);
  ps->SetMarkerStyle(20);
  return ps;
}

TEveStraightLineSet * AliHLTEveMuon::CreateTrackSet() {
  // See header file
  TEveStraightLineSet * lineset = new TEveStraightLineSet("MUON Tracks");
  lineset->SetMainColor(kRed);
  lineset->SetLineWidth(3);
  return lineset;
}

TEveTrackList * AliHLTEveMuon::CreateFullTrackList(){
  // See header file
  TEveTrackList * lineset = new TEveTrackList("MUON Full Tracks");
  lineset->SetMainColor(kBlue);
  return lineset;
}

void AliHLTEveMuon::ProcessHistogram(AliHLTHOMERBlockDesc * block ) {
  //See header file for documentation
  if(!fCanvas) {
    fCanvas = CreateCanvas("MUON QA", "MUON QA");
    fCanvas->Divide(3, 2);
  }
  AddHistogramsToCanvas(block, fCanvas, fHistoCount);
}

void AliHLTEveMuon::UpdateElements() {
  //See header file for documentation
  if(fCanvas) fCanvas->Update();
  if(fClusters) fClusters->ResetBBox();
  if(fTracks) fTracks->ElementChanged();
}

void AliHLTEveMuon::ResetElements(){
  //See header file for documentation
  fHistoCount = 0;
  
  if ( fClusters ) fClusters->Reset();
  if ( fTracks ){
    fTracks->Destroy();
    fTracks = NULL;
  }
  if ( fFullTrackList ){
    fFullTrackList->Destroy();
    fFullTrackList = NULL;
  }


}

void AliHLTEveMuon::ProcessClusters(AliHLTHOMERBlockDesc * block, TEvePointSet * clusters) {
    //See header file for documentation
  unsigned long size = block->GetSize();
  Int_t * buffer ;
  
  buffer = (Int_t *)block->GetData();
  //cout<<"block size : "<<size<<", buffer : "<<buffer<<", DataType : "<<block->GetDataType()<<endl;

  if(block->GetDataType().CompareTo("RECHITS") == 0){
    
    AliHLTMUONRecHitsBlockReader trackblock((char*)buffer, size);
    const AliHLTMUONRecHitStruct* hit = trackblock.GetArray();
    
    for(AliHLTUInt32_t ientry = 0; ientry < trackblock.Nentries(); ientry++){
      if(hit->fX!=0.0 && hit->fY!=0.0 && hit->fZ!=0.0)
	clusters->SetNextPoint(hit->fX,hit->fY,hit->fZ);
      hit++;
      
    }// track hit loop
  }
  
  else{// if rechits
    //     if(!strcmp((BlockType(ULong64_t(reader->GetBlockDataType(i)))).Data(),"TRIGRECS")){
    
    AliHLTMUONTriggerRecordsBlockReader trigblock(buffer, size);
    const AliHLTMUONTriggerRecordStruct* trigrec = trigblock.GetArray();
    for(AliHLTUInt32_t ientry = 0; ientry < trigblock.Nentries(); ientry++){
      
      const AliHLTMUONRecHitStruct* hit = &trigrec->fHit[0];
      for(AliHLTUInt32_t ch = 0; ch < 4; ch++)
	{
	  if(hit->fX!=0.0 && hit->fY!=0.0 && hit->fZ!=0.0)
	    clusters->SetNextPoint(hit->fX,hit->fY,hit->fZ);
	  hit++;
	}// trig chamber loop
      trigrec++;
    }//trig hit loop
  }//else trigger
  
}

void AliHLTEveMuon::ProcessTracks(AliHLTHOMERBlockDesc * block, TEveStraightLineSet * tracks) {
  //See header file for documentation  
  unsigned long size = block->GetSize();
  Int_t * buffer = (Int_t *)block->GetData();
  AliHLTMUONRecHitStruct hit1,hit2;
  hit1.fX = hit1.fY = hit1.fZ = hit2.fX = hit2.fY = hit2.fZ = 0;
  Int_t ch1=0, ch2=0;
  Float_t x0=0.0,y0=0.0,z0=0.0;
  Float_t x3=0.0,y3=0.0,z3=0.0;
  if(block->GetDataType().CompareTo("MANTRACK") == 0){  
    AliHLTMUONMansoTracksBlockReader mantrackblock(buffer, size);
    const AliHLTMUONMansoTrackStruct* mtrack = mantrackblock.GetArray();
    for(AliHLTUInt32_t ientry = 0; ientry < mantrackblock.Nentries(); ientry++){
      const AliHLTMUONRecHitStruct* hit = &mtrack->fHit[0];
      for(AliHLTUInt32_t ch = 0; ch < 4; ch++){
	// cout << setw(10) << left << ch + 7 << setw(0);
	// cout << setw(13) << left << hit->fX << setw(0);
	// cout << setw(13) << left << hit->fY << setw(0);
	// cout << hit->fZ << setw(0) << endl;
	if(hit->fZ != 0.0){
	  if(ch==0 || ch==1){
	    hit1 = *hit; ch1 = ch+6;
	  }else{
	    hit2 = *hit; ch2 = ch+6;
	  }
	}
	hit++;
      }// trig chamber loop
      // printf("ch : %d, (X,Y,Z) : (%f,%f,%f)\n",ch1,hit1.fX,hit1.fY,hit1.fZ);
      // printf("ch : %d, (X,Y,Z) : (%f,%f,%f)\n",ch2,hit2.fX,hit2.fY,hit2.fZ);
      // meminfo();
      z3 = AliMUONConstants::DefaultChamberZ(ch2+4);
      y3 =  hit1.fY - (hit1.fZ-z3)*(hit1.fY - hit2.fY)/(hit1.fZ - hit2.fZ) ;
      x3 =  hit1.fX - (hit1.fZ-z3)*(hit1.fX - hit2.fX)/(hit1.fZ - hit2.fZ) ;

      z0 = AliMUONConstants::DefaultChamberZ(ch1);
      y0 =  hit1.fY - (hit1.fZ-z0)*(hit1.fY - hit2.fY)/(hit1.fZ - hit2.fZ) ;
      x0 =  hit1.fX - (hit1.fZ-z0)*(hit1.fX - hit2.fX)/(hit1.fZ - hit2.fZ) ;
      

      tracks->AddLine(x0,y0,z0,x3,y3,z3);
      mtrack++;
    }
    //    cout<<"NofManso Tracks : "<<mantrackblock.Nentries()<<endl;
  }
}

int AliHLTEveMuon::MakeMUONTrack(AliMUONTrack *muonTrack, const AliHLTMUONTrackStruct *muonHLTTrack)
{
  // See header for documentation
  AliHLTUInt32_t clusterIndex = 0;  // for the cluster unique ID.			
  AliHLTMUONParticleSign sign;
  bool hitset[16];
  AliHLTMUONUtils::UnpackTrackFlags(
  				    muonHLTTrack->fFlags, sign, hitset
  				    );
  
  // add track parameters at vertex
  TVector3 mom(muonHLTTrack->fPx, muonHLTTrack->fPy, muonHLTTrack->fPz);
  AliMUONTrackParam paramAtVtx;
  if (mom.Mag() != 0)
    paramAtVtx.SetInverseBendingMomentum(muonHLTTrack->fInverseBendingMomentum);
  else
    paramAtVtx.SetInverseBendingMomentum(0.);
  paramAtVtx.SetNonBendingSlope(TMath::Tan(muonHLTTrack->fThetaX));
  paramAtVtx.SetBendingSlope(TMath::Tan(muonHLTTrack->fThetaY));
  paramAtVtx.SetZ(muonHLTTrack->fZ);
  paramAtVtx.SetBendingCoor(muonHLTTrack->fY);
  paramAtVtx.SetNonBendingCoor(muonHLTTrack->fX);
  muonTrack->SetTrackParamAtVertex(&paramAtVtx);

  //printf("(X,Y,Z) : (%8.3f,%8.3f,%8.3f)\n",muonHLTTrack->fX,muonHLTTrack->fY,muonHLTTrack->fZ);

  // add clusters
  Int_t nHits = 0;
  AliMUONVClusterStore* cStore = AliMUONESDInterface::NewClusterStore();
  if (!cStore) return -1;
  AliMUONVCluster* cluster = cStore->CreateCluster(0,0,0);
  AliMUONTrackParam trackParam;
  for (int i = 0; i < 16; i++)
    {
      if (not hitset[i]) continue;
				
      AliHLTUInt8_t chamber;
      AliHLTUInt16_t detElemId;
      AliHLTMUONUtils::UnpackRecHitFlags((muonHLTTrack->fHit[i]).fFlags, chamber, detElemId);
      
      cluster->SetUniqueID(AliMUONVCluster::BuildUniqueID(chamber, detElemId, clusterIndex++));
      cluster->SetXYZ((muonHLTTrack->fHit[i]).fX, (muonHLTTrack->fHit[i]).fY, (muonHLTTrack->fHit[i]).fZ);
      cluster->SetErrXY(    // Use nominal values.
			AliHLTMUONConstants::DefaultNonBendingReso(),
			AliHLTMUONConstants::DefaultBendingReso()
			);
      cluster->SetCharge(-1.);   // Indicate no total charge calculated.
      cluster->SetChi2(-1.);   // Indicate no fit made.
      trackParam.SetZ(cluster->GetZ());
      muonTrack->AddTrackParamAtCluster(trackParam, *cluster, kTRUE);
      nHits++;
    }
  
  // compute track parameters at each cluster
  if (nHits > 0) {
    AliMUONTrackParam *firstTrackParam = (AliMUONTrackParam*) muonTrack->GetTrackParamAtCluster()->First();
    trackParam = (*firstTrackParam);
    if (!AliMUONESDInterface::GetTracker()) AliMUONESDInterface::ResetTracker();
    if (!AliMUONESDInterface::GetTracker()->RefitTrack(*muonTrack, kFALSE) &&
	muonTrack->GetGlobalChi2() < AliMUONTrack::MaxChi2()) {
      *firstTrackParam = trackParam;
      muonTrack->UpdateCovTrackParamAtCluster();
    }
  }

  muonTrack->SetGlobalChi2(muonHLTTrack->fChi2);
  
  return 0;
}

Int_t AliHLTEveMuon::ProcessFullTracks(AliHLTHOMERBlockDesc * block, TEveTrackList * fullTracks) {

  // See header for documentation 

  Int_t iResult = 0;

  Double_t b[3], x[3];
  x[0] = 0.0 ; x[1] = 0.0 ; x[2] = -950.0;
  TGeoGlobalMagField::Instance()->Field(x,b);
  //" Field at (0.0, 0.0, -950.0) [at the middle of dipole magnet] 
  //should be (6.79, 0.03, -0.17) or similar value with change of sign"
  if(TMath::AreEqualAbs(b[0],0.0,1.0e-5) and TMath::AreEqualAbs(b[1],0.0,1.0e-5) and TMath::AreEqualAbs(b[2],0.0,1.0e-5)){
    printf("At (X,Y,Z) : (%6.2lf,%6.2lf,%6.2lf) Field (Bx,By,Bz) is (%6.2lf,%6.2lf,%6.2lf)\n",
	   x[0],x[1],x[2],b[0],b[1],b[2]);    
    cerr<<"Magnetic field is not properly set, MUON tracking will not possble"<<endl;
    return 1;
  }
  
  

  TEveRecTrack rt;
  
  unsigned long size = block->GetSize();
  Int_t * buffer = (Int_t *)block->GetData();

  AliHLTMUONTracksBlockReader muontrackblock(buffer, size);
  const AliHLTMUONTrackStruct* mtrack = muontrackblock.GetArray();
  //cout<<"NofTracks : "<<muontrackblock.Nentries()<<endl;
  for(AliHLTUInt32_t ientry = 0; ientry < muontrackblock.Nentries(); ientry++){
    
    AliMUONTrack *muonTrack = new AliMUONTrack();
    MakeMUONTrack(muonTrack,mtrack);
    if(muonTrack->GetNClusters()==0){
      delete muonTrack;
      continue;
    }
    
    rt.fLabel = ientry;
    AliEveMUONTrack* track = new AliEveMUONTrack(&rt, fullTracks->GetPropagator());
    track->MakeMUONTrack(muonTrack);
    //track->SetTitle(Form("HLT Track : %d, pt : %lf",ientry,TMath::Sqrt(((mtrack->fPx * mtrack->fPx) + (mtrack->fPy * mtrack->fPy)))));
    track->SetName(Form("HLT Track : %d, pt : %lf",ientry,TMath::Sqrt(((mtrack->fPx * mtrack->fPx) + (mtrack->fPy * mtrack->fPy)))));
    fullTracks->AddElement(track);
    
    mtrack++;
  }//track loop
  
  return iResult;

}
