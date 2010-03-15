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
#include "AliEveHOMERManager.h"
#include "TEveManager.h"
#include "AliMUONConstants.h"

ClassImp(AliHLTEveMuon);

AliHLTEveMuon::AliHLTEveMuon() : 
  AliHLTEveBase(),
  fTracks(NULL),
  fClusters(NULL)
{
  // Constructor.
}

AliHLTEveMuon::~AliHLTEveMuon()
{
  //Destructor
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
      fEventManager->GetEveManager()->AddElement(fClusters);
    }
    ProcessClusters( block, fClusters );
    
  }else if(block->GetDataType().CompareTo("MANTRACK") == 0){
    
    if ( !fTracks ) {
      fTracks = CreateTrackSet(); 
      fEventManager->GetEveManager()->AddElement(fTracks);
      gEve->AddElement(fTracks);
    }
    
    ProcessTracks( block, fTracks );
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
  TEveStraightLineSet * ls = new TEveStraightLineSet("MUON Tracks");
  ls->SetMainColor(kRed);
  ls->SetLineWidth(3);
  return ls;
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
