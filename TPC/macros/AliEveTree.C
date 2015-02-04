/// \file AliEveTree.C
///
/// ~~~{.cpp}
/// .L $ALICE_ROOT/TPC/macros/AliEveTree.C+
/// MakeEveTree();
/// 
/// TFile fesd("AliESDs.root");
/// TTree * treeESD = (TTree*)fesd.Get("esdTree");
///
/// TFile ftpc("TPCdebug.root");
/// TFile fits("ITSdebug.root");
/// TTree * treeTPC = (TTree*)ftpc.Get("Transform");
/// TTree * treeITS = (TTree*)fits.Get("Clusters");
///   
/// treeTPC->Draw("gx2/sqrt(gx0^2+gx1^2):sqrt(gx0^2+gx1^2)>>his(100,0,200,100,-1,1)","event==4","")
/// treeITS->Draw("gz/sqrt(gx^2+gy^2):sqrt(gx^2+gy^2)>>his(100,0,200,100,-1,1)","event==4","same")
/// 
/// treeTPC->Draw("atan2(gx1,gx0):sqrt(gx0^2+gx1^2)>>his(100,0,200,100,-1,1)","event==4","*")
/// treeITS->Draw("atan2(gy,gx):sqrt(gx^2+gy^2)>>his(100,0,200,100,-1,1)","event==4","same*");
/// ~~~

#include "TFile.h"
#include "TTree.h"
#include "TTreeStream.h"
#include "TPolyMarker3D.h"
#include "TVectorD.h"

#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDfriendTrack.h"
#include "AliTrackPointArray.h"


AliTrackPointArray* MakeArray(TPolyMarker3D *pol);
void MakeESDTree(AliESDEvent*esd, TTreeSRedirector &cstream);

void MakeEveTree(){
  TFile fesd("AliESDs.root");
  TTree * esdTree = (TTree*)fesd.Get("esdTree");
  AliESDEvent cesd; 
  esdTree->SetBranchStatus("*",kTRUE);
  cesd.ReadFromTree(esdTree);
  TTreeSRedirector cstream("eveTree.root");
  for (Int_t ievent=0; ievent<esdTree->GetEntries();ievent++) {
    //cout<<endl<<endl<<"********* Processing event number: "<<event<<"*******\n";
    esdTree->GetEvent(ievent);
    AliESDfriend *fESDfriend =  (AliESDfriend*)cesd.FindListObject("AliESDfriend");
    cesd.SetESDfriend(fESDfriend);
    MakeESDTree(&cesd,cstream);
  }
}


void MakeESDTree(AliESDEvent*esd, TTreeSRedirector &cstream){
  ///

  Float_t bz = esd->GetMagneticField();
  //AliTPCseed dummyTPC;
  //AliTPCseed dummyTRD;
  for (Int_t i=0;i<esd->GetNumberOfTracks(); i++){
    AliESDtrack * track = esd->GetTrack(i);
    if (!track) continue;
    TPolyMarker3D polA;
    TPolyMarker3D polV;
    TPolyMarker3D polI;
    TPolyMarker3D polO;
    //
    AliTrackPointArray * arrayA=0;
    AliTrackPointArray * arrayV=0;
    AliTrackPointArray * arrayI=0;
    AliTrackPointArray * arrayO=0;
    //
    if (track->GetInnerParam()) {
      ((AliExternalTrackParam*)track->GetInnerParam())->FillPolymarker(&polA,bz,0,350,1);}
    else{
      track->FillPolymarker(&polA,bz,0,300,5);
    }
    arrayA=MakeArray(&polA);
    track->FillPolymarker(&polV,bz,0,60,1);
    arrayV=MakeArray(&polV);
    if (track->GetInnerParam()) {
      ((AliExternalTrackParam*)track->GetInnerParam())->FillPolymarker(&polI,bz,60,170,1);
      arrayI = MakeArray(&polI);
    }
    if (track->GetOuterParam()) {
      ((AliExternalTrackParam*)track->GetOuterParam())->FillPolymarker(&polO,bz,170,350,1);    
      arrayO=MakeArray(&polO);
    }
    static  AliTrackPointArray cldummy(5);
    AliTrackPointArray *clarray= &cldummy;
    const AliESDfriendTrack *ftrack = track->GetFriendTrack();
    if (ftrack && ftrack->GetTrackPointArray()) {
      clarray=(AliTrackPointArray *)ftrack->GetTrackPointArray();
    }
    Int_t event = esd->GetEventNumberInFile();
    Int_t id = track->GetID();
    cstream<<"Tracks"<<
      "eventNr="<<event<<
      "trackNr="<<id<<
      "Tr.="<<track<<
      "Cl.="<<clarray<<
      //
      "pA.="<<&polA<<
      "pV.="<<&polV<<
      "pI.="<<&polI<<
      "pO.="<<&polO<<
      //
      "aA.="<<arrayA<<
      "aV.="<<arrayV<<
      "aI.="<<arrayI<<
      "aO.="<<arrayO<<
      "\n";
  }
}



AliTrackPointArray *MakeArray(TPolyMarker3D *pol){
  /// Make a aray of  points with errors

  Int_t entries = pol->GetN();
  AliTrackPointArray * array = new AliTrackPointArray(entries);
  for (Int_t i=0;i<entries;i++){
    Double_t xyz[3]={0,0,0};
    pol->GetPoint(i,xyz[0],xyz[1],xyz[2]);
    ((Float_t*)array->GetX())[i]=xyz[0];
    ((Float_t*)array->GetY())[i]=xyz[1];
    ((Float_t*)array->GetZ())[i]=xyz[2];
  }
  return array;
}

