
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
  // Load libraries

  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");
  
  
  .x ~/NimStyle.C
  gSystem->Load("libANALYSIS");
  gSystem->Load("libTPCcalib");

  // analyze results

  TFile f("CalibObjectsTrain2.root");
  AliTPCcalibMaterial *calibMaterial = (AliTPCcalibMaterial *)f->Get("alignMaterial");


*/

#include "Riostream.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THnSparse.h"
#include "TList.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TF1.h"
#include "TVectorD.h"
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TCanvas.h"

#include "AliTPCclusterMI.h"
#include "AliTPCseed.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDfriend.h"
#include "AliESDInputHandler.h"
#include "AliAnalysisManager.h"

#include "AliTracker.h"
#include "AliMagF.h"
#include "AliTPCCalROC.h"

#include "AliLog.h"

#include "AliTPCcalibMaterial.h"

#include "TTreeStream.h"
#include "AliTPCTracklet.h"
#include "TTimeStamp.h"
#include "AliTPCcalibDB.h"
#include "AliTPCcalibLaser.h"
#include "AliDCSSensorArray.h"
#include "AliDCSSensor.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"

#include "AliESDTZERO.h"

ClassImp(AliTPCcalibMaterial)

AliTPCcalibMaterial::AliTPCcalibMaterial():
  AliTPCcalibBase("calibMaterial","calibMaterial"),
  fHisMaterial(0),
  fHisMaterialRPhi(0)
{
  
}

AliTPCcalibMaterial::AliTPCcalibMaterial(const char * name, const char * title):
  AliTPCcalibBase(name,title),
  fHisMaterial(0),
  fHisMaterialRPhi(0)
{
  //
  //
  //
}

AliTPCcalibMaterial::~AliTPCcalibMaterial(){
  //
  // delete histograms
  // class is owner of all histograms
  //
  if (!fHisMaterial) return;
  delete fHisMaterial;
  delete fHisMaterialRPhi;
  fHisMaterial=0;
}


Long64_t AliTPCcalibMaterial::Merge(TCollection *li) {
  //
  // Merge histograms
  //
  TIterator* iter = li->MakeIterator();
  AliTPCcalibMaterial* cal = 0;

  while ((cal = (AliTPCcalibMaterial*)iter->Next())) {
    if (!cal->InheritsFrom(AliTPCcalibMaterial::Class())) {
      Error("Merge","Attempt to add object of class %s to a %s", cal->ClassName(), this->ClassName());
      return -1;
    }
    AliTPCcalibMaterial* calib= (AliTPCcalibMaterial*)(cal);
    //
    if (!fHisMaterial) fHisMaterial=MakeHisto();
    fHisMaterial->Add(calib->fHisMaterial);
    fHisMaterialRPhi->Add(calib->fHisMaterialRPhi);
  }
  return 0;
}



void AliTPCcalibMaterial::Process(AliESDEvent *event){
  //
  //
  //
  const Int_t kMinCl=40;
  const Float_t kMinRatio=0.7;
  const Float_t kMaxS=0.05;
  const Float_t kMinDist=5;
  const Double_t kStep=1.;
  if (!event) return;
  ProcessPairs(event);
  return;
  //  TTreeSRedirector * cstream =  GetDebugStreamer();
  //
  if (!fHisMaterial){
    MakeHisto();
  }
  
  //  
  // fill histogram of track prolongations
  Float_t dca[2];
  Int_t ntracks = event->GetNumberOfTracks();
  for (Int_t itrack=0; itrack<ntracks; itrack++){
    AliESDtrack *track=event->GetTrack(itrack);
    if (!track) continue;
    if (track->GetTPCNcls()<=kMinCl) continue;
    if ((1.+track->GetTPCNcls())/(1.+track->GetTPCNclsF())<=kMinRatio) continue;
    if ((1.+track->GetTPCnclsS())/(1.+track->GetTPCNcls())>kMaxS) continue;
    if (!track->GetInnerParam()) continue;
    if (track->GetKinkIndex(0)!=0) continue;
    //
    track->GetImpactParameters(dca[0],dca[1]);
    if (TMath::Abs(dca[0])<kMinDist && TMath::Abs(dca[1])<kMinDist) continue;
    AliExternalTrackParam param(*(track->GetInnerParam()));
    if (!AliTracker::PropagateTrackTo(&param,90,0.0005,10,kTRUE)) continue;
    Double_t x[5]={0,0,0,TMath::Sqrt(TMath::Abs(param.GetP()))*param.GetSign(),TMath::Sqrt(TMath::Abs(track->GetTPCsignal()))};
    //
    //
    for (Float_t radius=90; radius>0; radius-=kStep){
      if (!AliTracker::PropagateTrackTo(&param,radius,0.0005,kStep*0.5,kTRUE)) break;
      if (TMath::Abs(param.GetSnp())>0.8) break;
      param.GetXYZ(x);
      Double_t weight=1./TMath::Sqrt(1.+param.GetSnp()*param.GetSnp()+param.GetTgl()*param.GetTgl());
      fHisMaterial->Fill(x,weight);    
      Double_t r = TMath::Sqrt(x[0]*x[0]+x[1]*x[1]);
      Double_t phi = TMath::ATan2(x[1],x[0]);
      x[0]=r;
      x[1]=phi;
      fHisMaterialRPhi->Fill(x,weight);
    }
  }
}

THnSparse *AliTPCcalibMaterial::MakeHisto(){
  //
  // Make track prolongation histogram
  // 
  //
  //                    gX       gY     gz   p       dEdx
  Int_t    bins[5]   = {100,    100,   300,  40,   100};
  Double_t xmin[5]   = {-100,  -100,  -300,  -2,   5};
  Double_t xmax[5]   = {100,    100,   300,   2,   33};
  TString  axisName[5]={
    "gx",
    "gy",
    "gz",
    "p",
    "dedx"
  };
  TString  axisTitle[5]={
    "x    (cm)",
    "y    (cm)",
    "z    (cm)",
    "p    (GeV)",
    "dedx (a.u)"
  };

  Int_t    binsR[5]   = {30,    360,     300,  40,   100};
  Double_t xminR[5]   = { 0,    -3.14,  -300,  -2,   5};
  Double_t xmaxR[5]   = {30,    3.14,    300,   2,   33};
  TString  axisNameR[5]={
    "r",
    "rphi",
    "z",
    "p",
    "dedx"
  };
  TString  axisTitleR[5]={
    "r    (cm)",
    "rphi    (cm)",
    "z    (cm)",
    "p    (GeV)",
    "dedx (a.u)"
  };

  THnSparse *sparse = new THnSparseF("his_Material", "His Material", 5, bins, xmin, xmax);
  THnSparse *sparseR = new THnSparseF("his_MaterialRPhi", "His Material Rphi", 5, binsR, xminR, xmaxR);
  for (Int_t iaxis=0; iaxis<5; iaxis++){
    sparse->GetAxis(iaxis)->SetName(axisName[iaxis]);
    sparse->GetAxis(iaxis)->SetTitle(axisTitle[iaxis]);
    sparseR->GetAxis(iaxis)->SetName(axisNameR[iaxis]);
    sparseR->GetAxis(iaxis)->SetTitle(axisTitleR[iaxis]);
  }
  fHisMaterial=sparse;
  fHisMaterialRPhi=sparseR;
  return sparse;
}




void  AliTPCcalibMaterial::ProcessPairs(AliESDEvent *event){
  //
  // Process pairs of tracks to get a material budget map
  //
  //
  if (!event) return;
  if (event) AliKFParticle::SetField(event->GetMagneticField());  // set mean magnetic field for KF particles
  //
  // 1. Calculate total dEdx for all TPC tracks
  //
  const Int_t kMinCl=70;
  const Double_t kEpsilon=0.000001;
  const Float_t kMinRatio=0.7;
  const Float_t kMinDist=1.5;
  const Float_t kMinDistChi2=8;        // 
  const Float_t kMaxDistZ=280;         // max distanceZ
  const Float_t kMaxDistR=250;          // max distanceR
  const Double_t kMaxChi2   =36;     // maximal chi2 to define the vertex
  const Double_t kMaxDistVertexSec=3;     // maximal distance to secondary vertex

  if (!event) return;
  AliESDVertex *vertexSPD =  (AliESDVertex *)event->GetPrimaryVertexSPD();
  AliESDVertex * spdVertex   = (AliESDVertex *)event->GetPrimaryVertexSPD();
  AliESDVertex * trackVertex = (AliESDVertex *)event->GetPrimaryVertexTracks();
  AliESDVertex * tpcVertex   = (AliESDVertex *)event->GetPrimaryVertexTPC();
  AliESDTZERO  * tzero       = (AliESDTZERO  *)event->GetESDTZERO() ;
  //
  Double_t tpcSignalTotPrim=0; 
  Double_t tpcSignalTotSec=0; 
  Int_t ntracksTPC=0;
  Int_t nTPCPrim=0; 
  Int_t nTPCSec=0;  
  Int_t ntracks=event->GetNumberOfTracks();
  if ( ntracks<=2 ) return;
  
  //
  Float_t dca[2]={0};
  Float_t cov[3]={0};
  Float_t dca0[2]={0};
  Float_t dca1[2]={0};
  //
  //1. Calculate total dEdx for primary and secondary tracks
  //   and count primaries and secondaries
  Int_t *rejectTrack = new Int_t[ntracks];
 
  for (Int_t itrack=0; itrack<ntracks; itrack++){
    AliESDtrack *track=event->GetTrack(itrack);
    rejectTrack[itrack]=0;
    if (!track) continue;
    if (track->GetTPCNcls()<=kMinCl) continue;  // skip short tracks
    ntracksTPC++;
    if ((1.+track->GetTPCNcls())/(1.+track->GetTPCNclsF())<=kMinRatio) continue;
    if (!track->GetInnerParam()) continue;  // skip not TPC tracks
    if (track->GetKinkIndex(0)!=0)  {rejectTrack[itrack]+=16;continue;} // skip kinks
    track->GetImpactParameters(dca[0],dca[1]);
    if (TMath::Abs(dca[0])>kMaxDistR && TMath::Abs(dca[1])>kMaxDistZ) continue;
    // remove too dip secondaries
    if (TMath::Abs(dca[0])<kMinDist && TMath::Abs(dca[1])<kMinDist){
      tpcSignalTotPrim+=track->GetTPCsignal();
      nTPCPrim++;
    }else{
      tpcSignalTotSec+=track->GetTPCsignal();
      nTPCSec++;
    };
    if (cov[0]>kEpsilon &&TMath::Abs(dca[0])>kEpsilon &&TMath::Sqrt(dca[0]*dca[0]/(TMath::Abs(cov[0])))<kMinDistChi2) rejectTrack[itrack]+=1;  // primary
    if (cov[0]>kEpsilon &&TMath::Abs(dca[0])>kEpsilon &&TMath::Abs(dca[0])<kMinDist)  rejectTrack[itrack]+=1;  // primary
    if (track->GetTPCsignal()<40) rejectTrack[itrack]+=16;
    //
    if (CheckLooper(itrack, event))   rejectTrack[itrack]+=2;   // looper
    if (CheckV0(itrack,event))       rejectTrack[itrack]+=4;
  }
  
  
  //
  // 2. Find secondary vertices - double loop
  //    
  for (Int_t itrack0=0; itrack0<ntracks; itrack0++){
    AliESDtrack *track0=event->GetTrack(itrack0);
    if (!track0) continue;
    if (track0->GetTPCNcls()<=kMinCl) continue;  // skip short tracks
    if ((1.+track0->GetTPCNcls())/(1.+track0->GetTPCNclsF())<=kMinRatio) continue;
    if (!track0->GetInnerParam()) continue;  // skip not TPC tracks
    if (track0->GetKinkIndex(0)>0) continue; // skip kinks
    if (rejectTrack[itrack0]) continue;   // skip
    track0->GetImpactParameters(dca[0],dca[1]);
    track0->GetImpactParameters(dca0[0],dca0[1]);
    if (TMath::Abs(dca[0])>kMaxDistR && TMath::Abs(dca[1])>kMaxDistZ) continue;
    // remove too dip secondaries

    AliKFParticle part0(*track0,211);  //assuming pion mass
    if (track0->Charge()*part0.Q()<0) part0.Q()*=-1;  // change sign if opposite
    //
    for (Int_t itrack1=itrack0+1; itrack1<ntracks; itrack1++){
      AliESDtrack *track1=event->GetTrack(itrack1);
      if (!track1) continue;
      if (rejectTrack[itrack1]) continue;   // skip
      if (track1->GetTPCNcls()<=kMinCl) continue;  // skip short tracks
      if ((1.+track1->GetTPCNcls())/(1.+track1->GetTPCNclsF())<=kMinRatio) continue;
      if (!track1->GetInnerParam()) continue;  // skip not TPC tracks
      if (track1->GetKinkIndex(0)!=0) continue; // skip kinks
      track1->GetImpactParameters(dca1[0],dca1[1]);
      track1->GetImpactParameters(dca[0],dca[1]);
      if (TMath::Abs(dca[0])<kMinDist && TMath::Abs(dca[1])<kMinDist) continue;
      if (TMath::Abs(dca[0])>kMaxDistR && TMath::Abs(dca[1])>kMaxDistZ) continue;     
      AliKFParticle part1(*track1,211); // assuming pion mass
      if (track1->Charge()*part1.Q()<0) part1.Q()*=-1;  // change sign if opposite

      //
      //
      AliKFVertex vertex;
      vertex+=part0;
      vertex+=part1;
      if ((vertex.GetChi2()/vertex.GetNDF())> kMaxChi2) continue;
      if (TMath::Abs(vertex.GetX())>kMaxDistR) continue;
      if (TMath::Abs(vertex.GetY())>kMaxDistR) continue;
      if (TMath::Abs(vertex.GetZ())>kMaxDistZ) continue;
      Double_t errX2=vertex.GetErrX();
      Double_t errY2=vertex.GetErrY();
      Double_t errZ2=vertex.GetErrZ();
      //
      Double_t err3D=TMath::Sqrt(errX2*errX2+errY2*errY2+errZ2*errZ2/10.);  
      Double_t err2D=TMath::Sqrt(errX2*errX2+errY2*errY2);  
      if (err3D>kMaxDistVertexSec) continue;
      if (err3D*TMath::Sqrt(vertex.GetChi2()+0.00001)>kMaxDistVertexSec) continue;

      Double_t dvertex=0;
      dvertex += (vertexSPD->GetX()-vertex.GetX())*(vertexSPD->GetX()-vertex.GetX());
      dvertex += (vertexSPD->GetY()-vertex.GetY())*(vertexSPD->GetY()-vertex.GetY());
      dvertex += (vertexSPD->GetZ()-vertex.GetZ())*(vertexSPD->GetZ()-vertex.GetZ());
      dvertex=TMath::Sqrt(dvertex+0.00000001);
      if (err3D>0.2*dvertex) continue;    
      if (err3D*TMath::Sqrt(vertex.GetChi2()+0.000001)>0.1*dvertex) continue;
      Double_t radius = TMath::Sqrt((vertex.GetX()*vertex.GetX()+vertex.GetY()*vertex.GetY()));
      //
      AliKFVertex vertex2;
      vertex2+=part0;
      vertex2+=part1;
      //

      for (Int_t itrack2=0; itrack2<ntracks; itrack2++){  
	if (itrack2==itrack0) continue;
	if (itrack2==itrack1) continue;
	if (rejectTrack[itrack2]) continue;   // skip
	AliESDtrack *track2=event->GetTrack(itrack2);
	if (!track2) continue;
	if (track2->GetTPCNcls()<=kMinCl) continue;  // skip short tracks
	if ((1.+track2->GetTPCNcls())/(1.+track2->GetTPCNclsF())<=kMinRatio) continue;
	if (!track2->GetInnerParam()) continue;  // skip not TPC tracks
	if (track2->GetKinkIndex(0)>0) continue; // skip kinks
	track2->GetImpactParameters(dca[0],dca[1]);
	if (TMath::Abs(dca[0])<kMinDist && TMath::Abs(dca[1])<kMinDist) continue;
	if (TMath::Abs(dca[0])>kMaxDistR && TMath::Abs(dca[1])>kMaxDistZ) continue;     
	if (TMath::Abs(track2->GetD(vertex.GetX(), vertex.GetY(),event->GetMagneticField()))>kMaxDistVertexSec) continue;
	Double_t vtxx[3]={vertex2.GetX(),vertex2.GetY(),vertex2.GetZ()};
	Double_t svtxx[3]={vertex.GetErrX(),vertex.GetErrY(),vertex.GetErrZ()};
	AliESDVertex vtx(vtxx,svtxx);
	AliExternalTrackParam param=*track2;
	Double_t delta[2]={0,0};
	if (!param.PropagateToDCA(&vtx,event->GetMagneticField(),kMaxDistVertexSec,delta)) continue;	
	if (TMath::Abs(delta[0])>kMaxDistVertexSec) continue;
	if (TMath::Abs(delta[1])>kMaxDistVertexSec) continue;      
	if (TMath::Abs(delta[0])>6.*TMath::Sqrt(param.GetSigmaY2()+vertex2.GetErrY()*vertex2.GetErrY())+0.1) continue; 
	if (TMath::Abs(delta[1])>6.*TMath::Sqrt(param.GetSigmaZ2()+vertex2.GetErrZ()*vertex2.GetErrZ())+0.5) continue; 
	//
	AliKFParticle part2(param,211); // assuming pion mass
	if (track2->Charge()*part2.Q()<0) part2.Q()*=-1;  // change sign if opposite
	vertex2+=part2;
	rejectTrack[itrack0]+=10;  // do noit reuse the track
	rejectTrack[itrack1]+=10;  // do not reuse the track      
	rejectTrack[itrack2]+=10;
      }
      

      TTreeSRedirector *pcstream = GetDebugStreamer();      
      if (pcstream){
	//
	//
	//
	Float_t dedx0= track0->GetTPCsignal(); 
	Float_t dedx1= track1->GetTPCsignal();
	AliExternalTrackParam * p0= (AliExternalTrackParam *)track0->GetInnerParam();
	AliExternalTrackParam * p1= (AliExternalTrackParam *)track1->GetInnerParam();
	Double_t errX=vertex2.GetErrX();
	Double_t errY=vertex2.GetErrY();
	Double_t errZ=vertex2.GetErrZ();
	Double_t vx = vertex2.GetX();
	Double_t vy = vertex2.GetY();
	Double_t vz = vertex2.GetZ();
	(*pcstream)<<"mapTPC"<<
	  "run="<<fRun<<                     // run
	  "event="<<fEvent<<          //  event number
	  "time="<<fTime<<         // timeStamp
	  "trigger="<<fTrigger<<      //  trigger
	  "mag="<<fMagF<<             //  magnetic field
	  //
	  "spdV.="<<spdVertex<<            // spd vertex
	  "trackV.="<<trackVertex<<        // track vertex
	  "tpcV.="<<tpcVertex<<            // track vertex
	  "tzero.="<<tzero<<               // tzero info
	  //
	  "ntracks="<<ntracks<<
	  "ntracksTPC="<<ntracksTPC<<
	  "nPrim="<<nTPCPrim<<              // number of primaries
	  "nSec="<<nTPCSec<<                // number of secondaries
	  "sigPrim="<<tpcSignalTotPrim<<    // total dEdx in primaries
	  "sigSec="<<tpcSignalTotSec<<      // total dEdx in secondaries
	  "dedx0="<<dedx0<<                 // dedx part 0
	  "dedx1="<<dedx1<<                 // dedx part 1
	  "p0.="<<p0<<                      // part 0
	  "p1.="<<p1<<                       //part 1
	  "v.="<<&vertex<<                  // KF vertex
	  "v2.="<<&vertex2<<                // KF vertex all tracks
	  "z0="<<dca0[1]<<
	  "z1="<<dca1[1]<<
	  "rphi0="<<dca0[0]<<
	  "rphi1="<<dca1[0]<<
	  "dvertex="<<dvertex<<
	  "radius="<<radius<<
	  "vx="<<vx<<
	  "vy="<<vy<<
	  "vz="<<vz<<
	  "errX="<<errX<<
	  "errY="<<errY<<
	  "errZ="<<errZ<<
	  "err2D="<<err2D<<
	  "err3D="<<err3D<<	  
	  "\n";
	//
	if (vertex2.GetNDF()>2){
	  (*pcstream)<<"mapVertex"<<
	    "run="<<fRun<<                     // run
	    "event="<<fEvent<<          //  event number
	    "time="<<fTime<<         // timeStamp
	    "trigger="<<fTrigger<<      //  trigger
	    "mag="<<fMagF<<             //  magnetic field
	    //
	    "spdV.="<<spdVertex<<            // spd vertex
	    "trackV.="<<trackVertex<<        // track vertex
	    "tpcV.="<<tpcVertex<<            // track vertex
	    "tzero.="<<tzero<<               // tzero info
	    //
	    //
	    "ntracks="<<ntracks<<
	    "ntracksTPC="<<ntracksTPC<<
	    "nPrim="<<nTPCPrim<<              // number of primaries
	    "nSec="<<nTPCSec<<                // number of secondaries
	    "sigPrim="<<tpcSignalTotPrim<<    // total dEdx in primaries
	    "sigSec="<<tpcSignalTotSec<<      // total dEdx in secondaries
	    "dedx0="<<dedx0<<                 // dedx part 0
	    "dedx1="<<dedx1<<                 // dedx part 1
	    "p0.="<<p0<<                      // part 0
	    "p1.="<<p1<<                       //part 1
	    "v.="<<&vertex<<                  // KF vertex
	    "v2.="<<&vertex2<<                // KF vertex all tracks
	    "z0="<<dca0[1]<<
	    "z1="<<dca1[1]<<
	    "rphi0="<<dca0[0]<<
	    "rphi1="<<dca1[0]<<
	    "radius="<<radius<<
	    "vx="<<vx<<
	    "vy="<<vy<<
	    "vz="<<vz<<
	    "errX="<<errX<<
	    "errY="<<errY<<
	    "errZ="<<errZ<<
	    "err2D="<<err2D<<
	    "err3D="<<err3D<<	  
	    "dvertex="<<dvertex<<
	    "\n";	  
	}
      }
    }
  }
  TTreeSRedirector *pcstream = GetDebugStreamer();
  if (pcstream){    
    (*pcstream)<<"statTPC"<<
      "run="<<fRun<<                     // run
      "time="<<fTime<<         // timeStamp
      "trigger="<<fTrigger<<      //  trigger
      "mag="<<fMagF<<             //  magnetic field
      "ntracks="<<ntracks<<
      "ntracksTPC="<<ntracksTPC<<
      //
      "nPrim="<<nTPCPrim<<              // number of primaries
      "nSec="<<nTPCSec<<                // number of secondaries
      "sigPrim="<<tpcSignalTotPrim<<    // total dEdx in primaries
      "sigSec="<<tpcSignalTotSec<<      // total dEdx in secondaries
      //
      "spdV.="<<spdVertex<<            // spd vertex
      "trackV.="<<trackVertex<<        // track vertex
      "tpcV.="<<tpcVertex<<            // track vertex
      "tzero.="<<tzero<<               // tzero info
      "\n";
  }
  delete [] rejectTrack;
}


Bool_t AliTPCcalibMaterial::CheckLooper(Int_t index, AliESDEvent *event){
  //
  // check if given track is looper candidate
  // if looper return kTRUE
  // 
  Int_t ntracks=event->GetNumberOfTracks();
  Int_t index1=-1;
  const Double_t ktglCut=0.03;
  const Double_t kalphaCut=0.4;
  static Int_t counter=0;
  //
  AliESDtrack * track0 = event->GetTrack(index);
  AliESDtrack * track1P = 0;
  for (Int_t itrack1=0; itrack1<ntracks; itrack1++){
    if (itrack1==index) continue;
    AliESDtrack *track1=event->GetTrack(itrack1);
    if (!track1) continue;
    if (TMath::Abs(TMath::Abs(track1->GetTgl())-TMath::Abs(track0->GetTgl()))>ktglCut) continue;
    if (TMath::Abs(TMath::Abs(track1->GetAlpha())-TMath::Abs(track0->GetAlpha()))>kalphaCut) continue;
    index1=index;
    track1P=track1;
  }
  if (index1>=0){
    TTreeSRedirector *pcstream = GetDebugStreamer();      
    if (pcstream &&counter<100000){
      counter++;
      AliExternalTrackParam p0(*track0);
      AliExternalTrackParam p1(*track1P);
      (*pcstream)<<"checkLooper"<<
	"p0.="<<&p0<<
	"p1.="<<&p1<<
	"\n";
    }
    return kTRUE;
  }
  return kFALSE;
}

Bool_t AliTPCcalibMaterial::CheckV0(Int_t index, AliESDEvent *event){
  //
  // check if given track is V0 candidata
  // if looper return kTRUE
  // 
  return kFALSE;
  Int_t ntracks=event->GetNumberOfTracks();
  Int_t index1=-1;
  const Double_t kSigmaMass=0.001;
  const Int_t kChi2Cut=10;
  static Int_t counter=0;
  //
  AliESDtrack * track0 = event->GetTrack(index);
  AliExternalTrackParam pL(*track0);
  AliKFParticle part0El(*track0, 11);  //assuming  mass e
  AliKFParticle part0Pi(*track0, 211);  //assuming  mass pi0
  AliKFParticle part0P(*track0, 2212);  //assuming  mass proton
  if (track0->Charge()*part0El.Q()<0) {
    part0El.Q()*=-1;  // change sign if opposite
    part0Pi.Q()*=-1;  // change sign if opposite
    part0P.Q()*=-1;   // change sign if opposite
  }
  Bool_t isGamma=0;
  Bool_t isK0=0;
  
  for (Int_t itrack1=0; itrack1<ntracks; itrack1++){
    if (itrack1==index) continue;
    AliESDtrack *track1=event->GetTrack(itrack1);
    if (!track1) continue;
    if (track1->Charge()*track0->Charge()>0) continue;
    AliKFParticle part1El(*track1, 11);  //assuming  mass e
    AliKFParticle part1Pi(*track1, 211);  //assuming  mass e
    AliKFParticle part1P(*track1, 2212);  //assuming  mass e
    if (track1->Charge()*part1El.Q()<0) {
      part1El.Q()*=-1;  // change sign if opposite
      part1Pi.Q()*=-1;  // change sign if opposite
      part1P.Q()*=-1;   // change sign if opposite
    }
    //
    AliKFVertex vertexG;  // gamma conversion candidate
    vertexG+=part0El;
    vertexG+=part1El;
    AliKFVertex vertexGC;  // gamma conversion candidate
    vertexGC+=part0El;
    vertexGC+=part1El;
    vertexGC.SetMassConstraint(0,kSigmaMass);
    AliKFVertex vertexK0;  // gamma conversion candidate
    vertexK0+=part0Pi;
    vertexK0+=part1Pi;
    AliKFVertex vertexK0C;  // gamma conversion candidate
    vertexK0C+=part0Pi;
    vertexK0C+=part1Pi;
    vertexK0C.SetMassConstraint(0.497614,kSigmaMass);
    if (vertexGC.GetChi2()<kChi2Cut && vertexG.GetMass()<0.06)      isGamma=kTRUE;
    if (vertexK0C.GetChi2()<kChi2Cut&&TMath::Abs(vertexK0.GetMass()-0.5)<0.06)  isK0=kTRUE;
    if (isGamma||isK0) {
      index1=index;
      TTreeSRedirector *pcstream = GetDebugStreamer();      
      if (pcstream&&counter<2000){
	counter++;
	AliExternalTrackParam p0(*track0);
	AliExternalTrackParam p1(*track1);
	(*pcstream)<<"checkV0"<<
	  "p0.="<<&p0<<  //particle 0
	  "p1.="<<&p1<<  //particle 1
	  "isGamma="<<isGamma<<  // is gamma candidate
	  "isK0s="<<isK0<<      // is k0 candidate
	  "vG.="<<&vertexG<<     // is gamma candidate
	  "vGC.="<<&vertexGC<<   // is gamma candidate
	  "vK.="<<&vertexK0<<     // is K0s candidate
	  "vKC.="<<&vertexK0C<<   // is K0s candidate
	  "\n";
      }
      break;
    }
  }
  if (index1>0) return kTRUE;
  return kFALSE;
}

/*
  //AliXRDPROOFtoolkit::FilterList("mater.list","* mapTPC",1);
  AliXRDPROOFtoolkit toolkit;
  TChain *chain   = toolkit.MakeChainRandom("mater.list.Good","mapVertex",0,4000);
  TChain *chainTPC   = toolkit.MakeChainRandom("mater.list.Good","mapTPC",0,50000);
  
 
  TCut cutErr="sqrt(v.fChi2)*err3D<1.0";
  TCut cutOccu="sqrt(v.fChi2*(errX^2+errY^2))/min(sqrt(v.fP[0]^2+v.fP[1]^2+v.fP[2]^2/20),20)<0.2&&sqrt((errX^2+errY^2))/min(sqrt(v.fP[0]^2+v.fP[1]^2+v.fP[2]^2/20),20)<0.3&&sqrt(v.fChi2)*err3D<1.5";
  //
  chainTPC->Draw(">>listTPC",cutOccu,"entryList");
  TEntryList *elistTPC = (TEntryList*)gDirectory->Get("listTPC");
  chainTPC->SetEntryList(elistTPC);

  //
  //
  
  chainTPC.Draw("v.fP[1]:v.fP[0]","sqrt(v.fP[0]^2+v.fP[1]^2)<100","",10000000);

  TCut cutITS="abs(v.fP[2])-10<sqrt(v.fP[0]^2+v.fP[1]^2)";
  chainTPC.Draw("v.fP[1]:v.fP[0]",cutITS+"(sqrt(v.fP[0]^2+v.fP[1]^2)<50)&&err2D<1&&err3D/dvertex<0.05","",500000);


  chainTPC.Draw("v.fP[1]:v.fP[0]>>his(400,-100,100,400,-100,100)","err2D<1&&err2D/dvertex<0.2","colz",50000); 


  chainTPC.Draw("radius:vz>>his(400,-100,100,400,0,100)","err2D<1&&err2D/dvertex<0.05","colz",5000); 

  TTree * tree = chain->CloneTree();
  TFile f("material.root","recreate");
  tree->Write();
  f.Close();
  


  TCut cutChi2="v.fChi2<2&&sqrt(errX^2+errY^2)<2."; 
  TFile f("material.root")
  TTree * tree = (TTree*)f.Get("mapVertex");

  tree->Draw("v.fChi2",cutChi2,"",1000)

  


*/

 
