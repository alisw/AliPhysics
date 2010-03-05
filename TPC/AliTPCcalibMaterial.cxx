
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
    if(!calib) return 0;
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
