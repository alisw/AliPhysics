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

//====================================================================================================================================================
//
//      Parametric generator of primary pions and kaons
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TPDGCode.h"

#include "AliConst.h"
#include "AliRun.h"
#include "AliGenEventHeader.h"
#include "TDatabasePDG.h"
#include "AliPDG.h"
#include "TFile.h"
#include "TROOT.h"
#include "AliGenParamPionsKaons.h"
#include "TVector3.h"

ClassImp(AliGenParamPionsKaons)

//====================================================================================================================================================

AliGenParamPionsKaons::AliGenParamPionsKaons():
  AliGenerator(), 
  fGeneratePion(kTRUE),
  fGenerateKaon(kTRUE),
  fPtVsRapidityPrimaryPosPions(0x0),
  fPtVsRapidityPrimaryNegPions(0x0),
  fPtVsRapidityPrimaryPosKaons(0x0),
  fPtVsRapidityPrimaryNegKaons(0x0),
  fHistPdgCode(0x0) {

  // Default constructor

}

//====================================================================================================================================================

AliGenParamPionsKaons::AliGenParamPionsKaons(Int_t nPart, Char_t *inputFile):
  AliGenerator(nPart),
  fGeneratePion(kTRUE),
  fGenerateKaon(kTRUE),
  fPtVsRapidityPrimaryPosPions(0x0),
  fPtVsRapidityPrimaryNegPions(0x0),
  fPtVsRapidityPrimaryPosKaons(0x0),
  fPtVsRapidityPrimaryNegKaons(0x0),
  fHistPdgCode(0x0) {

  // Standard constructor

  fName  = "ParamPionsKaons";
  fTitle = "Parametric pions and kaons generator";

  LoadInputHistos(inputFile);

}

//====================================================================================================================================================

void AliGenParamPionsKaons::Generate() {

  // Generate one trigger
  
  Double_t polar[3]= {0,0,0};

  Double_t origin[3];
  Double_t p[3];

  Double_t mass=0., pt=0., rap=0., mom=0., energy=0, mt=0., phi=0., time=0.;
  Int_t nt;
  Int_t pdgCode;
  Double_t theta = 0.;

  for (Int_t j=0; j<3; j++) origin[j] = fOrigin[j];
  time = fTimeOrigin;
  if (fVertexSmear==kPerEvent) {
    Vertex();
    for (Int_t j=0; j<3; j++) origin[j] = fVertex[j];
    time = fTime;
  }

  Int_t nPartGenerated = 0;

  while (nPartGenerated<fNpart) {

    pdgCode = TMath::Nint(fHistPdgCode->GetRandom());
    if (TMath::Abs(pdgCode)==321 && !fGenerateKaon) continue;
    if (TMath::Abs(pdgCode)==211 && !fGeneratePion) continue;

    switch (pdgCode) {
    case 211:
      fPtVsRapidityPrimaryPosPions->GetRandom2(pt, rap);
      break;
    case -211:
      fPtVsRapidityPrimaryNegPions->GetRandom2(pt, rap);
      break;
    case 321:
      fPtVsRapidityPrimaryPosKaons->GetRandom2(pt, rap);
      break;
    case -321:
      fPtVsRapidityPrimaryNegKaons->GetRandom2(pt, rap);
      break;
    }

    mass = TDatabasePDG::Instance()->GetParticle(pdgCode)->Mass();

    mt = TMath::Sqrt(mass*mass + pt*pt);
    energy = mt * TMath::CosH(rap);
    mom = TMath::Sqrt(energy*energy - mass*mass);
    
    if (TestBit(kYRange))        if (rap<fYMin || rap>fYMax) continue;
    if (TestBit(kMomentumRange)) if (mom<fPMin || rap>fPMax) continue;
    if (TestBit(kPtRange))       if (pt<fPtMin || pt>fPtMax) continue;

    phi = fPhiMin + gRandom->Rndm()*(fPhiMax-fPhiMin);
    p[0] = pt*TMath::Cos(phi);
    p[1] = pt*TMath::Sin(phi);
    p[2] = mt*TMath::SinH(rap);
    
    TVector3 pv = TVector3(p);
    theta = pv.Theta();
    
    if (TestBit(kThetaRange)) if (theta<fThetaMin || theta>fThetaMax) continue;

    PushTrack(1, -1, Int_t(pdgCode), 
	      p[0],p[1],p[2],energy, 
	      origin[0],origin[1],origin[2],Double_t(time), 
	      polar[0],polar[1],polar[2], 
	      kPPrimary, nt, 1., 1);

    nPartGenerated++;

  }
  
  AliGenEventHeader* header = new AliGenEventHeader("ParamPionsKaons");
  header->SetPrimaryVertex(fVertex);
  header->SetNProduced(fNpart);
  header->SetInteractionTime(fTime);
  
  // Passes header either to the container or to gAlice
  if (fContainer) {
    fContainer->AddHeader(header);
  } 
  else {
    gAlice->SetGenEventHeader(header);	
  }

}

//====================================================================================================================================================

void AliGenParamPionsKaons::Init() {
  
  // Initialisation, check consistency of selected ranges
  if (TestBit(kPtRange) && TestBit(kMomentumRange)) 
    Fatal("Init","You should not set the momentum range and the pt range at the same time!\n");
  if ((!TestBit(kPtRange)) && (!TestBit(kMomentumRange))) 
    Fatal("Init","You should set either the momentum or the pt range!\n");
  if ((TestBit(kYRange) && TestBit(kThetaRange)) || (TestBit(kYRange) && TestBit(kEtaRange)) || (TestBit(kEtaRange) && TestBit(kThetaRange)))
    Fatal("Init","You should only set the range of one of these variables: y, eta or theta\n");
  if ((!TestBit(kYRange)) && (!TestBit(kEtaRange)) && (!TestBit(kThetaRange)))
    Fatal("Init","You should set the range of one of these variables: y, eta or theta\n");
  
  AliPDG::AddParticlesToPdgDataBase();
  
}

//====================================================================================================================================================

void AliGenParamPionsKaons::LoadInputHistos(Char_t *inputFile) {

  TFile *fileIn = new TFile(inputFile);

  TH2D *myPtVsRapidityPrimaryPosPions = (TH2D*) fileIn->Get("fPtVsRapidityPrimaryPosPions");
  TH2D *myPtVsRapidityPrimaryNegPions = (TH2D*) fileIn->Get("fPtVsRapidityPrimaryNegPions");
  TH2D *myPtVsRapidityPrimaryPosKaons = (TH2D*) fileIn->Get("fPtVsRapidityPrimaryPosKaons");
  TH2D *myPtVsRapidityPrimaryNegKaons = (TH2D*) fileIn->Get("fPtVsRapidityPrimaryNegKaons");
  TH1D *myHistPdgCode                 = (TH1D*) fileIn->Get("fHistPdgCode");

  myPtVsRapidityPrimaryPosPions -> SetName("myPtVsRapidityPrimaryPosPions");
  myPtVsRapidityPrimaryNegPions -> SetName("myPtVsRapidityPrimaryNegPions");
  myPtVsRapidityPrimaryPosKaons -> SetName("myPtVsRapidityPrimaryPosKaons");
  myPtVsRapidityPrimaryNegKaons -> SetName("myPtVsRapidityPrimaryNegKaons");
  myHistPdgCode                 -> SetName("myHistPdgCode");

  fPtVsRapidityPrimaryPosPions = new TH2D("fPtVsRapidityPrimaryPosPions","",
					  myPtVsRapidityPrimaryPosPions->GetXaxis()->GetNbins(),
					  myPtVsRapidityPrimaryPosPions->GetXaxis()->GetXmin(),
					  myPtVsRapidityPrimaryPosPions->GetXaxis()->GetXmax(),
					  myPtVsRapidityPrimaryPosPions->GetYaxis()->GetNbins(),
					  myPtVsRapidityPrimaryPosPions->GetYaxis()->GetXmin(),
					  myPtVsRapidityPrimaryPosPions->GetYaxis()->GetXmax());
  for (Int_t iBinX=0; iBinX<myPtVsRapidityPrimaryPosPions->GetXaxis()->GetNbins(); iBinX++) {
    for (Int_t iBinY=0; iBinY<myPtVsRapidityPrimaryPosPions->GetYaxis()->GetNbins(); iBinY++) {
      fPtVsRapidityPrimaryPosPions->SetBinContent(iBinX+1, iBinY+1, myPtVsRapidityPrimaryPosPions->GetBinContent(iBinX+1,iBinY+1));
    }
  }
					  
  fPtVsRapidityPrimaryNegPions = new TH2D("fPtVsRapidityPrimaryNegPions","",
					  myPtVsRapidityPrimaryNegPions->GetXaxis()->GetNbins(),
					  myPtVsRapidityPrimaryNegPions->GetXaxis()->GetXmin(),
					  myPtVsRapidityPrimaryNegPions->GetXaxis()->GetXmax(),
					  myPtVsRapidityPrimaryNegPions->GetYaxis()->GetNbins(),
					  myPtVsRapidityPrimaryNegPions->GetYaxis()->GetXmin(),
					  myPtVsRapidityPrimaryNegPions->GetYaxis()->GetXmax());
  for (Int_t iBinX=0; iBinX<myPtVsRapidityPrimaryNegPions->GetXaxis()->GetNbins(); iBinX++) {
    for (Int_t iBinY=0; iBinY<myPtVsRapidityPrimaryNegPions->GetYaxis()->GetNbins(); iBinY++) {
      fPtVsRapidityPrimaryNegPions->SetBinContent(iBinX+1, iBinY+1, myPtVsRapidityPrimaryNegPions->GetBinContent(iBinX+1,iBinY+1));
    }
  }
					  
  fPtVsRapidityPrimaryPosKaons = new TH2D("fPtVsRapidityPrimaryPosKaons","",
					  myPtVsRapidityPrimaryPosKaons->GetXaxis()->GetNbins(),
					  myPtVsRapidityPrimaryPosKaons->GetXaxis()->GetXmin(),
					  myPtVsRapidityPrimaryPosKaons->GetXaxis()->GetXmax(),
					  myPtVsRapidityPrimaryPosKaons->GetYaxis()->GetNbins(),
					  myPtVsRapidityPrimaryPosKaons->GetYaxis()->GetXmin(),
					  myPtVsRapidityPrimaryPosKaons->GetYaxis()->GetXmax());
  for (Int_t iBinX=0; iBinX<myPtVsRapidityPrimaryPosKaons->GetXaxis()->GetNbins(); iBinX++) {
    for (Int_t iBinY=0; iBinY<myPtVsRapidityPrimaryPosKaons->GetYaxis()->GetNbins(); iBinY++) {
      fPtVsRapidityPrimaryPosKaons->SetBinContent(iBinX+1, iBinY+1, myPtVsRapidityPrimaryPosKaons->GetBinContent(iBinX+1,iBinY+1));
    }
  }
					  
  fPtVsRapidityPrimaryNegKaons = new TH2D("fPtVsRapidityPrimaryNegKaons","",
					  myPtVsRapidityPrimaryNegKaons->GetXaxis()->GetNbins(),
					  myPtVsRapidityPrimaryNegKaons->GetXaxis()->GetXmin(),
					  myPtVsRapidityPrimaryNegKaons->GetXaxis()->GetXmax(),
					  myPtVsRapidityPrimaryNegKaons->GetYaxis()->GetNbins(),
					  myPtVsRapidityPrimaryNegKaons->GetYaxis()->GetXmin(),
					  myPtVsRapidityPrimaryNegKaons->GetYaxis()->GetXmax());
  for (Int_t iBinX=0; iBinX<myPtVsRapidityPrimaryNegKaons->GetXaxis()->GetNbins(); iBinX++) {
    for (Int_t iBinY=0; iBinY<myPtVsRapidityPrimaryNegKaons->GetYaxis()->GetNbins(); iBinY++) {
      fPtVsRapidityPrimaryNegKaons->SetBinContent(iBinX+1, iBinY+1, myPtVsRapidityPrimaryNegKaons->GetBinContent(iBinX+1,iBinY+1));
    }
  }
					  
  fHistPdgCode = new TH1D("fHistPdgCode","",
			  myHistPdgCode->GetXaxis()->GetNbins(),
			  myHistPdgCode->GetXaxis()->GetXmin(),
			  myHistPdgCode->GetXaxis()->GetXmax());
  for (Int_t iBinX=0; iBinX<myHistPdgCode->GetXaxis()->GetNbins(); iBinX++) {
    fHistPdgCode->SetBinContent(iBinX+1, myHistPdgCode->GetBinContent(iBinX+1));
  }

//   fHistPdgCode = new TH1D("fHistPdgCode", "fHistPdgCode", 10, 0., 10);
//   fHistPdgCode -> Fill(3.);

}

//====================================================================================================================================================
