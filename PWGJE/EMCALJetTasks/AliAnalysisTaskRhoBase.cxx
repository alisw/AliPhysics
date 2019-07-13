/************************************************************************************
 * Copyright (C) 2012, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#include <TFile.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TClonesArray.h>
#include <TGrid.h>

#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliVVZERO.h"

#include "AliAnalysisTaskRhoBase.h"

ClassImp(AliAnalysisTaskRhoBase)

AliAnalysisTaskRhoBase::AliAnalysisTaskRhoBase() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskRhoBase", kFALSE),
  fOutRhoName(),
  fOutRhoScaledName(),
  fCompareRhoName(),
  fCompareRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fInEventSigmaRho(35.83),
  fAttachToEvent(kTRUE),
  fIsPbPb(kTRUE),
  fOutRho(0),
  fOutRhoScaled(0),
  fCompareRho(0),
  fCompareRhoScaled(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistJetRhovsCent(0),
  fHistNjetvsCent(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fHistRhovsCent(0),
  fHistRhoScaledvsCent(0),
  fHistDeltaRhovsCent(0),
  fHistDeltaRhoScalevsCent(0),
  fHistRhovsNtrackvsV0Mult(0),
  fHistRhoScaledvsNtrackvsV0Mult(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0)
{
  for (Int_t i = 0; i < 4; i++) {
    fHistJetNconstVsPt[i] = 0;
    fHistJetRhovsEta[i] = 0;
  }
  for (Int_t i = 0; i < 12; i++) {
    fHistNjUEoverNjVsNj[i] = 0;
  }
}

AliAnalysisTaskRhoBase::AliAnalysisTaskRhoBase(const char *name, Bool_t histo) :
  AliAnalysisTaskEmcalJet(name, histo),
  fOutRhoName(),
  fOutRhoScaledName(),
  fCompareRhoName(),
  fCompareRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fInEventSigmaRho(35.83),
  fAttachToEvent(kTRUE),
  fIsPbPb(kTRUE),
  fOutRho(0),
  fOutRhoScaled(0),
  fCompareRho(0),
  fCompareRhoScaled(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistJetRhovsCent(0),
  fHistNjetvsCent(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fHistRhovsCent(0),
  fHistRhoScaledvsCent(0),
  fHistDeltaRhovsCent(0),
  fHistDeltaRhoScalevsCent(0),
  fHistRhovsNtrackvsV0Mult(0),
  fHistRhoScaledvsNtrackvsV0Mult(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0)
{
  for (Int_t i = 0; i < 4; i++) {
    fHistJetNconstVsPt[i] = 0;
    fHistJetRhovsEta[i] = 0;
  }
  for (Int_t i = 0; i < 12; i++) {
    fHistNjUEoverNjVsNj[i] = 0;
  }
  SetMakeGeneralHistograms(histo);
}

void AliAnalysisTaskRhoBase::UserCreateOutputObjects()
{
  if (!fCreateHisto)
    return;

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  //ranges for PbPb
  Float_t Ntrackrange[2] = {0, 6000};
  Float_t V0Mult[2] = {0.,25000.};
  if(!fIsPbPb){
     //set multiplicity related axes to a smaller max value
     Ntrackrange[1] = 200.;
     V0Mult[1] = 2000.;
  }
  
  fHistRhovsCent = new TH2F("fHistRhovsCent", "fHistRhovsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt*2);
  fHistRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhovsCent->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
  fOutput->Add(fHistRhovsCent);

  if (fParticleCollArray.GetEntriesFast()>0) {
    fHistRhovsNtrackvsV0Mult = new TH3F("fHistRhovsNtrackvsV0Mult", "fHistRhovsNtrackvsV0Mult", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt*2,100, V0Mult[0], V0Mult[1]);
    fHistRhovsNtrackvsV0Mult->GetXaxis()->SetTitle("No. of tracks");
    fHistRhovsNtrackvsV0Mult->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
    fHistRhovsNtrackvsV0Mult->GetZaxis()->SetTitle("V0 mult");
    fOutput->Add(fHistRhovsNtrackvsV0Mult);
  }

  if (fClusterCollArray.GetEntriesFast()>0) {
    fHistRhovsNcluster = new TH2F("fHistRhovsNcluster", "fHistRhovsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistRhovsNcluster->GetXaxis()->SetTitle("No. of clusters");
    fHistRhovsNcluster->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhovsNcluster);
  }

  if (fJetCollArray.GetEntriesFast()>0) {
    fHistJetPtvsCent = new TH2F("fHistJetPtvsCent", "fHistJetPtvsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt);
    fHistJetPtvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetPtvsCent->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
    fOutput->Add(fHistJetPtvsCent);

    fHistJetAreavsCent = new TH2F("fHistJetAreavsCent", "fHistJetAreavsCent", 101, -1, 100, 100, 0, 1);
    fHistJetAreavsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetAreavsCent->GetYaxis()->SetTitle("Jet area");
    fOutput->Add(fHistJetAreavsCent);

    fHistJetRhovsCent = new TH2F("fHistJetRhovsCent", "fHistJetRhovsCent", 101, -1, 100, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistJetRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetRhovsCent->GetYaxis()->SetTitle("Jet #rho (GeV/c)");
    fOutput->Add(fHistJetRhovsCent);

    fHistNjetvsCent = new TH2F("fHistNjetvsCent",  "fHistNjetvsCent", 101, -1, 100, 150, -0.5, 149.5);
    fHistNjetvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistNjetvsCent->GetYaxis()->SetTitle("No. of jets");
    fOutput->Add(fHistNjetvsCent);


    if (fParticleCollArray.GetEntriesFast()>0) {
      fHistJetPtvsNtrack = new TH2F("fHistJetPtvsNtrack", "fHistJetPtvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt);
      fHistJetPtvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistJetPtvsNtrack->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
      fOutput->Add(fHistJetPtvsNtrack);

      fHistJetAreavsNtrack = new TH2F("fHistJetAreavsNtrack", "fHistJetAreavsNtrack", 150, Ntrackrange[0], Ntrackrange[1], 100, 0, 1);
      fHistJetAreavsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistJetAreavsNtrack->GetYaxis()->SetTitle("Jet area");
      fOutput->Add(fHistJetAreavsNtrack);

      fHistNjetvsNtrack = new TH2F("fHistNjetvsNtrack", "fHistNjetvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], 150, -0.5, 149.5);
      fHistNjetvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistNjetvsNtrack->GetYaxis()->SetTitle("No. of jets");
      fOutput->Add(fHistNjetvsNtrack);
    }


    TString name;
    for (Int_t i = 0; i < 4; i++) {
      name = Form("fHistJetNconstVsPt_%d",i);
      fHistJetNconstVsPt[i] = new TH2F(name, name, 150, -0.5, 149.5, fNbins, fMinBinPt, fMaxBinPt);
      fHistJetNconstVsPt[i]->GetXaxis()->SetTitle("No. of constituents");
      fHistJetNconstVsPt[i]->GetYaxis()->SetTitle("p_{T,jet} (GeV/c)");
      fOutput->Add(fHistJetNconstVsPt[i]);

      name = Form("fHistJetRhovsEta_%d",i);
      fHistJetRhovsEta[i] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt*2, 16, -0.8, 0.8);
      fHistJetRhovsEta[i]->GetXaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
      fHistJetRhovsEta[i]->GetYaxis()->SetTitle("#eta");
      fOutput->Add(fHistJetRhovsEta[i]);

      for (Int_t j = 0; j < 3; j++) {
        name = Form("NjUEoverNjVsNj_%d_%d",i,j+1);
        fHistNjUEoverNjVsNj[i*3+j] = new TH2F(name, name, 150, -0.5, 149.5, 120, 0.01, 1.21);
        fHistNjUEoverNjVsNj[i*3+j]->GetXaxis()->SetTitle("N_{jet}");
        fHistNjUEoverNjVsNj[i*3+j]->GetYaxis()->SetTitle("N_{jet_{UE}} / N_{jet}");
        fOutput->Add(fHistNjUEoverNjVsNj[i*3+j]);
      }
    }
  }
  
  if (!fCompareRhoName.IsNull()) {
    fHistDeltaRhovsCent = new TH2F("fHistDeltaRhovsCent", "fHistDeltaRhovsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistDeltaRhovsCent->GetYaxis()->SetTitle("#Delta#rho (GeV/c * rad^{-1})");
    fOutput->Add(fHistDeltaRhovsCent);

    if (fParticleCollArray.GetEntriesFast()>0) {
      fHistDeltaRhovsNtrack = new TH2F("fHistDeltaRhovsNtrack", "fHistDeltaRhovsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, -fMaxBinPt, fMaxBinPt);
      fHistDeltaRhovsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistDeltaRhovsNtrack->GetYaxis()->SetTitle("#Delta#rho (GeV/c * rad^{-1})");
      fOutput->Add(fHistDeltaRhovsNtrack);
    }
  }

  if (fScaleFunction) {
    fHistRhoScaledvsCent = new TH2F("fHistRhoScaledvsCent", "fHistRhoScaledvsCent", 101, -1, 100, fNbins, fMinBinPt , fMaxBinPt*2);
    fHistRhoScaledvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoScaledvsCent->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhoScaledvsCent);

    if (fParticleCollArray.GetEntriesFast()>0) {
      fHistRhoScaledvsNtrackvsV0Mult = new TH3F("fHistRhoScaledvsNtrackvsV0Mult", "fHistRhoScaledvsNtrackvsV0Mult", 150, Ntrackrange[0], Ntrackrange[1], fNbins, fMinBinPt, fMaxBinPt*2,100, V0Mult[0], V0Mult[1]);
      fHistRhoScaledvsNtrackvsV0Mult->GetXaxis()->SetTitle("No. of tracks");
      fHistRhoScaledvsNtrackvsV0Mult->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
      fHistRhoScaledvsNtrackvsV0Mult->GetZaxis()->SetTitle("V0 mult");
      fOutput->Add(fHistRhoScaledvsNtrackvsV0Mult);
    }

    if (fClusterCollArray.GetEntriesFast()>0) {
      fHistRhoScaledvsNcluster = new TH2F("fHistRhoScaledvsNcluster", "fHistRhoScaledvsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt*2);
      fHistRhoScaledvsNcluster->GetXaxis()->SetTitle("No. of clusters");
      fHistRhoScaledvsNcluster->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistRhoScaledvsNcluster);
    }

    if (!fCompareRhoScaledName.IsNull()) {
      fHistDeltaRhoScalevsCent = new TH2F("fHistDeltaRhoScalevsCent", "fHistDeltaRhoScalevsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
      fHistDeltaRhoScalevsCent->GetXaxis()->SetTitle("Centrality (%)");
      fHistDeltaRhoScalevsCent->GetYaxis()->SetTitle("#Delta#rho_{scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistDeltaRhoScalevsCent);
      
      if (fParticleCollArray.GetEntriesFast()>0) {
        fHistDeltaRhoScalevsNtrack = new TH2F("fHistDeltaRhoScalevsNtrack", "fHistDeltaRhoScalevsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, -fMaxBinPt, fMaxBinPt);
        fHistDeltaRhoScalevsNtrack->GetXaxis()->SetTitle("No. of tracks");
        fHistDeltaRhoScalevsNtrack->GetYaxis()->SetTitle("#Delta#rho_{scaled} (GeV/c * rad^{-1})");
        fOutput->Add(fHistDeltaRhoScalevsNtrack);
      }
    }
  }
}

Bool_t AliAnalysisTaskRhoBase::Run() 
{
  Double_t rho = GetRhoFactor(fCent);
  fOutRho->SetVal(rho);

  if (fScaleFunction) {
    Double_t rhoScaled = rho * GetScaleFactor(fCent);
    fOutRhoScaled->SetVal(rhoScaled);
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskRhoBase::FillHistograms() 
{
  Int_t Ntracks   = 0;
  Int_t Nclusters = 0;

  AliVVZERO* vV0 = InputEvent()->GetVZEROData();
  Float_t multV0A = -1.;
  Float_t multV0C = -1.;
  if(vV0) {  
    multV0A = vV0->GetMTotV0A();
    multV0C = vV0->GetMTotV0C();
  }

  // Loop over all possible contianers
  AliParticleContainer * partCont = 0;
  TIter nextPartCont(&fParticleCollArray);
  while ((partCont = static_cast<AliParticleContainer*>(nextPartCont()))) {
    Ntracks += partCont->GetNAcceptedParticles();
  }
  AliClusterContainer * clusCont = 0;
  TIter nextClusCont(&fClusterCollArray);
  while ((clusCont = static_cast<AliClusterContainer*>(nextClusCont()))) {
    Nclusters += clusCont->GetNAcceptedClusters();
  }

  if (fJets) {
    Int_t    Njets         = fJets->GetEntries();
    Int_t    NjetAcc       = 0;
    Int_t    NjetUE1Sigma  = 0;
    Int_t    NjetUE2Sigma  = 0;
    Int_t    NjetUE3Sigma  = 0;
    Double_t rhoPlus1Sigma = fOutRho->GetVal() + fInEventSigmaRho;
    Double_t rhoPlus2Sigma = fOutRho->GetVal() + 2*fInEventSigmaRho;
    Double_t rhoPlus3Sigma = fOutRho->GetVal() + 3*fInEventSigmaRho;

    for (Int_t i = 0; i < Njets; ++i) {

      AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(i));
      if (!jet) {
        AliError(Form("%s: Could not receive jet %d", GetName(), i));
        continue;
      } 

      if (!AcceptJet(jet))
        continue;

      fHistJetPtvsCent->Fill(fCent, jet->Pt());
      fHistJetAreavsCent->Fill(fCent, jet->Area());
      fHistJetRhovsCent->Fill(fCent, jet->Pt() / jet->Area());
      fHistJetRhovsEta[fCentBin]->Fill(jet->Pt() / jet->Area(), jet->Eta());

      if (fTracks) {
        fHistJetPtvsNtrack->Fill(Ntracks, jet->Pt());
        fHistJetAreavsNtrack->Fill(Ntracks, jet->Area());
      }

      fHistJetNconstVsPt[fCentBin]->Fill(jet->GetNumberOfConstituents(), jet->Pt());

      if (jet->Pt() < rhoPlus1Sigma * jet->Area())
        NjetUE1Sigma++;

      if (jet->Pt() < rhoPlus2Sigma * jet->Area())
        NjetUE2Sigma++;

      if (jet->Pt() < rhoPlus3Sigma * jet->Area())
        NjetUE3Sigma++;

      NjetAcc++;
    }

    if (NjetAcc>0) {
      fHistNjUEoverNjVsNj[fCentBin*3  ]->Fill(NjetAcc,1.*NjetUE1Sigma/NjetAcc);
      fHistNjUEoverNjVsNj[fCentBin*3+1]->Fill(NjetAcc,1.*NjetUE2Sigma/NjetAcc);
      fHistNjUEoverNjVsNj[fCentBin*3+2]->Fill(NjetAcc,1.*NjetUE3Sigma/NjetAcc);
    }

    fHistNjetvsCent->Fill(fCent, NjetAcc);
    if (fTracks)
      fHistNjetvsNtrack->Fill(Ntracks, NjetAcc);
  }

  fHistRhovsCent->Fill(fCent, fOutRho->GetVal());

  if (fTracks)
    fHistRhovsNtrackvsV0Mult->Fill(Ntracks, fOutRho->GetVal(),multV0A+multV0C);
  if (fCaloClusters)
    fHistRhovsNcluster->Fill(Nclusters, fOutRho->GetVal());
  if (fCompareRho) {
    fHistDeltaRhovsCent->Fill(fCent, fOutRho->GetVal() - fCompareRho->GetVal());
    if (fTracks)
      fHistDeltaRhovsNtrack->Fill(Ntracks, fOutRho->GetVal() - fCompareRho->GetVal());
  }

  if (fOutRhoScaled) {
    fHistRhoScaledvsCent->Fill(fCent, fOutRhoScaled->GetVal());
    if (fTracks)
      fHistRhoScaledvsNtrackvsV0Mult->Fill(Ntracks, fOutRhoScaled->GetVal(),multV0A+multV0C);
    if (fCaloClusters)
      fHistRhoScaledvsNcluster->Fill(Nclusters,  fOutRhoScaled->GetVal());
    if (fCompareRhoScaled) {
      fHistDeltaRhoScalevsCent->Fill(fCent, fOutRhoScaled->GetVal() - fCompareRhoScaled->GetVal());
      if (fTracks)
        fHistDeltaRhoScalevsNtrack->Fill(Ntracks, fOutRhoScaled->GetVal() - fCompareRhoScaled->GetVal());
    }
  }

  return kTRUE;
}      


void AliAnalysisTaskRhoBase::ExecOnce() 
{
  if (!fOutRho) {
    fOutRho = new AliRhoParameter(fOutRhoName, 0);

    if (fAttachToEvent) {
      if (!(InputEvent()->FindListObject(fOutRhoName))) {
        InputEvent()->AddObject(fOutRho);
      } else {
        AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fOutRhoName.Data()));
        return;
      }
    }
  }

  if (fScaleFunction && !fOutRhoScaled) {
    fOutRhoScaled = new AliRhoParameter(fOutRhoScaledName, 0);

    if (fAttachToEvent) {
      if (!(InputEvent()->FindListObject(fOutRhoScaledName))) {
        InputEvent()->AddObject(fOutRhoScaled);
      } else {
        AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fOutRhoScaledName.Data()));
        return;
      }
    }
  }

  if (!fCompareRhoName.IsNull() && !fCompareRho) {
    fCompareRho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fCompareRhoName));
    if (!fCompareRho) {
      AliWarning(Form("%s: Could not retrieve rho %s!", GetName(), fCompareRhoName.Data()));
    }
  }

  if (!fCompareRhoScaledName.IsNull() && !fCompareRhoScaled) {
    fCompareRhoScaled = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fCompareRhoScaledName));
    if (!fCompareRhoScaled) {
      AliWarning(Form("%s: Could not retrieve rho %s!", GetName(), fCompareRhoScaledName.Data()));
    }
  }

  AliAnalysisTaskEmcalJet::ExecOnce();
}

Double_t AliAnalysisTaskRhoBase::GetRhoFactor(Double_t cent)
{
  Double_t rho = 0;
  if (fRhoFunction)
    rho = fRhoFunction->Eval(cent);
  return rho;
}

Double_t AliAnalysisTaskRhoBase::GetScaleFactor(Double_t cent)
{
  Double_t scale = 1;
  if (fScaleFunction)
    scale = fScaleFunction->Eval(cent);
  return scale;
}

TF1* AliAnalysisTaskRhoBase::LoadRhoFunction(const char* path, const char* name)
{
  TString fname(path);
  if (fname.BeginsWith("alien://")) {
    TGrid::Connect("alien://");
  }

  TFile* file = TFile::Open(path);

  if (!file || file->IsZombie()) {
    ::Error("AddTaskRho", "Could not open scale function file");
    return 0;
  }

  TF1* sfunc = dynamic_cast<TF1*>(file->Get(name));

  if (sfunc) {
    ::Info("AliAnalysisTaskRhoBase::LoadRhoFunction", "Scale function %s loaded from file %s.", name, path);
  }
  else {
    ::Error("AliAnalysisTaskRhoBase::LoadRhoFunction", "Scale function %s not found in file %s.", name, path);
    return 0;
  }

  fScaleFunction = static_cast<TF1*>(sfunc->Clone());

  file->Close();
  delete file;

  return fScaleFunction;
}
