/************************************************************************************
 * Copyright (C) 2014, Copyright Holders of the ALICE Collaboration                 *
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
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliEmcalJet.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"

#include "AliAnalysisTaskRhoMassBase.h"

ClassImp(AliAnalysisTaskRhoMassBase)

AliAnalysisTaskRhoMassBase::AliAnalysisTaskRhoMassBase() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskRhoMassBase", kFALSE),
  fOutRhoMassName(),
  fOutRhoMassScaledName(),
  fCompareRhoMassName(),
  fCompareRhoMassScaledName(),
  fRhoMassFunction(0),
  fScaleFunction(0),
  fAttachToEvent(kTRUE),
  fOutRhoMass(0),
  fOutRhoMassScaled(0),
  fCompareRhoMass(0),
  fCompareRhoMassScaled(0),
  fHistJetMassvsCent(0),
  fHistRhoMassvsCent(0),
  fHistRhoMassScaledvsCent(0),
  fHistDeltaRhoMassvsCent(0),
  fHistDeltaRhoMassScalevsCent(0),
  fHistRhoMassvsNtrack(0),
  fHistRhoMassScaledvsNtrack(0),
  fHistDeltaRhoMassvsNtrack(0),
  fHistDeltaRhoMassScalevsNtrack(0),
  fHistRhoMassvsNcluster(0),
  fHistRhoMassScaledvsNcluster(0),
  fHistGammaVsNtrack(0)
{
}

AliAnalysisTaskRhoMassBase::AliAnalysisTaskRhoMassBase(const char *name, Bool_t histo) :
  AliAnalysisTaskEmcalJet(name, histo),
  fOutRhoMassName(),
  fOutRhoMassScaledName(),
  fCompareRhoMassName(),
  fCompareRhoMassScaledName(),
  fRhoMassFunction(0),
  fScaleFunction(0),
  fAttachToEvent(kTRUE),
  fOutRhoMass(0),
  fOutRhoMassScaled(0),
  fCompareRhoMass(0),
  fCompareRhoMassScaled(0),
  fHistJetMassvsCent(0),
  fHistRhoMassvsCent(0),
  fHistRhoMassScaledvsCent(0),
  fHistDeltaRhoMassvsCent(0),
  fHistDeltaRhoMassScalevsCent(0),
  fHistRhoMassvsNtrack(0),
  fHistRhoMassScaledvsNtrack(0),
  fHistDeltaRhoMassvsNtrack(0),
  fHistDeltaRhoMassScalevsNtrack(0),
  fHistRhoMassvsNcluster(0),
  fHistRhoMassScaledvsNcluster(0),
  fHistGammaVsNtrack(0)
{
  SetMakeGeneralHistograms(histo);
}

void AliAnalysisTaskRhoMassBase::UserCreateOutputObjects()
{
  if (!fCreateHisto)
    return;

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  //ranges for PbPb
  Float_t Ntrackrange[2] = {0, 6000};

  Int_t nBinsRhom = 200;
  Double_t minRhom = 0.;
  Double_t maxRhom = 20.;
  
  Int_t nBinsM  = 100;
  Double_t minM = -20.;
  Double_t maxM = 80.;

  if(!fIsPbPb){
     //set multiplicity related axes to a smaller max value
     Ntrackrange[1] = 200.;
     maxRhom = 0.25;
  }
  
  fHistRhoMassvsCent = new TH2F("fHistRhoMassvsCent", "fHistRhoMassvsCent", 101, -1,  100, nBinsRhom,minRhom,maxRhom);
  fHistRhoMassvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhoMassvsCent->GetYaxis()->SetTitle("#rho_{m} (GeV/c * rad^{-1})");
  fOutput->Add(fHistRhoMassvsCent);

  if (fParticleCollArray.GetEntriesFast()>0) {
    fHistRhoMassvsNtrack = new TH2F("fHistRhoMassvsNtrack", "fHistRhoMassvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], nBinsRhom,minRhom,maxRhom);
    fHistRhoMassvsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistRhoMassvsNtrack->GetYaxis()->SetTitle("#rho_{m} (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhoMassvsNtrack);

    //fHistGammaVsNtrack
    fHistGammaVsNtrack = new TH2F("fHistGammaVsNtrack", "fHistGammaVsNtrack", 150, Ntrackrange[0], Ntrackrange[1], 100,0.,10.);
    fHistGammaVsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistGammaVsNtrack->GetYaxis()->SetTitle("#gamma = #LT E #GT / #LT M #GT");
    fOutput->Add(fHistGammaVsNtrack);
  }

  if (fClusterCollArray.GetEntriesFast()>0) {
    fHistRhoMassvsNcluster = new TH2F("fHistRhoMassvsNcluster", "fHistRhoMassvsNcluster", 50, 0, 1500, nBinsRhom,minRhom,maxRhom);
    fHistRhoMassvsNcluster->GetXaxis()->SetTitle("No. of tracks");
    fHistRhoMassvsNcluster->GetYaxis()->SetTitle("#rho_{m} (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhoMassvsNcluster);
  }

  if (fJetCollArray.GetEntriesFast()>0) {
    fHistJetMassvsCent = new TH2F("fHistJetMassvsCent", "fHistJetMassvsCent", 101, -1,  100, nBinsM,minM,maxM);
    fHistJetMassvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetMassvsCent->GetYaxis()->SetTitle("#it{M}_{jet} (GeV/c)");
    fOutput->Add(fHistJetMassvsCent);
  }
  
  if (!fCompareRhoMassName.IsNull()) {
    fHistDeltaRhoMassvsCent = new TH2F("fHistDeltaRhoMassvsCent", "fHistDeltaRhoMassvsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaRhoMassvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistDeltaRhoMassvsCent->GetYaxis()->SetTitle("#Delta#rho (GeV/c * rad^{-1})");
    fOutput->Add(fHistDeltaRhoMassvsCent);

    if (fParticleCollArray.GetEntriesFast()>0) {
      fHistDeltaRhoMassvsNtrack = new TH2F("fHistDeltaRhoMassvsNtrack", "fHistDeltaRhoMassvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, -fMaxBinPt, fMaxBinPt);
      fHistDeltaRhoMassvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistDeltaRhoMassvsNtrack->GetYaxis()->SetTitle("#Delta#rho (GeV/c * rad^{-1})");
      fOutput->Add(fHistDeltaRhoMassvsNtrack);
    }
  }

  if (fScaleFunction) {
    fHistRhoMassScaledvsCent = new TH2F("fHistRhoMassScaledvsCent", "fHistRhoMassScaledvsCent", 101, -1, 100, nBinsRhom,minRhom,maxRhom);
    fHistRhoMassScaledvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoMassScaledvsCent->GetYaxis()->SetTitle("#rho_{m,scaled} (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhoMassScaledvsCent);

    if (fParticleCollArray.GetEntriesFast()>0) {
      fHistRhoMassScaledvsNtrack = new TH2F("fHistRhoMassScaledvsNtrack", "fHistRhoMassScaledvsNtrack", 150, Ntrackrange[0], Ntrackrange[1], nBinsRhom,minRhom,maxRhom);
      fHistRhoMassScaledvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistRhoMassScaledvsNtrack->GetYaxis()->SetTitle("#rho_{m,scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistRhoMassScaledvsNtrack);
    }

    if (fClusterCollArray.GetEntriesFast()>0) {
      fHistRhoMassScaledvsNcluster = new TH2F("fHistRhoMassScaledvsNcluster", "fHistRhoMassScaledvsNcluster", 50, 0, 1500, nBinsRhom,minRhom,maxRhom);
      fHistRhoMassScaledvsNcluster->GetXaxis()->SetTitle("No. of clusters");
      fHistRhoMassScaledvsNcluster->GetYaxis()->SetTitle("#rho_{m,scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistRhoMassScaledvsNcluster);
    }

    if (!fCompareRhoMassScaledName.IsNull()) {
      fHistDeltaRhoMassScalevsCent = new TH2F("fHistDeltaRhoMassScalevsCent", "fHistDeltaRhoMassScalevsCent", 101, -1, 100, nBinsRhom,minRhom,maxRhom);
      fHistDeltaRhoMassScalevsCent->GetXaxis()->SetTitle("Centrality (%)");
      fHistDeltaRhoMassScalevsCent->GetYaxis()->SetTitle("#Delta#rho_{m,scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistDeltaRhoMassScalevsCent);
      
      if (fParticleCollArray.GetEntriesFast()>0) {
        fHistDeltaRhoMassScalevsNtrack = new TH2F("fHistDeltaRhoMassScalevsNtrack", "fHistDeltaRhoMassScalevsNtrack", 150, Ntrackrange[0], Ntrackrange[1], fNbins, -fMaxBinPt, fMaxBinPt);
        fHistDeltaRhoMassScalevsNtrack->GetXaxis()->SetTitle("No. of tracks");
        fHistDeltaRhoMassScalevsNtrack->GetYaxis()->SetTitle("#Delta#rho_{m,scaled} (GeV/c * rad^{-1})");
        fOutput->Add(fHistDeltaRhoMassScalevsNtrack);
      }
    }
  }
}

Bool_t AliAnalysisTaskRhoMassBase::Run() 
{
  Double_t rhom = GetRhoMassFactor(fCent);
  fOutRhoMass->SetVal(rhom);

  if (fScaleFunction) {
    Double_t rhomScaled = rhom * GetScaleFactor(fCent);
    fOutRhoMassScaled->SetVal(rhomScaled);
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskRhoMassBase::FillHistograms() 
{
  Int_t Ntracks   = 0;
  Int_t Nclusters = 0;

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

    for (Int_t i = 0; i < Njets; ++i) {
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(i));
      if (!jet) {
        AliError(Form("%s: Could not receive jet %d", GetName(), i));
        continue;
      } 

      if (!AcceptJet(jet))
        continue;

      fHistJetMassvsCent->Fill(fCent, jet->M());
      NjetAcc++;
    }
  }

  fHistRhoMassvsCent->Fill(fCent, fOutRhoMass->GetVal());

  if (fTracks)
    fHistRhoMassvsNtrack->Fill(Ntracks, fOutRhoMass->GetVal());
  if (fCaloClusters)
    fHistRhoMassvsNcluster->Fill(Nclusters, fOutRhoMass->GetVal());
  if (fCompareRhoMass) {
    fHistDeltaRhoMassvsCent->Fill(fCent, fOutRhoMass->GetVal() - fCompareRhoMass->GetVal());
    if (fTracks)
      fHistDeltaRhoMassvsNtrack->Fill(Ntracks, fOutRhoMass->GetVal() - fCompareRhoMass->GetVal());
  }

  if (fOutRhoMassScaled) {
    fHistRhoMassScaledvsCent->Fill(fCent, fOutRhoMassScaled->GetVal());
    if (fTracks)
      fHistRhoMassScaledvsNtrack->Fill(Ntracks, fOutRhoMassScaled->GetVal());
    if (fCaloClusters)
      fHistRhoMassScaledvsNcluster->Fill(Nclusters,  fOutRhoMassScaled->GetVal());
    if (fCompareRhoMassScaled) {
      fHistDeltaRhoMassScalevsCent->Fill(fCent, fOutRhoMassScaled->GetVal() - fCompareRhoMassScaled->GetVal());
      if (fTracks)
        fHistDeltaRhoMassScalevsNtrack->Fill(Ntracks, fOutRhoMassScaled->GetVal() - fCompareRhoMassScaled->GetVal());
    }
  }

  return kTRUE;
}      


void AliAnalysisTaskRhoMassBase::ExecOnce() 
{
  if (!fOutRhoMass) {
    fOutRhoMass = new AliRhoParameter(fOutRhoMassName, 0);

    if (fAttachToEvent) {
      if (!(InputEvent()->FindListObject(fOutRhoMassName))) {
        InputEvent()->AddObject(fOutRhoMass);
      } else {
        AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fOutRhoMassName.Data()));
        return;
      }
    }
  }

  if (fScaleFunction && !fOutRhoMassScaled) {
    fOutRhoMassScaled = new AliRhoParameter(fOutRhoMassScaledName, 0);

    if (fAttachToEvent) {
      if (!(InputEvent()->FindListObject(fOutRhoMassScaledName))) {
        InputEvent()->AddObject(fOutRhoMassScaled);
      } else {
        AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fOutRhoMassScaledName.Data()));
        return;
      }
    }
  }

  if (!fCompareRhoMassName.IsNull() && !fCompareRhoMass) {
    fCompareRhoMass = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fCompareRhoMassName));
    if (!fCompareRhoMass) {
      AliWarning(Form("%s: Could not retrieve rho %s!", GetName(), fCompareRhoMassName.Data()));
    }
  }

  if (!fCompareRhoMassScaledName.IsNull() && !fCompareRhoMassScaled) {
    fCompareRhoMassScaled = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fCompareRhoMassScaledName));
    if (!fCompareRhoMassScaled) {
      AliWarning(Form("%s: Could not retrieve rho %s!", GetName(), fCompareRhoMassScaledName.Data()));
    }
  }

  AliAnalysisTaskEmcalJet::ExecOnce();
}

Double_t AliAnalysisTaskRhoMassBase::GetRhoMassFactor(Double_t cent)
{
  Double_t rho = 0;
  if (fRhoMassFunction)
    rho = fRhoMassFunction->Eval(cent);
  return rho;
}

Double_t AliAnalysisTaskRhoMassBase::GetScaleFactor(Double_t cent)
{
  Double_t scale = 1;
  if (fScaleFunction)
    scale = fScaleFunction->Eval(cent);
  return scale;
}
