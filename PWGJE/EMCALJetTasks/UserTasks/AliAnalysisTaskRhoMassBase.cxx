// $Id$
//
// Base class for rho mass calculation.
// Calculates parameterized rho mass for given centrality independent of input.
//
// Author: M. Verweij. Similar to AliAnalysisTaskRhoBase

#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliEmcalJet.h"

#include "AliAnalysisTaskRhoMassBase.h"

ClassImp(AliAnalysisTaskRhoMassBase)

//________________________________________________________________________
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
  fHistJetRhoMassvsCent(0),
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
  // Constructor.

  for (Int_t i = 0; i < 4; i++)
    fHistJetRhoMassvsEta[i] = 0;
}

//________________________________________________________________________
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
  fHistJetRhoMassvsCent(0),
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
  // Constructor.

  for (Int_t i = 0; i < 4; i++)
    fHistJetRhoMassvsEta[i] = 0;

  SetMakeGeneralHistograms(histo);
}

//________________________________________________________________________
void AliAnalysisTaskRhoMassBase::UserCreateOutputObjects()
{
  // User create output objects, called at the beginning of the analysis.

  if (!fCreateHisto)
    return;

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fHistRhoMassvsCent = new TH2F("fHistRhoMassvsCent", "fHistRhoMassvsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt/2.);
  fHistRhoMassvsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhoMassvsCent->GetYaxis()->SetTitle("#rho_{m} (GeV/c * rad^{-1})");
  fOutput->Add(fHistRhoMassvsCent);

  if (fParticleCollArray.GetEntriesFast()>0) {
    fHistRhoMassvsNtrack = new TH2F("fHistRhoMassvsNtrack", "fHistRhoMassvsNtrack", 150, 0, 6000, fNbins, fMinBinPt, fMaxBinPt/2);
    fHistRhoMassvsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistRhoMassvsNtrack->GetYaxis()->SetTitle("#rho_{m} (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhoMassvsNtrack);

    //fHistGammaVsNtrack
    fHistGammaVsNtrack = new TH2F("fHistGammaVsNtrack", "fHistGammaVsNtrack", 150, 0, 6000, 100,0.,10.);
    fHistGammaVsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistGammaVsNtrack->GetYaxis()->SetTitle("#gamma = #LT E #GT / #LT M #GT");
    fOutput->Add(fHistGammaVsNtrack);
  }

  if (fClusterCollArray.GetEntriesFast()>0) {
    fHistRhoMassvsNcluster = new TH2F("fHistRhoMassvsNcluster", "fHistRhoMassvsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt/2);
    fHistRhoMassvsNcluster->GetXaxis()->SetTitle("No. of tracks");
    fHistRhoMassvsNcluster->GetYaxis()->SetTitle("#rho_{m} (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhoMassvsNcluster);
  }

  if (fJetCollArray.GetEntriesFast()>0) {
    fHistJetMassvsCent = new TH2F("fHistJetMassvsCent", "fHistJetMassvsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt);
    fHistJetMassvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetMassvsCent->GetYaxis()->SetTitle("#it{M}_{jet} (GeV/c)");
    fOutput->Add(fHistJetMassvsCent);

    fHistJetRhoMassvsCent = new TH2F("fHistJetRhoMassvsCent", "fHistJetRhoMassvsCent", 101, -1, 100, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistJetRhoMassvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetRhoMassvsCent->GetYaxis()->SetTitle("Jet #rho_{m} (GeV/c)");
    fOutput->Add(fHistJetRhoMassvsCent);

 
    TString name;
    for (Int_t i = 0; i < 4; i++) {
      name = Form("fHistJetRhoMassvsEta_%d",i);
      fHistJetRhoMassvsEta[i] = new TH2F(name, name, fNbins, fMinBinPt, fMaxBinPt*2, 16, -0.8, 0.8);
      fHistJetRhoMassvsEta[i]->GetXaxis()->SetTitle("#rho_{m} (GeV/c)");
      fHistJetRhoMassvsEta[i]->GetYaxis()->SetTitle("#eta");
      fOutput->Add(fHistJetRhoMassvsEta[i]);
    }
  }
  
  if (!fCompareRhoMassName.IsNull()) {
    fHistDeltaRhoMassvsCent = new TH2F("fHistDeltaRhoMassvsCent", "fHistDeltaRhoMassvsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaRhoMassvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistDeltaRhoMassvsCent->GetYaxis()->SetTitle("#Delta#rho (GeV/c * rad^{-1})");
    fOutput->Add(fHistDeltaRhoMassvsCent);

    if (fParticleCollArray.GetEntriesFast()>0) {
      fHistDeltaRhoMassvsNtrack = new TH2F("fHistDeltaRhoMassvsNtrack", "fHistDeltaRhoMassvsNtrack", 150, 0, 6000, fNbins, -fMaxBinPt, fMaxBinPt);
      fHistDeltaRhoMassvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistDeltaRhoMassvsNtrack->GetYaxis()->SetTitle("#Delta#rho (GeV/c * rad^{-1})");
      fOutput->Add(fHistDeltaRhoMassvsNtrack);
    }
  }

  if (fScaleFunction) {
    fHistRhoMassScaledvsCent = new TH2F("fHistRhoMassScaledvsCent", "fHistRhoMassScaledvsCent", 101, -1, 100, fNbins, fMinBinPt , fMaxBinPt*2);
    fHistRhoMassScaledvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistRhoMassScaledvsCent->GetYaxis()->SetTitle("#rho_{m,scaled} (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhoMassScaledvsCent);

    if (fParticleCollArray.GetEntriesFast()>0) {
      fHistRhoMassScaledvsNtrack = new TH2F("fHistRhoMassScaledvsNtrack", "fHistRhoMassScaledvsNtrack", 150, 0, 6000, fNbins, fMinBinPt, fMaxBinPt*2);
      fHistRhoMassScaledvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistRhoMassScaledvsNtrack->GetYaxis()->SetTitle("#rho_{m,scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistRhoMassScaledvsNtrack);
    }

    if (fClusterCollArray.GetEntriesFast()>0) {
      fHistRhoMassScaledvsNcluster = new TH2F("fHistRhoMassScaledvsNcluster", "fHistRhoMassScaledvsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt*2);
      fHistRhoMassScaledvsNcluster->GetXaxis()->SetTitle("No. of clusters");
      fHistRhoMassScaledvsNcluster->GetYaxis()->SetTitle("#rho_{m,scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistRhoMassScaledvsNcluster);
    }

    if (!fCompareRhoMassScaledName.IsNull()) {
      fHistDeltaRhoMassScalevsCent = new TH2F("fHistDeltaRhoMassScalevsCent", "fHistDeltaRhoMassScalevsCent", 101, -1, 100, fNbins, -fMaxBinPt, fMaxBinPt);
      fHistDeltaRhoMassScalevsCent->GetXaxis()->SetTitle("Centrality (%)");
      fHistDeltaRhoMassScalevsCent->GetYaxis()->SetTitle("#Delta#rho_{m,scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistDeltaRhoMassScalevsCent);
      
      if (fParticleCollArray.GetEntriesFast()>0) {
	fHistDeltaRhoMassScalevsNtrack = new TH2F("fHistDeltaRhoMassScalevsNtrack", "fHistDeltaRhoMassScalevsNtrack", 150, 0, 6000, fNbins, -fMaxBinPt, fMaxBinPt);
	fHistDeltaRhoMassScalevsNtrack->GetXaxis()->SetTitle("No. of tracks");
	fHistDeltaRhoMassScalevsNtrack->GetYaxis()->SetTitle("#Delta#rho_{m,scaled} (GeV/c * rad^{-1})");
	fOutput->Add(fHistDeltaRhoMassScalevsNtrack);
      }
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoMassBase::Run() 
{
  // Run the analysis.

  Double_t rhom = GetRhoMassFactor(fCent);
  fOutRhoMass->SetVal(rhom);

  if (fScaleFunction) {
    Double_t rhomScaled = rhom * GetScaleFactor(fCent);
    fOutRhoMassScaled->SetVal(rhomScaled);
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoMassBase::FillHistograms() 
{
  // Fill histograms.

  Int_t Ntracks   = 0;
  Int_t Nclusters = 0;

  if (fTracks)
    Ntracks = fTracks->GetEntries();
  if (fCaloClusters)
    Nclusters = fCaloClusters->GetEntries();

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
      fHistJetRhoMassvsCent->Fill(fCent, jet->M() / jet->Area());
      fHistJetRhoMassvsEta[fCentBin]->Fill(jet->M() / jet->Area(), jet->Eta());

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


//________________________________________________________________________
void AliAnalysisTaskRhoMassBase::ExecOnce() 
{
  // Init the analysis.

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

//________________________________________________________________________
Double_t AliAnalysisTaskRhoMassBase::GetRhoMassFactor(Double_t cent)
{
  // Return rho per centrality.

  Double_t rho = 0;
  if (fRhoMassFunction)
    rho = fRhoMassFunction->Eval(cent);
  return rho;
}

//________________________________________________________________________
Double_t AliAnalysisTaskRhoMassBase::GetScaleFactor(Double_t cent)
{
  // Get scale factor.

  Double_t scale = 1;
  if (fScaleFunction)
    scale = fScaleFunction->Eval(cent);
  return scale;
}
