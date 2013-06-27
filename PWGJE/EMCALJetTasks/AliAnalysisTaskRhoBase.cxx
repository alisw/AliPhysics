// $Id$
//
// Base class for rho calculation.
// Calculates parameterized rho for given centrality independent of input.
//
// Author: S.Aiola

#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TClonesArray.h>

#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliEmcalJet.h"

#include "AliAnalysisTaskRhoBase.h"

ClassImp(AliAnalysisTaskRhoBase)

//________________________________________________________________________
AliAnalysisTaskRhoBase::AliAnalysisTaskRhoBase() : 
  AliAnalysisTaskEmcalJet("AliAnalysisTaskRhoBase", kFALSE),
  fRhoScaledName(),
  fCompareRhoName(),
  fCompareRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fInEventSigmaRho(35.83),
  fAttachToEvent(kTRUE),
  fRhoScaled(0),
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
  fHistRhovsNtrack(0),
  fHistRhoScaledvsNtrack(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0)
{
  // Constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistJetNconstVsPt[i] = 0;
    fHistJetRhovsEta[i] = 0;
  }
  for (Int_t i = 0; i < 12; i++) {
    fHistNjUEoverNjVsNj[i] = 0;
  }
}

//________________________________________________________________________
AliAnalysisTaskRhoBase::AliAnalysisTaskRhoBase(const char *name, Bool_t histo) :
  AliAnalysisTaskEmcalJet(name, histo),
  fRhoScaledName(),
  fCompareRhoName(),
  fCompareRhoScaledName(),
  fRhoFunction(0),
  fScaleFunction(0),
  fInEventSigmaRho(35.83),
  fAttachToEvent(kTRUE),
  fRhoScaled(0),
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
  fHistRhovsNtrack(0),
  fHistRhoScaledvsNtrack(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistRhovsNcluster(0),
  fHistRhoScaledvsNcluster(0)
{
  // Constructor.

  for (Int_t i = 0; i < 4; i++) {
    fHistJetNconstVsPt[i] = 0;
    fHistJetRhovsEta[i] = 0;
  }
  for (Int_t i = 0; i < 12; i++) {
    fHistNjUEoverNjVsNj[i] = 0;
  }
  SetMakeGeneralHistograms(histo);
}

//________________________________________________________________________
void AliAnalysisTaskRhoBase::UserCreateOutputObjects()
{
  // User create output objects, called at the beginning of the analysis.

  if (!fCreateHisto)
    return;

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  fHistRhovsCent = new TH2F("fHistRhovsCent", "fHistRhovsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt*2);
  fHistRhovsCent->GetXaxis()->SetTitle("Centrality (%)");
  fHistRhovsCent->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
  fOutput->Add(fHistRhovsCent);

  if (!fTracksName.IsNull()) {
    fHistRhovsNtrack = new TH2F("fHistRhovsNtrack", "fHistRhovsNtrack", 150, 0, 6000, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistRhovsNtrack->GetXaxis()->SetTitle("No. of tracks");
    fHistRhovsNtrack->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhovsNtrack);
  }

  if (!fCaloName.IsNull()) {
    fHistRhovsNcluster = new TH2F("fHistRhovsNcluster", "fHistRhovsNcluster", 50, 0, 1500, fNbins, fMinBinPt, fMaxBinPt*2);
    fHistRhovsNcluster->GetXaxis()->SetTitle("No. of tracks");
    fHistRhovsNcluster->GetYaxis()->SetTitle("#rho (GeV/c * rad^{-1})");
    fOutput->Add(fHistRhovsNcluster);
  }

  if (!fJetsName.IsNull()) {
    fHistJetPtvsCent = new TH2F("fHistJetPtvsCent", "fHistJetPtvsCent", 101, -1,  100, fNbins, fMinBinPt, fMaxBinPt);
    fHistJetPtvsCent->GetXaxis()->SetTitle("Centrality (%)");
    fHistJetPtvsCent->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
    fOutput->Add(fHistJetPtvsCent);

    fHistJetAreavsCent = new TH2F("fHistJetAreavsCent", "fHistJetAreavsCent", 101, -1, 100, 30, 0, fJetRadius * fJetRadius * TMath::Pi() * 3);
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


    if (!fTracksName.IsNull()) {
      fHistJetPtvsNtrack = new TH2F("fHistJetPtvsNtrack", "fHistJetPtvsNtrack", 150, 0, 6000, fNbins, fMinBinPt, fMaxBinPt);
      fHistJetPtvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistJetPtvsNtrack->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/c)");
      fOutput->Add(fHistJetPtvsNtrack);

      fHistJetAreavsNtrack = new TH2F("fHistJetAreavsNtrack", "fHistJetAreavsNtrack", 150, 0, 6000, 30, 0, fJetRadius * fJetRadius * TMath::Pi() * 3);
      fHistJetAreavsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistJetAreavsNtrack->GetYaxis()->SetTitle("Jet area");
      fOutput->Add(fHistJetAreavsNtrack);

      fHistNjetvsNtrack = new TH2F("fHistNjetvsNtrack", "fHistNjetvsNtrack", 150, 0, 6000, 150, -0.5, 149.5);
      fHistJetAreavsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistJetAreavsNtrack->GetYaxis()->SetTitle("Jet area");
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

    if (!fTracksName.IsNull()) {
      fHistDeltaRhovsNtrack = new TH2F("fHistDeltaRhovsNtrack", "fHistDeltaRhovsNtrack", 150, 0, 6000, fNbins, -fMaxBinPt, fMaxBinPt);
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

    if (!fTracksName.IsNull()) {
      fHistRhoScaledvsNtrack = new TH2F("fHistRhoScaledvsNtrack", "fHistRhoScaledvsNtrack", 150, 0, 6000, fNbins, fMinBinPt, fMaxBinPt*2);
      fHistRhoScaledvsNtrack->GetXaxis()->SetTitle("No. of tracks");
      fHistRhoScaledvsNtrack->GetYaxis()->SetTitle("#rho_{scaled} (GeV/c * rad^{-1})");
      fOutput->Add(fHistRhoScaledvsNtrack);
    }

    if (!fCaloName.IsNull()) {
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
      
      if (!fTracksName.IsNull()) {
	fHistDeltaRhoScalevsNtrack = new TH2F("fHistDeltaRhoScalevsNtrack", "fHistDeltaRhoScalevsNtrack", 150, 0, 6000, fNbins, -fMaxBinPt, fMaxBinPt);
	fHistDeltaRhoScalevsNtrack->GetXaxis()->SetTitle("No. of tracks");
	fHistDeltaRhoScalevsNtrack->GetYaxis()->SetTitle("#Delta#rho_{scaled} (GeV/c * rad^{-1})");
	fOutput->Add(fHistDeltaRhoScalevsNtrack);
      }
    }
  }
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoBase::Run() 
{
  // Run the analysis.

  Double_t rho = GetRhoFactor(fCent);
  fRho->SetVal(rho);

  if (fScaleFunction) {
    Double_t rhoScaled = rho * GetScaleFactor(fCent);
    fRhoScaled->SetVal(rhoScaled);
  }

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskRhoBase::FillHistograms() 
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
    Int_t    NjetUE1Sigma  = 0;
    Int_t    NjetUE2Sigma  = 0;
    Int_t    NjetUE3Sigma  = 0;
    Double_t rhoPlus1Sigma = fRho->GetVal() + fInEventSigmaRho;
    Double_t rhoPlus2Sigma = fRho->GetVal() + 2*fInEventSigmaRho;
    Double_t rhoPlus3Sigma = fRho->GetVal() + 3*fInEventSigmaRho;

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
  
  fHistRhovsCent->Fill(fCent, fRho->GetVal());

  if (fTracks)
    fHistRhovsNtrack->Fill(Ntracks, fRho->GetVal());
  if (fCaloClusters)
    fHistRhovsNcluster->Fill(Nclusters, fRho->GetVal());
  if (fCompareRho) {
    fHistDeltaRhovsCent->Fill(fCent, fRho->GetVal() - fCompareRho->GetVal());
    if (fTracks)
      fHistDeltaRhovsNtrack->Fill(Ntracks, fRho->GetVal() - fCompareRho->GetVal());
  }

  if (fRhoScaled) {
    fHistRhoScaledvsCent->Fill(fCent, fRhoScaled->GetVal());
    if (fTracks)
      fHistRhoScaledvsNtrack->Fill(Ntracks, fRhoScaled->GetVal());
    if (fCaloClusters)
      fHistRhoScaledvsNcluster->Fill(Nclusters,  fRhoScaled->GetVal());
    if (fCompareRhoScaled) {
      fHistDeltaRhoScalevsCent->Fill(fCent, fRhoScaled->GetVal() - fCompareRhoScaled->GetVal());
      if (fTracks)
	fHistDeltaRhoScalevsNtrack->Fill(Ntracks, fRhoScaled->GetVal() - fCompareRhoScaled->GetVal());
    }
  }

  return kTRUE;
}      


//________________________________________________________________________
void AliAnalysisTaskRhoBase::ExecOnce() 
{
  // Init the analysis.

  if (!fRho) {
    fRho = new AliRhoParameter(fRhoName, 0);

    if (fAttachToEvent) {
      if (!(InputEvent()->FindListObject(fRhoName))) {
	InputEvent()->AddObject(fRho);
      } else {
	AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fRhoName.Data()));
	return;
      }
    }
    
    if (fScaleFunction && !fRhoScaled) {
      fRhoScaled = new AliRhoParameter(fRhoScaledName, 0);
      if (fAttachToEvent) {
	if (!(InputEvent()->FindListObject(fRhoScaledName))) {
	  InputEvent()->AddObject(fRhoScaled);
	} else {
	  AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fRhoScaledName.Data()));
	  return;
	}
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

//________________________________________________________________________
Double_t AliAnalysisTaskRhoBase::GetRhoFactor(Double_t cent)
{
  // Return rho per centrality.

  Double_t rho = 0;
  if (fRhoFunction)
    rho = fRhoFunction->Eval(cent);
  return rho;
}

//________________________________________________________________________
Double_t AliAnalysisTaskRhoBase::GetScaleFactor(Double_t cent)
{
  // Get scale factor.

  Double_t scale = 1;
  if (fScaleFunction)
    scale = fScaleFunction->Eval(cent);
  return scale;
}
