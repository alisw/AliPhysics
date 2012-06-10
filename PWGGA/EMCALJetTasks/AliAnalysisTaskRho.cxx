// $Id$
//
// Calculation of rho from a collection of jets.
// If scale function is given the scaled rho will be exported
// with the name as "fRhoName".apppend("_Scaled").
//
// Authors: R.Reed, S.Aiola

#include "AliAnalysisTaskRho.h"

#include <TClonesArray.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TList.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliEmcalJet.h"
#include "AliLog.h"
#include "AliRhoParameter.h"
#include "AliVCluster.h"
#include "AliVEventHandler.h"

ClassImp(AliAnalysisTaskRho)

//________________________________________________________________________
AliAnalysisTaskRho::AliAnalysisTaskRho() : 
  AliAnalysisTaskRhoBase(),
  fTracksName("tracks"),
  fJetsName("KtJets"),
  fRhoScaledName(""),
  fPhiMin(0),
  fPhiMax(0),
  fEtaMin(0),
  fEtaMax(0),
  fAreaCut(0),
  fNExclLeadJets(0),
  fScaleFunction(0),
  fCreateHisto(kFALSE),
  fTracks(0),
  fJets(0),
  fOutputList(0),
  fHistCentrality(0),
  fHistJetPt(0),
  fHistJetArea(0),
  fHistRhovsCent(0),
  fHistDeltaRhovsCent(0),
  fHistDeltaRhoScalevsCent(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistNjetvsCent(0),
  fHistRhovsNtrack(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fRhoScaled(0)
{
  // Constructor
}

//________________________________________________________________________
AliAnalysisTaskRho::AliAnalysisTaskRho(const char *name, Bool_t histo) :
  AliAnalysisTaskRhoBase(name),
  fTracksName("tracks"),
  fJetsName("KtJets"),
  fRhoScaledName(""),
  fPhiMin(0),
  fPhiMax(TMath::TwoPi()),
  fEtaMin(-0.5),
  fEtaMax(+0.5),
  fAreaCut(0.01),
  fNExclLeadJets(1),
  fScaleFunction(0),
  fCreateHisto(histo),
  fTracks(0),
  fJets(0),
  fOutputList(0),
  fHistCentrality(0),
  fHistJetPt(0),
  fHistJetArea(0),
  fHistRhovsCent(0),
  fHistDeltaRhovsCent(0),
  fHistDeltaRhoScalevsCent(0),
  fHistJetPtvsCent(0),
  fHistJetAreavsCent(0),
  fHistNjetvsCent(0),
  fHistRhovsNtrack(0),
  fHistDeltaRhovsNtrack(0),
  fHistDeltaRhoScalevsNtrack(0),
  fHistJetPtvsNtrack(0),
  fHistJetAreavsNtrack(0),
  fHistNjetvsNtrack(0),
  fRhoScaled(0)
{
  // Constructor

  if (fCreateHisto)
    DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AliAnalysisTaskRho::UserCreateOutputObjects()
{
  // User create output objects, called at the beginning of the analysis.

  AliAnalysisTaskRhoBase::UserCreateOutputObjects();

  if (fScaleFunction) {
    fRhoScaledName = fRhoName;
    fRhoScaledName += "_Scaled";
    fRhoScaled = new AliRhoParameter(fRhoScaledName, 0);  
  }

  if (!fCreateHisto)
    return;

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();

  fHistCentrality             = new TH1F("Centrality",            "Centrality",            101, -1,  100);
  fHistRhovsCent              = new TH2F("RhovsCent",             "RhovsCent",             101, -1,  100,   500,  0,   500);
  fHistDeltaRhovsCent         = new TH2F("DeltaRhovsCent",        "DetlaRhovsCent",        101, -1,  100,   500, -250, 250);
  fHistDeltaRhoScalevsCent    = new TH2F("DeltaRhoScalevsCent",   "DeltaRhoScalevsCent",   101, -1,  100,   500, -250, 250);
  fHistJetPtvsCent            = new TH2F("JetPtvsCent",           "JetPtvsCent",           101, -1,  100,   200,  0,   500);
  fHistJetAreavsCent          = new TH2F("JetAreavsCent",         "JetAreavsCent",         101, -1,  100,   100,  0,   1.0);
  fHistNjetvsCent             = new TH2F("NjetvsCent",            "NjetvsCent",            101, -1,  100,   100,  0,   100);

  fHistRhovsNtrack            = new TH2F("RhovsNtrack",           "RhovsNtrack",           500,  0,  2500,  500,  0,   500);
  fHistDeltaRhovsNtrack       = new TH2F("DeltaRhovsNtrack",      "DeltaRhovsNtrack",      500,  0,  2500,  500, -250, 250);
  fHistDeltaRhoScalevsNtrack  = new TH2F("DeltaRhoScalevsNtrack", "DeltaRhoScalevsNtrack", 500,  0,  2500,  500, -250, 250);
  fHistJetPtvsNtrack          = new TH2F("JetPtvsNtrack",         "JetPtvsNtrack",         500,  0,  2500,  200,  0,   500);
  fHistJetAreavsNtrack        = new TH2F("JetAreavsNtrack",       "JetAreavsNtrack",       500,  0,  2500,  100,  0,   1.0);
  fHistNjetvsNtrack           = new TH2F("NjetvsNtrack",          "rNjetvsNtrack",         500,  0,  2500,  100,  0,   100);

  fHistJetPt                  = new TH1F("JetPt",                  "Jet Pt",               100,   0,    250);
  fHistJetArea                = new TH1F("JetArea",                "Jet Area",             100, 0.0,    1.0);
  
  fOutputList->Add(fHistCentrality);
  fOutputList->Add(fHistRhovsCent);
  fOutputList->Add(fHistJetPt);
  fOutputList->Add(fHistJetArea);
  fOutputList->Add(fHistDeltaRhovsCent);
  fOutputList->Add(fHistDeltaRhoScalevsCent);
  fOutputList->Add(fHistJetPtvsCent);
  fOutputList->Add(fHistJetAreavsCent);
  fOutputList->Add(fHistNjetvsCent);
  
  fOutputList->Add(fHistRhovsNtrack);
  fOutputList->Add(fHistDeltaRhovsNtrack);
  fOutputList->Add(fHistDeltaRhoScalevsNtrack);
  fOutputList->Add(fHistJetPtvsNtrack);
  fOutputList->Add(fHistJetAreavsNtrack);
  fOutputList->Add(fHistNjetvsNtrack);
  
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskRho::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  
  if (!fIsInit) {
    ExecOnce();
    fIsInit = 1;
  }

  fRho->SetVal(-1);
  if (fRhoScaled)
    fRhoScaled->SetVal(-1);

  DetermineCent();

  const Int_t Njets   = fJets->GetEntries();

  Int_t maxJetIds[]   = {-1, -1};
  Float_t maxJetPts[] = { 0,  0};

  if (fNExclLeadJets > 0) {
    for (Int_t ij = 0; ij < Njets; ++ij) {
      AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(ij));
      if (!jet) {
	AliError(Form("%s: Could not receive jet %d", GetName(), ij));
	continue;
      } 
      if (jet->Pt() > maxJetPts[0]) {
	maxJetPts[1] = maxJetPts[0];
	maxJetIds[1] = maxJetIds[0];
	maxJetPts[0] = jet->Pt();
	maxJetIds[0] = ij;
      } else if (jet->Pt() > maxJetPts[1]) {
	maxJetPts[1] = jet->Pt();
	maxJetIds[1] = ij;
      }
    }
    if (fNExclLeadJets < 2) {
      maxJetIds[1] = -1;
      maxJetPts[1] = -1;
    }
  }

  static Double_t rhovec[999];
  Int_t NjetAcc = 0;

  // push all jets within selected acceptance into stack
  for (Int_t iJets = 0; iJets < Njets; ++iJets) {

    // exlcuding lead jets
    if (iJets == maxJetIds[0] || iJets == maxJetIds[1])
      continue;

    AliEmcalJet *jet = static_cast<AliEmcalJet*>(fJets->At(iJets));
    if (!jet)
      continue; 

    // applying some other cuts
    if (jet->Area() < fAreaCut)
      continue;
    if ((jet->Phi() < fPhiMin) || (jet->Phi() > fPhiMax))
      continue;
    if ((jet->Eta() < fEtaMin) || (jet->Eta() > fEtaMax))
      continue;

    rhovec[NjetAcc] = jet->Pt() / jet->Area();

    ++NjetAcc;

    if (fCreateHisto) {
      // filling histograms
      const Int_t Ntracks = fTracks->GetEntriesFast();
      fHistJetPt->Fill(jet->Pt());
      fHistJetArea->Fill(jet->Area());
      fHistJetPtvsCent->Fill(fCent, jet->Pt());
      fHistJetPtvsNtrack->Fill(Ntracks, jet->Pt());
      fHistJetAreavsCent->Fill(fCent, jet->Area());
      fHistJetAreavsNtrack->Fill(Ntracks, jet->Area());
    }
  }

  if (NjetAcc > 0) {
    //find median value
    Double_t rho0      = TMath::Median(NjetAcc, rhovec);
    fRho->SetVal(rho0);

    Double_t rhoScaled = rho0;
    if (fRhoScaled) {
      rhoScaled *= GetScaleFactor(fCent);
      fRhoScaled->SetVal(rhoScaled);
    }

    if (fCreateHisto) {
      // filling other histograms
      Double_t rho        = GetRhoFactor(fCent);
      const Int_t Ntracks = fTracks->GetEntriesFast();
      fHistRhovsCent->Fill(fCent, rho0);
      fHistDeltaRhovsCent->Fill(fCent, rho0 - rho);
      fHistDeltaRhoScalevsCent->Fill(fCent, rhoScaled - rho);
      fHistRhovsNtrack->Fill(Ntracks, rho0);
      fHistDeltaRhovsNtrack->Fill(Ntracks, rho0 - rho);
      fHistDeltaRhoScalevsNtrack->Fill(Ntracks, rhoScaled - rho);
    }
  }

  if (fCreateHisto) {
    const Int_t Ntracks = fTracks->GetEntriesFast();
    fHistCentrality->Fill(fCent);
    fHistNjetvsCent->Fill(fCent, NjetAcc);
    fHistNjetvsNtrack->Fill(Ntracks, NjetAcc);
    PostData(1, fOutputList);
  }
}      

//________________________________________________________________________
void AliAnalysisTaskRho::ExecOnce() 
{
  // Initialize some settings that need to be determined in UserExec.

  AliAnalysisTaskRhoBase::ExecOnce();

  if (fRhoScaled) {
    // add rho to event if not yet there
    if (!(InputEvent()->FindListObject(fRhoScaledName))) {
      InputEvent()->AddObject(fRhoScaled);
    } else {
      AliFatal(Form("%s: Container with same name %s already present. Aborting", GetName(), fRhoScaledName.Data()));
      return;
    }
  }

  if (fCreateHisto) {
    fTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksName));
    if (!fTracks) {
      AliError(Form("%s: Pointer to tracks %s == 0", GetName(), fTracksName.Data() ));
      return;
    }
  }

  fJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJetsName));
  if (!fJets) {
    AliError(Form("%s: Pointer to jets %s == 0", GetName(), fJetsName.Data() ));
    return;
  }
}

//________________________________________________________________________
Double_t AliAnalysisTaskRho::GetScaleFactor(Double_t cent)
{
  // Get scale factor.

  Double_t scale = 1;
  if (fScaleFunction)
    scale = fScaleFunction->Eval(cent);
  return scale;
}
