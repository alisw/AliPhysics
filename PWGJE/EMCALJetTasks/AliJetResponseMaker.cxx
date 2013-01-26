// $Id$
//
// Emcal jet response matrix maker task.
//
// Author: S. Aiola

#include "AliJetResponseMaker.h"

#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TLorentzVector.h>

#include "AliVCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliGenPythiaEventHeader.h"
#include "AliLog.h"
#include "AliMCEvent.h"

ClassImp(AliJetResponseMaker)

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker() : 
  AliAnalysisTaskEmcalJet("AliJetResponseMaker", kTRUE),
  fMCTracksName("MCParticles"),
  fMCJetsName("MCJets"),
  fMaxDistance(0.25),
  fDoWeighting(kFALSE),
  fEventWeightHist(kFALSE),
  fMCFiducial(kFALSE),
  fMCminEta(0),
  fMCmaxEta(0),
  fMCminPhi(0),
  fMCmaxPhi(0),
  fSelectPtHardBin(-999),
  fDoMatching(kTRUE),
  fPythiaHeader(0),
  fEventWeight(0),
  fPtHardBin(0),
  fNTrials(0),
  fMCTracks(0),
  fMCJets(0),
  fHistNTrials(0),
  fHistEvents(0),
  fHistMCJetsPhiEta(0),
  fHistMCJetsPtArea(0),
  fHistMCJetsPhiEtaFiducial(0),
  fHistMCJetsPtAreaFiducial(0),
  fHistMCJetsNEFvsPt(0),
  fHistMCJetsZvsPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsCorrPtArea(0),
  fHistJetsNEFvsPt(0),
  fHistJetsZvsPt(0),
  fHistMatchingLevelMCPt(0),
  fHistClosestDeltaEtaPhiMCPt(0),
  fHistClosestDeltaPtMCPt(0),
  fHistClosestDeltaCorrPtMCPt(0),
  fHistNonMatchedMCJetsPtArea(0),
  fHistNonMatchedJetsPtArea(0),
  fHistNonMatchedJetsCorrPtArea(0),
  fHistPartvsDetecPt(0),
  fHistPartvsDetecCorrPt(0),
  fHistMissedMCJetsPtArea(0)
{
  // Default constructor.

  for (Int_t i = 0; i < 11; i++) {
    fHistEventWeight[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fMCTracksName("MCParticles"),
  fMCJetsName("MCJets"),
  fMaxDistance(0.25),
  fDoWeighting(kFALSE),
  fEventWeightHist(kFALSE),
  fMCFiducial(kFALSE),
  fMCminEta(0),
  fMCmaxEta(0),
  fMCminPhi(0),
  fMCmaxPhi(0),
  fSelectPtHardBin(-999),
  fDoMatching(kTRUE),
  fPythiaHeader(0),
  fEventWeight(0),
  fPtHardBin(0),
  fNTrials(0),
  fMCTracks(0),
  fMCJets(0),
  fHistNTrials(0),
  fHistEvents(0),
  fHistMCJetsPhiEta(0),
  fHistMCJetsPtArea(0),
  fHistMCJetsPhiEtaFiducial(0),
  fHistMCJetsPtAreaFiducial(0),
  fHistMCJetsNEFvsPt(0),
  fHistMCJetsZvsPt(0),
  fHistJetsPhiEta(0),
  fHistJetsPtArea(0),
  fHistJetsCorrPtArea(0),
  fHistJetsNEFvsPt(0),
  fHistJetsZvsPt(0),
  fHistMatchingLevelMCPt(0),
  fHistClosestDeltaEtaPhiMCPt(0),
  fHistClosestDeltaPtMCPt(0),
  fHistClosestDeltaCorrPtMCPt(0),
  fHistNonMatchedMCJetsPtArea(0),
  fHistNonMatchedJetsPtArea(0),
  fHistNonMatchedJetsCorrPtArea(0),
  fHistPartvsDetecPt(0),
  fHistPartvsDetecCorrPt(0),
  fHistMissedMCJetsPtArea(0)
{
  // Standard constructor.

  for (Int_t i = 0; i < 11; i++) {
    fHistEventWeight[i] = 0;
  }

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetResponseMaker::~AliJetResponseMaker()
{
  // Destructor
}

//________________________________________________________________________
void AliJetResponseMaker::UserCreateOutputObjects()
{
  // Create user objects.

  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  fHistNTrials = new TH1F("fHistNTrials", "fHistNTrials", 11, 0, 11);
  fHistNTrials->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistNTrials->GetYaxis()->SetTitle("trials");
  fOutput->Add(fHistNTrials);

  fHistEvents = new TH1F("fHistEvents", "fHistEvents", 11, 0, 11);
  fHistEvents->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistEvents->GetYaxis()->SetTitle("total events");
  fOutput->Add(fHistEvents);

  if (fEventWeightHist) {
    for (Int_t i = 0; i < 11; i++) {
      TString name(Form("fHistEventWeight_%d", i+1));
      fHistEventWeight[i] = new TH1F(name, name, 10, 0, 10);
      fOutput->Add(fHistEventWeight[i]);
    }
  }

  for (Int_t i = 1; i < 12; i++) {
    fHistNTrials->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
    fHistEvents->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
  }

  fHistJetsPhiEta = new TH2F("fHistJetsPhiEta", "fHistJetsPhiEta", 20, -2, 2, 32, 0, 6.4);
  fHistJetsPhiEta->GetXaxis()->SetTitle("#eta");
  fHistJetsPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJetsPhiEta);
  
  fHistJetsPtArea = new TH2F("fHistJetsPtArea", "fHistJetsPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistJetsPtArea->GetXaxis()->SetTitle("area");
  fHistJetsPtArea->GetYaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
  fHistJetsPtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsPtArea);

  fHistJetsCorrPtArea = new TH2F("fHistJetsCorrPtArea", "fHistJetsCorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, (Int_t)(1.5*fNbins), -fMaxBinPt/2, fMaxBinPt);
  fHistJetsCorrPtArea->GetXaxis()->SetTitle("area");
  fHistJetsCorrPtArea->GetYaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
  fHistJetsCorrPtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJetsCorrPtArea);
  
  fHistJetsZvsPt = new TH2F("fHistJetsZvsPt", "fHistJetsZvsPt", fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJetsZvsPt->GetXaxis()->SetTitle("Z");
  fHistJetsZvsPt->GetYaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
  fOutput->Add(fHistJetsZvsPt);
  
  fHistJetsNEFvsPt = new TH2F("fHistJetsNEFvsPt", "fHistJetsNEFvsPt", fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistJetsNEFvsPt->GetXaxis()->SetTitle("NEF");
  fHistJetsNEFvsPt->GetYaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
  fOutput->Add(fHistJetsNEFvsPt);

  fHistMCJetsPhiEta = new TH2F("fHistMCJetsPhiEta", "fHistMCJetsPhiEta", 20, -2, 2, 32, 0, 6.4);
  fHistMCJetsPhiEta->GetXaxis()->SetTitle("#eta");
  fHistMCJetsPhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistMCJetsPhiEta);
  
  fHistMCJetsPtArea = new TH2F("fHistMCJetsPtArea", "fHistMCJetsPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistMCJetsPtArea->GetXaxis()->SetTitle("area");
  fHistMCJetsPtArea->GetYaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fHistMCJetsPtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistMCJetsPtArea);

  fHistMCJetsPhiEtaFiducial = new TH2F("fHistMCJetsPhiEtaFiducial", "fHistMCJetsPhiEtaFiducial", 20, -2, 2, 32, 0, 6.4);
  fHistMCJetsPhiEtaFiducial->GetXaxis()->SetTitle("#eta");
  fHistMCJetsPhiEtaFiducial->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistMCJetsPhiEtaFiducial);
  
  fHistMCJetsPtAreaFiducial = new TH2F("fHistMCJetsPtAreaFiducial", "fHistMCJetsPtAreaFiducial", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistMCJetsPtAreaFiducial->GetXaxis()->SetTitle("area");
  fHistMCJetsPtAreaFiducial->GetYaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fHistMCJetsPtAreaFiducial->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistMCJetsPtAreaFiducial);
  
  fHistMCJetsZvsPt = new TH2F("fHistMCJetsZvsPt", "fHistMCJetsZvsPt", fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistMCJetsZvsPt->GetXaxis()->SetTitle("Z");
  fHistMCJetsZvsPt->GetYaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fOutput->Add(fHistMCJetsZvsPt);

  fHistMCJetsNEFvsPt = new TH2F("fHistMCJetsNEFvsPt", "fHistMCJetsNEFvsPt", fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistMCJetsNEFvsPt->GetXaxis()->SetTitle("NEF");
  fHistMCJetsNEFvsPt->GetYaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fOutput->Add(fHistMCJetsNEFvsPt);

  fHistMatchingLevelMCPt = new TH2F("fHistMatchingLevelMCPt", "fHistMatchingLevelMCPt", fNbins, 0, 1.2, fNbins, fMinBinPt, fMaxBinPt);
  fHistMatchingLevelMCPt->GetXaxis()->SetTitle("Matching level");
  fHistMatchingLevelMCPt->GetYaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fOutput->Add(fHistMatchingLevelMCPt);

  fHistClosestDeltaEtaPhiMCPt = new TH3F("fHistClosestDeltaEtaPhiMCPt", "fHistClosestDeltaEtaPhiMCPt", TMath::CeilNint(fJetMaxEta - fJetMinEta) * 20, fJetMinEta * 2, fJetMaxEta * 2, 64, -1.6, 4.8, fNbins, fMinBinPt, fMaxBinPt);
  fHistClosestDeltaEtaPhiMCPt->GetXaxis()->SetTitle("#Delta#eta");
  fHistClosestDeltaEtaPhiMCPt->GetYaxis()->SetTitle("#Delta#phi");
  fHistClosestDeltaEtaPhiMCPt->GetZaxis()->SetTitle("p_{T}^{gen}");
  fOutput->Add(fHistClosestDeltaEtaPhiMCPt);

  fHistClosestDeltaPtMCPt = new TH2F("fHistClosestDeltaPtMCPt", "fHistClosestDeltaPtMCPt", fNbins, fMinBinPt, fMaxBinPt, (Int_t)(fNbins*1.5), -fMaxBinPt / 2, fMaxBinPt);
  fHistClosestDeltaPtMCPt->GetXaxis()->SetTitle("p_{T}^{gen}");  
  fHistClosestDeltaPtMCPt->GetYaxis()->SetTitle("#Deltap_{T}^{rec} (GeV/c)");
  fHistClosestDeltaPtMCPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistClosestDeltaPtMCPt);

  fHistClosestDeltaCorrPtMCPt = new TH2F("fHistClosestDeltaCorrPtMCPt", "fHistClosestDeltaCorrPtMCPt", fNbins, fMinBinPt, fMaxBinPt, (Int_t)(fNbins*1.5), -fMaxBinPt / 2, fMaxBinPt);
  fHistClosestDeltaCorrPtMCPt->GetXaxis()->SetTitle("p_{T}^{gen}");  
  fHistClosestDeltaCorrPtMCPt->GetYaxis()->SetTitle("#Deltap_{T}^{corr} (GeV/c)");
  fHistClosestDeltaCorrPtMCPt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistClosestDeltaCorrPtMCPt);

  fHistNonMatchedMCJetsPtArea = new TH2F("fHistNonMatchedMCJetsPtArea", "fHistNonMatchedMCJetsPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistNonMatchedMCJetsPtArea->GetXaxis()->SetTitle("area");
  fHistNonMatchedMCJetsPtArea->GetYaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fHistNonMatchedMCJetsPtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistNonMatchedMCJetsPtArea);

  fHistNonMatchedJetsPtArea = new TH2F("fHistNonMatchedJetsPtArea", "fHistNonMatchedJetsPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistNonMatchedJetsPtArea->GetXaxis()->SetTitle("area");
  fHistNonMatchedJetsPtArea->GetYaxis()->SetTitle("p_{T}^{rec} (GeV/c)");
  fHistNonMatchedJetsPtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistNonMatchedJetsPtArea);

  fHistNonMatchedJetsCorrPtArea = new TH2F("fHistNonMatchedJetsCorrPtArea", "fHistNonMatchedJetsCorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, (Int_t)(1.5*fNbins), -fMaxBinPt/2, fMaxBinPt);
  fHistNonMatchedJetsCorrPtArea->GetXaxis()->SetTitle("area");
  fHistNonMatchedJetsCorrPtArea->GetYaxis()->SetTitle("p_{T}^{corr} (GeV/c)");
  fHistNonMatchedJetsCorrPtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistNonMatchedJetsCorrPtArea);

  fHistPartvsDetecPt = new TH2F("fHistPartvsDetecPt", "fHistPartvsDetecPt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistPartvsDetecPt->GetXaxis()->SetTitle("p_{T}^{rec}");
  fHistPartvsDetecPt->GetYaxis()->SetTitle("p_{T}^{gen}");
  fOutput->Add(fHistPartvsDetecPt);

  fHistPartvsDetecCorrPt = new TH2F("fHistPartvsDetecCorrPt", "fHistPartvsDetecCorrPt", (Int_t)(1.5*fNbins), -fMaxBinPt/2, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistPartvsDetecCorrPt->GetXaxis()->SetTitle("p_{T}^{corr}");
  fHistPartvsDetecCorrPt->GetYaxis()->SetTitle("p_{T}^{gen}");
  fOutput->Add(fHistPartvsDetecCorrPt);

  fHistMissedMCJetsPtArea = new TH2F("fHistMissedMCJetsPtArea", "fHistMissedMCJetsPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistMissedMCJetsPtArea->GetXaxis()->SetTitle("area");  
  fHistMissedMCJetsPtArea->GetYaxis()->SetTitle("p_{T} (GeV/c)");
  fHistMissedMCJetsPtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistMissedMCJetsPtArea);

  PostData(1, fOutput); // Post data for ALL output slots > 0 here, to get at least an empty histogram
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::AcceptJet(AliEmcalJet *jet) const
{   
  // Return true if jet is accepted.

  if (jet->Pt() <= fJetPtCut)
    return kFALSE;
  if (jet->Area() <= fJetAreaCut)
    return kFALSE;

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::ExecOnce()
{
  // Execute once.

  if (!fMCJetsName.IsNull() && !fMCJets) {
    fMCJets = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCJetsName));
    if (!fMCJets) {
      AliError(Form("%s: Could not retrieve mc jets %s!", GetName(), fMCJetsName.Data()));
      return;
    }
    else if (!fMCJets->GetClass()->GetBaseClass("AliEmcalJet")) {
      AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fMCJetsName.Data())); 
      fMCJets = 0;
      return;
    }
  }

  if (!fMCTracksName.IsNull() && !fMCTracks) {
    fMCTracks = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fMCTracksName));
    if (!fMCTracks) {
      AliError(Form("%s: Could not retrieve mc tracks %s!", GetName(), fMCTracksName.Data())); 
      return;
    }
    else {
      TClass *cl = fMCTracks->GetClass();
      if (!cl->GetBaseClass("AliVParticle") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fMCTracksName.Data())); 
	fMCTracks = 0;
	return;
      }
    }
  }

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fMCFiducial) {
    fMCminEta = fJetMinEta;
    fMCmaxEta = fJetMaxEta;
    fMCminPhi = fJetMinPhi;
    fMCmaxPhi = fJetMaxPhi;
  }
  else {
    fMCminEta = fJetMinEta - fJetRadius;
    fMCmaxEta = fJetMaxEta + fJetRadius;
    fMCminPhi = fJetMinPhi - fJetRadius;
    fMCmaxPhi = fJetMaxPhi + fJetRadius;
  }
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::IsEventSelected()
{
  // Check if event is selected

  if (fSelectPtHardBin != -999 && fSelectPtHardBin != fPtHardBin) 
    return kFALSE;

  return AliAnalysisTaskEmcalJet::IsEventSelected();
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::RetrieveEventObjects()
{
  // Retrieve event objects.

  if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
    return kFALSE;
  
  fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());

  if (!fPythiaHeader)
    return kFALSE;

  if (fDoWeighting)
    fEventWeight = fPythiaHeader->EventWeight();
  else
    fEventWeight = 1;

  const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};

  Double_t pthard = fPythiaHeader->GetPtHard();

  for (fPtHardBin = 0; fPtHardBin < 11; fPtHardBin++) {
    if (pthard >= ptHardLo[fPtHardBin] && pthard < ptHardHi[fPtHardBin])
      break;
  }

  fNTrials = fPythiaHeader->Trials();

  return kTRUE;
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::Run()
{
  // Find the closest jets

  if (!fDoMatching)
    return kTRUE;

  DoJetLoop(fJets, fMCJets, kFALSE);
  DoJetLoop(fMCJets, fJets, kTRUE);

  const Int_t nMCJets = fMCJets->GetEntriesFast();

  for (Int_t i = 0; i < nMCJets; i++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fMCJets->At(i));

    if (!jet) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    if (jet->Eta() < fMCminEta || jet->Eta() > fMCmaxEta || jet->Phi() < fMCminPhi || jet->Phi() > fMCmaxPhi)
      continue;

    if (jet->Pt() > fMaxBinPt)
      continue;

    if (jet->ClosestJet() && jet->ClosestJet()->ClosestJet() == jet && 
        jet->ClosestJetDistance() < fMaxDistance) {    // Matched jet found
      jet->SetMatchedToClosest();
      jet->ClosestJet()->SetMatchedToClosest();
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::DoJetLoop(TClonesArray *jets1, TClonesArray *jets2, Bool_t mc)
{
  // Do the jet loop.

  Int_t nJets1 = jets1->GetEntriesFast();
  Int_t nJets2 = jets2->GetEntriesFast();

  for (Int_t i = 0; i < nJets1; i++) {

    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(jets1->At(i));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet1))
      continue;

    if (!mc) {
      if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
	continue;
    }
    else {
      if (jet1->Eta() < fMCminEta || jet1->Eta() > fMCmaxEta || jet1->Phi() < fMCminPhi || jet1->Phi() > fMCmaxPhi)
	continue;
    }

    for (Int_t j = 0; j < nJets2; j++) {
      
      AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(jets2->At(j));
      
      if (!jet2) {
	AliError(Form("Could not receive jet %d", j));
	continue;
      }  
      
      if (!AcceptJet(jet2))
	continue;

      if (mc) {
	if (jet2->Eta() < fJetMinEta || jet2->Eta() > fJetMaxEta || jet2->Phi() < fJetMinPhi || jet2->Phi() > fJetMaxPhi)
	  continue;
      }
      else {
	if (jet1->Eta() < fMCminEta || jet1->Eta() > fMCmaxEta || jet1->Phi() < fMCminPhi || jet1->Phi() > fMCmaxPhi)
	  continue;
      }
      
      Double_t deta = jet2->Eta() - jet1->Eta();
      Double_t dphi = jet2->Phi() - jet1->Phi();
      Double_t d = TMath::Sqrt(deta * deta + dphi * dphi);

      if (d < jet1->ClosestJetDistance()) {
	jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
	jet1->SetClosestJet(jet2, d);
      }
      else if (d < jet1->SecondClosestJetDistance()) {
	jet1->SetSecondClosestJet(jet2, d);
      }
    }
  }
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::FillHistograms()
{
  // Fill histograms.

  fHistEvents->SetBinContent(fPtHardBin + 1, fHistEvents->GetBinContent(fPtHardBin + 1) + 1);
  fHistNTrials->SetBinContent(fPtHardBin + 1, fHistNTrials->GetBinContent(fPtHardBin + 1) + fNTrials);
  if (fEventWeightHist)
    fHistEventWeight[fPtHardBin]->Fill(fPythiaHeader->EventWeight());

  const Int_t nMCJets = fMCJets->GetEntriesFast();

  for (Int_t i = 0; i < nMCJets; i++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fMCJets->At(i));

    if (!jet) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet))
      continue;

    if (jet->Eta() < fMCminEta || jet->Eta() > fMCmaxEta || jet->Phi() < fMCminPhi || jet->Phi() > fMCmaxPhi)
      continue;

    if (jet->Pt() > fMaxBinPt)
      continue;

    if (jet->MatchedJet()) {

      if (!AcceptBiasJet(jet->MatchedJet()) || 
	  jet->MatchedJet()->MaxTrackPt() > fMaxTrackPt || jet->MatchedJet()->MaxClusterPt() > fMaxClusterPt ||
	  jet->MatchedJet()->Pt() > fMaxBinPt) {
	fHistMissedMCJetsPtArea->Fill(jet->Area(), jet->Pt(), fEventWeight);
      }
      else {
	fHistMatchingLevelMCPt->Fill(jet->ClosestJetDistance(), jet->Pt(), fEventWeight);

	Double_t deta = jet->MatchedJet()->Eta() - jet->Eta();
	Double_t dphi = jet->MatchedJet()->Phi() - jet->Phi();
	fHistClosestDeltaEtaPhiMCPt->Fill(deta, dphi, jet->Pt(), fEventWeight);

	Double_t dpt = jet->MatchedJet()->Pt() - jet->Pt();
	fHistClosestDeltaPtMCPt->Fill(jet->Pt(), dpt, fEventWeight);
	fHistClosestDeltaCorrPtMCPt->Fill(jet->Pt(), dpt - fRhoVal * jet->MatchedJet()->Area(), fEventWeight);

	fHistPartvsDetecPt->Fill(jet->MatchedJet()->Pt(), jet->Pt(), fEventWeight);
	fHistPartvsDetecCorrPt->Fill(jet->MatchedJet()->Pt() - fRhoVal * jet->MatchedJet()->Area(), jet->Pt(), fEventWeight);
      }
    }
    else {
      fHistNonMatchedMCJetsPtArea->Fill(jet->Area(), jet->Pt(), fEventWeight);
      fHistMissedMCJetsPtArea->Fill(jet->Area(), jet->Pt(), fEventWeight);
    }

    fHistMCJetsPtArea->Fill(jet->Area(), jet->Pt(), fEventWeight);
    fHistMCJetsPhiEta->Fill(jet->Eta(), jet->Phi(), fEventWeight);

    fHistMCJetsNEFvsPt->Fill(jet->NEF(), jet->Pt(), fEventWeight);
      
    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVParticle *track = jet->TrackAt(it, fMCTracks);
      if (track)
	fHistMCJetsZvsPt->Fill(track->Pt() / jet->Pt(), jet->Pt(), fEventWeight);
    }

    if (!AcceptBiasJet(jet))
      continue;
    if (jet->Eta() < fJetMinEta || jet->Eta() > fJetMaxEta || jet->Phi() < fJetMinPhi || jet->Phi() > fJetMaxPhi)
      continue;
    
    fHistMCJetsPtAreaFiducial->Fill(jet->Area(), jet->Pt(), fEventWeight);
    fHistMCJetsPhiEtaFiducial->Fill(jet->Eta(), jet->Phi(), fEventWeight);
  }

  const Int_t nJets = fJets->GetEntriesFast();

  for (Int_t i = 0; i < nJets; i++) {

    AliEmcalJet* jet = static_cast<AliEmcalJet*>(fJets->At(i));

    if (!jet) {
      AliError(Form("Could not receive mc jet %d", i));
      continue;
    }  
    
    if (!AcceptJet(jet))
      continue;
    if (!AcceptBiasJet(jet))
      continue;
    if (jet->MaxTrackPt() > fMaxTrackPt || jet->MaxClusterPt() > fMaxClusterPt)
      continue;
    if (jet->Eta() < fJetMinEta || jet->Eta() > fJetMaxEta || jet->Phi() < fJetMinPhi || jet->Phi() > fJetMaxPhi)
      continue;

    if (!jet->MatchedJet()) {
      fHistNonMatchedJetsPtArea->Fill(jet->Area(), jet->Pt(), fEventWeight);
      fHistNonMatchedJetsCorrPtArea->Fill(jet->Area(), jet->Pt() - fRhoVal * jet->Area(), fEventWeight);
    }

    fHistJetsPtArea->Fill(jet->Area(), jet->Pt(), fEventWeight);
    fHistJetsCorrPtArea->Fill(jet->Area(), jet->Pt() - fRhoVal * jet->Area(), fEventWeight);

    fHistJetsPhiEta->Fill(jet->Eta(), jet->Phi(), fEventWeight);

    fHistJetsNEFvsPt->Fill(jet->NEF(), jet->Pt(), fEventWeight);

    for (Int_t it = 0; it < jet->GetNumberOfTracks(); it++) {
      AliVParticle *track = jet->TrackAt(it, fTracks);
      if (track)
	fHistJetsZvsPt->Fill(track->Pt() / jet->Pt(), jet->Pt(), fEventWeight);
    }

    for (Int_t ic = 0; ic < jet->GetNumberOfClusters(); ic++) {
      AliVCluster *cluster = jet->ClusterAt(ic, fCaloClusters);
      if (cluster) {
	TLorentzVector nP;
	cluster->GetMomentum(nP, fVertex);
	fHistJetsZvsPt->Fill(nP.Pt() / jet->Pt(), jet->Pt(), fEventWeight);
      }
    }
  }

  return kTRUE;
}
