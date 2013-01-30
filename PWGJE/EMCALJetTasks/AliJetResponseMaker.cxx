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
#include "AliRhoParameter.h"

ClassImp(AliJetResponseMaker)

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker() : 
  AliAnalysisTaskEmcalJet("AliJetResponseMaker", kTRUE),
  fTracks2Name(""),
  fCalo2Name(""),
  fJets2Name(""),
  fRho2Name(""),
  fPtBiasJet2Track(0),
  fPtBiasJet2Clus(0),
  fAreCollections1MC(kFALSE),  
  fAreCollections2MC(kTRUE),
  fMatching(kNoMatching),
  fMatchingPar(0),
  fJet2MinEta(-999),
  fJet2MaxEta(-999),
  fJet2MinPhi(-999),
  fJet2MaxPhi(-999),
  fSelectPtHardBin(-999),
  fPythiaHeader(0),
  fPtHardBin(0),
  fNTrials(0),
  fTracks2(0),
  fCaloClusters2(0),
  fJets2(0),
  fRho2(0),
  fRho2Val(0),
  fHistNTrials(0),
  fHistEvents(0),
  fHistJets1PhiEta(0),
  fHistJets1PtArea(0),
  fHistJets1CorrPtArea(0),
  fHistJets2PhiEta(0),
  fHistJets2PtArea(0),
  fHistJets2CorrPtArea(0),
  fHistMatchingLevelvsJet2Pt(0),
  fHistClosestDeltaEtaPhivsJet2Pt(0),
  fHistClosestDeltaPtvsJet2Pt(0),
  fHistClosestDeltaCorrPtvsJet2Pt(0),
  fHistNonMatchedJets1PtArea(0),
  fHistNonMatchedJets2PtArea(0),
  fHistNonMatchedJets1CorrPtArea(0),
  fHistNonMatchedJets2CorrPtArea(0),
  fHistJet1PtvsJet2Pt(0),
  fHistJet1CorrPtvsJet2CorrPt(0),
  fHistMissedJets2PtArea(0)
{
  // Default constructor.

  SetMakeGeneralHistograms(kTRUE);
}

//________________________________________________________________________
AliJetResponseMaker::AliJetResponseMaker(const char *name) : 
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fTracks2Name("MCParticles"),
  fCalo2Name(""),
  fJets2Name("MCJets"),
  fRho2Name(""),
  fPtBiasJet2Track(0),
  fPtBiasJet2Clus(0),
  fAreCollections1MC(kFALSE),  
  fAreCollections2MC(kTRUE),
  fMatching(kNoMatching),
  fMatchingPar(0.25),
  fJet2MinEta(-999),
  fJet2MaxEta(-999),
  fJet2MinPhi(-999),
  fJet2MaxPhi(-999),
  fSelectPtHardBin(-999),
  fPythiaHeader(0),
  fPtHardBin(0),
  fNTrials(0),
  fTracks2(0),
  fCaloClusters2(0),
  fJets2(0),
  fRho2(0),
  fRho2Val(0),
  fHistNTrials(0),
  fHistEvents(0),
  fHistJets1PhiEta(0),
  fHistJets1PtArea(0),
  fHistJets1CorrPtArea(0),
  fHistJets2PhiEta(0),
  fHistJets2PtArea(0),
  fHistJets2CorrPtArea(0),
  fHistMatchingLevelvsJet2Pt(0),
  fHistClosestDeltaEtaPhivsJet2Pt(0),
  fHistClosestDeltaPtvsJet2Pt(0),
  fHistClosestDeltaCorrPtvsJet2Pt(0),
  fHistNonMatchedJets1PtArea(0),
  fHistNonMatchedJets2PtArea(0),
  fHistNonMatchedJets1CorrPtArea(0),
  fHistNonMatchedJets2CorrPtArea(0),
  fHistJet1PtvsJet2Pt(0),
  fHistJet1CorrPtvsJet2CorrPt(0),
  fHistMissedJets2PtArea(0)
{
  // Standard constructor.

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

  for (Int_t i = 1; i < 12; i++) {
    fHistNTrials->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
    fHistEvents->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
  }

  fHistJets1PhiEta = new TH2F("fHistJets1PhiEta", "fHistJets1PhiEta", 20, -2, 2, 32, 0, 6.4);
  fHistJets1PhiEta->GetXaxis()->SetTitle("#eta");
  fHistJets1PhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets1PhiEta);
  
  fHistJets1PtArea = new TH2F("fHistJets1PtArea", "fHistJets1PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets1PtArea->GetXaxis()->SetTitle("area");
  fHistJets1PtArea->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistJets1PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets1PtArea);

  if (!fRhoName.IsNull()) {
    fHistJets1CorrPtArea = new TH2F("fHistJets1CorrPtArea", "fHistJets1CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets1CorrPtArea->GetXaxis()->SetTitle("area");
    fHistJets1CorrPtArea->GetYaxis()->SetTitle("p_{T,1}^{corr} (GeV/c)");
    fHistJets1CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets1CorrPtArea);
  }

  fHistJets2PhiEta = new TH2F("fHistJets2PhiEta", "fHistJets2PhiEta", 20, -2, 2, 32, 0, 6.4);
  fHistJets2PhiEta->GetXaxis()->SetTitle("#eta");
  fHistJets2PhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets2PhiEta);
  
  fHistJets2PtArea = new TH2F("fHistJets2PtArea", "fHistJets2PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2PtArea->GetXaxis()->SetTitle("area");
  fHistJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2PtArea);

  if (!fRho2Name.IsNull()) {
    fHistJets2CorrPtArea = new TH2F("fHistJets2CorrPtArea", "fHistJets2CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets2CorrPtArea->GetXaxis()->SetTitle("area");
    fHistJets2CorrPtArea->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistJets2CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets2CorrPtArea);
  }

  fHistMatchingLevelvsJet2Pt = new TH2F("fHistMatchingLevelvsJet2Pt", "fHistMatchingLevelvsJet2Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistMatchingLevelvsJet2Pt->GetXaxis()->SetTitle("Matching level");
  fHistMatchingLevelvsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");  
  fHistMatchingLevelvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistMatchingLevelvsJet2Pt);

  fHistClosestDeltaEtaPhivsJet2Pt = new TH3F("fHistClosestDeltaEtaPhivsJet2Pt", "fHistClosestDeltaEtaPhivsJet2Pt", 40, -1, 1, 128, -1.6, 4.8, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistClosestDeltaEtaPhivsJet2Pt->GetXaxis()->SetTitle("#Delta#eta");
  fHistClosestDeltaEtaPhivsJet2Pt->GetYaxis()->SetTitle("#Delta#phi");
  fHistClosestDeltaEtaPhivsJet2Pt->GetZaxis()->SetTitle("p_{T,2}");
  fOutput->Add(fHistClosestDeltaEtaPhivsJet2Pt);

  fHistClosestDeltaPtvsJet2Pt = new TH2F("fHistClosestDeltaPtvsJet2Pt", "fHistClosestDeltaPtvsJet2Pt", fNbins/2, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistClosestDeltaPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
  fHistClosestDeltaPtvsJet2Pt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistClosestDeltaPtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistClosestDeltaPtvsJet2Pt);

  if (!fRhoName.IsNull() || !fRho2Name.IsNull()) {  
    fHistClosestDeltaCorrPtvsJet2Pt = new TH2F("fHistClosestDeltaCorrPtvsJet2Pt", "fHistClosestDeltaCorrPtvsJet2Pt", fNbins/2, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistClosestDeltaCorrPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
    fHistClosestDeltaCorrPtvsJet2Pt->GetYaxis()->SetTitle("#Deltap_{T}^{corr} (GeV/c)");
    fHistClosestDeltaCorrPtvsJet2Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistClosestDeltaCorrPtvsJet2Pt);
  }

  fHistNonMatchedJets1PtArea = new TH2F("fHistNonMatchedJets1PtArea", "fHistNonMatchedJets1PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistNonMatchedJets1PtArea->GetXaxis()->SetTitle("area");
  fHistNonMatchedJets1PtArea->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
  fHistNonMatchedJets1PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistNonMatchedJets1PtArea);

  fHistNonMatchedJets2PtArea = new TH2F("fHistNonMatchedJets2PtArea", "fHistNonMatchedJets2PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistNonMatchedJets2PtArea->GetXaxis()->SetTitle("area");
  fHistNonMatchedJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistNonMatchedJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistNonMatchedJets2PtArea);

  if (!fRhoName.IsNull()) {  
    fHistNonMatchedJets1CorrPtArea = new TH2F("fHistNonMatchedJets1CorrPtArea", "fHistNonMatchedJets1CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistNonMatchedJets1CorrPtArea->GetXaxis()->SetTitle("area");
    fHistNonMatchedJets1CorrPtArea->GetYaxis()->SetTitle("p_{T,1} (GeV/c)");
    fHistNonMatchedJets1CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistNonMatchedJets1CorrPtArea);
  }

  if (!fRho2Name.IsNull()) {  
    fHistNonMatchedJets2CorrPtArea = new TH2F("fHistNonMatchedJets2CorrPtArea", "fHistNonMatchedJets2CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistNonMatchedJets2CorrPtArea->GetXaxis()->SetTitle("area");
    fHistNonMatchedJets2CorrPtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
    fHistNonMatchedJets2CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistNonMatchedJets2CorrPtArea);
  }

  fHistJet1PtvsJet2Pt = new TH2F("fHistJet1PtvsJet2Pt", "fHistJet1PtvsJet2Pt", fNbins, fMinBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
  fHistJet1PtvsJet2Pt->GetXaxis()->SetTitle("p_{T,1}");
  fHistJet1PtvsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");
  fHistJet1PtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJet1PtvsJet2Pt);

  if (!fRhoName.IsNull() || !fRho2Name.IsNull()) {
    if (fRhoName.IsNull()) 
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", fNbins, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    else if (fRho2Name.IsNull()) 
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", 2*fNbins, -fMaxBinPt, fMaxBinPt, fNbins, fMinBinPt, fMaxBinPt);
    else
      fHistJet1CorrPtvsJet2CorrPt = new TH2F("fHistJet1CorrPtvsJet2CorrPt", "fHistJet1CorrPtvsJet2CorrPt", 2*fNbins, -fMaxBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJet1CorrPtvsJet2CorrPt->GetXaxis()->SetTitle("p_{T,1}^{corr}");
    fHistJet1CorrPtvsJet2CorrPt->GetYaxis()->SetTitle("p_{T,2}^{corr}");
    fHistJet1CorrPtvsJet2CorrPt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJet1CorrPtvsJet2CorrPt);
  }

  fHistMissedJets2PtArea = new TH2F("fHistMissedJets2PtArea", "fHistMissedJets2PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistMissedJets2PtArea->GetXaxis()->SetTitle("area");  
  fHistMissedJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistMissedJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistMissedJets2PtArea);

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
Bool_t AliJetResponseMaker::AcceptBiasJet2(AliEmcalJet *jet) const
{ 
  // Accept jet with a bias.

  if (fLeadingHadronType == 0) {
    if (jet->MaxTrackPt() < fPtBiasJet2Track) return kFALSE;
  }
  else if (fLeadingHadronType == 1) {
    if (jet->MaxClusterPt() < fPtBiasJet2Clus) return kFALSE;
  }
  else {
    if (jet->MaxTrackPt() < fPtBiasJet2Track && jet->MaxClusterPt() < fPtBiasJet2Clus) return kFALSE;
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::ExecOnce()
{
  // Execute once.

  if (!fJets2Name.IsNull() && !fJets2) {
    fJets2 = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fJets2Name));
    if (!fJets2) {
      AliError(Form("%s: Could not retrieve jets2 %s!", GetName(), fJets2Name.Data()));
      return;
    }
    else if (!fJets2->GetClass()->GetBaseClass("AliEmcalJet")) {
      AliError(Form("%s: Collection %s does not contain AliEmcalJet objects!", GetName(), fJets2Name.Data())); 
      fJets2 = 0;
      return;
    }
  }

  if (!fTracks2Name.IsNull() && !fTracks2) {
    fTracks2 = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracks2Name));
    if (!fTracks2) {
      AliError(Form("%s: Could not retrieve tracks2 %s!", GetName(), fTracks2Name.Data())); 
      return;
    }
    else {
      TClass *cl = fTracks2->GetClass();
      if (!cl->GetBaseClass("AliVParticle") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVParticle nor AliEmcalParticle objects!", GetName(), fTracks2Name.Data())); 
	fTracks2 = 0;
	return;
      }
    }
  }

  if (!fCalo2Name.IsNull() && !fCaloClusters2) {
    fCaloClusters2 =  dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fCalo2Name));
    if (!fCaloClusters2) {
      AliError(Form("%s: Could not retrieve clusters %s!", GetName(), fCalo2Name.Data())); 
      return;
    } else {
      TClass *cl = fCaloClusters2->GetClass();
      if (!cl->GetBaseClass("AliVCluster") && !cl->GetBaseClass("AliEmcalParticle")) {
	AliError(Form("%s: Collection %s does not contain AliVCluster nor AliEmcalParticle objects!", GetName(), fCalo2Name.Data())); 
	fCaloClusters2 = 0;
	return;
      }
    }
  }

  if (!fRho2Name.IsNull() && !fRho2) {
    fRho2 = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(fRho2Name));
    if (!fRho2) {
      AliError(Form("%s: Could not retrieve rho %s!", GetName(), fRho2Name.Data()));
      fInitialized = kFALSE;
      return;
    }
  }

  if (fJet2MinEta == -999)
    fJet2MinEta = fJetMinEta - fJetRadius;
  if (fJet2MaxEta == -999)
    fJet2MaxEta = fJetMaxEta + fJetRadius;
  if (fJet2MinPhi == -999)
    fJet2MinPhi = fJetMinPhi - fJetRadius;
  if (fJet2MaxPhi == -999)
    fJet2MaxPhi = fJetMaxPhi + fJetRadius;

  AliAnalysisTaskEmcalJet::ExecOnce();
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

  if (fRho2)
    fRho2Val = fRho2->GetVal();
  
  fPythiaHeader = dynamic_cast<AliGenPythiaEventHeader*>(MCEvent()->GenEventHeader());

  if (!fPythiaHeader)
    return kFALSE;

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
  
  if (fMatching == kNoMatching)
    return kTRUE;

  DoJetLoop(fJets, fJets2, kFALSE);
  DoJetLoop(fJets2, fJets, kTRUE);

  const Int_t nJets2 = fJets2->GetEntriesFast();

  for (Int_t i = 0; i < nJets2; i++) {

    AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(fJets2->At(i));

    if (!jet2) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet2))
      continue;

    if (jet2->Eta() < fJet2MinEta || jet2->Eta() > fJet2MaxEta || jet2->Phi() < fJet2MinPhi || jet2->Phi() > fJet2MaxPhi)
      continue;

    if (jet2->Pt() > fMaxBinPt)
      continue;

    if (jet2->ClosestJet() && jet2->ClosestJet()->ClosestJet() == jet2 && 
        jet2->ClosestJetDistance() < fMatchingPar) {    // Matched jet found
      jet2->SetMatchedToClosest();
      jet2->ClosestJet()->SetMatchedToClosest();
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::DoJetLoop(TClonesArray *jets1, TClonesArray *jets2, Bool_t first)
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

    if (first) {
     if (jet1->Eta() < fJet2MinEta || jet1->Eta() > fJet2MaxEta || jet1->Phi() < fJet2MinPhi || jet1->Phi() > fJet2MaxPhi)
	continue;
    }
    else {
      if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
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

      if (first) {
	if (jet2->Eta() < fJetMinEta || jet2->Eta() > fJetMaxEta || jet2->Phi() < fJetMinPhi || jet2->Phi() > fJetMaxPhi)
	  continue;
      }
      else {
	if (jet1->Eta() < fJet2MinEta || jet1->Eta() > fJet2MaxEta || jet1->Phi() < fJet2MinPhi || jet1->Phi() > fJet2MaxPhi)
	  continue;
      }

      Double_t d = GetMatchingLevel(jet1, jet2);

      if (d < 0)
	continue;

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
Double_t AliJetResponseMaker::GetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2) const
{
  Double_t r = -1;

  switch (fMatching) {
  case kGeometrical:
    {
      Double_t deta = jet2->Eta() - jet1->Eta();
      Double_t dphi = jet2->Phi() - jet1->Phi();
      r = TMath::Sqrt(deta * deta + dphi * dphi);
    }
    break;
  case kMCLabel:
    AliError("MC label matching not implemented!");
    break;
  default:
    ;
  }
  
  return r;
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::FillHistograms()
{
  // Fill histograms.

  fHistEvents->SetBinContent(fPtHardBin + 1, fHistEvents->GetBinContent(fPtHardBin + 1) + 1);
  fHistNTrials->SetBinContent(fPtHardBin + 1, fHistNTrials->GetBinContent(fPtHardBin + 1) + fNTrials);

  const Int_t nJets2 = fJets2->GetEntriesFast();

  for (Int_t i = 0; i < nJets2; i++) {

    AliEmcalJet* jet2 = static_cast<AliEmcalJet*>(fJets2->At(i));

    if (!jet2) {
      AliError(Form("Could not receive jet2 %d", i));
      continue;
    }  

    if (!AcceptJet(jet2))
      continue;

    if (!AcceptBiasJet2(jet2))
      continue;

    if (jet2->Eta() < fJet2MinEta || jet2->Eta() > fJet2MaxEta || jet2->Phi() < fJet2MinPhi || jet2->Phi() > fJet2MaxPhi)
      continue;

    if (jet2->Pt() > fMaxBinPt)
      continue;

    if (jet2->MatchedJet()) {

      if (!AcceptBiasJet(jet2->MatchedJet()) || 
	  jet2->MatchedJet()->MaxTrackPt() > fMaxTrackPt || jet2->MatchedJet()->MaxClusterPt() > fMaxClusterPt ||
	  jet2->MatchedJet()->Pt() > fMaxBinPt) {
	fHistMissedJets2PtArea->Fill(jet2->Area(), jet2->Pt());
      }
      else {
	fHistMatchingLevelvsJet2Pt->Fill(jet2->ClosestJetDistance(), jet2->Pt());

	Double_t deta = jet2->MatchedJet()->Eta() - jet2->Eta();
	Double_t dphi = jet2->MatchedJet()->Phi() - jet2->Phi();
	fHistClosestDeltaEtaPhivsJet2Pt->Fill(deta, dphi, jet2->Pt());

	Double_t dpt = jet2->MatchedJet()->Pt() - jet2->Pt();
	fHistClosestDeltaPtvsJet2Pt->Fill(jet2->Pt(), dpt);

	fHistJet1PtvsJet2Pt->Fill(jet2->MatchedJet()->Pt(), jet2->Pt());
	
	if (!fRhoName.IsNull() || !fRho2Name.IsNull()) {
	  dpt -= fRhoVal * jet2->MatchedJet()->Area() - fRho2Val * jet2->Area();
	  fHistClosestDeltaCorrPtvsJet2Pt->Fill(jet2->Pt(), dpt);
	  fHistJet1CorrPtvsJet2CorrPt->Fill(jet2->MatchedJet()->Pt() - fRhoVal * jet2->MatchedJet()->Area(), jet2->Pt() - fRho2Val * jet2->Area());
	}
      }
    }
    else {
      fHistNonMatchedJets2PtArea->Fill(jet2->Area(), jet2->Pt());
      fHistMissedJets2PtArea->Fill(jet2->Area(), jet2->Pt());

      if (!fRho2Name.IsNull())
	fHistNonMatchedJets2CorrPtArea->Fill(jet2->Area(), jet2->Pt() - fRhoVal * jet2->Area());
    }

    fHistJets2PtArea->Fill(jet2->Area(), jet2->Pt());
    fHistJets2PhiEta->Fill(jet2->Eta(), jet2->Phi());

    if (!fRho2Name.IsNull())
      fHistJets2CorrPtArea->Fill(jet2->Area(), jet2->Pt() - fRho2Val * jet2->Area());
  }

  const Int_t nJets1 = fJets->GetEntriesFast();

  for (Int_t i = 0; i < nJets1; i++) {

    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(fJets->At(i));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  
    
    if (!AcceptJet(jet1))
      continue;

    if (!AcceptBiasJet(jet1))
      continue;

    if (jet1->MaxTrackPt() > fMaxTrackPt || jet1->MaxClusterPt() > fMaxClusterPt)
      continue;

    if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
      continue;

    if (!jet1->MatchedJet()) {
      fHistNonMatchedJets1PtArea->Fill(jet1->Area(), jet1->Pt());
      if (!fRhoName.IsNull())
	fHistNonMatchedJets1CorrPtArea->Fill(jet1->Area(), jet1->Pt() - fRhoVal * jet1->Area());
    }

    fHistJets1PtArea->Fill(jet1->Area(), jet1->Pt());
    fHistJets1PhiEta->Fill(jet1->Eta(), jet1->Phi());

    if (!fRhoName.IsNull())
      fHistJets1CorrPtArea->Fill(jet1->Area(), jet1->Pt() - fRhoVal * jet1->Area());
  }

  return kTRUE;
}
