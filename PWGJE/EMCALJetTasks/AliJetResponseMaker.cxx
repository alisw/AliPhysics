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
  fMatchingPar1(0),
  fMatchingPar2(0),
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
  fTracks2Map(0),
  fHistNTrials(0),
  fHistEvents(0),
  fHistJets1PhiEta(0),
  fHistJets1PtArea(0),
  fHistJets1CorrPtArea(0),
  fHistJets2PhiEta(0),
  fHistJets2PtArea(0),
  fHistJets2CorrPtArea(0),
  fHistJets2PhiEtaAcceptance(0),
  fHistJets2PtAreaAcceptance(0),
  fHistJets2CorrPtAreaAcceptance(0),
  fHistMatchingLevelvsJet2Pt(0),
  fHistDistancevsCommonEnergy(0),
  fHistDeltaEtaPhivsJet2Pt(0),
  fHistDeltaPtvsJet2Pt(0),
  fHistDeltaPtvsMatchingLevel(0),
  fHistDeltaCorrPtvsJet2Pt(0),
  fHistDeltaCorrPtvsMatchingLevel(0),
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
  fMatchingPar1(0),
  fMatchingPar2(0),
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
  fTracks2Map(0),
  fHistNTrials(0),
  fHistEvents(0),
  fHistJets1PhiEta(0),
  fHistJets1PtArea(0),
  fHistJets1CorrPtArea(0),
  fHistJets2PhiEta(0),
  fHistJets2PtArea(0),
  fHistJets2CorrPtArea(0),
  fHistJets2PhiEtaAcceptance(0),
  fHistJets2PtAreaAcceptance(0),
  fHistJets2CorrPtAreaAcceptance(0),
  fHistMatchingLevelvsJet2Pt(0),
  fHistDistancevsCommonEnergy(0),
  fHistDeltaEtaPhivsJet2Pt(0),
  fHistDeltaPtvsJet2Pt(0),
  fHistDeltaPtvsMatchingLevel(0),
  fHistDeltaCorrPtvsJet2Pt(0),
  fHistDeltaCorrPtvsMatchingLevel(0),
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

  fHistJets1PhiEta = new TH2F("fHistJets1PhiEta", "fHistJets1PhiEta", 40, -1, 1, 40, 0, TMath::Pi()*2);
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

  fHistJets2PhiEta = new TH2F("fHistJets2PhiEta", "fHistJets2PhiEta", 40, -1, 1, 40, 0, TMath::Pi()*2);
  fHistJets2PhiEta->GetXaxis()->SetTitle("#eta");
  fHistJets2PhiEta->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets2PhiEta);

  fHistJets2PhiEtaAcceptance = new TH2F("fHistJets2PhiEtaAcceptance", "fHistJets2PhiEtaAcceptance", 40, -1, 1, 40, 0, TMath::Pi()*2);
  fHistJets2PhiEtaAcceptance->GetXaxis()->SetTitle("#eta");
  fHistJets2PhiEtaAcceptance->GetYaxis()->SetTitle("#phi");
  fOutput->Add(fHistJets2PhiEtaAcceptance);
  
  fHistJets2PtArea = new TH2F("fHistJets2PtArea", "fHistJets2PtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2PtArea->GetXaxis()->SetTitle("area");
  fHistJets2PtArea->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2PtArea->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2PtArea);

  fHistJets2PtAreaAcceptance = new TH2F("fHistJets2PtAreaAcceptance", "fHistJets2PtAreaAcceptance", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, fNbins, fMinBinPt, fMaxBinPt);
  fHistJets2PtAreaAcceptance->GetXaxis()->SetTitle("area");
  fHistJets2PtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2} (GeV/c)");
  fHistJets2PtAreaAcceptance->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistJets2PtAreaAcceptance);

  if (!fRho2Name.IsNull()) {
    fHistJets2CorrPtArea = new TH2F("fHistJets2CorrPtArea", "fHistJets2CorrPtArea", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets2CorrPtArea->GetXaxis()->SetTitle("area");
    fHistJets2CorrPtArea->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistJets2CorrPtArea->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets2CorrPtArea);

    fHistJets2CorrPtAreaAcceptance = new TH2F("fHistJets2CorrPtAreaAcceptance", "fHistJets2CorrPtAreaAcceptance", 40, 0, fJetRadius * fJetRadius * TMath::Pi() * 3, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistJets2CorrPtAreaAcceptance->GetXaxis()->SetTitle("area");
    fHistJets2CorrPtAreaAcceptance->GetYaxis()->SetTitle("p_{T,2}^{corr} (GeV/c)");
    fHistJets2CorrPtAreaAcceptance->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistJets2CorrPtAreaAcceptance);
  }

  fHistMatchingLevelvsJet2Pt = new TH2F("fHistMatchingLevelvsJet2Pt", "fHistMatchingLevelvsJet2Pt", fNbins/2, 0, 1.2, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistMatchingLevelvsJet2Pt->GetXaxis()->SetTitle("Matching level");
  fHistMatchingLevelvsJet2Pt->GetYaxis()->SetTitle("p_{T,2}");  
  fHistMatchingLevelvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistMatchingLevelvsJet2Pt);

  fHistDistancevsCommonEnergy = new TH2F("fHistDistancevsCommonEnergy", "fHistDistancevsCommonEnergy", fNbins/2, 0, 1.2, fNbins/2, 0, 1.2);
  fHistDistancevsCommonEnergy->GetXaxis()->SetTitle("Distance");
  fHistDistancevsCommonEnergy->GetYaxis()->SetTitle("Common energy (%)");  
  fHistDistancevsCommonEnergy->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDistancevsCommonEnergy);

  fHistDeltaEtaPhivsJet2Pt = new TH3F("fHistDeltaEtaPhivsJet2Pt", "fHistDeltaEtaPhivsJet2Pt", 40, -1, 1, 128, -1.6, 4.8, fNbins/2, fMinBinPt, fMaxBinPt);
  fHistDeltaEtaPhivsJet2Pt->GetXaxis()->SetTitle("#Delta#eta");
  fHistDeltaEtaPhivsJet2Pt->GetYaxis()->SetTitle("#Delta#phi");
  fHistDeltaEtaPhivsJet2Pt->GetZaxis()->SetTitle("p_{T,2}");
  fOutput->Add(fHistDeltaEtaPhivsJet2Pt);

  fHistDeltaPtvsJet2Pt = new TH2F("fHistDeltaPtvsJet2Pt", "fHistDeltaPtvsJet2Pt", fNbins/2, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
  fHistDeltaPtvsJet2Pt->GetYaxis()->SetTitle("#deltap_{T} (GeV/c)");
  fHistDeltaPtvsJet2Pt->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsJet2Pt);

  fHistDeltaPtvsMatchingLevel = new TH2F("fHistDeltaPtvsMatchingLevel", "fHistDeltaPtvsMatchingLevel", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
  fHistDeltaPtvsMatchingLevel->GetXaxis()->SetTitle("Matching level");  
  fHistDeltaPtvsMatchingLevel->GetYaxis()->SetTitle("#Deltap_{T} (GeV/c)");
  fHistDeltaPtvsMatchingLevel->GetZaxis()->SetTitle("counts");
  fOutput->Add(fHistDeltaPtvsMatchingLevel);

  if (!fRhoName.IsNull() || !fRho2Name.IsNull()) {  
    fHistDeltaCorrPtvsJet2Pt = new TH2F("fHistDeltaCorrPtvsJet2Pt", "fHistDeltaCorrPtvsJet2Pt", fNbins/2, fMinBinPt, fMaxBinPt, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsJet2Pt->GetXaxis()->SetTitle("p_{T,2}");  
    fHistDeltaCorrPtvsJet2Pt->GetYaxis()->SetTitle("#Deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsJet2Pt->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsJet2Pt);

    fHistDeltaCorrPtvsMatchingLevel = new TH2F("fHistDeltaCorrPtvsMatchingLevel", "fHistDeltaCorrPtvsMatchingLevel", fNbins/2, 0, 1.2, 2*fNbins, -fMaxBinPt, fMaxBinPt);
    fHistDeltaCorrPtvsMatchingLevel->GetXaxis()->SetTitle("Matching level");  
    fHistDeltaCorrPtvsMatchingLevel->GetYaxis()->SetTitle("#Deltap_{T}^{corr} (GeV/c)");
    fHistDeltaCorrPtvsMatchingLevel->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistDeltaCorrPtvsJet2Pt);
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

    if (fAreCollections2MC) {
      fTracks2Map = dynamic_cast<TH1*>(InputEvent()->FindListObject(fTracks2Name + "_Map"));
      // this is needed to map the MC labels with the indexes of the MC particle collection
      // if teh map is not given, the MC labels are assumed to be consistent with the indexes (which is not the case if AliEmcalMCTrackSelector is used)
      if (!fTracks2Map) {
	AliWarning(Form("%s: Could not retrieve map for tracks2 %s! Will assume MC labels consistent with indexes...", GetName(), fTracks2Name.Data())); 
	fTracks2Map = new TH1I("tracksMap","tracksMap",9999,0,1);
	for (Int_t i = 0; i < 9999; i++) {
	  fTracks2Map->SetBinContent(i,i);
	}
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
  else
    return DoJetMatching();
}

//________________________________________________________________________
Bool_t AliJetResponseMaker::DoJetMatching()
{
  DoJetLoop(kFALSE);

  const Int_t nJets = fJets->GetEntriesFast();

  for (Int_t i = 0; i < nJets; i++) {

    AliEmcalJet* jet1 = static_cast<AliEmcalJet*>(fJets->At(i));

    if (!jet1) {
      AliError(Form("Could not receive jet %d", i));
      continue;
    }  

    if (!AcceptJet(jet1))
      continue;

    if (jet1->Eta() < fJetMinEta || jet1->Eta() > fJetMaxEta || jet1->Phi() < fJetMinPhi || jet1->Phi() > fJetMaxPhi)
      continue;

    if (jet1->Pt() > fMaxBinPt)
      continue;

    if (jet1->ClosestJet() && jet1->ClosestJet()->ClosestJet() == jet1 && 
        jet1->ClosestJetDistance() < fMatchingPar1 && jet1->ClosestJet()->ClosestJetDistance() < fMatchingPar2) {    // Matched jet found
      jet1->SetMatchedToClosest(fMatching);
      jet1->ClosestJet()->SetMatchedToClosest(fMatching);
    }
  }

  return kTRUE;
}

//________________________________________________________________________
void AliJetResponseMaker::DoJetLoop(Bool_t order)
{
  // Do the jet loop.

  TClonesArray *jets1 = 0;
  TClonesArray *jets2 = 0;

  if (order) {
    jets1 = fJets2;
    jets2 = fJets;
  }
  else {
    jets1 = fJets;
    jets2 = fJets2;
  }

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

    if (order) {
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

      if (order) {
	if (jet2->Eta() < fJetMinEta || jet2->Eta() > fJetMaxEta || jet2->Phi() < fJetMinPhi || jet2->Phi() > fJetMaxPhi)
	  continue;
      }
      else {
	if (jet1->Eta() < fJet2MinEta || jet1->Eta() > fJet2MaxEta || jet1->Phi() < fJet2MinPhi || jet1->Phi() > fJet2MaxPhi)
	  continue;
      }

      SetMatchingLevel(jet1, jet2, fMatching);
    }
  }
}

//________________________________________________________________________
void AliJetResponseMaker::GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const
{
  Double_t deta = jet2->Eta() - jet1->Eta();
  Double_t dphi = jet2->Phi() - jet1->Phi();
  d = TMath::Sqrt(deta * deta + dphi * dphi);
}

//________________________________________________________________________
void AliJetResponseMaker::GetMCLabelMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const
{ 
  // d1 and d2 represent the matching level: 0 = maximum level of matching, 1 = the two jets are completely unrelated
  d1 = jet1->Pt();
  d2 = jet2->Pt();
  Double_t totalPt1 = d1; // the total pt of the reconstructed jet will be cleaned from the background

  for (Int_t iTrack2 = 0; iTrack2 < jet2->GetNumberOfTracks(); iTrack2++) {
    Bool_t track2Found = kFALSE;
    Int_t index2 = jet2->TrackAt(iTrack2);
    for (Int_t iTrack = 0; iTrack < jet1->GetNumberOfTracks(); iTrack++) {
      AliVParticle *track = jet1->TrackAt(iTrack,fTracks);
      if (!track) {
	AliWarning(Form("Could not find track %d!", iTrack));
	continue;
      }
      Int_t MClabel = track->GetLabel();
      Int_t index = -1;
	  
      if (MClabel <= 0) {// this is not a MC particle; remove it completely
	AliDebug(3,Form("Track %d (pT = %f) is not a MC particle (MClabel = %d)!",iTrack,track->Pt(),MClabel));
	totalPt1 -= track->Pt();
	d1 -= track->Pt();
	continue;
      }
      else if (MClabel < fTracks2Map->GetNbinsX()-2) {
	index = fTracks2Map->GetBinContent(MClabel);
      }
	  
      if (index < 0) {
	AliDebug(2,Form("Track %d (pT = %f) does not have an associated MC particle (MClabel = %d)!",iTrack,track->Pt(),MClabel));
	continue;
      }

      if (index2 == index) { // found common particle
	track2Found = kTRUE;
	d1 -= track->Pt();
	AliVParticle *MCpart = static_cast<AliVParticle*>(fTracks2->At(index2));
	AliDebug(3,Form("Track %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
			iTrack,track->Pt(),track->Eta(),track->Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));
	d2 -= MCpart->Pt();
	break;
      }
    }
    for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
      AliVCluster *clus = jet1->ClusterAt(iClus,fCaloClusters);
      if (!clus) {
	AliWarning(Form("Could not find cluster %d!", iClus));
	continue;
      }
      TLorentzVector part;
      clus->GetMomentum(part, const_cast<Double_t*>(fVertex));
	  
      if (fCaloCells) { // if the cell colection is available, look for cells with a matched MC particle
	for (Int_t iCell = 0; iCell < clus->GetNCells(); iCell++) {
	  Int_t cellId = clus->GetCellAbsId(iCell);
	  Double_t cellFrac = clus->GetCellAmplitudeFraction(iCell);

	  Int_t MClabel = fCaloCells->GetCellMCLabel(cellId);
	  Int_t index = -1;
	  
	  if (MClabel <= 0) {// this is not a MC particle; remove it completely
	    AliDebug(3,Form("Cell %d (frac = %f) is not a MC particle (MClabel = %d)!",iCell,cellFrac,MClabel));
	    totalPt1 -= part.Pt() * cellFrac;
	    d1 -= part.Pt() * cellFrac;
	    continue;
	  }
	  else if (MClabel < fTracks2Map->GetNbinsX()-2) {
	    index = fTracks2Map->GetBinContent(MClabel);
	  }

	  if (index < 0) {
	    AliDebug(3,Form("Cell %d (frac = %f) does not have an associated MC particle (MClabel = %d)!",iCell,cellFrac,MClabel));
	    continue;
	  }
	  if (index2 == index) { // found common particle
	    d1 -= part.Pt() * cellFrac;
		
	    if (!track2Found) {// only if it is not already found among charged tracks (charged particles are most likely already found)
	      AliVParticle *MCpart = static_cast<AliVParticle*>(fTracks2->At(index2));
	      AliDebug(3,Form("Cell %d belonging to cluster %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
			      iCell,iClus,part.Pt(),part.Eta(),part.Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));		  
	      d2 -= MCpart->Pt() * cellFrac;
	    }
	    break;
	  }
	}
      }
      else { //otherwise look for the first contributor to the cluster, and if matched to a MC label remove it
	Int_t MClabel = clus->GetLabel();
	Int_t index = -1;
	    
	if (MClabel <= 0) {// this is not a MC particle; remove it completely
	  AliDebug(3,Form("Cluster %d (pT = %f) is not a MC particle (MClabel = %d)!",iClus,part.Pt(),MClabel));
	  totalPt1 -= part.Pt();
	  d1 -= part.Pt();
	  continue;
	}
	else if (MClabel < fTracks2Map->GetNbinsX()-2) {
	  index = fTracks2Map->GetBinContent(MClabel);
	}
	 
	if (index < 0) {
	  AliDebug(3,Form("Cluster %d (pT = %f) does not have an associated MC particle (MClabel = %d)!",iClus,part.Pt(),MClabel));
	  continue;
	}
	if (index2 == index) { // found common particle
	  d1 -= part.Pt();

	  if (!track2Found) {// only if it is not already found among charged tracks (charged particles are most likely already found)
	    AliVParticle *MCpart = static_cast<AliVParticle*>(fTracks2->At(index2));
	    AliDebug(3,Form("Cluster %d (pT = %f, eta = %f, phi = %f) is associated with the MC particle %d (pT = %f, eta = %f, phi = %f)!",
			    iClus,part.Pt(),part.Eta(),part.Phi(),MClabel,MCpart->Pt(),MCpart->Eta(),MCpart->Phi()));
		
	    d2 -= MCpart->Pt();
	  }
	  break;
	}
      }
    }
  }
  if (d1 <= 0 || totalPt1 < 1)
    d1 = 0;
  else
    d1 /= totalPt1;

  if (jet2->Pt() > 0 && d2 > 0)
    d2 /= jet2->Pt();
  else
    d2 = 0;
}

//________________________________________________________________________
void AliJetResponseMaker::GetSameCollectionsMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const
{ 
  // d1 and d2 represent the matching level: 0 = maximum level of matching, 1 = the two jets are completely unrelated
  d1 = jet1->Pt();
  d2 = jet2->Pt();

  if (fTracks && fTracks2) {

    for (Int_t iTrack2 = 0; iTrack2 < jet2->GetNumberOfTracks(); iTrack2++) {
      Int_t index2 = jet2->TrackAt(iTrack2);
      for (Int_t iTrack = 0; iTrack < jet1->GetNumberOfTracks(); iTrack++) {
	Int_t index = jet1->TrackAt(iTrack);
	if (index2 == index) { // found common particle
	  AliVParticle *part = static_cast<AliVParticle*>(fTracks->At(index));
	  if (!part) {
	    AliWarning(Form("Could not find track %d!", index));
	    continue;
	  }
	  AliVParticle *part2 = static_cast<AliVParticle*>(fTracks2->At(index2));
	  if (!part2) {
	    AliWarning(Form("Could not find track %d!", index2));
	    continue;
	  }

	  d1 -= part->Pt();
	  d2 -= part2->Pt();
	  break;
	}
      }
    }

  }

  if (fCaloClusters && fCaloClusters2) {

    for (Int_t iClus2 = 0; iClus2 < jet2->GetNumberOfClusters(); iClus2++) {
      Int_t index2 = jet2->ClusterAt(iClus2);
      for (Int_t iClus = 0; iClus < jet1->GetNumberOfClusters(); iClus++) {
	Int_t index = jet1->ClusterAt(iClus);
	AliVCluster *clus =  static_cast<AliVCluster*>(fCaloClusters->At(index));
	if (!clus) {
	  AliWarning(Form("Could not find cluster %d!", index));
	  continue;
	}
	AliVCluster *clus2 =  static_cast<AliVCluster*>(fCaloClusters2->At(index2));
	if (!clus2) {
	  AliWarning(Form("Could not find cluster %d!", index2));
	  continue;
	}
	TLorentzVector part, part2;
	clus->GetMomentum(part, const_cast<Double_t*>(fVertex));
	clus2->GetMomentum(part2, const_cast<Double_t*>(fVertex));

	d1 -= part.Pt();
	d2 -= part2.Pt();
	break;
      }
    }

  }

  if (jet1->Pt() > 0 && d1 > 0)
    d1 /= jet1->Pt();
  else
    d1 = 0;

  if (jet2->Pt() > 0 && d2 > 0)
    d2 /= jet2->Pt();
  else
    d2 = 0;
}

//________________________________________________________________________
void AliJetResponseMaker::SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, MatchingType matching) 
{
  Double_t d1 = -1;
  Double_t d2 = -1;

  switch (matching) {
  case kGeometrical:
    GetGeometricalMatchingLevel(jet1,jet2,d1);
    d2 = d1;
    break;
  case kMCLabel: // jet1 = detector level and jet2 = particle level!
    GetMCLabelMatchingLevel(jet1,jet2,d1,d2);
    break;
  case kSameCollections:
    GetSameCollectionsMatchingLevel(jet1,jet2,d1,d2);
    break;
  default:
    ;
  }

  if (d1 > 0) {

    if (d1 < jet1->ClosestJetDistance()) {
      jet1->SetSecondClosestJet(jet1->ClosestJet(), jet1->ClosestJetDistance());
      jet1->SetClosestJet(jet2, d1);
    }
    else if (d1 < jet1->SecondClosestJetDistance()) {
      jet1->SetSecondClosestJet(jet2, d1);
    }
  }
  
  if (d2 > 0) {
    
    if (d2 < jet2->ClosestJetDistance()) {
      jet2->SetSecondClosestJet(jet2->ClosestJet(), jet2->ClosestJetDistance());
      jet2->SetClosestJet(jet1, d2);
    }
    else if (d2 < jet2->SecondClosestJetDistance()) {
      jet2->SetSecondClosestJet(jet1, d2);
    }
  }
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

    if (jet2->Pt() > fMaxBinPt)
      continue;

    if (!AcceptJet(jet2))
      continue;

    if (AcceptBiasJet(jet2) &&
	(jet2->Eta() > fJetMinEta && jet2->Eta() < fJetMaxEta && jet2->Phi() > fJetMinPhi && jet2->Phi() < fJetMaxPhi)) {
      
      fHistJets2PtAreaAcceptance->Fill(jet2->Area(), jet2->Pt());
      fHistJets2PhiEtaAcceptance->Fill(jet2->Eta(), jet2->Phi());
      
      if (!fRho2Name.IsNull())
	fHistJets2CorrPtAreaAcceptance->Fill(jet2->Area(), jet2->Pt() - fRho2Val * jet2->Area());
    }

    if (!AcceptBiasJet2(jet2))
      continue;

    if (jet2->Eta() < fJet2MinEta || jet2->Eta() > fJet2MaxEta || jet2->Phi() < fJet2MinPhi || jet2->Phi() > fJet2MaxPhi)
      continue;

    fHistJets2PtArea->Fill(jet2->Area(), jet2->Pt());
    fHistJets2PhiEta->Fill(jet2->Eta(), jet2->Phi());

    if (!fRho2Name.IsNull())
      fHistJets2CorrPtArea->Fill(jet2->Area(), jet2->Pt() - fRho2Val * jet2->Area());

    if (jet2->MatchedJet()) {

      if (!AcceptBiasJet(jet2->MatchedJet()) || 
	  jet2->MatchedJet()->MaxTrackPt() > fMaxTrackPt || jet2->MatchedJet()->MaxClusterPt() > fMaxClusterPt ||
	  jet2->MatchedJet()->Pt() > fMaxBinPt) {
	fHistMissedJets2PtArea->Fill(jet2->Area(), jet2->Pt());
      }
      else {
	Double_t d1=-1, d2=-1;
	if (jet2->GetMatchingType() == kGeometrical) {

	  if (fAreCollections2MC && !fAreCollections1MC)
	    GetMCLabelMatchingLevel(jet2->MatchedJet(), jet2, d1, d2);
	  else if ((fAreCollections1MC && fAreCollections2MC) || (!fAreCollections1MC && !fAreCollections2MC))
	    GetSameCollectionsMatchingLevel(jet2->MatchedJet(), jet2, d1, d2);

	  fHistDistancevsCommonEnergy->Fill(jet2->ClosestJetDistance(), d2);
	}
	else if (jet2->GetMatchingType() == kMCLabel || jet2->GetMatchingType() == kSameCollections) {
	  GetGeometricalMatchingLevel(jet2->MatchedJet(), jet2, d1);
	  fHistDistancevsCommonEnergy->Fill(d1, jet2->ClosestJetDistance());
	}
	  
	fHistMatchingLevelvsJet2Pt->Fill(jet2->ClosestJetDistance(), jet2->Pt());

	Double_t deta = jet2->MatchedJet()->Eta() - jet2->Eta();
	Double_t dphi = jet2->MatchedJet()->Phi() - jet2->Phi();
	fHistDeltaEtaPhivsJet2Pt->Fill(deta, dphi, jet2->Pt());

	Double_t dpt = jet2->MatchedJet()->Pt() - jet2->Pt();
	fHistDeltaPtvsJet2Pt->Fill(jet2->Pt(), dpt);
	fHistDeltaPtvsMatchingLevel->Fill(jet2->ClosestJetDistance(), dpt);

	fHistJet1PtvsJet2Pt->Fill(jet2->MatchedJet()->Pt(), jet2->Pt());
	
	if (!fRhoName.IsNull() || !fRho2Name.IsNull()) {
	  dpt -= fRhoVal * jet2->MatchedJet()->Area() - fRho2Val * jet2->Area();
	  fHistDeltaCorrPtvsJet2Pt->Fill(jet2->Pt(), dpt);
	  fHistDeltaCorrPtvsMatchingLevel->Fill(jet2->ClosestJetDistance(), dpt);
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
