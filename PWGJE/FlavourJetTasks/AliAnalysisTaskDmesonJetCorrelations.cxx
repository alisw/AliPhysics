/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//
//  Analysis Taks for D meson - jet correlation analysis
//
//-----------------------------------------------------------------------
// Author:
// Salvatore Aiola, Yale University, salvatore.aiola@cern.ch
//-----------------------------------------------------------------------

// Root
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <THnSparse.h>
#include <TParticle.h>
#include <TMath.h>

// Aliroot general
#include "AliLog.h"

// Aliroot HF
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliRDHFCutsDStartoKpipi.h"

// Aliroot EMCal jet framework
#include "AliEmcalJet.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliAnalysisTaskSEDmesonsFilterCJ.h"
#include "AliEmcalParticle.h"

#include "AliAnalysisTaskDmesonJetCorrelations.h"

ClassImp(AliAnalysisTaskDmesonJetCorrelations)

//_______________________________________________________________________________
AliAnalysisTaskDmesonJetCorrelations::AliAnalysisTaskDmesonJetCorrelations() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskDmesonJetCorrelations", kTRUE),
  fCuts(0),
  fCandidateType(kDstartoKpipi),
  fMinMass(0.),
  fMaxMass(1.),
  fNBinsMass(80),
  fMaxR(0.2),
  fShowPositionD(kTRUE),
  fShowInvMass(kFALSE),
  fShow2ProngInvMass(kFALSE),
  fShowDeltaInvMass(kTRUE),
  fShowSoftPionPt(kFALSE),
  fShowDmesonZ(kTRUE),
  fShowDeltaR(kTRUE),
  fShowDeltaEta(kFALSE),
  fShowDeltaPhi(kFALSE),
  fShowPositionJet(kTRUE),
  fShowLeadingPt(kFALSE),
  fShowJetArea(kFALSE),
  fShowJetConstituents(kFALSE),
  fShowMatchingLevel(kFALSE),
  fShowDaughterDistance(0),
  fInhibitTask(kFALSE),
  fMatchingType(kGeometricalMatching),
  fOnlyAcceptedJets(kTRUE),
  fOnlySingleMatches(kFALSE),
  fCheckTrackColl(kTRUE),
  fParticleLevel(kFALSE),
  fUseExchangeContainer(kTRUE),
  fAliEmcalParticleMode(kFALSE),
  fUseMCInfo(kFALSE),
  fAodEvent(0),
  fCandidateArray(0),
  fHistRejectionReason(0),
  fHistTracksNotInJetsPt(0),
  fHistTracksNotInJetsEtaPhi(0),
  fDmesons(0)
{
  // Default ctor

}

//_______________________________________________________________________________
AliAnalysisTaskDmesonJetCorrelations::AliAnalysisTaskDmesonJetCorrelations(const char* name, AliRDHFCuts* cuts, ECandidateType cand, Bool_t useExchCont) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fCuts(cuts),
  fCandidateType(cand),
  fMinMass(0.),
  fMaxMass(1.),
  fNBinsMass(80),
  fMaxR(0.2),
  fShowPositionD(kTRUE),
  fShowInvMass(kFALSE),
  fShow2ProngInvMass(kFALSE),
  fShowDeltaInvMass(kTRUE),
  fShowSoftPionPt(kFALSE),
  fShowDmesonZ(kTRUE),
  fShowDeltaR(kTRUE),
  fShowDeltaEta(kFALSE),
  fShowDeltaPhi(kFALSE),
  fShowPositionJet(kTRUE),
  fShowLeadingPt(kFALSE),
  fShowJetArea(kFALSE),
  fShowJetConstituents(kFALSE),
  fShowMatchingLevel(kFALSE),
  fShowDaughterDistance(0),
  fInhibitTask(kFALSE),
  fMatchingType(kGeometricalMatching),
  fOnlyAcceptedJets(kTRUE),
  fOnlySingleMatches(kFALSE),
  fCheckTrackColl(kTRUE),
  fParticleLevel(kFALSE),
  fUseExchangeContainer(useExchCont),
  fAliEmcalParticleMode(kFALSE),
  fUseMCInfo(kFALSE),
  fAodEvent(0),
  fCandidateArray(0),
  fHistRejectionReason(0),
  fHistTracksNotInJetsPt(0),
  fHistTracksNotInJetsEtaPhi(0),
  fDmesons(0)
{
  // Constructor.
  
  if (fUseExchangeContainer) DefineInput(1, TClonesArray::Class());
  DefineOutput(2, AliRDHFCuts::Class()); // my cuts

  if (fCandidateType == kDstartoKpipi) {
    SetMassLimits(0.50, 413);   // Set mass limits of the D*
    SetShowInvMass(kFALSE);
    SetShowDeltaInvMass(kTRUE);
    SetShowSoftPionPt(kTRUE);
  }
  else if (fCandidateType == kD0toKpi) {
    SetMassLimits(0.50, 421);   // Set mass limits of the D0
    SetShowInvMass(kTRUE);
    SetShowDeltaInvMass(kFALSE);
    SetShowSoftPionPt(kFALSE);
  }

  SetMakeGeneralHistograms(kTRUE);
}

//_______________________________________________________________________________
AliAnalysisTaskDmesonJetCorrelations::~AliAnalysisTaskDmesonJetCorrelations()
{
  // Destructor.
  
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::Init()
{
  // Initialization.
   
  Info("AliAnalysisTaskDmesonJetCorrelations::Init()", "Entering method");
   
  switch (fCandidateType) {
  case kD0toKpi: 
    {
      AliRDHFCutsD0toKpi* copyfCutsDzero = new AliRDHFCutsD0toKpi(*(static_cast<AliRDHFCutsD0toKpi*>(fCuts)));
      copyfCutsDzero->SetName("AnalysisCutsDzero");
      PostData(2, copyfCutsDzero);  // Post the data
    } break;
  case kDstartoKpipi: 
    {
      AliRDHFCutsDStartoKpipi* copyfCutsDstar = new AliRDHFCutsDStartoKpipi(*(static_cast<AliRDHFCutsDStartoKpipi*>(fCuts)));
      copyfCutsDstar->SetName("AnalysisCutsDStar");
      PostData(2, copyfCutsDstar); // Post the cuts
    } break;
  default:
    return;
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::UserCreateOutputObjects()
{ 
  // Creates the output containers.
  
  Info("UserCreateOutputObjects", "CreateOutputObjects of task %s", GetName());
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
   
  // define histograms
  // the TList fOutput is already defined in  AliAnalysisTaskEmcal::UserCreateOutputObjects()

  if (fParticleLevel) {
    if (fMatchingType == kDaughterConstituentMatching) {
      AliWarning("Daughter constituent matching not available at particle level. Switching to geometrical matching.");
      fMatchingType = kGeometricalMatching;
    }
    if (fShowDaughterDistance) {
      AliWarning("Daughter distance not implemented for particle level. Switching it off.");
      fShowDaughterDistance = kFALSE;
    }
  }

  TString histname;
  
  histname = "fHistRejectionReason";
  fHistRejectionReason = new TH2F(histname, histname, 32, 0, 32, 100, 0, 250);
  fHistRejectionReason->GetXaxis()->SetTitle("Rejection reason");
  fHistRejectionReason->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/#it{c})");
  fHistRejectionReason->GetZaxis()->SetTitle("counts");
  SetRejectionReasonLabels(fHistRejectionReason->GetXaxis());
  fOutput->Add(fHistRejectionReason);

  if (fCheckTrackColl) {
    histname = "fHistTracksNotInJetsPt";
    fHistTracksNotInJetsPt = new TH1F(histname, histname, 100, 0, 50);
    fHistTracksNotInJetsPt->GetXaxis()->SetTitle("#it{p}_{T,track} (GeV/#it{c})");
    fHistTracksNotInJetsPt->GetYaxis()->SetTitle("counts");
    fOutput->Add(fHistTracksNotInJetsPt);

    histname = "fHistTracksNotInJetsEtaPhi";
    fHistTracksNotInJetsEtaPhi = new TH2F(histname, histname, 50, -1, 1, 150, 0, TMath::TwoPi());
    fHistTracksNotInJetsEtaPhi->GetXaxis()->SetTitle("#eta_{track}");
    fHistTracksNotInJetsEtaPhi->GetYaxis()->SetTitle("#phi_{track}");
    fHistTracksNotInJetsEtaPhi->GetZaxis()->SetTitle("counts");
    fOutput->Add(fHistTracksNotInJetsEtaPhi);
  }

  AllocateTHnSparse();
   
  PostData(1, fOutput);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDmesonJetCorrelations::IsEventSelected()
{
  // Return true if the event is selected.

  if (!fAodEvent) return kTRUE;
  
  Bool_t iseventselected = kFALSE;

  if (fCuts) iseventselected = fCuts->IsEventSelected(fAodEvent);
  
  if (!iseventselected) {
    if (fGeneralHistograms) fHistEventRejection->Fill("HFevSel",1);
    return kFALSE;
  }
  
  return AliAnalysisTaskEmcalJet::IsEventSelected();
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskDmesonJetCorrelations::FillHistograms()
{
  // Fill the histograms.

  if (fMatchingType != kJetLoop) return kTRUE;

  AliJetContainer* jets = GetJetContainer(0);
  
  Int_t ftag = AliEmcalJet::kD0;
  if (fCandidateType == kDstartoKpipi) ftag = AliEmcalJet::kDStar;
  
  jets->ResetCurrentID();
  AliEmcalJet* jet = 0;
  while ((jet = jets->GetNextJet())) {
    if (!jet->TestFlavourTag(ftag)) continue;
    AliVParticle* HFcand = 0;
    Int_t itrk = 0;
    while ((HFcand = jet->GetFlavourTrack(itrk))) {
      FillHistograms(HFcand, jet, kSingleMatch, 0.);
      itrk++;
    }
  }

  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskDmesonJetCorrelations::ExtractHFcandAttributes(AliVParticle* HFcand, TLorentzVector& Dvector, Double_t& invMassD, Double_t& softPionPtD, Double_t& invMass2prong, UInt_t i)
{
  if (fParticleLevel) {
    AliAODMCParticle* part = static_cast<AliAODMCParticle*>(HFcand);
    return ExtractParticleLevelHFAttributes(part, Dvector, invMassD, softPionPtD, invMass2prong, i);
  }
  else {
    AliAODRecoDecayHF2Prong* Dcand = static_cast<AliAODRecoDecayHF2Prong*>(HFcand);
    return ExtractRecoDecayAttributes(Dcand, Dvector, invMassD, softPionPtD, invMass2prong, i);
  }
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskDmesonJetCorrelations::ExtractParticleLevelHFAttributes(AliAODMCParticle* part, TLorentzVector& Dvector, Double_t& invMassD, Double_t& /*softPionPtD*/, Double_t& /*invMass2prong*/, UInt_t i)
{
  if (i > 0) return kFALSE;
  Dvector.SetPtEtaPhiM(part->Pt(), part->Eta(), part->Phi(), part->M());
  invMassD = part->M();
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskDmesonJetCorrelations::ExtractRecoDecayAttributes(AliAODRecoDecayHF2Prong* Dcand, TLorentzVector& Dvector, Double_t& invMassD, Double_t& softPionPtD, Double_t& invMass2prong, UInt_t i)
{
  if (fCandidateType == kD0toKpi) {
    Int_t MCtruthPdgCode = 0;

    if (fUseMCInfo) {
      Int_t pdg = 0;
      if (fCandidateType == kD0toKpi) {
        pdg = 421;
      }
      else {
        pdg = 413;
      }

      AliParticleContainer* mcpartCont = GetParticleContainer(1);
      TClonesArray* mcpartArray = mcpartCont->GetArray();
      Int_t mcLab = Dcand->MatchToMC(pdg, mcpartArray);
      AliDebug(2, Form("MC label is %d", mcLab));
      if (mcLab >= 0) {
        AliAODMCParticle* aodMcPart = static_cast<AliAODMCParticle*>(mcpartCont->GetParticle(mcLab));
        MCtruthPdgCode = aodMcPart->PdgCode();
        AliDebug(2, Form("MC truth pdg code is %d",MCtruthPdgCode));
      }
    }

    AliDebug(2,"Checking if D0 meson is selected");
    Int_t isSelected = fCuts->IsSelected(Dcand, AliRDHFCuts::kAll, fAodEvent);
    if (isSelected == 1 && (MCtruthPdgCode == 0 || MCtruthPdgCode == 421)) {  // selected as a D0
      if (i > 0) return kFALSE;
      AliDebug(2,"Selected as D0");
      invMassD = Dcand->InvMassD0();
    }
    else if (isSelected == 2 && (MCtruthPdgCode == 0 || MCtruthPdgCode == -421)) { // selected as a D0bar
      if (i > 0) return kFALSE;
      AliDebug(2,"Selected as D0bar");
      invMassD = Dcand->InvMassD0bar();
    }
    else if (isSelected == 3) { // selected as either a D0bar or a D0 (PID on K and pi undecisive)
      AliDebug(2,"Selected as either D0 or D0bar");

      if (MCtruthPdgCode == 0) {
        // Select D0 or D0bar depending on the i-parameter
        if (i == 0) {
          AliDebug(2, "Returning invariant mass with D0 hypothesis");
          invMassD = Dcand->InvMassD0();
        }
        else if (i == 1) {
          AliDebug(2, "Returning invariant mass with D0bar hypothesis");
          invMassD = Dcand->InvMassD0bar();
        }
        else {
          return kFALSE;
        }
      }
      else {
        if (i > 0) return kFALSE;
        if (MCtruthPdgCode == 421) {
          AliDebug(2, "MC truth is D0");
          invMassD = Dcand->InvMassD0();
        }
        else if (MCtruthPdgCode == -421) {
          AliDebug(2, "MC truth is D0bar");
          invMassD = Dcand->InvMassD0bar();
        }
        else {
          return kFALSE;
        }
      }
    }
  }
  else if (fCandidateType == kDstartoKpipi) {
    if (i > 0) return kFALSE;
    AliAODRecoCascadeHF* DstarCand = static_cast<AliAODRecoCascadeHF*>(Dcand);
    invMassD = DstarCand->InvMassDstarKpipi();
    softPionPtD = DstarCand->GetBachelor()->Pt();
    invMass2prong = DstarCand->InvMassD0();
  }

  Dvector.SetPtEtaPhiM(Dcand->Pt(), Dcand->Eta(), Dcand->Phi(), invMassD);
  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskDmesonJetCorrelations::Run()
{
  // Run the analysis.
  
  if (fMatchingType == kJetLoop) {
    DoJetLoop();
  }
  else {
    DoDmesonLoop();
  }

  return kTRUE;
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::DoJetLoop()
{
  // Run the loop over the jets.

  AliJetContainer* jets = GetJetContainer(0);
  AliParticleContainer* particles = jets->GetParticleContainer();
  
  Int_t ftag = AliEmcalJet::kD0;
  Int_t pdgCode = 421;
  TString recoDecayClassName("AliAODRecoDecayHF2Prong");
  if (fCandidateType == kDstartoKpipi) {
    ftag = AliEmcalJet::kDStar;
    pdgCode = 413;
    recoDecayClassName = "AliAODRecoCascadeHF";
  }

  TClass* recoDecayClass = TClass::GetClass(recoDecayClassName);
  
  jets->ResetCurrentID();
  AliEmcalJet* jet = 0;
  while ((jet = jets->GetNextJet())) {
    Int_t ntrk = jet->GetNumberOfTracks();
    AliDebug(2, Form("Jet n %d, pT = %.2f, eta = %.2f, phi = %.2f, nTrk = %d", jets->GetCurrentID(), jet->Pt(), jet->Eta(), jet->Phi(), ntrk));
    for (Int_t itrk = 0; itrk < ntrk; itrk++) {
      AliVParticle* part = jet->TrackAt(itrk, particles->GetArray());
      if (!part) continue;
      Bool_t matched = kFALSE;

      if (fParticleLevel) {
        Int_t partPdg = TMath::Abs(part->PdgCode());
        if (partPdg == pdgCode) matched = kTRUE;
      }
      else {
        AliEmcalParticle* emcpart = dynamic_cast<AliEmcalParticle*>(part);
        if (emcpart) part = emcpart->GetTrack();
        if (part->InheritsFrom(recoDecayClass)) matched = kTRUE;
      }
      
      if (matched) {
        jet->AddFlavourTag(ftag);
        jet->AddFlavourTrack(part);
      }

      AliDebug(2, Form("%s: pT = %.2f, eta = %.2f, phi = %.2f", part->ClassName(), part->Pt(), part->Eta(), part->Phi()));
    }
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::DoDmesonLoop()
{
  // Run the loop over the D meson candidates.

  const Int_t nDcand = fCandidateArray->GetEntriesFast();

  Int_t pdgCode = 421;
  TString recoDecayClassName("AliAODRecoDecayHF2Prong");
  if (fCandidateType == kDstartoKpipi) {
    pdgCode = 413;
    recoDecayClassName = "AliAODRecoCascadeHF";
  }

  TArrayD matchingLevel(5);
  TList matchedJets;
  matchedJets.SetOwner(kFALSE);

  Int_t ftag = AliEmcalJet::kD0;
  if (fCandidateType == kDstartoKpipi) ftag = AliEmcalJet::kDStar;

  AliDebug(2,"Starting D meson candidate loop");
  for (Int_t icand = 0; icand < nDcand; icand++) {
    AliDebug(2,Form("D meson candidate %d", icand));
    AliVParticle* vpart = static_cast<AliVParticle*>(fCandidateArray->At(icand));
    if (!vpart) continue;

    AliVParticle* HFcand = 0;
    if (fAliEmcalParticleMode) {
      AliEmcalParticle* emcpart = static_cast<AliEmcalParticle*>(vpart);
      HFcand = emcpart->GetTrack();
    }
    else {
      HFcand = vpart;
    }
    if (!HFcand) continue;

    if (!fUseExchangeContainer) {
      if (fParticleLevel) {
        AliAODMCParticle* part = static_cast<AliAODMCParticle*>(HFcand);
        if (TMath::Abs(part->PdgCode()) != pdgCode) continue;
      }
      else {
        if (!HFcand->InheritsFrom(recoDecayClassName)) continue;
      }
    }

    // Look for D-jet correlation
    Int_t n = FindMatchedJet(fMatchingType, HFcand, matchingLevel, matchedJets);

    AliEmcalJet* jet = 0;
    Int_t matchingStatus = kNotMatched;
 
    if (n == 1) {
      matchingStatus = kSingleMatch;
    }
    else if (n > 1) {
      matchingStatus = kMultipleMatches;
    }

    if (matchingStatus == kSingleMatch || (matchingStatus == kMultipleMatches && !fOnlySingleMatches)) {
      jet = static_cast<AliEmcalJet*>(matchedJets.At(0));
      jet->AddFlavourTrack(HFcand);
      jet->AddFlavourTag(ftag);
    }
    
    FillHistograms(HFcand, jet, matchingStatus, matchingLevel[0]);
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::FillHistograms(AliVParticle* HFcand, AliEmcalJet* jet, Int_t matchingStatus, Double_t matchingLevel)
{
  // Fill the histograms.

  // Extracting D meson attributes
  TLorentzVector Dvector;
  Double_t invMassD = 0;
  Double_t softPionPtD = 0;
  Double_t invMass2prong = 0;

  TLorentzVector jetVector;
  Double_t leadPtJet = 0;
  Double_t areaJet = 0;
  Int_t constJet = 0;

  Double_t daughterDist[5] = {0.};

  // Extracting jet attributes (if a jet was found)
  if (jet) {
    AliJetContainer* jets = GetJetContainer(0);
    
    // Check if the jet passes the cuts
    Bool_t acceptedJet = jets->AcceptJet(jet);
    if (!acceptedJet) {
      fHistRejectionReason->Fill(jets->GetRejectionReasonBitPosition(), jet->Pt());
      matchingStatus = kJetNotAccepted;
    }
    if (acceptedJet || !fOnlyAcceptedJets) {
      jetVector.SetPtEtaPhiM(jet->Pt(), jet->Eta(), jet->Phi(), 0);
        
      leadPtJet = jet->MaxPartPt();
      areaJet = jet->Area();
      constJet = jet->N();

      if (fShowDaughterDistance > 0 && !fParticleLevel) {
        AliAODRecoDecayHF2Prong* Dcand = static_cast<AliAODRecoDecayHF2Prong*>(HFcand);
        TObjArray daughters(5);
        AliAnalysisTaskSEDmesonsFilterCJ::AddDaughters(Dcand, daughters);
        for (Int_t i = 0; i < daughters.GetEntriesFast(); i++) {
          AliVTrack* track = static_cast<AliVTrack*>(daughters.At(i));
          if (!track) continue;
          daughterDist[i] = jet->DeltaR(track);
        }
      }
    }
  }

  if (ExtractHFcandAttributes(HFcand, Dvector, invMassD, softPionPtD, invMass2prong, 0)) {
    FillTHnSparse(Dvector, softPionPtD, invMass2prong, jetVector, leadPtJet, areaJet, constJet, matchingStatus, matchingLevel, daughterDist);
  }

  if (ExtractHFcandAttributes(HFcand, Dvector, invMassD, softPionPtD, invMass2prong, 1)) {
    FillTHnSparse(Dvector, softPionPtD, invMass2prong, jetVector, leadPtJet, areaJet, constJet, matchingStatus, matchingLevel, daughterDist);
  }
}

//_______________________________________________________________________________
Int_t AliAnalysisTaskDmesonJetCorrelations::FindMatchedJet(EMatchingType matchType, AliVParticle* cand, TArrayD& matchingLevel, TList& matchedJets)
{
  // Find jet matched to a reconstructed decay candidate, using a constituent-based algorithm.
  // The returned value is the number of jets found that share the decay products.

  Int_t nj = 0;

  matchingLevel.Reset(1);
  matchedJets.Clear();
  
  AliJetContainer* jets = GetJetContainer(0);
  AliEmcalJet* jet = 0;

  AliDebug(2, Form("D candidate pt = %.3f eta = %.3f phi = %.3f", cand->Pt(), cand->Eta(), cand->Phi()));

  jets->ResetCurrentID();
  Bool_t reset = kTRUE;
  while ((jet = jets->GetNextJet())) {
    Double_t m = CalculateMatchingLevel(matchType, cand, jet, reset);
    reset = kFALSE;
    
    if (m >= 1) continue; // jet is not matched

    if (m < 0) break; // this is the signal from CalculateMatchingLevel that the search is completed

    if (nj >= matchingLevel.GetSize()) {
      matchingLevel.Set(nj*2+1);
    }

    Int_t pos = nj;
    
    while (pos >= 0) {
      if (pos == 0 || matchingLevel[pos-1] < m) {
        matchingLevel[pos] = m;
        break;
      }
      matchingLevel[pos] = matchingLevel[pos-1];
      pos--;
    }
    
    matchedJets.AddAt(jet, pos);
    nj++;
  }

  AliDebug(2, Form("Matching levels\n%.3f\n%.3f\n%.3f", matchingLevel[0], matchingLevel[1], matchingLevel[2]));
  AliDebug(2, Form("Final number of matched jets %d",nj));
    
  return nj;
}

//_______________________________________________________________________________
Double_t AliAnalysisTaskDmesonJetCorrelations::CalculateMatchingLevel(EMatchingType matchType, AliVParticle* cand, AliEmcalJet* jet, Bool_t reset)
{
  // Calculate the matchin level between cand and jet.
  
  Double_t m = 1;
  
  if (matchType == kGeometricalMatching) {    
    m = CalculateGeometricalMatchingLevel(cand, jet);
  }
  else if (matchType == kDaughterConstituentMatching) {
    AliAODRecoDecay* recoDecay = static_cast<AliAODRecoDecay*>(cand);
    m = CalculateDaughterConstituentMatchingLevel(recoDecay, jet, reset);
  }
  else if (matchType == kCandidateConstituentMatching) {
    m = CalculateCandidateConstituentMatchingLevel(cand, jet);
  }
  else {
    AliWarning(Form("Matching algorithm type %d not implemented!", matchType));  
  }

  return m;
}

//_______________________________________________________________________________
Double_t AliAnalysisTaskDmesonJetCorrelations::CalculateGeometricalMatchingLevel(AliVParticle* cand, AliEmcalJet* jet)
{
  // Calculate the matching level using a geometrical algorithm.
  // The matching level is the inverse of the distance.

  Double_t deltaR = jet->DeltaR(cand);
  Double_t m = deltaR / fMaxR;
  if (m < 0) m = 1;
  
  return m;
}

//_______________________________________________________________________________
Double_t AliAnalysisTaskDmesonJetCorrelations::CalculateDaughterConstituentMatchingLevel(AliAODRecoDecay* cand, AliEmcalJet* jet, Bool_t reset)
{
  // Calculate the matching level using a constituent-based algorithm.
  // If the jet is not provided but reset==kTRUE then it only initializes the daughters list.
  // The matching level is the sum of the pt of the decay products contained in the jet.

  if (!cand) return 1;

  static AliAODRecoDecay* prevCand = 0;
  static TObjArray daughters(5);
  static Double_t ptTot = 0;

  if (prevCand != cand || reset) {
    daughters.Clear();
    ptTot = AliAnalysisTaskSEDmesonsFilterCJ::AddDaughters(cand, daughters);
    prevCand = cand;

    if (fCheckTrackColl) {
      AliParticleContainer* tracks = GetParticleContainer(0);
      TClonesArray* trackArray = 0;
      if (tracks) trackArray = tracks->GetArray();

      if (trackArray) {
        AliAODTrack* daughter = 0;
        TIter next(&daughters);
        while ((daughter = static_cast<AliAODTrack*>(next()))) {
          if (!trackArray->Contains(daughter) || !daughter->IsHybridGlobalConstrainedGlobal()) {
            fHistTracksNotInJetsPt->Fill(daughter->Pt());
            fHistTracksNotInJetsEtaPhi->Fill(daughter->Eta(), daughter->Phi());
          }
        }
      }
    }
  }

  if (!jet) return 1;

  Int_t nt = daughters.GetEntriesFast();
  if (nt == 0 || ptTot == 0) return -1;  // this will end the search for matched jets

  AliParticleContainer* tracks = GetParticleContainer(0);
  if (!tracks) return 1;

  // To save computation time, if the jet is too far from the candidate, let's skip it
  Double_t mlg = CalculateGeometricalMatchingLevel(cand, jet);
  if (mlg >= 1.) return 1.;

  Double_t pt = 0;

  for (Int_t it = 0; it < nt; it++) {
    AliVTrack* track = static_cast<AliVTrack*>(daughters.At(it));
    if (!track) continue;
    // Check if the jet contains the track
    if (jet->ContainsTrack(track, tracks->GetArray()) >= 0) {
      pt += track->Pt();
      AliDebug(2, Form("Jet pt = %.3f eta = %.3f phi = %.3f contains daughter pT = %.3f eta = %.3f phi = %.3f", jet->Pt(), jet->Eta(), jet->Phi(), track->Pt(), track->Eta(), track->Phi()));
      daughters.RemoveAt(it);
    }
  }

  Double_t m = 1 - pt / ptTot;
  
  return m;
}

//_______________________________________________________________________________
Double_t AliAnalysisTaskDmesonJetCorrelations::CalculateCandidateConstituentMatchingLevel(AliVParticle* cand, AliEmcalJet* jet)
{
  // Calculate the matching level using a constituent-based algorithm.
  // If the jet is not provided but reset==kTRUE then it only initializes the daughters list.
  // The matching level is the sum of the pt of the decay products contained in the jet.
  
  if (!cand) return 1.;
  if (!jet) return 1.;

  Bool_t emcalPart = kFALSE;

  AliParticleContainer* tracks = GetParticleContainer(0);
  if (!tracks) return 1;

  TClonesArray* trkArray = tracks->GetArray();
  if (!trkArray) return 1;

  if (strcmp(trkArray->GetClass()->GetName(), "AliEmcalParticle") == 0) emcalPart = kTRUE;

  // To save computation time, if the jet is too far from the candidate, let's skip it
  Double_t mlg = CalculateGeometricalMatchingLevel(cand, jet);
  if (mlg >= 1.) return 1.;
  
  Int_t nt = jet->GetNumberOfTracks();

  for (Int_t it = 0; it < nt; it++) {
    AliVParticle* part = jet->TrackAt(it, trkArray);
    if (emcalPart) {
      AliEmcalParticle* emcPart = static_cast<AliEmcalParticle*>(part);
      part = emcPart->GetTrack();
    }

    if (part == cand) return 0.;
  }

  return 1.;
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::ExecOnce()
{
  // Execute only once (for the first event).

  if (fInhibitTask) return;

  if (!fParticleLevel) {
    if (!fCuts) {
      AliError(Form("%s: Cuts not provided. Task will not run.", GetName()));
      fInhibitTask = kTRUE;
      return;
    }

    // Load the event
    fAodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);

    if (!fAodEvent) {
      if (AODEvent() && IsStandardAOD()) {

        // In case there is an AOD handler writing a standard AOD, use the AOD
        // event in memory rather than the input (ESD) event.
        fAodEvent = dynamic_cast<AliAODEvent*>(AODEvent());
      }

    }

    if (!fAodEvent) {
      AliError(Form("This task need an AOD event! Task '%s' will be disabled!", GetName()));
      fInhibitTask = kTRUE;
      return;
    }
  }

  if (fAliEmcalParticleMode && fUseExchangeContainer) {
    AliError(Form("%s: AliEmcalParticle mode is incompatible with using the exchange container.", GetName()));
    fInhibitTask = kTRUE;
    return;
  }

  if (fParticleLevel && fUseMCInfo) {
    fUseMCInfo = kFALSE; // no point in specifing "use MC info" if the analysis is performed at particle level
  }

  if (fUseMCInfo) {
    AliParticleContainer* mcPartCont = GetParticleContainer(1);
    if (!mcPartCont) {
      AliError(Form("%s: Could not find the MC particle container at position 1. The MC info won't be used in this analysis!", GetName()));
      fUseMCInfo = kFALSE;
    }
    else {
      mcPartCont->SetClassName("AliAODMCParticle");
    }
  }

  TString className;
  if (fAliEmcalParticleMode) {
    className = "AliEmcalParticle";
  }
  else {
    if (fParticleLevel) {
      className = "AliAODMCParticle";
    }
    else {
      if (fCandidateType == kD0toKpi) {
        className = "AliAODRecoDecayHF2Prong";
      }
      else if (fCandidateType == kDstartoKpipi) {
        className = "AliAODRecoCascadeHF";
      }
      else {
        AliError(Form("%s: Candidate type %d not recognized!",
                      GetName(), (Int_t)fCandidateType));
        fCandidateArray = 0;
        fInhibitTask = kTRUE;
        return;
      }
    }
  }

  AliAnalysisTaskEmcalJet::ExecOnce();

  if (fUseExchangeContainer) {
    fCandidateArray = dynamic_cast<TClonesArray*>(GetInputData(1));
  }
  else {
    AliParticleContainer* partCont = GetParticleContainer(0);
    if (!partCont) {
      AliError(Form("%s: Unable to find the candidate array in input slot 1.", GetName()));
      fInhibitTask = kTRUE;
      return;
    }
    fCandidateArray = partCont->GetArray();
  }

  if (fCandidateArray) {
    TString objname(fCandidateArray->GetClass()->GetName());
    TClass cls(objname);

    if (!cls.InheritsFrom(className)) {
      AliError(Form("%s: Objects of type %s in %s are not inherited from %s! Task will not run.",
                    GetName(), cls.GetName(), fCandidateArray->GetName(), className.Data()));
      fCandidateArray = 0;
      fInhibitTask = kTRUE;
      return;
    }
  }
  else {
    AliError(Form("%s: Unable to find the candidate array.", GetName()));
    fInhibitTask = kTRUE;
    return;
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::Terminate(Option_t*)
{    
  // The Terminate() function is the last function to be called during
  // a query. It always runs on the client, it can be used to present
  // the results graphically or save the results to file.
   
  Info("Terminate"," terminate");
  AliAnalysisTaskEmcalJet::Terminate();
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::AllocateTHnSparse()
{
  // Allocate the THnSparse histogram.

  TString title[20]= {""};
  Int_t nbins[20]  = {0};
  Double_t min[20] = {0.};
  Double_t max[20] = {0.};
  Int_t dim = 0;

  title[dim] = "#it{p}_{T,D} (GeV/#it{c})";
  nbins[dim] = fNbins;
  min[dim] = 0;
  max[dim] = 100;
  dim++;

  if (fShowPositionD) {
    title[dim] = "#eta_{D}";
    nbins[dim] = 50;
    min[dim] = -1;
    max[dim] = 1;
    dim++;

    title[dim] = "#phi_{D} (rad)";
    nbins[dim] = 150;
    min[dim] = 0;
    max[dim] = TMath::TwoPi();
    dim++;
  }

  if (fShowInvMass) {
    title[dim] = "#it{M}_{D} (GeV/#it{c}^{2})";
    nbins[dim] = fNBinsMass;
    min[dim] = fMinMass;
    max[dim] = fMaxMass;
    dim++;
  }

  if (fShow2ProngInvMass) {
    title[dim] = "#it{M}_{2-prong} (GeV/#it{c}^{2})";
    nbins[dim] = fNBinsMass;
    CalculateMassLimits(fMaxMass - fMinMass, 421, fNBinsMass, min[dim], max[dim]);
    dim++;
  }

  if (fShowDeltaInvMass) {
    title[dim] = "#it{M}_{D*} - #it{M}_{D_{0}} (GeV/#it{c}^{2})";
    nbins[dim] = fNBinsMass*6;
    CalculateMassLimits(0.20, 413, fNBinsMass, min[dim], max[dim]);

    // subtract mass of D0
    Double_t D0mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(421))->Mass();
    min[dim] -= D0mass;  
    max[dim] -= D0mass;
    
    dim++;
  }

  if (fShowSoftPionPt) {
    title[dim] = "#it{p}_{T,#pi} (GeV/#it{c})";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 25;
    dim++;
  }

  if (fShowDmesonZ) {
    title[dim] = "#it{z}_{D}";
    nbins[dim] = 200;
    min[dim] = 0;
    max[dim] = 2.0;
    dim++;
  }

  if (fShowDeltaR) {
    title[dim] = "#Delta R_{D-jet}";
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = fMaxR * 1.1;
    dim++;
  }

  if (fShowDeltaEta) {
    title[dim] = "#eta_{D} - #eta_{jet}";
    nbins[dim] = 100;
    min[dim] = -fMaxR * 1.1;
    max[dim] = fMaxR * 1.1;
    dim++;
  }
  
  if (fShowDeltaPhi) {
    title[dim] = "#phi_{D} - #phi_{jet} (rad)";
    nbins[dim] = 100;
    min[dim] = -fMaxR * 1.1;
    max[dim] = fMaxR * 1.1;
    dim++;
  }

  title[dim] = "#it{p}_{T,jet} (GeV/#it{c})";
  nbins[dim] = fNbins;
  min[dim] = fMinBinPt;
  max[dim] = fMaxBinPt;
  dim++;

  if (fShowPositionJet) {
    title[dim] = "#eta_{jet}";
    nbins[dim] = 50;
    min[dim] = -1;
    max[dim] = 1;
    dim++;

    title[dim] = "#phi_{jet} (rad)";
    nbins[dim] = 150;
    min[dim] = 0;
    max[dim] = TMath::TwoPi();
    dim++;
  }

  if (fShowLeadingPt) {
    title[dim] = "#it{p}_{T,particle}^{leading} (GeV/#it{c})";
    nbins[dim] = 120;
    min[dim] = 0;
    max[dim] = 120;
    dim++;
  }

  if (fShowJetArea) {
    title[dim] = "#it{A}_{jet}";
    nbins[dim] = 150;
    min[dim] = 0;
    max[dim] = 1.5;
    dim++;
  }

  if (fShowJetConstituents) {
    title[dim] = "No. of constituents";
    nbins[dim] = 250;
    min[dim] = -0.5;
    max[dim] = 249.5;
    dim++;
  }

  title[dim] = "Matching status";
  nbins[dim] = 4;
  min[dim] = 0;
  max[dim] = 4;
  dim++;

  if (fShowMatchingLevel) {
    title[dim] = "Matching level";
    nbins[dim] = 50;
    min[dim] = 0;
    max[dim] = 1;
    dim++;
  }

  for (Int_t i = 0; i < fShowDaughterDistance; i++) {
    title[dim] = Form("#Delta R_{d%d-jet}", i);
    nbins[dim] = 100;
    min[dim] = 0;
    max[dim] = 4;
    dim++;
  }
  
  fDmesons = new THnSparseD("fDmesons","fDmesons",dim,nbins,min,max);
  fOutput->Add(fDmesons);
  for (Int_t i = 0; i < dim; i++) {
    fDmesons->GetAxis(i)->SetTitle(title[i]);
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::FillTHnSparse(TLorentzVector D, Double_t softPionPtD, Double_t invMass2prong,
                                                         TLorentzVector jet, Double_t leadPtJet, Double_t areaJet, Int_t constJet, Int_t matchingStatus, Double_t matchingLevel, Double_t daughterDist[5])
{
  // Fill the THnSparse histogram.

  Double_t contents[20] = {0.};

  Double_t z = 1.;
  Double_t deltaR = 1.;
  Double_t deltaPhi = 1.;
  Double_t deltaEta = 1.;

  if (jet.Pt() > 0) {
    TVector3 dvect = D.Vect();
    TVector3 jvect = jet.Vect();
    z = (dvect * jvect) / (jvect * jvect);
    
    if (z == 1 || (z > 1 && z - 1 < 1e-3)) z = 0.999; // so that it will contribute to the bin 0.9-1 rather than 1-1.1
    
    deltaPhi = TVector2::Phi_mpi_pi(D.Phi() - jet.Phi());
    deltaEta = D.Eta() - jet.Eta();
    
    deltaR = TMath::Sqrt(deltaPhi*deltaPhi + deltaEta*deltaEta);
  }
  
  for (Int_t i = 0; i < fDmesons->GetNdimensions(); i++) {
    TString title(fDmesons->GetAxis(i)->GetTitle());
    if      (title=="#it{p}_{T,D} (GeV/#it{c})")                     contents[i] = D.Pt();
    else if (title=="#eta_{D}")                                      contents[i] = D.Eta();
    else if (title=="#phi_{D} (rad)")                                contents[i] = D.Phi() < 0 ? D.Phi()+TMath::TwoPi() : D.Phi();
    else if (title=="#it{M}_{D} (GeV/#it{c}^{2})")                   contents[i] = D.M();
    else if (title=="#it{M}_{2-prong} (GeV/#it{c}^{2})")             contents[i] = invMass2prong;
    else if (title=="#it{M}_{D*} - #it{M}_{D_{0}} (GeV/#it{c}^{2})") contents[i] = D.M() - invMass2prong;
    else if (title=="#it{p}_{T,#pi} (GeV/#it{c})")                   contents[i] = softPionPtD;
    else if (title=="#it{z}_{D}")                                    contents[i] = z;
    else if (title=="#Delta R_{D-jet}")                              contents[i] = deltaR;
    else if (title=="#eta_{D} - #eta_{jet}")                         contents[i] = deltaEta;
    else if (title=="#phi_{D} - #phi_{jet} (rad)")                   contents[i] = deltaPhi;
    else if (title=="#it{p}_{T,jet} (GeV/#it{c})")                   contents[i] = jet.Pt();
    else if (title=="#eta_{jet}")                                    contents[i] = jet.Eta();
    else if (title=="#phi_{jet} (rad)")                              contents[i] = jet.Phi() < 0 ? jet.Phi()+TMath::TwoPi() : jet.Phi();
    else if (title=="#it{p}_{T,particle}^{leading} (GeV/#it{c})")    contents[i] = leadPtJet;
    else if (title=="#it{A}_{jet}")                                  contents[i] = areaJet;
    else if (title=="No. of constituents")                           contents[i] = constJet;
    else if (title=="Matching status")                               contents[i] = matchingStatus;
    else if (title=="Matching level")                                contents[i] = matchingLevel;
    else if (title=="#Delta R_{d0-jet}")                             contents[i] = daughterDist[0];
    else if (title=="#Delta R_{d1-jet}")                             contents[i] = daughterDist[1];
    else if (title=="#Delta R_{d2-jet}")                             contents[i] = daughterDist[2];
    else AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }

  fDmesons->Fill(contents);
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::CalculateMassLimits(Double_t range, Int_t pdg, Int_t nbins, Double_t& minMass, Double_t& maxMass)
{
  // Set the mass limits for the histograms using information from TDatabasePDG.

  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg));
  
  Float_t mass = part->Mass();

  // To make sure that the PDG mass value is not at the edge of a bin
  if (nbins % 2 == 0) {
    minMass = mass - range / 2 - range / nbins / 2;
    maxMass = mass + range / 2 - range / nbins / 2;
  }
  else {
    minMass = mass - range / 2;
    maxMass = mass + range / 2;
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::SetMassLimits(Double_t range, Int_t pdg)
{
  // Set the mass limits for the histograms using information from TDatabasePDG.

  CalculateMassLimits(range, pdg, fNBinsMass, fMinMass, fMaxMass);

  AliInfo(Form("Mass limits set for particle %d.", pdg));
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::SetMassLimits(Double_t lowlimit, Double_t uplimit)
{
  // Set the mass limits for the histograms.
  
  if (lowlimit < 0.) {
    lowlimit = 0.;
  }

  if (uplimit <= lowlimit) {
    AliError(Form("%s: wrong mass limits!", GetName()));
    return;
  }

  AliInfo(Form("Setting mass limits to %f, %f", lowlimit, uplimit));
  
  fMinMass = lowlimit;
  fMaxMass = uplimit;
}
