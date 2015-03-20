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

#include "AliAnalysisTaskDmesonJetCorrelations.h"

ClassImp(AliAnalysisTaskDmesonJetCorrelations)

//_______________________________________________________________________________
AliAnalysisTaskDmesonJetCorrelations::AliAnalysisTaskDmesonJetCorrelations() :
  AliAnalysisTaskEmcalJet("AliAnalysisTaskDmesonJetCorrelations", kTRUE),
  fCuts(0),
  fCandidateType(kDstartoKpipi),
  fMinMass(0.),
  fMaxMass(1.),
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
  fMinDeltaPhiHisto(-TMath::Pi()),
  fInhibitTask(kFALSE),
  fAodEvent(0),
  fCandidateArray(0),
  fHistRejectionReason(0),
  fDmesons(0)
{
  // Default ctor

}

//_______________________________________________________________________________
AliAnalysisTaskDmesonJetCorrelations::AliAnalysisTaskDmesonJetCorrelations(const char* name, AliRDHFCuts* cuts, ECandidateType cand) :
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fCuts(cuts),
  fCandidateType(cand),
  fMinMass(0.),
  fMaxMass(1.),
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
  fMinDeltaPhiHisto(-TMath::Pi()),
  fInhibitTask(kFALSE),
  fAodEvent(0),
  fCandidateArray(0),
  fHistRejectionReason(0),
  fDmesons(0)
{
  // Constructor.
  
  DefineInput(1, TClonesArray::Class());
  DefineInput(2, TClonesArray::Class());
  DefineOutput(2, AliRDHFCuts::Class()); // my cuts

  if (fCandidateType == kDstartoKpipi) {
    SetMassLimits(0.04, 413);   // Set mass limits of the D*
    SetShowInvMass(kFALSE);
    SetShowDeltaInvMass(kTRUE);
    SetShowSoftPionPt(kTRUE);
  }
  else if (fCandidateType == kD0toKpi) {
    SetMassLimits(0.04, 421);   // Set mass limits of the D0
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

  TString histname("fHistRejectionReason");
  fHistRejectionReason = new TH2F(histname, histname, 32, 0, 32, 100, 0, 250);
  fHistRejectionReason->GetXaxis()->SetTitle("Rejection reason");
  fHistRejectionReason->GetYaxis()->SetTitle("#it{p}_{T,jet} (GeV/#it{c})");
  fHistRejectionReason->GetZaxis()->SetTitle("counts");
  SetRejectionReasonLabels(fHistRejectionReason->GetXaxis());
  fOutput->Add(fHistRejectionReason);

  AllocateTHnSparse();
   
  PostData(1, fOutput);
}

//________________________________________________________________________
Bool_t AliAnalysisTaskDmesonJetCorrelations::IsEventSelected()
{
  // Return true if the event is selected.
  
  Bool_t iseventselected = kFALSE;

  if (fCuts) iseventselected = fCuts->IsEventSelected(fAodEvent);
  
  if (!iseventselected) {
    if (fGeneralHistograms) fHistEventRejection->Fill("HFevSel",1);
    return kFALSE;
  }
  
  return AliAnalysisTaskEmcalJet::IsEventSelected();
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskDmesonJetCorrelations::Run()
{
  // Run the analysis.

  return kTRUE;
}

//_______________________________________________________________________________
Bool_t AliAnalysisTaskDmesonJetCorrelations::FillHistograms()
{
  // Fill the histograms.

  AliJetContainer *jets = GetJetContainer(0);

  const Int_t nDcand = fCandidateArray->GetEntriesFast();

  Bool_t fillRejReason = kTRUE;

  AliDebug(2,"Starting D meson candidate loop");
  for (Int_t icand = 0; icand < nDcand; icand++) {
    AliDebug(2,Form("D meson candidate %d", icand));
    AliAODRecoDecayHF2Prong* Dcand = static_cast<AliAODRecoDecayHF2Prong*>(fCandidateArray->At(icand));
    if (!Dcand) continue;

    TLorentzVector Dvector;
    Double_t invMassD = 0;
    Double_t softPionPtD = 0;
    Double_t invMass2prong = 0;

    TLorentzVector jetVector;
    Double_t leadPtJet = 0;
    Double_t areaJet = 0;
    Int_t constJet = 0;

    if (fCandidateType == kD0toKpi) {
      AliDebug(2,"Checking if D0 meson is selected");
      Int_t isSelected = fCuts->IsSelected(Dcand, AliRDHFCuts::kAll, fAodEvent);
      if (isSelected == 0) {  // Not selected as D0/D0bar
        AliDebug(2,"Not selected");
        continue;
      }
      else if (isSelected == 1) {  // selected as a D0
        AliDebug(2,"Selected as D0");
        invMassD = Dcand->InvMassD0();
      }
      else if (isSelected == 2) { // selected as a D0bar
        AliDebug(2,"Selected as D0bar");
        invMassD = Dcand->InvMassD0bar();
      }
      else if (isSelected == 3) { // selected as a D0bar/D0 (PID on K and pi undecisive)
        AliDebug(2,"Selected as either D0 or D0bar");
        Double_t massD0 = Dcand->InvMassD0();
        Double_t massD0bar = Dcand->InvMassD0bar();

        TParticlePDG* D0part = TDatabasePDG::Instance()->GetParticle(TMath::Abs(421));
        Float_t pdgMass = D0part->Mass();

        AliDebug(2,Form("D0 inv mass = %.3f, D0bar inv mass = %.3f, PDG mass = %.3f", massD0, massD0bar, pdgMass));
        
        // Select D0 or D0bar depending on which one gives a mass closest to the PDG value
        if (TMath::Abs(massD0 - pdgMass) < TMath::Abs(massD0bar - pdgMass)) {
          AliDebug(2, "Mass closer to D0");
          invMassD = massD0;
        }
        else {
          AliDebug(2, "Mass closer to D0bar");
          invMassD = massD0bar;
        }
      }
    }
    else if (fCandidateType == kDstartoKpipi) {
      AliAODRecoCascadeHF* DstarCand = static_cast<AliAODRecoCascadeHF*>(Dcand);
      invMassD = DstarCand->InvMassDstarKpipi();
      softPionPtD = DstarCand->GetBachelor()->Pt();
      invMass2prong = DstarCand->InvMassD0();
    }

    Dvector.SetPtEtaPhiM(Dcand->Pt(), Dcand->Eta(), Dcand->Phi(), invMassD);

    // Look for D-jet correlation
    AliEmcalJet* jet = 0;
    Double_t deltaR = fMaxR;
    if (jets) {
      AliEmcalJet* jet_new = 0;
      jets->ResetCurrentID();
      AliDebug(2,"Starting jet loop");
      while ((jet_new = jets->GetNextJet())) {
        AliDebug(2,Form("Jet %d", jets->GetCurrentID()));
        
        if (!jets->AcceptJet(jet_new)) {
          if (fillRejReason) fHistRejectionReason->Fill(jets->GetRejectionReasonBitPosition(), jet_new->Pt());
          AliDebug(2, Form("Rejecting jet %d, pT = %.3f, eta = %.3f, phi = %.3f, reason = %d", jets->GetCurrentID(), jet_new->Pt(), jet_new->Eta(), jet_new->Phi(), jets->GetRejectionReasonBitPosition()));
          continue;
        }
        
        Double_t deltaR_new = jet_new->DeltaR(Dcand);
        if (deltaR_new < deltaR) {
          jet = jet_new;
          deltaR = deltaR_new;
        }
      }
      fillRejReason = kFALSE;
    }

    if (jet) {
      jetVector.SetPtEtaPhiM(jet->Pt(), jet->Eta(), jet->Phi(), 0);
      
      leadPtJet = jet->MaxPartPt();
      areaJet = jet->Area();
      constJet = jet->N();
    }

    AliDebug(2,"Filling THnSparse");
    FillTHnSparse(Dvector, softPionPtD, invMass2prong, jetVector, leadPtJet, areaJet, constJet);
  }
  
  return kTRUE;
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::ExecOnce()
{
  // Execute only once (for the first event).

  if (fInhibitTask) return;

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
  
  fCandidateArray = dynamic_cast<TClonesArray*>(GetInputData(1));
  if (fCandidateArray) {
    TString objname(fCandidateArray->GetClass()->GetName());
    TClass cls(objname);
    TString className;
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
    if (!cls.InheritsFrom(className)) {
      AliError(Form("%s: Objects of type %s in %s are not inherited from %s! Task will not run.", 
                    GetName(), cls.GetName(), fCandidateArray->GetName(), className.Data())); 
      fCandidateArray = 0;
      fInhibitTask = kTRUE;
      return;
    }
  }
  else {
    AliError(Form("%s: Unable to find the candidate array in input slot 1.", GetName()));
    return;
  }
  
  AliAnalysisTaskEmcalJet::ExecOnce();
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
    nbins[dim] = 101;
    min[dim] = 0;
    max[dim] = 2*TMath::Pi()*nbins[dim]/(nbins[dim]-1);
    dim++;
  }

  if (fShowInvMass) {
    title[dim] = "#it{M}_{D} (GeV/#it{c}^{2})";
    nbins[dim] = 200;
    min[dim] = fMinMass;
    max[dim] = fMaxMass;
    dim++;
  }

  if (fShow2ProngInvMass) {
    title[dim] = "#it{M}_{2-prong} (GeV/#it{c}^{2})";
    nbins[dim] = 200;
    
    Double_t D0mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(421))->Mass();
    
    min[dim] = D0mass - 0.04;
    max[dim] = D0mass + 0.04;
    dim++;
  }

  if (fShowDeltaInvMass) {
    title[dim] = "#it{M}_{D*} - #it{M}_{D_{0}} (GeV/#it{c}^{2})";
    nbins[dim] = 200;
    
    Double_t D0mass = TDatabasePDG::Instance()->GetParticle(TMath::Abs(421))->Mass();
    
    min[dim] = fMinMass - D0mass;  // subtract mass of D0
    max[dim] = fMaxMass - D0mass;
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
    max[dim] = fMaxR*1.1;
    dim++;
  }

  if (fShowDeltaEta) {
    title[dim] = "#eta_{D} - #eta_{jet}";
    nbins[dim] = 100;
    min[dim] = -1;
    max[dim] = 1;
    dim++;
  }
  
  if (fShowDeltaPhi) {
    title[dim] = "#phi_{D} - #phi_{jet} (rad)";
    nbins[dim] = 201;
    min[dim] = fMinDeltaPhiHisto;
    max[dim] = TMath::TwoPi()*nbins[dim]/(nbins[dim]-1) + fMinDeltaPhiHisto;
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
    nbins[dim] = 101;
    min[dim] = 0;
    max[dim] = 2*TMath::Pi()*nbins[dim]/(nbins[dim]-1);
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

  fDmesons = new THnSparseD("fDmesons","fDmesons",dim,nbins,min,max);
  fOutput->Add(fDmesons);
  for (Int_t i = 0; i < dim; i++) {
    fDmesons->GetAxis(i)->SetTitle(title[i]);
  }
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::FillTHnSparse(TLorentzVector D, Double_t softPionPtD, Double_t invMass2prong,
                                                         TLorentzVector jet, Double_t leadPtJet, Double_t areaJet, Int_t constJet)
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
    
    deltaPhi = D.Phi() - jet.Phi();
    if (deltaPhi < 0)               deltaPhi += TMath::TwoPi();
    if (deltaPhi > TMath::TwoPi())  deltaPhi -= TMath::TwoPi();

    Double_t deltaPhi_min = deltaPhi;
    if (deltaPhi_min > TMath::Pi()) deltaPhi_min -= TMath::TwoPi();

    while (deltaPhi < fMinDeltaPhiHisto) deltaPhi += TMath::TwoPi();
    while (deltaPhi > fMinDeltaPhiHisto + TMath::TwoPi()) deltaPhi -= TMath::TwoPi();
    deltaEta = D.Eta() - jet.Eta();
    deltaR = TMath::Sqrt(deltaPhi_min*deltaPhi_min + deltaEta*deltaEta);
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
    else AliWarning(Form("Unable to fill dimension %s!",title.Data()));
  }

  fDmesons->Fill(contents);
}

//_______________________________________________________________________________
void AliAnalysisTaskDmesonJetCorrelations::SetMassLimits(Double_t range, Int_t pdg)
{
  // Set the mass limits for the histograms using information from TDatabasePDG.

  TParticlePDG* part = TDatabasePDG::Instance()->GetParticle(TMath::Abs(pdg));

  Info("AliAnalysisTaskDmesonJetCorrelations::SetMassLimits", "Setting mass limit to particle %d '%s'.", pdg, part->GetName());
  
  Float_t mass = part->Mass();

  SetMassLimits(mass - range, mass + range);
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
