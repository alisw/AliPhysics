#ifndef ALIANALYSISTASKFLAVOURJETCORRELATIONS_H
#define ALIANALYSISTASKFLAVOURJETCORRELATIONS_H
// $Id$
/**************************************************************************
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

//-----------------------------------------------------------------------
// Author : A. Grelli,  Utrecht University
//          C. Bianchin, Utrecht University
//          X. Zhang, LBNL
//-----------------------------------------------------------------------


#include <TH2F.h>
#include "AliAODEvent.h"
#include "AliPicoTrack.h"
#include "AliAnalysisTaskEmcalJet.h"

class TH3F;
class TParticle ;
class TClonesArray ;
class AliMCParticle;
class AliAODMCParticle;
class AliRDHFCuts;
class AliEmcalJet;
class AliAODRecoDecayHF;
class AliAODRecoCascadeHF;
class AliAODEvent;

class AliAnalysisTaskFlavourJetCorrelations : public AliAnalysisTaskEmcalJet 
{

 public:

  enum ECandidateType{ kD0toKpi, kDstartoKpipi };

  AliAnalysisTaskFlavourJetCorrelations();
  AliAnalysisTaskFlavourJetCorrelations(const Char_t* name,AliRDHFCuts* cuts, ECandidateType candtype, TString jetArrName);
  virtual ~AliAnalysisTaskFlavourJetCorrelations();

  virtual void     UserCreateOutputObjects();
  virtual void     UserExec(Option_t *option);
  virtual void     Terminate(Option_t *);
  virtual void     Init();
  virtual void     LocalInit() {Init();}

  // inizializations
  Bool_t DefineHistoForAnalysis();  

  // set MC usage
  void   SetMC(Bool_t theMCon) {fUseMCInfo = theMCon;}
  Bool_t GetMC() const {return fUseMCInfo;}
  
  void SetMassLimits(Double_t range, Int_t pdg);
  void SetMassLimits(Double_t lowlimit, Double_t uplimit);

  //jet reconstruction algorithm
  void SetJetArrayName(TString jetArrName) {fJetArrName=jetArrName;};
  TString GetJetArrayName() const {return fJetArrName;};

  // trigger on jet events
  void SetTriggerOnLeadingJet(Bool_t triggerOnLeadingJet) {fLeadingJetOnly=triggerOnLeadingJet;};
  Bool_t GetTriggerOnLeadingJet() const {return fLeadingJetOnly;}


  // Array of D0 width for the Dstar
  Bool_t SetD0WidthForDStar(Int_t nptbins,Float_t* width);

  //Bool_t   FillMCDJetInfo(AliPicoTrack *jetTrk,AliEmcalJet* jet, TClonesArray *mcArray,Double_t ptjet);
  void FillHistogramsRecoJetCorr(AliAODRecoDecayHF* candidate, AliEmcalJet *jet);
  void FillHistogramsD0JetCorr(AliAODRecoDecayHF* candidate, Double_t dPhi, Double_t z, Double_t ptD, Double_t deltaR, AliAODEvent* aodEvent);

  void FillHistogramsDstarJetCorr(AliAODRecoCascadeHF* dstar, Double_t dPhi, Double_t z, Double_t ptD, Double_t deltaR);

  void FillHistogramsD(Double_t mass,Double_t dphi, Double_t z,Double_t ptD, Double_t deltaR);

 private:
  
  AliAnalysisTaskFlavourJetCorrelations(const AliAnalysisTaskFlavourJetCorrelations &source);
  AliAnalysisTaskFlavourJetCorrelations& operator=(const AliAnalysisTaskFlavourJetCorrelations& source); 

  Double_t Z(AliVParticle* part,AliEmcalJet* jet) const;
  Float_t DeltaR(AliVParticle *p1, AliVParticle *p2) const;

  /*
  Bool_t IsD(Int_t pdg) const;
  Bool_t IsD(Int_t pdg,Int_t abspdgD) const;
  Bool_t PartFromC(AliMCParticle* mother) const;
  Int_t GetFirstMother(Int_t lab,TClonesArray* mcarr) const; 
  Int_t FindPDGInFamily(Int_t labpart,Int_t pdgcode, TClonesArray *mcArray) const;
  */

  Bool_t fUseMCInfo;             //  Use MC info
  Int_t  fCandidateType;         // Dstar or D0
  Int_t  fPDGmother;             // PDG code of D meson
  Int_t  fNProngs;               // number of prong of the decay channel  
  Int_t  fPDGdaughters[4];       // PDG codes of daughters
  Float_t fSigmaD0[30];            //
  TString fBranchName;           // AOD branch name
  TList *fOutput;                  //! user output
  AliRDHFCuts *fCuts;  // Cuts 

  Double_t fMinMass;             // mass lower limit histogram
  Double_t fMaxMass;             // mass upper limit histogram

  TString  fJetArrName;          // name of the jet array, taken from the task running the jet finder
  TString fCandArrName;          // string which correspond to the candidate type
  Bool_t fLeadingJetOnly;        // use only the leading jet in the event to make the correlations

  ClassDef(AliAnalysisTaskFlavourJetCorrelations,1); // class for charm-jet correlations
};

#endif
