#ifndef ALIANALYSISTASKDMESONJETCORRELATIONS_H
#define ALIANALYSISTASKDMESONJETCORRELATIONS_H
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
// Author : Salvatore Aiola, Yale University, salvatore.aiola@cern.ch
//-----------------------------------------------------------------------

class TClonesArray;
class TList;
class THnSparse;
class TH2;
class TLorentzVector;
class AliRDHFCuts;
class AliAODEvent;
class AliAODRecoDecay;
class AliAODRecoDecayHF2Prong;
class AliVParticle;
class AliAODMCParticle;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskDmesonJetCorrelations : public AliAnalysisTaskEmcalJet 
{
 public:

  enum ECandidateType  { kD0toKpi, kDstartoKpipi };
  enum EMatchingType   { kGeometricalMatching, kConstituentMatching, kJetLoop };
  enum EMatchingStatus { kSingleMatch = 0, kMultipleMatches = 1, kJetNotAccepted = 2, kNotMatched = 3 };
   
  AliAnalysisTaskDmesonJetCorrelations();
  AliAnalysisTaskDmesonJetCorrelations(const char* name, AliRDHFCuts* cuts, ECandidateType cand, Bool_t useExchCont=kFALSE);
  virtual ~AliAnalysisTaskDmesonJetCorrelations();

  void SetShowPositionD(Bool_t b = kTRUE)       { fShowPositionD       = b   ; }
  void SetShowInvMass(Bool_t b = kTRUE)         { fShowInvMass         = b   ; }
  void SetShow2ProngInvMass(Bool_t b = kTRUE)   { fShow2ProngInvMass   = b   ; }
  void SetShowDeltaInvMass(Bool_t b = kTRUE)    { fShowDeltaInvMass    = b   ; }
  void SetShowSoftPionPt(Bool_t b = kTRUE)      { fShowSoftPionPt      = b   ; }
  void SetShowDmesonZ(Bool_t b = kTRUE)         { fShowDmesonZ         = b   ; }
  void SetShowDeltaR(Bool_t b = kTRUE)          { fShowDeltaR          = b   ; }
  void SetShowDeltaEta(Bool_t b = kTRUE)        { fShowDeltaEta        = b   ; }
  void SetShowDeltaPhi(Bool_t b = kTRUE)        { fShowDeltaPhi        = b   ; }
  void SetShowPositionJet(Bool_t b = kTRUE)     { fShowPositionJet     = b   ; }
  void SetShowLeadingPt(Bool_t b = kTRUE)       { fShowLeadingPt       = b   ; }
  void SetShowJetArea(Bool_t b = kTRUE)         { fShowJetArea         = b   ; }
  void SetShowJetConstituents(Bool_t b = kTRUE) { fShowJetConstituents = b   ; }
  void SetShowMatchingLevel(Bool_t b = kTRUE)   { fShowMatchingLevel   = b   ; }
  void SetShowDaughterDistance(Int_t n = 2)     { fShowDaughterDistance= n   ; }
  void SetCandidateType(ECandidateType c)       { fCandidateType       = c   ; }
  void SetMaxR(Double_t r)                      { fMaxR       = r > 0 ? r : 1; }
  void SetMatchingType(EMatchingType m)         { fMatchingType        = m   ; }
  void SetPlotOnlyAcceptedJets(Bool_t b)        { fOnlyAcceptedJets    = b   ; }
  void SetPlotOnlySingleMatches(Bool_t b)       { fOnlySingleMatches   = b   ; }
  void SetCheckTrackColl(Bool_t b)              { fCheckTrackColl      = b   ; }
  void SetNBinsMass(Int_t n)                    { fNBinsMass           = n   ; }
  void SetParticleLevel(Bool_t s)               { fParticleLevel       = s   ; }
  void SetAliEmcalParticleMode(Bool_t m)        { fAliEmcalParticleMode= m   ; }
  void SetMassLimits(Double_t range, Int_t pdg);
  void SetMassLimits(Double_t lowlimit, Double_t uplimit);

  virtual void     UserCreateOutputObjects();
  virtual void     ExecOnce();
  virtual Bool_t   IsEventSelected();
  virtual Bool_t   Run();
  virtual Bool_t   FillHistograms();
  virtual void     Init();
  virtual void     LocalInit() { Init(); }
  virtual void     Terminate(Option_t *);

 protected:

  void                 AllocateTHnSparse();
  void                 FillTHnSparse(TLorentzVector D, Double_t softPionPtD, Double_t invMass2prong,
                                     TLorentzVector jet, Double_t leadPtJet, Double_t areaJet, Int_t constJet,
                                     Int_t matchingStatus, Double_t matchingLevel, Double_t daughterDist[5]);
  Int_t                FindMatchedJet(EMatchingType matchType, AliVParticle* cand, TArrayD& matchingLevel, TList& matchedJets);
  Double_t             CalculateMatchingLevel(EMatchingType matchType, AliVParticle* cand, AliEmcalJet* jet, Bool_t reset=kFALSE);
  Double_t             CalculateGeometricalMatchingLevel(AliVParticle* cand, AliEmcalJet* jet);
  Double_t             CalculateConstituentMatchingLevel(AliAODRecoDecay* cand, AliEmcalJet* jet, Bool_t reset=kFALSE);

  Bool_t                 ExtractHFcandAttributes(AliVParticle* HFcand, TLorentzVector& Dvector, Double_t& invMassD, Double_t& softPionPtD, Double_t& invMass2prong, UInt_t i);
  Bool_t                 ExtractParticleLevelHFAttributes(AliAODMCParticle* part, TLorentzVector& Dvector, Double_t& invMassD, Double_t& softPionPtD, Double_t& invMass2prong, UInt_t i);
  Bool_t                 ExtractRecoDecayAttributes(AliAODRecoDecayHF2Prong* Dcand, TLorentzVector& Dvector, Double_t& invMassD, Double_t& softPionPtD, Double_t& invMass2prong, UInt_t i);

  void                 DoJetLoop();
  void                 DoDmesonLoop();
  void                 FillHistograms(AliVParticle* HFcand, AliEmcalJet* jet, Int_t matchingStatus, Double_t matchingLevel);
  
  static void          CalculateMassLimits(Double_t range, Int_t pdg, Int_t nbins, Double_t& minMass, Double_t& maxMass);

  AliRDHFCuts     *fCuts                  ; //  Analysis cuts
  ECandidateType   fCandidateType         ; //  Candidate type, D0 or D*
  Double_t         fMinMass               ; //  Min mass in histogram axis
  Double_t         fMaxMass               ; //  Max mass in histogram axis
  Int_t            fNBinsMass             ; //  Number of bins in mass axis
  Double_t         fMaxR                  ; //  Maximum distance for a jet to be matched with a D candidate
  Bool_t           fShowPositionD         ; //  Add the D meson eta/phi axis in the THnSparse
  Bool_t           fShowInvMass           ; //  Add the invariant mass axis in the THnSparse
  Bool_t           fShow2ProngInvMass     ; //  Add the 2 prong invariant mass axis in the THnSparse (for D* this is the inv mass of the D0)
  Bool_t           fShowDeltaInvMass      ; //  Add the delta invariant mass (mass(D*) - mass(D0)) axis in the THnSparse (for D*)
  Bool_t           fShowSoftPionPt        ; //  Add the soft pion pt axis in the THnSparse (for D*)
  Bool_t           fShowDmesonZ           ; //  Add the z of the D meson axis in the THnSparse
  Bool_t           fShowDeltaR            ; //  Add the delta R axis in the THnSparse
  Bool_t           fShowDeltaEta          ; //  Add the delta eta axis in the THnSparse
  Bool_t           fShowDeltaPhi          ; //  Add the delta phi axis in the THnSparse
  Bool_t           fShowPositionJet       ; //  Add the jet eta/phi axis in the THnSparse
  Bool_t           fShowLeadingPt         ; //  Add the leading pt axis in the THnSparse
  Bool_t           fShowJetArea           ; //  Add the jet area axis in the THnSparse
  Bool_t           fShowJetConstituents   ; //  Add the jet constituent axis in the THnSparse
  Bool_t           fShowMatchingLevel     ; //  Add the matching level axis in the THnSparse
  Int_t            fShowDaughterDistance  ; //  Add n axis in THnSparse with the daughters' distances
  Bool_t           fInhibitTask           ; //  Inhibit the task
  EMatchingType    fMatchingType          ; //  Dmeson-jet matching type
  Bool_t           fOnlyAcceptedJets      ; //  Only matched accepted jets are plotted (in any case jets are rejected AFTER matching)
  Bool_t           fOnlySingleMatches     ; //  Only unambiguosly matched jets are plotted
  Bool_t           fCheckTrackColl        ; //  if true, makes histograms with tracks found in D meson candidates, that are NOT found in the track collection used for jet finding
  Bool_t           fParticleLevel         ; //  set particle level analysis
  Bool_t           fUseExchangeContainer  ; //  use exchange container for the D candidate list
  Bool_t           fAliEmcalParticleMode  ; //  use AliEmcalParticle objects in fCandidateArray

  AliAODEvent     *fAodEvent                  ; //! AOD event
  TClonesArray    *fCandidateArray            ; //! D meson candidate array
  TH2             *fHistRejectionReason       ; //! Rejection reason vs. jet pt
  TH1             *fHistTracksNotInJetsPt     ; //! Pt of tracks found in D mesons but not found in collection used for jet finding
  TH2             *fHistTracksNotInJetsEtaPhi ; //! Eta/phi of tracks found in D mesons but not found in collection used for jet finding
  THnSparse       *fDmesons                   ; //! D meson histogram
   
 private:
   
  AliAnalysisTaskDmesonJetCorrelations(const AliAnalysisTaskDmesonJetCorrelations &source);
  AliAnalysisTaskDmesonJetCorrelations& operator=(const AliAnalysisTaskDmesonJetCorrelations& source); 

  ClassDef(AliAnalysisTaskDmesonJetCorrelations, 3); // class for D meson - jet correlations
};

#endif
