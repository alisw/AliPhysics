/// \class AliAnalysisTaskDmesonJets
/// \brief Analysis task for D meson jets
///
/// This task selects D meson candidates according to predefined cuts,
/// then runs a jet finder to reconstruct the jets that contain
/// the D meson candidates.
///
/// The main output is stored in a THnSparse histogram
///
/// \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
/// \date Feb 6, 2016

#ifndef ALIANALYSISTASKDMESONJETS_H
#define ALIANALYSISTASKDMESONJETS_H

/**************************************************************************
* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
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

class TClonesArray;
class AliRDHFCuts;
class AliAODEvent;
class AliAODRecoDecay;
class AliAODRecoDecayHF2Prong;
class AliAODRecoCascadeHF;
class AliVParticle;
class AliAODMCParticle;
class AliHFAODMCParticleContainer;
class AliHFTrackContainer;
class AliParticleContainer;
class AliClusterContainer;
class THnSparse;
class AliFJWrapper;
class THashList;

#include <list>
#include <vector>

#include "AliTLorentzVector.h"
#include "THistManager.h"

#include "AliAnalysisTaskEmcalJet.h"

using std::vector;
using std::list;

class AliAnalysisTaskDmesonJets : public AliAnalysisTaskEmcal
{
 public:

  typedef AliJetContainer::EJetType_t EJetType_t;
  typedef AliJetContainer::EJetAlgo_t EJetAlgo_t;
  typedef AliJetContainer::ERecoScheme_t ERecoScheme_t;

  enum ECandidateType_t  { kD0toKpi, kDstartoKpipi };
  enum EMCMode_t { kNoMC, kSignalOnly, kBackgroundOnly, kMCTruth };
  enum EMesonOrigin_t {
    kUnknownQuark = BIT(0),
    kFromBottom   = BIT(1),
    kFromCharm    = BIT(2)
  };

  enum EMesonDecayChannel_t {
    kDecayOther          = BIT(0),
    kDecayD0toKpi        = BIT(1),
    kDecayDStartoKpipi   = BIT(2),
    kAnyDecay            = kDecayOther | kDecayD0toKpi | kDecayDStartoKpipi
  };

  enum EAxis_t {
   kPositionD         = BIT(0) , // Add the D meson eta/phi axis in the THnSparse
   kInvMass           = BIT(1) , // Add the invariant mass axis in the THnSparse
   k2ProngInvMass     = BIT(2) , // Add the 2 prong invariant mass axis in the THnSparse (for D* this is the inv mass of the D0)
   kSoftPionPt        = BIT(3) , // Add the soft pion pt axis in the THnSparse (for D*)
   kDeltaR            = BIT(4) , // Add the delta R axis in the THnSparse
   kDeltaEta          = BIT(5) , // Add the delta eta axis in the THnSparse
   kDeltaPhi          = BIT(6) , // Add the delta phi axis in the THnSparse
   kPositionJet       = BIT(7) , // Add the jet eta/phi axis in the THnSparse
   kLeadingPt         = BIT(8), // Add the leading pt axis in the THnSparse
   kJetConstituents   = BIT(9), // Add the jet constituent axis in the THnSparse
   kDaughterDistances = BIT(10), // Add n axis in THnSparse with the daughters' distances
  };

  /// \class AliEmcalDmesonJetInfo
  /// \brief Struct that encapsulates D meson jets
  ///
  /// This struct encapsulates D meson jet
  /// information that can be easily passed to a function.
  struct AliDmesonJetInfo {
    AliDmesonJetInfo() : fD(), fSoftPionPt(0), fInvMass2Prong(0), fJet(), fJetLeadingPt(0), fJetNConstituents(0), fDaughterDistances(3) {}

    AliTLorentzVector fD                   ; ///< 4-momentum of the D meson candidate
    Double_t          fSoftPionPt          ; ///< Transverse momentum of the soft pion of the D* candidate
    Double_t          fInvMass2Prong       ; ///< 2-prong mass of the D* candidate (w/o the soft pion)
    AliTLorentzVector fJet                 ; ///< 4-momentum of the jet
    Double_t          fJetLeadingPt        ; ///< Transverse momentum of the leading particle of the jet
    Int_t             fJetNConstituents    ; ///< Number of constituents of the jet
    TArrayD           fDaughterDistances   ; ///< Distance of the D meson daughters from the jet axis

    void Reset();

    void Print() const;
  };

  class AliJetDefinition : public TObject {
  public:
    AliJetDefinition();
    AliJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco);
    AliJetDefinition(const AliJetDefinition &source);

    AliJetDefinition& operator=(const AliJetDefinition& source);

    const char* GetName() const;

    friend bool        operator< (const AliJetDefinition& lhs, const AliJetDefinition& rhs);
    friend inline bool operator> (const AliJetDefinition& lhs, const AliJetDefinition& rhs){ return rhs < lhs    ; }
    friend inline bool operator<=(const AliJetDefinition& lhs, const AliJetDefinition& rhs){ return !(lhs > rhs) ; }
    friend inline bool operator>=(const AliJetDefinition& lhs, const AliJetDefinition& rhs){ return !(lhs < rhs) ; }

    friend bool        operator==(const AliJetDefinition& lhs, const AliJetDefinition& rhs);
    friend inline bool operator!=(const AliJetDefinition& lhs, const AliJetDefinition& rhs){ return !(lhs == rhs); }

  protected:
    friend class AliAnalysisTaskDmesonJets;
    friend class AnalysisEngine;

    EJetType_t                fJetType       ; ///<  Jet type (charged, full, neutral)
    Double_t                  fRadius        ; ///<  Jet radius
    EJetAlgo_t                fJetAlgo       ; ///<  Jet algorithm (kt, anti-kt,...)
    ERecoScheme_t             fRecoScheme    ; ///<  Jet recombination scheme (pt scheme, E scheme, ...)
    vector<AliDmesonJetInfo>  fDmesonJets    ; //!<! Array containing the D meson jets

  private:
    /// \cond CLASSIMP
    ClassDef(AliJetDefinition, 1);
    /// \endcond
  };

  /// \class AnalysisEngine
  /// \brief Struct that encapsulates analysis parameters
  ///
  /// This struct encapsulates analysis parameters
  /// for the D meson jet analysis.
  class AnalysisEngine : public TObject {
  public:
    static EMesonOrigin_t CheckOrigin(AliAODMCParticle* part, TClonesArray* mcArray);
    static EMesonDecayChannel_t CheckDecayChannel(AliAODMCParticle* part, TClonesArray* mcArray);

    AnalysisEngine();
    AnalysisEngine(ECandidateType_t type, EMCMode_t MCmode, AliRDHFCuts* cuts = 0, Int_t nBins=80, Double_t range = 0.50);
    AnalysisEngine(const AnalysisEngine &source);
    AnalysisEngine& operator=(const AnalysisEngine& source);

    virtual ~AnalysisEngine();

    void SetCandidateType(ECandidateType_t t)       { fCandidateType  = t    ; }
    void SetMCMode(EMCMode_t m)                     { fMCMode         = m    ; }
    void SetNMassBins(Int_t n)                      { fNMassBins      = n    ; }
    void SetMassRange(Double_t min, Double_t max)   { fMinMass        = min  ; fMaxMass        = max  ; }
    void AdoptRDHFCuts(AliRDHFCuts* cuts);
    void SetRDHFCuts(AliRDHFCuts* cuts);
    void SetRejectedOriginMap(UInt_t m)             { fRejectedOrigin = m    ; }
    void SetAcceptedDecayMap(UInt_t m)              { fAcceptedDecay  = m    ; }

    const char* GetName() const;
    const char* GetName(const AliJetDefinition& jetDef) const;

    AliJetDefinition* AddJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco);
    AliJetDefinition* AddJetDefinition(const AliJetDefinition& def);
    std::list<AliJetDefinition>::iterator FindJetDefinition(const AliJetDefinition& eng);

    friend bool        operator< (const AnalysisEngine& lhs, const AnalysisEngine& rhs);
    friend inline bool operator> (const AnalysisEngine& lhs, const AnalysisEngine& rhs){ return rhs < lhs    ; }
    friend inline bool operator<=(const AnalysisEngine& lhs, const AnalysisEngine& rhs){ return !(lhs > rhs) ; }
    friend inline bool operator>=(const AnalysisEngine& lhs, const AnalysisEngine& rhs){ return !(lhs < rhs) ; }

    friend bool        operator==(const AnalysisEngine& lhs, const AnalysisEngine& rhs);
    friend inline bool operator!=(const AnalysisEngine& lhs, const AnalysisEngine& rhs){ return !(lhs == rhs); }

  protected:
    void RunAnalysis();

    ECandidateType_t                   fCandidateType         ; ///<  Candidate type
    TString                            fCandidateName         ; ///<  Candidate name
    UInt_t                             fCandidatePDG          ; ///<  Candidate PDG
    UChar_t                            fNDaughters            ; ///<  Number of daughters
    TArrayI                            fPDGdaughters          ; ///<  List of the PDG code of the daughters
    TString                            fBranchName            ; ///<  AOD branch where the D meson candidate are found
    EMCMode_t                          fMCMode                ; ///<  MC mode: No MC (data and MC detector level), background-only (MC), signal-only (MC), MC truth (particle level)
    Int_t                              fNMassBins             ; ///<  Mass number of bins
    Double_t                           fMinMass               ; ///<  Min mass in histogram axis
    Double_t                           fMaxMass               ; ///<  Max mass in histogram axis
    AliRDHFCuts                       *fRDHFCuts              ; ///<  D meson candidates cuts
    UInt_t                             fRejectedOrigin        ; ///<  Bit mask with D meson origins that are rejected
    UInt_t                             fAcceptedDecay         ; ///<  Bit mask with D meson decays that are accepted
    Bool_t                             fInhibit               ; ///<  Inhibit the task
    list<AliJetDefinition>             fJetDefinitions        ; ///<  Jet definitions
    TClonesArray                      *fCandidateArray        ; //!<! D meson candidate array
    AliHFAODMCParticleContainer       *fMCContainer           ; //!<! MC particle container
    AliHFTrackContainer               *fTrackContainer        ; //!<! Track container
    AliClusterContainer               *fClusterContainer      ; //!<! Cluster container
    AliAODEvent                       *fAodEvent              ; //!<! AOD event
    AliFJWrapper                      *fFastJetWrapper        ; //!<! Fastjet wrapper
    THistManager                      *fHistManager           ; //!<! Histogram manager

    friend class AliAnalysisTaskDmesonJets;

  private:

    void                AddInputVectors(AliEmcalContainer* cont, Int_t offset, TH2* rejectHist);
    void                SetCandidateProperties(Double_t range);
    AliAODMCParticle*   MatchToMC() const;
    void                RunDetectorLevelAnalysis();
    void                RunParticleLevelAnalysis();

    Bool_t              ExtractParticleLevelHFAttributes(const AliAODMCParticle* part, AliDmesonJetInfo& DmesonJet);
    Bool_t              ExtractRecoDecayAttributes(const AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, UInt_t i);
    Bool_t              ExtractD0Attributes(const AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, UInt_t i);
    Bool_t              ExtractDstarAttributes(const AliAODRecoCascadeHF* DstarCand, AliDmesonJetInfo& DmesonJet, UInt_t i);
    Bool_t              FindJet(AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, AliJetDefinition& jetDef);

    /// \cond CLASSIMP
    ClassDef(AnalysisEngine, 2);
    /// \endcond
  };

  AliAnalysisTaskDmesonJets();
  AliAnalysisTaskDmesonJets(const char* name);
  virtual ~AliAnalysisTaskDmesonJets();

  AnalysisEngine* AddAnalysisEngine(ECandidateType_t type, EMCMode_t bkgMode, EJetType_t jettype, Double_t jetradius, TString cutfname = "");
  AnalysisEngine* AddAnalysisEngine(ECandidateType_t type, EMCMode_t bkgMode, const AliJetDefinition& jetDef, TString cutfname = "");
  std::list<AnalysisEngine>::iterator FindAnalysisEngine(const AnalysisEngine& eng);

  void SetShowPositionD(Bool_t b = kTRUE)         { fEnabledAxis = b ?  fEnabledAxis | kPositionD         : fEnabledAxis & ~kPositionD         ; }
  void SetShowInvMass(Bool_t b = kTRUE)           { fEnabledAxis = b ?  fEnabledAxis | kInvMass           : fEnabledAxis & ~kInvMass           ; }
  void SetShow2ProngInvMass(Bool_t b = kTRUE)     { fEnabledAxis = b ?  fEnabledAxis | k2ProngInvMass     : fEnabledAxis & ~k2ProngInvMass     ; }
  void SetShowSoftPionPt(Bool_t b = kTRUE)        { fEnabledAxis = b ?  fEnabledAxis | kSoftPionPt        : fEnabledAxis & ~kSoftPionPt        ; }
  void SetShowDeltaR(Bool_t b = kTRUE)            { fEnabledAxis = b ?  fEnabledAxis | kDeltaR            : fEnabledAxis & ~kDeltaR            ; }
  void SetShowDeltaEta(Bool_t b = kTRUE)          { fEnabledAxis = b ?  fEnabledAxis | kDeltaEta          : fEnabledAxis & ~kDeltaEta          ; }
  void SetShowDeltaPhi(Bool_t b = kTRUE)          { fEnabledAxis = b ?  fEnabledAxis | kDeltaPhi          : fEnabledAxis & ~kDeltaPhi          ; }
  void SetShowPositionJet(Bool_t b = kTRUE)       { fEnabledAxis = b ?  fEnabledAxis | kPositionJet       : fEnabledAxis & ~kPositionJet       ; }
  void SetShowLeadingPt(Bool_t b = kTRUE)         { fEnabledAxis = b ?  fEnabledAxis | kLeadingPt         : fEnabledAxis & ~kLeadingPt         ; }
  void SetShowJetConstituents(Bool_t b = kTRUE)   { fEnabledAxis = b ?  fEnabledAxis | kJetConstituents   : fEnabledAxis & ~kJetConstituents   ; }
  void SetShowDaughterDistances(Bool_t b = kTRUE) { fEnabledAxis = b ?  fEnabledAxis | kDaughterDistances : fEnabledAxis & ~kDaughterDistances ; }

  virtual void         UserCreateOutputObjects();
  virtual void         ExecOnce();
  virtual Bool_t       Run();
  virtual Bool_t       FillHistograms();

 protected:

  AliRDHFCuts*         LoadDMesonCutsFromFile(TString cutfname, TString cutsname);
  void                 AllocateTHnSparse(const AnalysisEngine& param);
  void                 FillTHnSparse(THnSparse* h, const AliDmesonJetInfo& DmesonJet);
  
  static const char*   GetHFEventRejectionReasonLabel(UInt_t& bitmap);
  static void          CalculateMassLimits(Double_t range, Int_t pdg, Int_t nbins, Double_t& minMass, Double_t& maxMass);

  list<AnalysisEngine> fAnalysisEngines           ; ///<  Array of analysis parameters
  UInt_t               fEnabledAxis               ; ///<  Use bit defined in EAxis_t to enable axis in the THnSparse
  THistManager         fHistManager               ; ///<  Histogram manager
  AliAODEvent         *fAodEvent                  ; //!<! AOD event
  AliFJWrapper        *fFastJetWrapper            ; //!<! Fastjet wrapper

 private:
   
  AliAnalysisTaskDmesonJets(const AliAnalysisTaskDmesonJets &source);
  AliAnalysisTaskDmesonJets& operator=(const AliAnalysisTaskDmesonJets& source);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskDmesonJets, 1);
  /// \endcond
};

#endif
