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
class TTree;
class AliEMCALGeometry;

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
  typedef AliJetContainer::JetAcceptanceType EJetAcceptanceType_t;

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
   kJetConstituents   = BIT(8)   // Add the jet constituent axis in the THnSparse
  };

  /// \class AliJetInfo
  /// \brief Class that encapsulates jets
  ///
  /// This class encapsulates jet
  /// information that can be easily passed to a function.
  class AliJetInfo {
  public:
    AliJetInfo() : fMomentum(), fNConstituents(0), fNEF(0), fMaxChargedPt(0), fMaxNeutralPt(0) {}

    Double_t Pt()        const { return fMomentum.Pt()       ; }
    Double_t Eta()       const { return fMomentum.Eta()      ; }
    Double_t Phi()       const { return fMomentum.Phi()      ; }
    Double_t Phi_0_2pi() const { return fMomentum.Phi_0_2pi(); }

    AliTLorentzVector fMomentum             ; ///< 4-momentum of the jet
    Int_t             fNConstituents        ; ///< Number of constituents of the jet
    Double_t          fNEF                  ; ///< Neutral Energy Fraction of the jet
    Double_t          fMaxChargedPt         ; ///< Transverse momentum of the leading charged particle (or track)
    Double_t          fMaxNeutralPt         ; ///< Transverse momentum of the leading neutral particle (or cluster)
  };

  /// \class AliEmcalDmesonJetInfo
  /// \brief Class that encapsulates D meson jets
  ///
  /// This class encapsulates D meson jet
  /// information that can be easily passed to a function.
  class AliDmesonJetInfo {
  public:
    AliDmesonJetInfo() : fD(), fSoftPionPt(0), fInvMass2Prong(0), fJets() {}

    AliTLorentzVector  fD                       ; //!<! 4-momentum of the D meson candidate
    Double_t           fSoftPionPt              ; //!<! Transverse momentum of the soft pion of the D* candidate
    Double_t           fInvMass2Prong           ; //!<! 2-prong mass of the D* candidate (w/o the soft pion)
    std::map<std::string, AliJetInfo>
                       fJets                    ; //!<! list of jets

    void Reset();
    Double_t GetZ(std::string n) const;
    Double_t GetDistance(std::string n, Double_t& deta, Double_t& dphi) const;
    Double_t GetDistance(std::string n) const;
    void Print() const;
  };

  /// \class AliJetInfoSummary
  /// \brief Lightweight class that encapsulates D meson jets
  ///
  /// This class encapsulates D meson jet
  /// information in a very compact data structure (49 bits)
  class AliJetInfoSummary {
  public:
    AliJetInfoSummary() : fPt(0), fEta(0), fPhi(0), fR(0), fZ(0) {;}
    AliJetInfoSummary(const AliDmesonJetInfo& source, std::string n);

    virtual void Set(const AliDmesonJetInfo& source, std::string n);

    /// Transverse momentum of the jet in GeV/c
    Double32_t  fPt        ; //[0,200,12]
    /// Eta of the jet
    Double32_t  fEta       ; //[-2,2,10]
    /// Phi of the jet
    Double32_t  fPhi       ; //[0,2*pi,10]
    /// Distance between D meson and jet axis
    Double32_t  fR         ; //[0,2,7]
    /// Z of the D meson
    Double32_t  fZ         ; //[0,1,10]

    /// \cond CLASSIMP
    ClassDef(AliJetInfoSummary, 1);
    /// \endcond
  };

  /// \class AliDmesonInfoSummary
  /// \brief Lightweight class that encapsulates D meson jets
  ///
  /// This class encapsulates D meson
  /// information in a very compact data structure (30 bits)
  class AliDmesonInfoSummary {
  public:
    AliDmesonInfoSummary() : fPt(0), fEta(0), fPhi(0) {;}
    AliDmesonInfoSummary(const AliDmesonJetInfo& source);

    virtual void Set(const AliDmesonJetInfo& source);

    /// Transverse momentum of the D meson in GeV/c
    Double32_t   fPt     ; //[0,200,12]
    /// Eta of the jet
    Double32_t   fEta    ; //[-2,2,9]
    /// Phi of the jet
    Double32_t   fPhi    ; //[0,2*pi,9]

    /// \cond CLASSIMP
    ClassDef(AliDmesonInfoSummary, 1);
    /// \endcond
  };

  /// \class AliD0InfoSummary
  /// \brief Lightweight class that encapsulates D0
  ///
  /// This class encapsulates D0 jet
  /// information in a very compact data structure (42 bits)
  class AliD0InfoSummary : public AliDmesonInfoSummary {
  public:
    AliD0InfoSummary() : AliDmesonInfoSummary(), fInvMass(0) {}
    AliD0InfoSummary(const AliDmesonJetInfo& source);

    virtual void Set(const AliDmesonJetInfo& source);

    /// Invariant mass of the D0 meson candidate in GeV/c2
    Double32_t   fInvMass   ; //[0,5,12]

    /// \cond CLASSIMP
    ClassDef(AliD0InfoSummary, 1);
    /// \endcond
  };

  /// \class AliDStarJetInfoSummary
  /// \brief Lightweight class that encapsulates D*
  ///
  /// This class encapsulates D*
  /// information in a very compact data structure (54 bits)
  class AliDStarInfoSummary : public AliDmesonInfoSummary {
  public:
    AliDStarInfoSummary() : AliDmesonInfoSummary(), f2ProngInvMass(0), fDeltaInvMass(0) {}
    AliDStarInfoSummary(const AliDmesonJetInfo& source);

    virtual void Set(const AliDmesonJetInfo& source);

    ///< Invariant mass of the D0 meson candidate in GeV/c2
    Double32_t   f2ProngInvMass   ; //[0,5,12]
    ///< Difference between the Kpipi and the Kpi invariant masses in GeV/c2
    Double32_t   fDeltaInvMass    ; //[0,1,12]

    /// \cond CLASSIMP
    ClassDef(AliDStarInfoSummary, 1);
    /// \endcond
  };

  class AliHFJetDefinition : public TObject {
  public:
    AliHFJetDefinition();
    AliHFJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco);
    AliHFJetDefinition(const AliHFJetDefinition &source);

    AliHFJetDefinition& operator=(const AliHFJetDefinition& source);

    const char* GetName() const;

    void SetJetPhiRange(Double_t min, Double_t max)       { fMinJetPhi    = min; fMaxJetPhi    = max; }
    void SetJetEtaRange(Double_t min, Double_t max)       { fMinJetEta    = min; fMaxJetEta    = max; }
    void SetJetPtMin(Double_t min)                        { fMinJetPt     = min;                      }
    void SetChargedPtRange(Double_t min, Double_t max)    { fMinChargedPt = min; fMaxChargedPt = max; }
    void SetNeutralPtRange(Double_t min, Double_t max)    { fMinNeutralPt = min; fMaxNeutralPt = max; }
    void SetAcceptanceType(EJetAcceptanceType_t a)        { fAcceptance   = a  ;                      }

    Bool_t IsJetInAcceptance(const AliJetInfo& jet) const;

    friend bool        operator< (const AliHFJetDefinition& lhs, const AliHFJetDefinition& rhs);
    friend inline bool operator> (const AliHFJetDefinition& lhs, const AliHFJetDefinition& rhs){ return rhs < lhs    ; }
    friend inline bool operator<=(const AliHFJetDefinition& lhs, const AliHFJetDefinition& rhs){ return !(lhs > rhs) ; }
    friend inline bool operator>=(const AliHFJetDefinition& lhs, const AliHFJetDefinition& rhs){ return !(lhs < rhs) ; }

    friend bool        operator==(const AliHFJetDefinition& lhs, const AliHFJetDefinition& rhs);
    friend inline bool operator!=(const AliHFJetDefinition& lhs, const AliHFJetDefinition& rhs){ return !(lhs == rhs); }

  protected:
    friend class AliAnalysisTaskDmesonJets;
    friend class AnalysisEngine;

    void                      SetDetectorJetEtaPhiRange(const AliEMCALGeometry* const geom, Int_t run);

    EJetType_t                fJetType       ; ///<  Jet type (charged, full, neutral)
    Double_t                  fRadius        ; ///<  Jet radius
    EJetAlgo_t                fJetAlgo       ; ///<  Jet algorithm (kt, anti-kt,...)
    ERecoScheme_t             fRecoScheme    ; ///<  Jet recombination scheme (pt scheme, E scheme, ...)
    EJetAcceptanceType_t      fAcceptance    ; ///<  Jet acceptance
    Double_t                  fMinJetPt      ; ///<  Minimum jet pT
    Double_t                  fMinJetPhi     ; ///<  Minimum jet phi
    Double_t                  fMaxJetPhi     ; ///<  Maximum jet phi
    Double_t                  fMinJetEta     ; ///<  Minimum jet eta
    Double_t                  fMaxJetEta     ; ///<  Maximum jet eta
    Double_t                  fMinChargedPt  ; ///<  Minimum pt of the leading charged particle (or track)
    Double_t                  fMaxChargedPt  ; ///<  Maximum pt of the leading charged particle (or track)
    Double_t                  fMinNeutralPt  ; ///<  Minimum pt of the leading neutral particle (or cluster)
    Double_t                  fMaxNeutralPt  ; ///<  Maximum pt of the leading neutral particle (or cluster)

  private:
    /// \cond CLASSIMP
    ClassDef(AliHFJetDefinition, 2);
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
    const char* GetName(const AliHFJetDefinition& jetDef) const;

    AliHFJetDefinition* AddJetDefinition(EJetType_t type, Double_t r, EJetAlgo_t algo, ERecoScheme_t reco);
    AliHFJetDefinition* AddJetDefinition(const AliHFJetDefinition& def);
    std::vector<AliHFJetDefinition>::iterator FindJetDefinition(const AliHFJetDefinition& eng);
    std::vector<AliAnalysisTaskDmesonJets::AliHFJetDefinition>& GetJetDefinitions() { return fJetDefinitions; }
    Bool_t IsAnyJetInAcceptance(const AliDmesonJetInfo& dMesonJet) const;

    void Init(const AliEMCALGeometry* const geom, Int_t runNumber);
    TTree* BuildTree();
    TTree* GetTree() { return fTree; }
    Bool_t FillTree();

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
    vector<AliHFJetDefinition>         fJetDefinitions        ; ///<  Jet definitions
    TTree                             *fTree                  ; //!<! Output tree
    AliDmesonInfoSummary              *fCurrentDmesonJetInfo  ; //!<! Current D meson jet info
    AliJetInfoSummary                **fCurrentJetInfo        ; //!<! Current jet info
    vector<AliDmesonJetInfo>           fDmesonJets            ; //!<! Array containing the D meson jets
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
    Bool_t              FindJet(AliAODRecoDecayHF2Prong* Dcand, AliDmesonJetInfo& DmesonJet, AliHFJetDefinition& jetDef);

    /// \cond CLASSIMP
    ClassDef(AnalysisEngine, 2);
    /// \endcond
  };

  AliAnalysisTaskDmesonJets();
  AliAnalysisTaskDmesonJets(const char* name);
  virtual ~AliAnalysisTaskDmesonJets();

  AnalysisEngine* AddAnalysisEngine(ECandidateType_t type, EMCMode_t bkgMode, EJetType_t jettype, Double_t jetradius, TString cutfname = "");
  AnalysisEngine* AddAnalysisEngine(ECandidateType_t type, EMCMode_t bkgMode, const AliHFJetDefinition& jetDef, TString cutfname = "");
  std::list<AnalysisEngine>::iterator FindAnalysisEngine(const AnalysisEngine& eng);

  void SetShowPositionD(Bool_t b = kTRUE)         { fEnabledAxis = b ?  fEnabledAxis | kPositionD         : fEnabledAxis & ~kPositionD         ; }
  void SetShowInvMass(Bool_t b = kTRUE)           { fEnabledAxis = b ?  fEnabledAxis | kInvMass           : fEnabledAxis & ~kInvMass           ; }
  void SetShow2ProngInvMass(Bool_t b = kTRUE)     { fEnabledAxis = b ?  fEnabledAxis | k2ProngInvMass     : fEnabledAxis & ~k2ProngInvMass     ; }
  void SetShowSoftPionPt(Bool_t b = kTRUE)        { fEnabledAxis = b ?  fEnabledAxis | kSoftPionPt        : fEnabledAxis & ~kSoftPionPt        ; }
  void SetShowDeltaR(Bool_t b = kTRUE)            { fEnabledAxis = b ?  fEnabledAxis | kDeltaR            : fEnabledAxis & ~kDeltaR            ; }
  void SetShowDeltaEta(Bool_t b = kTRUE)          { fEnabledAxis = b ?  fEnabledAxis | kDeltaEta          : fEnabledAxis & ~kDeltaEta          ; }
  void SetShowDeltaPhi(Bool_t b = kTRUE)          { fEnabledAxis = b ?  fEnabledAxis | kDeltaPhi          : fEnabledAxis & ~kDeltaPhi          ; }
  void SetShowPositionJet(Bool_t b = kTRUE)       { fEnabledAxis = b ?  fEnabledAxis | kPositionJet       : fEnabledAxis & ~kPositionJet       ; }
  void SetShowJetConstituents(Bool_t b = kTRUE)   { fEnabledAxis = b ?  fEnabledAxis | kJetConstituents   : fEnabledAxis & ~kJetConstituents   ; }

  void SetTreeOutput(Bool_t b)                    { fTreeOutput = b; }

  virtual void         UserCreateOutputObjects();
  virtual void         ExecOnce();
  virtual Bool_t       Run();
  virtual Bool_t       FillHistograms();

 protected:

  AliRDHFCuts*         LoadDMesonCutsFromFile(TString cutfname, TString cutsname);
  void                 AllocateTHnSparse(const AnalysisEngine& param);
  void                 FillTHnSparse(THnSparse* h, const AliDmesonJetInfo& DmesonJet, std::string n);
  
  static const char*   GetHFEventRejectionReasonLabel(UInt_t& bitmap);
  static void          CalculateMassLimits(Double_t range, Int_t pdg, Int_t nbins, Double_t& minMass, Double_t& maxMass);

  list<AnalysisEngine> fAnalysisEngines           ; ///<  Array of analysis parameters
  UInt_t               fEnabledAxis               ; ///<  Use bit defined in EAxis_t to enable axis in the THnSparse
  Bool_t               fTreeOutput                ; ///<  If true, output will be posted in a TTree rather than a THnSparse
  THistManager         fHistManager               ; ///<  Histogram manager
  AliAODEvent         *fAodEvent                  ; //!<! AOD event
  AliFJWrapper        *fFastJetWrapper            ; //!<! Fastjet wrapper

 private:
   
  AliAnalysisTaskDmesonJets(const AliAnalysisTaskDmesonJets &source);
  AliAnalysisTaskDmesonJets& operator=(const AliAnalysisTaskDmesonJets& source);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskDmesonJets, 2);
  /// \endcond
};

#endif
