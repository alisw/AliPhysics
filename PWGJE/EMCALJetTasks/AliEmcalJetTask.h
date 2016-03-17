#ifndef ALIEMCALJETTASK_H
#define ALIEMCALJETTASK_H

class TClonesArray;
class TObjArray;
class AliVEvent;
class AliEMCALGeometry;
class AliEmcalJetUtility;

#include "AliLog.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliFJWrapper.h"
#include "FJ_includes.h"
#include "AliAODMCParticle.h"
#include "AliEmcalJet.h"
#include "AliJetContainer.h"

namespace fastjet {
  class PseudoJet;
}

class AliEmcalJetTask : public AliAnalysisTaskEmcal {
 public:

  typedef AliJetContainer::EJetType_t EJetType_t;
  typedef AliJetContainer::EJetAlgo_t EJetAlgo_t;
  typedef AliJetContainer::ERecoScheme_t ERecoScheme_t;

#if !defined(__CINT__) && !defined(__MAKECINT__)
  typedef fastjet::JetAlgorithm FJJetAlgo;
  typedef fastjet::RecombinationScheme FJRecoScheme;
#endif

  AliEmcalJetTask();
  AliEmcalJetTask(const char *name);
  virtual ~AliEmcalJetTask();

  Bool_t Run();

  void                   SetGhostArea(Double_t gharea)              { if (IsLocked()) return; fGhostArea        = gharea; }
  void                   SetJetsName(const char *n)                 { if (IsLocked()) return; fJetsTag          = n     ; }
  void                   SetJetEtaRange(Double_t emi, Double_t ema) { if (IsLocked()) return; fJetEtaMin        = emi   ; fJetEtaMax = ema; }
  void                   SetJetPhiRange(Double_t pmi, Double_t pma) { if (IsLocked()) return; fJetPhiMin        = pmi   ; fJetPhiMax = pma; }
  void                   SetJetAlgo(EJetAlgo_t a)                   { if (IsLocked()) return; fJetAlgo          = a     ; }
  void                   SetJetType(EJetType_t t)                   { if (IsLocked()) return; fJetType          = t     ; }
  void                   SetLocked()                                { fLocked = kTRUE;}
  void                   SetMinJetArea(Double_t a)                  { if (IsLocked()) return; fMinJetArea       = a     ; }
  void                   SetMinJetPt(Double_t j)                    { if (IsLocked()) return; fMinJetPt         = j     ; }
  void                   SetRecombScheme(ERecoScheme_t scheme)      { if (IsLocked()) return; fRecombScheme     = scheme; }
  void                   SetTrackEfficiency(Double_t t)             { if (IsLocked()) return; fTrackEfficiency  = t     ; }
  void                   SetLegacyMode(Bool_t mode)                 { if (IsLocked()) return; fLegacyMode       = mode  ; }
  void                   SetFillGhost(Bool_t b=kTRUE)               { if (IsLocked()) return; fFillGhost        = b     ; }
  void                   SetRadius(Double_t r)                      { if (IsLocked()) return; fRadius           = r     ; }

  void                   SetEtaRange(Double_t emi, Double_t ema);
  void                   SetMinJetClusPt(Double_t min);
  void                   SetMinJetClusE(Double_t min);
  void                   SetMinJetTrackPt(Double_t min);
  void                   SetPhiRange(Double_t pmi, Double_t pma);

  AliEmcalJetUtility*    AddUtility(AliEmcalJetUtility* utility);

  Double_t               GetGhostArea()                   { return fGhostArea         ; }
  const char*            GetJetsName()                    { return fJetsName.Data()   ; }
  Double_t               GetJetEtaMin()                   { return fJetEtaMin         ; }
  Double_t               GetJetEtaMax()                   { return fJetEtaMax         ; }
  Double_t               GetJetPhiMin()                   { return fJetPhiMin         ; }
  Double_t               GetJetPhiMax()                   { return fJetPhiMax         ; }
  UInt_t                 GetJetType()                     { return fJetType           ; }
  UInt_t                 GetJetAlgo()                     { return fJetAlgo           ; }
  Bool_t                 GetLegacyMode()                  { return fLegacyMode        ; }
  Double_t               GetMinJetArea()                  { return fMinJetArea        ; }
  Double_t               GetMinJetPt()                    { return fMinJetPt          ; }
  Int_t                  GetMinMCLabel()                  { return fMinMCLabel        ; }
  Double_t               GetRadius()                      { return fRadius            ; }
  Int_t                  GetRecombScheme()                { return fRecombScheme      ; }
  Double_t               GetTrackEfficiency()             { return fTrackEfficiency   ; }

  TClonesArray*          GetJets()                        { return fJets              ; }
  TObjArray*             GetUtilities()                   { return fUtilities         ; }

  void                   FillJetConstituents(AliEmcalJet *jet, std::vector<fastjet::PseudoJet>& constituents,
                                             std::vector<fastjet::PseudoJet>& constituents_sub, Int_t flag = 0, TClonesArray *particles_sub = 0);

  Int_t                  GetIndexSub(Double_t phi_sub, std::vector<fastjet::PseudoJet>& constituents_unsub);

  Bool_t                 IsLocked() const;
  void                   SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kMB);
  void                   SetType(Int_t t);

#if !defined(__CINT__) && !defined(__MAKECINT__)
  static FJJetAlgo       ConvertToFJAlgo(EJetAlgo_t algo);
  static FJRecoScheme    ConvertToFJRecoScheme(ERecoScheme_t reco);
#endif

 protected:

  Int_t                  FindJets();
  void                   FillJetBranch();
  void                   ExecOnce();
  void                   InitUtilities();
  void                   PrepareUtilities();
  void                   ExecuteUtilities(AliEmcalJet* jet, Int_t ij);
  void                   TerminateUtilities();
  Bool_t                 GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;

  TString                fJetsTag;                // tag of jet collection (usually = "Jets")

  EJetType_t             fJetType;                // jet type (full, charged, neutral)
  EJetAlgo_t             fJetAlgo;                // jet algorithm (kt, akt, etc)
  ERecoScheme_t          fRecombScheme;           // recombination scheme used by fastjet
  Double_t               fRadius;                 // jet radius
  Double_t               fMinJetArea;             // min area to keep jet in output
  Double_t               fMinJetPt;               // min jet pt to keep jet in output
  Double_t               fJetPhiMin;              // minimum phi to keep jet in output
  Double_t               fJetPhiMax;              // maximum phi to keep jet in output
  Double_t               fJetEtaMin;              // minimum eta to keep jet in output
  Double_t               fJetEtaMax;              // maximum eta to keep jet in output
  Double_t               fGhostArea;              // ghost area
  Double_t               fTrackEfficiency;        // artificial tracking inefficiency (0...1)
  TObjArray             *fUtilities;              // jet utilities (gen subtractor, constituent subtractor etc.)
  Bool_t                 fLocked;                 // true if lock is set

  TString                fJetsName;               //!name of jet collection
  Bool_t                 fIsInit;                 //!=true if already initialized
  Bool_t                 fIsPSelSet;              //!=true if physics selection was set
  Bool_t                 fIsEmcPart;              //!=true if emcal particles are given as input (for clusters)
  Bool_t                 fLegacyMode;             //!=true to enable FJ 2.x behavior
  Bool_t                 fFillGhost;              //!=true ghost particles will be filled in AliEmcalJet obj

  TClonesArray          *fJets;                   //!jet collection
  AliFJWrapper           fFastJetWrapper;         //!fastjet wrapper

  static const Int_t     fgkConstIndexShift;      //!contituent index shift

 private:
  AliEmcalJetTask(const AliEmcalJetTask&);            // not implemented
  AliEmcalJetTask &operator=(const AliEmcalJetTask&); // not implemented

  ClassDef(AliEmcalJetTask, 23) // Jet producing task
};
#endif
