#ifndef ALIEMCALJETTASK_H
#define ALIEMCALJETTASK_H

class TClonesArray;
class TObjArray;
class AliVEvent;
class AliEMCALGeometry;
class AliEmcalJetUtility;
class AliFJWrapper;

namespace fastjet {
  class PseudoJet;
}

#include "AliLog.h"
#include "AliAnalysisTaskEmcal.h"
#include "AliFJWrapper.h"
#include "AliAODMCParticle.h"
#include "AliEmcalJet.h"

class AliEmcalJetTask : public AliAnalysisTaskEmcal {
 public:

  enum JetType {
    kNone=0,
    kKT=1<<0,
    kAKT=1<<1,
    kFullJet=1<<2,
    kChargedJet=1<<3,
    kNeutralJet=1<<4,
    kR020Jet=1<<5,
    kR030Jet=1<<6,
    kR040Jet=1<<7,
    kRX1Jet=1<<8,  // user defined radii, use SetRadius(Double_t)
    kRX2Jet=1<<9,
    kRX3Jet=1<<10
  };

  AliEmcalJetTask();
  AliEmcalJetTask(const char *name, Int_t useExchangeCont=0);
  virtual ~AliEmcalJetTask();

  Bool_t Run();

  void                   SetAlgo(Int_t a)                           { if (IsLocked()) return; if (a==0) fJetType |= kKT; else fJetType |= kAKT; }  // for backward compatibility only
  void                   SetClusLabelRange(Int_t min, Int_t max)    { if (IsLocked()) return; fMinLabelClusters = min   ; fMaxLabelClusters = max; }

  void                   SetGhostArea(Double_t gharea)              { if (IsLocked()) return; fGhostArea        = gharea; }
  void                   SetJetsName(const char *n)                 { if (IsLocked()) return; fJetsName         = n     ; }
  void                   SetJetEtaRange(Double_t emi, Double_t ema) { if (IsLocked()) return; fJetEtaMin        = emi   ; fJetEtaMax = ema; }
  void                   SetJetPhiRange(Double_t pmi, Double_t pma) { if (IsLocked()) return; fJetPhiMin        = pmi   ; fJetPhiMax = pma; }
  void                   SetJetType(UInt_t t)                       { if (IsLocked()) return; fJetType          = t     ; }
  void                   SetLocked()                                { fLocked = kTRUE;}
  void                   SetMinJetArea(Double_t a)                  { if (IsLocked()) return; fMinJetArea       = a     ; }
  void                   SetMinJetPt(Double_t j)                    { if (IsLocked()) return; fMinJetPt         = j     ; }
  void                   SetRecombScheme(Int_t scheme)              { if (IsLocked()) return; fRecombScheme     = scheme; }
  void                   SetTrackEfficiency(Double_t t)             { if (IsLocked()) return; fTrackEfficiency  = t     ; }
  void                   SetTrackLabelRange(Int_t min, Int_t max)   { if (IsLocked()) return; fMinLabelTracks   = min   ; fMaxLabelTracks = max; }
  void                   SetLegacyMode(Bool_t mode)                 { if (IsLocked()) return; fLegacyMode       = mode  ; }
  void                   SetFillGhost(Bool_t b=kTRUE)               { if (IsLocked()) return; fFillGhost          = b   ; }

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
  void                   SetRadius(Double_t r);
  void                   SetType(Int_t t);

 protected:

  Int_t                  FindJets();
  void                   FillJetBranch();
  void                   ExecOnce();
  void                   InitUtilities();
  void                   PrepareUtilities();
  void                   ExecuteUtilities(AliEmcalJet* jet, Int_t ij);
  void                   TerminateUtilities();
  Bool_t                 GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;

  TString                fJetsName;               // name of jet collection
  UInt_t                 fJetType;                // jet type (algorithm, radius, constituents)
  Int_t                  fMinLabelTracks;         // select track constituents with a minimum label index (track->GetLabel())
  Int_t                  fMaxLabelTracks;         // select track constituents with a maximum label index (track->GetLabel())
  Int_t                  fMinLabelClusters;       // select cluster constituents with a minimum label index (cluster->GetLabel())
  Int_t                  fMaxLabelClusters;       // select cluster constituents with a maximum label index (cluster->GetLabel())
  Int_t                  fMinMCLabel;             // minimum MC label value for the tracks/clusters being considered MC particles
  Double_t               fRadius;                 // jet radius
  Double_t               fMinJetArea;             // min area to keep jet in output
  Double_t               fMinJetPt;               // min jet pt to keep jet in output
  Double_t               fJetPhiMin;              // minimum phi to keep jet in output
  Double_t               fJetPhiMax;              // maximum phi to keep jet in output
  Double_t               fJetEtaMin;              // minimum eta to keep jet in output
  Double_t               fJetEtaMax;              // maximum eta to keep jet in output
  Double_t               fGhostArea;              // ghost area
  Int_t                  fRecombScheme;           // recombination scheme used by fastjet
  Double_t               fTrackEfficiency;        // artificial tracking inefficiency (0...1)
  TObjArray             *fUtilities;              // jet utilities (gen subtractor, constituent subtractor etc.)
  Int_t                  fUseExchangeCont;        // use exchange containers as input
  Bool_t                 fLocked;                 // true if lock is set

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

  ClassDef(AliEmcalJetTask, 22) // Jet producing task
};
#endif
