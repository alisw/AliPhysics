#ifndef ALIEMCALJETTASK_H
#define ALIEMCALJETTASK_H

class TClonesArray;
class AliVEvent;
class AliRhoParameter;

namespace fastjet {
  class PseudoJet;
}

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliFJWrapper.h"

class AliEmcalJetTask : public AliAnalysisTaskSE {
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
  AliEmcalJetTask(const char *name);
  virtual ~AliEmcalJetTask();

  void                   UserCreateOutputObjects();
  void                   UserExec(Option_t *option);
  void                   Terminate(Option_t *option);

  void                   SelectConstituents(UInt_t constSel, UInt_t MCconstSel)  { fConstSel = constSel; fMCConstSel = MCconstSel; };
  void                   SetAlgo(Int_t a)                 { if (a==0) fJetType |= kKT; else fJetType |= kAKT; }  // for backward compatibility only
  void                   SetClusName(const char *n)       { fCaloName      = n     ; }
  void                   SetJetsName(const char *n)       { fJetsName      = n     ; }
  void                   SetJetsSubName(const char *n)    { fJetsSubName   = n     ; }
  void                   SetJetType(UInt_t t)             { fJetType       = t     ; }
  void                   SetMarkConstituents(UInt_t m)    { fMarkConst     = m     ; }
  void                   SetMinJetArea(Double_t a)        { fMinJetArea    = a     ; }
  void                   SetMinJetClusPt(Double_t min)    { fMinJetClusPt  = min   ; }
  void                   SetMinJetPt(Double_t j)          { fMinJetPt      = j     ; }
  void                   SetMinJetTrackPt(Double_t min)   { fMinJetTrackPt = min   ; }
  void                   SetRadius(Double_t r)            { fRadius        = r     ; if ((fJetType & (kRX1Jet|kRX2Jet|kRX3Jet)) == 0) AliWarning("Radius value will be ignored if jet type is not set to a user defined radius (kRX1Jet,kRX2Jet,kRX3Jet)."); }
  void                   SetTrackEfficiency(Double_t t)   { fTrackEfficiency = t   ; }
  void                   SetTracksName(const char *n)     { fTracksName    = n     ; }
  void                   SetType(Int_t t)                 { if (t==0) fJetType |= kFullJet; 
                                                            else if (t==1) fJetType |= kChargedJet; 
                                                            else if (t==2) fJetType |= kNeutralJet; } // for backward compatibility only
  void                   SetEtaRange(Double_t emi, Double_t ema) {fEtaMin = emi; fEtaMax = ema; }
  void                   SetPhiRange(Double_t pmi, Double_t pma) {fPhiMin = pmi; fPhiMax = pma; }
  void                   SetJetEtaRange(Double_t emi, Double_t ema) {fJetEtaMin = emi; fJetEtaMax = ema; }
  void                   SetJetPhiRange(Double_t pmi, Double_t pma) {fJetPhiMin = pmi; fJetPhiMax = pma; }
  void                   SetGhostArea(Double_t gharea)    { fGhostArea      = gharea;  }
  void                   SetMinMCLabel(Int_t s)           { fMinMCLabel     = s     ;  }
  void                   SetRecombScheme(Int_t scheme)    { fRecombScheme   = scheme;  }
  void                   SelectCollisionCandidates(UInt_t offlineTriggerMask = AliVEvent::kMB)
  {
    if(!fIsPSelSet)
    {
      fIsPSelSet = kTRUE;
      fOfflineTriggerMask = offlineTriggerMask;
    }
    else
    {
      AliWarning("Manually setting the event selection for jet finders is not allowed! Using trigger=old_trigger | your_trigger");  
      fOfflineTriggerMask = fOfflineTriggerMask | offlineTriggerMask;
    }
  }
  void                   SetLegacyMode(Bool_t mode)       { fLegacyMode = mode; }
  void                   SetCodeDebug(Bool_t val)         { fCodeDebug = val; }

  void                   SetRhoName(const char *n)              { fRhoName      = n            ; }
  void                   SetRhomName(const char *n)             { fRhomName     = n            ; }
  void                   SetGenericSubtraction(Bool_t b)        { fDoGenericSubtraction     = b; }
  void                   SetConstituentSubtraction(Bool_t b)    { fDoConstituentSubtraction = b; }
  void                   SetUseExternalBkg(Bool_t b, Double_t rho, Double_t rhom) { fUseExternalBkg = b; fRho = rho; fRhom = rhom;}

  UInt_t                 GetJetType()                     { return fJetType; }
  Bool_t                 GetLegacyMode()                  { return fLegacyMode; }
  TString                GetJetsSubName()                 { return fJetsSubName; }

  const char*            GetJetsName()                    { return fJetsName.Data(); }
  Double_t               GetRadius()                      { return fRadius; }
  const char*            GetTracksName()                  { return fTracksName.Data(); }
  const char*            GetClusName()                    { return fCaloName.Data(); }
  UInt_t                 GetMarkConstituents()            { return fMarkConst; }
  Double_t               GetMinJetArea()                  { return fMinJetArea; }
  Double_t               GetMinJetClusPt()                { return fMinJetClusPt; }
  Double_t               GetMinJetPt()                    { return fMinJetPt; }
  Double_t               GetMinJetTrackPt()               { return fMinJetTrackPt; }
  Double_t               GetTrackEfficiency()             { return fTrackEfficiency; }
  Double_t               GetGhostArea()                   { return fGhostArea; }
  Int_t                  GetMinMCLabel()                  { return fMinMCLabel; }
  Int_t                  GetRecombScheme()                { return fRecombScheme; }
  Double_t               GetEtaMin()                      { return fEtaMin; }
  Double_t               GetEtaMax()                      { return fEtaMax; }
  Double_t               GetPhiMin()                      { return fPhiMin; }
  Double_t               GetPhiMax()                      { return fPhiMax; }
  Double_t               GetJetEtaMin()                   { return fJetEtaMin; }
  Double_t               GetJetEtaMax()                   { return fJetEtaMax; }
  Double_t               GetJetPhiMin()                   { return fJetPhiMin; }
  Double_t               GetJetPhiMax()                   { return fJetPhiMax; }

  AliVEvent*             GetEvent()                       { return fEvent;}
  TClonesArray*          GetClusters()                    { return fClus;}
  TClonesArray*          GetTracks()                      { return fTracks;}
  TClonesArray*          GetJets()                        { return fJets;}

 protected:
  void                   FindJets();
  Bool_t                 DoInit();
  Bool_t                 GetSortedArray(Int_t indexes[], std::vector<fastjet::PseudoJet> array) const;

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fJetsName;               // name of jet collection
  TString                fJetsSubName;            // name of subtracted jet collection
  UInt_t                 fJetType;                // jet type (algorithm, radius, constituents)
  UInt_t                 fConstSel;               // select constituents from a previous jet finding
  UInt_t                 fMCConstSel;             // select MC constituents (label!=0) from a previous jet finding
  UInt_t                 fMarkConst;              // constituents are marked (via TObject::SetBit) as belonging to the # leading jets
  Double_t               fRadius;                 // jet radius
  Double_t               fMinJetTrackPt;          // min jet track momentum   (applied before clustering)
  Double_t               fMinJetClusPt;           // min jet cluster momentum (applied before clustering)
  Double_t               fPhiMin;                 // minimum phi for constituents (applied before clustering)
  Double_t               fPhiMax;                 // maximum phi for constituents (applied before clustering)
  Double_t               fEtaMin;                 // minimum eta for constituents (applied before clustering)
  Double_t               fEtaMax;                 // maximum eta for constituents (applied before clustering)
  Double_t               fMinJetArea;             // min area to keep jet in output
  Double_t               fMinJetPt;               // min jet pt to keep jet in output
  Double_t               fJetPhiMin;              // minimum phi to keep jet in output
  Double_t               fJetPhiMax;              // maximum phi to keep jet in output
  Double_t               fJetEtaMin;              // minimum eta to keep jet in output
  Double_t               fJetEtaMax;              // maximum eta to keep jet in output
  Double_t               fGhostArea;              // ghost area
  Int_t                  fMinMCLabel;             // minimum MC label value for the tracks/clusters being considered MC particles
  Int_t                  fRecombScheme;           // recombination scheme used by fastjet
  Double_t               fTrackEfficiency;        // artificial tracking inefficiency (0...1)
  Bool_t                 fIsInit;                 //!=true if already initialized
  Bool_t                 fIsPSelSet;              //!=true if physics selection was set
  Bool_t                 fIsMcPart;               //!=true if MC particles are given as input
  Bool_t                 fIsEmcPart;              //!=true if emcal particles are given as input (for clusters)
  Bool_t                 fLegacyMode;             //! if true, enable FJ 2.x behavior
  Bool_t                 fCodeDebug;              // use nontested code changes 

  Bool_t                 fDoGenericSubtraction;   // calculate generic subtraction
  Bool_t                 fDoConstituentSubtraction; // calculate constituent subtraction
  Bool_t                 fUseExternalBkg;         // use external background for generic subtractor
  TString                fRhoName;                // name of rho
  TString                fRhomName;               // name of rhom
  Double_t               fRho;                    // pT background density
  Double_t               fRhom;                   // mT background density

  TClonesArray          *fJets;                   //!jet collection
  TClonesArray          *fJetsSub;                //!subtracted jet collection
  AliVEvent             *fEvent;                  //!current event
  TClonesArray          *fTracks;                 //!tracks collection
  TClonesArray          *fClus;                   //!cluster collection
  AliRhoParameter       *fRhoParam;               //!event rho
  AliRhoParameter       *fRhomParam;              //!event rhom

 private:
  AliEmcalJetTask(const AliEmcalJetTask&);            // not implemented
  AliEmcalJetTask &operator=(const AliEmcalJetTask&); // not implemented

  ClassDef(AliEmcalJetTask, 12) // Jet producing task
};
#endif
