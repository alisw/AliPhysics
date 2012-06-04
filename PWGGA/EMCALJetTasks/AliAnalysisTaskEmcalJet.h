#ifndef ALIANALYSISTASKEMCALJET_H
#define ALIANALYSISTASKEMCALJET_H

// $Id: $

class TClonesArray;
class TString;
class TList;
class AliVParticle;
class AliVCluster;
class AliEmcalJet;

#include "AliAnalysisTaskEmcal.h"

class AliAnalysisTaskEmcalJet : public AliAnalysisTaskEmcal {
 public:

  AliAnalysisTaskEmcalJet();
  AliAnalysisTaskEmcalJet(const char *name);
  AliAnalysisTaskEmcalJet(const char *name, Bool_t histo); 
  virtual ~AliAnalysisTaskEmcalJet();

  virtual void                UserExec(Option_t *option);
  virtual void                Init();

  void                        SetJetAreaCut(Float_t cut)                           { fJetAreaCut     = cut        ; }
  void                        SetJetPtCut(Float_t cut)                             { fJetPtCut       = cut        ; }
  void                        SetJetRadius(Float_t r)                              { fJetRadius      = r          ; } 
  void                        SetJetsName(const char *n)                           { fJetsName       = n          ; }
  void                        SetPtBiasJetClus(Float_t b)                          { fPtBiasJetClus  = b          ; }
  void                        SetPtBiasJetTrack(Float_t b)                         { fPtBiasJetTrack = b          ; }
  void                        SetEtaLimits(Float_t min, Float_t max)               { fMinEta = min, fMaxEta = max ; }
  void                        SetPhiLimits(Float_t min, Float_t max)               { fMinPhi = min, fMaxPhi = max ; }
  void                        SetInitialized(Bool_t ini = kTRUE)                   { fInitialized    = ini        ; }

 protected:

  Bool_t                      AcceptJet(AliEmcalJet* jet, Bool_t bias = kTRUE)                     const;
  Bool_t                      IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted = kTRUE)    const;
  Bool_t                      IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted = kTRUE)   const;

  virtual Bool_t              RetrieveEventObjects();

  Float_t                     fJetRadius;                  // jet radius
  TString                     fJetsName;                   // name of jet collection
  Float_t                     fPtBiasJetTrack;             // select jets with a minimum pt track
  Float_t                     fPtBiasJetClus;              // select jets with a minimum pt cluster
  Float_t                     fJetPtCut;                   // cut on jet pt
  Float_t                     fJetAreaCut;                 // cut on jet area
  Float_t                     fMinEta;                     // minimum eta jet accepatance
  Float_t                     fMaxEta;                     // maximum eta jet accepatance
  Float_t                     fMinPhi;                     // minimum phi jet accepatance
  Float_t                     fMaxPhi;                     // maximum phi jet accepatance  

  TClonesArray               *fJets;                       //!jets

 private:
  AliAnalysisTaskEmcalJet(const AliAnalysisTaskEmcalJet&);            // not implemented
  AliAnalysisTaskEmcalJet &operator=(const AliAnalysisTaskEmcalJet&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJet, 2) // EMCAL Jet base analysis task
};
#endif
