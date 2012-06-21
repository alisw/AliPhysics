#ifndef ALIANALYSISTASKEMCALJET_H
#define ALIANALYSISTASKEMCALJET_H

// $Id$

class TClonesArray;
class TList;
class TString;
class AliEmcalJet;
class AliRhoParameter;
class AliVCluster;
class AliVParticle;

#include "AliAnalysisTaskEmcal.h"

class AliAnalysisTaskEmcalJet : public AliAnalysisTaskEmcal {
 public:
  AliAnalysisTaskEmcalJet();
  AliAnalysisTaskEmcalJet(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliAnalysisTaskEmcalJet();

  void                        SetEtaLimits(Float_t min, Float_t max)               { fMinEta = min, fMaxEta = max ; }
  void                        SetJetAreaCut(Float_t cut)                           { fJetAreaCut     = cut        ; }
  void                        SetJetPtCut(Float_t cut)                             { fJetPtCut       = cut        ; }
  void                        SetJetRadius(Float_t r)                              { fJetRadius      = r          ; } 
  void                        SetJetsName(const char *n)                           { fJetsName       = n          ; }
  void                        SetMaxClusterPt(Float_t b)                           { fMaxClusterPt  = b           ; }
  void                        SetMaxTrackPt(Float_t b)                             { fMaxTrackPt = b              ; }
  void                        SetPhiLimits(Float_t min, Float_t max)               { fMinPhi = min, fMaxPhi = max ; }
  void                        SetPtBiasJetClus(Float_t b)                          { fPtBiasJetClus  = b          ; }
  void                        SetPtBiasJetTrack(Float_t b)                         { fPtBiasJetTrack = b          ; }
 
 protected:
  Bool_t                      AcceptJet(AliEmcalJet* jet, Bool_t bias = kTRUE, Bool_t upCut = kTRUE)   const;
  Bool_t                      AcceptBiasJet(AliEmcalJet* jet)                                          const;
  void                        ExecOnce()                                                                    ;
  AliRhoParameter            *GetRhoFromEvent(const char *name);
  Bool_t                      IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted = kTRUE)        const;
  Bool_t                      IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted = kTRUE)       const;
  Bool_t                      RetrieveEventObjects()                                                        ;

  Float_t                     fJetRadius;                  // jet radius
  TString                     fJetsName;                   // name of jet collection
  Float_t                     fPtBiasJetTrack;             // select jets with a minimum pt track
  Float_t                     fPtBiasJetClus;              // select jets with a minimum pt cluster
  Float_t                     fJetPtCut;                   // cut on jet pt
  Float_t                     fJetAreaCut;                 // cut on jet area
  Float_t                     fMinEta;                     // minimum eta jet acceptance
  Float_t                     fMaxEta;                     // maximum eta jet acceptance
  Float_t                     fMinPhi;                     // minimum phi jet acceptance
  Float_t                     fMaxPhi;                     // maximum phi jet acceptance  
  Float_t                     fMaxClusterPt;               // maximum cluster constituent pt to accept the jet
  Float_t                     fMaxTrackPt;                 // maximum track constituent pt to accept the jet
  TClonesArray               *fJets;                       //!jets

 private:
  AliAnalysisTaskEmcalJet(const AliAnalysisTaskEmcalJet&);            // not implemented
  AliAnalysisTaskEmcalJet &operator=(const AliAnalysisTaskEmcalJet&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJet, 3) // EMCAL Jet base analysis task
};
#endif
