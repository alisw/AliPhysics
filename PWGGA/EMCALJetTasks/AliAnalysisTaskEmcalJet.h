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

  void                        SetJetAreaCut(Float_t cut)                           { fJetAreaCut     = cut        ; }
  void                        SetJetPtCut(Float_t cut)                             { fJetPtCut       = cut        ; }
  void                        SetJetRadius(Float_t r)                              { fJetRadius      = r          ; } 
  void                        SetJetsName(const char *n)                           { fJetsName       = n          ; }
  void                        SetPtBiasJetClus(Float_t b)                          { fPtBiasJetClus  = b          ; }
  void                        SetPtBiasJetTrack(Float_t b)                         { fPtBiasJetTrack = b          ; }

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

  TClonesArray               *fJets;                       //!jets

 private:
  AliAnalysisTaskEmcalJet(const AliAnalysisTaskEmcalJet&);            // not implemented
  AliAnalysisTaskEmcalJet &operator=(const AliAnalysisTaskEmcalJet&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJet, 1) // EMCAL Jet base analysis task
};
#endif
