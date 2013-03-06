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

  void                        SetJetEtaLimits(Float_t min, Float_t max)            { fJetMinEta = min, fJetMaxEta = max ; }
  void                        SetJetPhiLimits(Float_t min, Float_t max)            { fJetMinPhi = min, fJetMaxPhi = max ; }
  void                        SetJetAreaCut(Float_t cut)                           { fJetAreaCut     = cut              ; }
  void                        SetPercAreaCut(Float_t p)                            { fPercAreaCut    = p                ; }
  void                        SetAreaEmcCut(Double_t a = 0.99)                     { fAreaEmcCut     = a                ; }
  void                        SetJetPtCut(Float_t cut)                             { fJetPtCut       = cut              ; }
  void                        SetJetRadius(Float_t r)                              { fJetRadius      = r                ; } 
  void                        SetJetsName(const char *n)                           { fJetsName       = n                ; }
  virtual void                SetRhoName(const char *n)                            { fRhoName        = n                ; }
  void                        SetMaxClusterPt(Float_t b)                           { fMaxClusterPt   = b                ; }
  void                        SetMaxTrackPt(Float_t b)                             { fMaxTrackPt     = b                ; }
  void                        SetPtBiasJetClus(Float_t b)                          { fPtBiasJetClus  = b                ; }
  void                        SetPtBiasJetTrack(Float_t b)                         { fPtBiasJetTrack = b                ; }
  void                        SetLeadingHadronType(Int_t t)                        { fLeadingHadronType = t             ; }
  void                        SetNLeadingJets(Int_t t)                             { fNLeadingJets   = t                ; }
  void                        SetJetBitMap(UInt_t m)                               { fJetBitMap      = m                ; }
 
 protected:
  Float_t*                    GenerateFixedBinArray(Int_t n, Float_t min, Float_t max) const;
  virtual Bool_t              AcceptJet(AliEmcalJet* jet)                                              const;
  Bool_t                      AcceptBiasJet(AliEmcalJet* jet)                                          const;
  Double_t                    GetLeadingHadronPt(AliEmcalJet* jet)                                     const;
  void                        ExecOnce()                                                                    ;
  AliRhoParameter            *GetRhoFromEvent(const char *name)                                             ;
  Bool_t                      GetSortedArray(Int_t indexes[], TClonesArray *array, Double_t rho=0)     const;
  Bool_t                      IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted = kTRUE)        const;
  Bool_t                      IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted = kTRUE)       const;
  Bool_t                      RetrieveEventObjects()                                                        ;

  Float_t                     fJetRadius;                  // jet radius
  TString                     fJetsName;                   // name of jet collection
  TString                     fRhoName;                    // Name of rho object
  Float_t                     fPtBiasJetTrack;             // select jets with a minimum pt track
  Float_t                     fPtBiasJetClus;              // select jets with a minimum pt cluster
  Float_t                     fJetPtCut;                   // cut on jet pt
  Float_t                     fJetAreaCut;                 // cut on jet area
  Float_t                     fPercAreaCut;                // cut on jet area as a percentage of average jet area
  Float_t                     fAreaEmcCut;                 // minimum cut on jet emcal area
  Float_t                     fJetMinEta;                  // minimum eta jet acceptance
  Float_t                     fJetMaxEta;                  // maximum eta jet acceptance
  Float_t                     fJetMinPhi;                  // minimum phi jet acceptance
  Float_t                     fJetMaxPhi;                  // maximum phi jet acceptance  
  Float_t                     fMaxClusterPt;               // maximum cluster constituent pt to accept the jet
  Float_t                     fMaxTrackPt;                 // maximum track constituent pt to accept the jet
  Int_t                       fLeadingHadronType;          // 0 = charged, 1 = neutral, 2 = both
  Int_t                       fNLeadingJets;               // how many jets are to be considered the leading jet(s)
  UInt_t                      fJetBitMap;                  // bit map of accepted jets
  TClonesArray               *fJets;                       //!jets
  AliRhoParameter            *fRho;                        //!event rho
  Double_t                    fRhoVal;                     //!event rho value

 private:
  AliAnalysisTaskEmcalJet(const AliAnalysisTaskEmcalJet&);            // not implemented
  AliAnalysisTaskEmcalJet &operator=(const AliAnalysisTaskEmcalJet&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJet, 7) // EMCAL Jet base analysis task
};
#endif
