#ifndef ALIANALYSISTASKEMCALJET_H
#define ALIANALYSISTASKEMCALJET_H

// $Id$

class TClonesArray;
class TList;
class TString;
class AliEmcalJet;
class AliRhoParameter;
class AliLocalRhoParameter;
class AliVCluster;
class AliVParticle;
class AliJetContainer;

#include "AliAnalysisTaskEmcal.h"

class AliAnalysisTaskEmcalJet : public AliAnalysisTaskEmcal {
 public:
  AliAnalysisTaskEmcalJet();
  AliAnalysisTaskEmcalJet(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliAnalysisTaskEmcalJet();

  //these should all point to the jet container
  void                        SetAnaType(UInt_t t, Int_t c = 0) { SetJetAcceptanceType(t,c); }
  void                        SetJetAcceptanceType(UInt_t t, Int_t c = 0);
  void                        SetJetAcceptanceType(TString cutType, Int_t c = 0);
  void                        SetJetEtaLimits(Float_t min, Float_t max, Int_t c = 0);
  void                        SetJetPhiLimits(Float_t min, Float_t max, Int_t c = 0);
  void                        SetJetAreaCut(Float_t cut, Int_t c = 0);
  void                        SetPercAreaCut(Float_t p, Int_t c = 0);
  void                        SetZLeadingCut(Float_t zemc, Float_t zch, Int_t c = 0);
  void                        SetNEFCut(Float_t min, Float_t max, Int_t c = 0);
  void                        SetAreaEmcCut(Double_t a = 0.99, Int_t c = 0);
  void                        SetJetPtCut(Float_t cut, Int_t c = 0);
  void                        SetJetRadius(Float_t r, Int_t c = 0);
  void                        SetMaxClusterPt(Float_t b, Int_t c = 0);
  void                        SetMaxTrackPt(Float_t b, Int_t c = 0);
  void                        SetPtBiasJetClus(Float_t b, Int_t c = 0);
  void                        SetPtBiasJetTrack(Float_t b, Int_t c = 0);
  void                        SetLeadingHadronType(Int_t t, Int_t c = 0);
  void                        SetNLeadingJets(Int_t t, Int_t c = 0);
  void                        SetJetBitMap(UInt_t m, Int_t c = 0);
  void                        SetJetTrigger(UInt_t t, Int_t c = 0);
  void                        SetIsParticleLevel(Bool_t b, Int_t c = 0);
  void                        SetJetsName(const char *n)                   { AddJetContainer(n); }
  virtual void                SetRhoName(const char *n, Int_t c = 0);
  virtual void                SetLocalRhoName(const char *n)               { fLocalRhoName   = n; }
  const TString&              GetRhoName(Int_t c = 0) const;
  AliJetContainer            *AddJetContainer(const char *n, TString defaultCutType = "", Float_t jetRadius = 0.4);
  void                        RemoveJetContainer(Int_t i)                        { fJetCollArray.RemoveAt(i);} 
  AliJetContainer            *GetJetContainer(Int_t i=0)                                               const;
  AliJetContainer            *GetJetContainer(const char* name)                                        const;

 protected:
  Float_t*                    GenerateFixedBinArray(Int_t n, Float_t min, Float_t max) const;
  virtual Bool_t              AcceptJet(AliEmcalJet* jet, Int_t c =0);
  Bool_t                      AcceptBiasJet(AliEmcalJet* jet, Int_t c =0);
  Double_t                    GetLeadingHadronPt(AliEmcalJet* jet, Int_t c =0);
  void                        ExecOnce()                                                                    ;

  AliRhoParameter            *GetRhoFromEvent(const char *name)                                             ;
  AliLocalRhoParameter       *GetLocalRhoFromEvent(const char *name)                                        ;
  Bool_t                      IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted = kFALSE)       const;
  Bool_t                      IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted = kFALSE)      const;
  Bool_t                      RetrieveEventObjects()                                                        ;
  Double_t                    GetJetRadius(Int_t i=0)                                                  const;
  TClonesArray               *GetJetArray(Int_t i=0)                                                   const;
  AliEmcalJet                *GetJetFromArray(Int_t j, Int_t c=0)                                      const;
  AliEmcalJet                *GetAcceptJetFromArray(Int_t j, Int_t c=0)                                const;
  Int_t                       GetNJets(Int_t i=0)                                                      const;
  Double_t                    GetRhoVal(Int_t i=0)                                                     const;

  TString                     fRhoName;                    // rho name
  TString                     fLocalRhoName;               // name for local rho
  TObjArray                   fJetCollArray;               // jet collection array

  TClonesArray               *fJets;                       //!jets
  AliRhoParameter            *fRho;                        //!event rho
  AliLocalRhoParameter       *fLocalRho;                   //!local event rho
  Double_t                    fRhoVal;                     //!event rho value, same for local rho

 private:
  AliAnalysisTaskEmcalJet(const AliAnalysisTaskEmcalJet&);            // not implemented
  AliAnalysisTaskEmcalJet &operator=(const AliAnalysisTaskEmcalJet&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJet, 3) // EMCAL Jet base analysis task
};
#endif
