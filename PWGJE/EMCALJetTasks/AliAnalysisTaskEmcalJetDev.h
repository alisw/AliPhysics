#ifndef ALIANALYSISTASKEMCALJETDEV_H
#define ALIANALYSISTASKEMCALJETDEV_H

// $Id$

class TClonesArray;
class TList;
class TString;
class AliEmcalJet;
class AliRhoParameter;
class AliVCluster;
class AliVParticle;
class AliJetContainer;

#include "AliAnalysisTaskEmcalDev.h"

class AliAnalysisTaskEmcalJetDev : public AliAnalysisTaskEmcalDev {
 public:
  AliAnalysisTaskEmcalJetDev();
  AliAnalysisTaskEmcalJetDev(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliAnalysisTaskEmcalJetDev();

  //these should all point to the jet container
  void                SetJetEtaLimits(Float_t min, Float_t max, Int_t c = 0);
  void                SetJetPhiLimits(Float_t min, Float_t max, Int_t c = 0);
  void                SetJetAreaCut(Float_t cut, Int_t c = 0);
  void                SetPercAreaCut(Float_t p, Int_t c = 0);
  void                SetAreaEmcCut(Double_t a = 0.99, Int_t c = 0);
  void                SetJetPtCut(Float_t cut, Int_t c = 0);
  void                SetJetRadius(Float_t r, Int_t c = 0);
  void                SetMaxClusterPt(Float_t b, Int_t c = 0);
  void                SetMaxTrackPt(Float_t b, Int_t c = 0);
  void                SetPtBiasJetClus(Float_t b, Int_t c = 0);
  void                SetPtBiasJetTrack(Float_t b, Int_t c = 0);
  void                SetLeadingHadronType(Int_t t, Int_t c = 0);
  void                SetNLeadingJets(Int_t t, Int_t c = 0);
  void                SetJetBitMap(UInt_t m, Int_t c = 0);

  void                SetJetsName(const char *n)                   { fJetsName       = n; AddJetContainer(n); }
  virtual void                SetRhoName(const char *n, Int_t c = 0);

  void                        AddJetContainer(const char *n, TString defaultCutType = "");
  void                        RemoveJetContainer(Int_t i)          { fJetCollArray.RemoveAt(i);} 

  AliJetContainer            *GetJetContainer(Int_t i)                const;
  TClonesArray               *GetJetArray(Int_t i)                    const;
  AliEmcalJet                *GetAcceptJetFromArray(Int_t j, Int_t c) const;
  Int_t                       GetNJets(Int_t i)                       const;
  Double_t                    GetRhoVal(Int_t i)                      const;

 protected:
  Float_t*                    GenerateFixedBinArray(Int_t n, Float_t min, Float_t max) const;
  virtual Bool_t              AcceptJet(AliEmcalJet* jet, Int_t c =0);
  Bool_t                      AcceptBiasJet(AliEmcalJet* jet, Int_t c =0);
  Double_t                    GetLeadingHadronPt(AliEmcalJet* jet, Int_t c =0);
  void                        ExecOnce()                                                                    ;

  AliRhoParameter            *GetRhoFromEvent(const char *name)                                             ;
  Bool_t                      GetSortedArray(Int_t indexes[], TClonesArray *array, Double_t rho=0, Int_t c = 0)     ;
  Bool_t                      IsJetTrack(AliEmcalJet* jet, Int_t itrack, Bool_t sorted = kFALSE)       const;
  Bool_t                      IsJetCluster(AliEmcalJet* jet, Int_t iclus, Bool_t sorted = kFALSE)      const;
  Bool_t                      RetrieveEventObjects()                                                        ;

  TString                     fJetsName;                   // name of jet collection
  TString                     fRhoName;                    // Name of rho object

  TClonesArray               *fJets;                       //!jets
  AliRhoParameter            *fRho;                        //!event rho
  Double_t                    fRhoVal;                     //!event rho value

  TObjArray                   fJetCollArray;               // jet collection array

 private:
  AliAnalysisTaskEmcalJetDev(const AliAnalysisTaskEmcalJetDev&);            // not implemented
  AliAnalysisTaskEmcalJetDev &operator=(const AliAnalysisTaskEmcalJetDev&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetDev, 1) // EMCAL Jet base analysis task
};
#endif
