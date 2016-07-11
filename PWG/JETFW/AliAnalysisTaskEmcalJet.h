#ifndef ALIANALYSISTASKEMCALJET_H
#define ALIANALYSISTASKEMCALJET_H
/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

class TClonesArray;
class TList;
class TString;
class AliEmcalJet;
class AliRhoParameter;
class AliLocalRhoParameter;
class AliVCluster;
class AliVParticle;

#include "AliJetContainer.h"
#include "AliAnalysisTaskEmcal.h"

/**
 * @class AliAnalysisTaskEmcalJet
 * @brief Base task in the EMCAL jet framework
 * @author Marta Verweij
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 *
 * This class is the base class for Analysis Tasks using the
 * core EMCAL jet framework. User tasks should inherit from it.
 * This class derives from AliAnalysisTaskEmcal and adds
 * some additional features, useful for jet analysis. The key feature
 * is the handling of the jet containers.
 */
class AliAnalysisTaskEmcalJet : public AliAnalysisTaskEmcal {
 public:
  typedef AliJetContainer::EJetType_t EJetType_t;
  typedef AliJetContainer::EJetAlgo_t EJetAlgo_t;
  typedef AliJetContainer::ERecoScheme_t ERecoScheme_t;

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
  virtual void                SetRhoName(const char *n, Int_t c = 0);
  virtual void                SetLocalRhoName(const char *n)               { fLocalRhoName   = n; }
  const TString&              GetRhoName(Int_t c = 0) const;
  AliJetContainer            *AddJetContainer(const char *n, TString defaultCutType, Float_t jetRadius = 0.4);
  AliJetContainer            *AddJetContainer(const char *n, UInt_t accType = AliJetContainer::kUser, Float_t jetRadius = 0.4);
  AliJetContainer            *AddJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius,
      UInt_t accType, AliParticleContainer* partCont, AliClusterContainer* clusCont, TString tag = "Jet");
  AliJetContainer            *AddJetContainer(EJetType_t jetType, EJetAlgo_t jetAlgo, ERecoScheme_t recoScheme, Double_t radius,
      UInt_t accType, TString tag = "Jet");
  void                        AdoptJetContainer(AliJetContainer* cont)           { fJetCollArray.Add(cont)  ;}

  void                        RemoveJetContainer(Int_t i)                        { fJetCollArray.RemoveAt(i);} 
  AliJetContainer            *GetJetContainer(Int_t i=0)                                               const;
  AliJetContainer            *GetJetContainer(const char* name)                                        const;

 protected:
  virtual Bool_t              AcceptJet(AliEmcalJet* jet, Int_t c =0);
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

  TString                     fRhoName;                    ///<  rho name
  TString                     fLocalRhoName;               ///<  name for local rho
  TObjArray                   fJetCollArray;               ///<  jet collection array

  TClonesArray               *fJets;                       //!<! jets
  AliRhoParameter            *fRho;                        //!<! event rho
  AliLocalRhoParameter       *fLocalRho;                   //!<! local event rho
  Double_t                    fRhoVal;                     //!<! event rho value, same for local rho

 private:
  AliAnalysisTaskEmcalJet(const AliAnalysisTaskEmcalJet&);            // not implemented
  AliAnalysisTaskEmcalJet &operator=(const AliAnalysisTaskEmcalJet&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJet, 3);
  /// \endcond
};
#endif
