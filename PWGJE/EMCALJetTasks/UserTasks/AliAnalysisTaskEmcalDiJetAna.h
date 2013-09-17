#ifndef ALIANALYSISTASKEMCALDIJETANA_H
#define ALIANALYSISTASKEMCALDIJETANA_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisUtils;
class AliAnalysisManager;
class AliGenPythiaEventHeader;

#include "AliJetContainer.h"

#include "AliAnalysisTaskEmcalDiJetBase.h"

class AliAnalysisTaskEmcalDiJetAna : public AliAnalysisTaskEmcalDiJetBase {
 public:
  AliAnalysisTaskEmcalDiJetAna();
  AliAnalysisTaskEmcalDiJetAna(const char *name);
  virtual ~AliAnalysisTaskEmcalDiJetAna();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  //Setters
  void                        SetMatchFullCharged(Bool_t b)        { fDoMatchFullCharged = b;}


  //Getters
  Int_t                       GetPtTriggerBin(Double_t pt);

 protected:
  Bool_t                      Run()              ;
  void                        CorrelateJets(const Int_t type);
  void                        CorrelateAllJets(const Int_t type);
  void                        CorrelateTwoJets(const Int_t type);
  void                        CorrelateLeadingSubleadingJets(const Int_t type);
  AliEmcalJet                *GetLeadingJet(const Int_t type);
  AliEmcalJet                *GetLeadingAssociatedJet(const Int_t type, AliEmcalJet *jetTrig);

  Bool_t                      FillHistograms()   ;
  void                        FillDiJetHistos(const AliEmcalJet *jet1 = 0, const AliEmcalJet *jet2 = 0, const Int_t mode = 0);
  Bool_t                      RetrieveEventObjects();

  void                        FillMatchFullChargedHistos(Int_t cFull,Int_t cCharged);
  Int_t                       MatchFullAndChargedJets(Int_t cFull, Int_t cCharged);

 private:
  Bool_t            fDoMatchFullCharged;                  //  do full-charged matching histos
  TH2F             *fh2CentRhoCh;                         //! cent vs rho charged
  TH2F             *fh2CentRhoScaled;                     //! cent vs rho scaled
  TH3F             *fh3PtEtaPhiJetFull;                   //! pt,eta,phi of full jets
  TH3F             *fh3PtEtaPhiJetCharged;                //! pt,eta,phi of charged jets

  THnSparse        *fhnDiJetVarsFull;                     //! sparse with di-jet properties (full-full)
  THnSparse        *fhnDiJetVarsCh;                       //! sparse with di-jet properties (charged-charged)
  THnSparse        *fhnDiJetVarsFullCharged;              //! sparse with di-jet properties (full-charged)
  THnSparse        *fhnMatchingFullCharged;               //! sparse comparing full with matched charged jet
  TH3F             *fh3JetPtFullFractionDR;               //! full jet pt vs highest shared charged fraction vs DeltaR

  TH3F             *fh3DiJetKtNEFPtAssoc[4];              //! dijet kt vs NEF vs pTassoc for 4 trigger intervals

  TH3F             *fCentCorrPtAssocCh[4];                //! default(V0A) vs ZNA centrality vs pT trigger assoc
  TH3F             *fCentCorrPtAssocFuCh[4];              //! default(V0A) vs ZNA centrality vs pT trigger assoc

  TH3F             *fAjPtAssocCentCh[4];                  //! Aj vs pT trigger assoc vs centrality
  TH3F             *fAjPtAssocCentFuCh[4];                //! Aj vs pT trigger assoc vs centrality


  AliAnalysisTaskEmcalDiJetAna(const AliAnalysisTaskEmcalDiJetAna&);            // not implemented
  AliAnalysisTaskEmcalDiJetAna &operator=(const AliAnalysisTaskEmcalDiJetAna&); // not implemented

  ClassDef(AliAnalysisTaskEmcalDiJetAna, 7) // dijet analysis task
};
#endif
