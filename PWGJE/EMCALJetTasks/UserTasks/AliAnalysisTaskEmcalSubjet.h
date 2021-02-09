#ifndef ALIANALYSISTASKEMCALSUBJET_H
#define ALIANALYSISTASKEMCALSUBJET_H

#include "AliAnalysisTaskEmcalJet.h"

class THistManager;
class AliEmcalJetFinder;
//=============================================================================

class AliAnalysisTaskEmcalSubjet : public AliAnalysisTaskEmcalJet {

 public :

  AliAnalysisTaskEmcalSubjet();
  AliAnalysisTaskEmcalSubjet(const char *name, const Bool_t bHistos=kTRUE);
  virtual ~AliAnalysisTaskEmcalSubjet();

  virtual void UserCreateOutputObjects();
//virtual void Terminate(Option_t *opt);
//=============================================================================

//AliEmcalJetFinder *GetSubjetFinder() { return fSubjetFinder; }

  void SetSubjetR(Double_t d)       { fSubjetRadius    = d; }
  void SetSubjetAlgorithm(UInt_t w) { fSubjetAlgorithm = w; }
  static AliAnalysisTaskEmcalSubjet *AddTask(const TString sTrks  = "usedefault",
                                             const TString sClus  = "usedefault",
                                             const TString sCells = "usedefault");
//=============================================================================

 protected :

//virtual void UserExecOnce();
//virtual Bool_t RetrieveEventObjects();
//virtual Bool_t IsEventSelected();
//virtual Bool_t FillHistograms();
  virtual Bool_t Run();
//=============================================================================

 private :

  AliAnalysisTaskEmcalSubjet(const AliAnalysisTaskEmcalSubjet &);
  AliAnalysisTaskEmcalSubjet& operator=(const AliAnalysisTaskEmcalSubjet &);
//=============================================================================

  void CreateHistoJets();
  void CreateHistoSubjets();
  void CreateHistoJetConstis();

  void LoopJets      (AliJetContainer const *pc);
  void LoopSubjets   (AliEmcalJet const *pj, AliJetContainer const *pc);
  void LoopJetConstis(AliEmcalJet const *pj, AliJetContainer const *pc);
//=============================================================================

  Double_t           fSubjetRadius; //
  UInt_t             fSubjetAlgorithm; // subjet algo (0: antikt, 1: kt)

  AliEmcalJetFinder *fSubjetFinder; //!
  THistManager      *fHistMgr; //!
//=============================================================================

  ClassDef(AliAnalysisTaskEmcalSubjet, 3)
};

#endif
