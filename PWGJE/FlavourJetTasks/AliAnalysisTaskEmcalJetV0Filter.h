#ifndef ALIANALYSISTASKEMCALJETV0FILTER_H
#define ALIANALYSISTASKEMCALJETV0FILTER_H

#include <TMap.h>

#include "AliAnalysisTaskEmcalJet.h"
//=============================================================================

class TList;
class TString;
class TClonesArray;

class AliPicoHeaderJet;
//=============================================================================

class AliAnalysisTaskEmcalJetV0Filter : public AliAnalysisTaskEmcalJet {

 public :

  AliAnalysisTaskEmcalJetV0Filter();
  AliAnalysisTaskEmcalJetV0Filter(const char *name, const Bool_t bHistos=kTRUE);
  virtual ~AliAnalysisTaskEmcalJetV0Filter();

  virtual void Init();
  virtual void LocalInit() { Init(); }
  virtual void UserCreateOutputObjects();
  virtual void Terminate(Option_t *opt);
//=============================================================================

  void AddMultEsti(const TString s) {
    if (fMult.IsNull()) {
      fMult = s;
    } else {
      fMult += Form(":%s",s.Data());
    }
    return;
  }

  void SetPicoV0(const TString s, const Bool_t b=kFALSE) {
    fV0sName = s;
    fIsMC   = b;
    return;
  }
//=============================================================================

 protected :

  virtual void   ExecOnce();
  virtual Bool_t FillGeneralHistograms();
  virtual Bool_t FillHistograms();
  virtual Bool_t IsEventSelected();
  virtual Bool_t RetrieveEventObjects();
  virtual Bool_t Run();
//=============================================================================

 private :

  AliAnalysisTaskEmcalJetV0Filter(const AliAnalysisTaskEmcalJetV0Filter &);
  AliAnalysisTaskEmcalJetV0Filter& operator=(const AliAnalysisTaskEmcalJetV0Filter &);
//=============================================================================

  TString fMult; //

  Bool_t fIsMC;     //
  TString fV0sName; //

  TClonesArray *fV0s; //!
//=============================================================================

  AliPicoHeaderJet *fPicoHeader; //!

  TMap fMapJets; //!
  TClonesArray *fPicoV0sClArr;  //!
//=============================================================================

  TList *fListUserOutputs;  //!
//=============================================================================

  ClassDef(AliAnalysisTaskEmcalJetV0Filter, 3);
};

#endif
