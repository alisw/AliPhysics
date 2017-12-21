#ifndef ALIANALYSISTASKSEPICOV0FILTER_H
#define ALIANALYSISTASKSEPICOV0FILTER_H
//=============================================================================

#include <TString.h>

#include "AliAnalysisTaskSE.h"

class TClonesArray;

class AliPicoHeaderV0;
//=============================================================================

class AliAnalysisTaskSEPicoV0Filter : public AliAnalysisTaskSE {

 public :

  AliAnalysisTaskSEPicoV0Filter();
  AliAnalysisTaskSEPicoV0Filter(const char *name);
  virtual ~AliAnalysisTaskSEPicoV0Filter();

  virtual void Init();
  virtual void LocalInit() { Init(); }
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);
  virtual void Terminate(Option_t *opt);

  void SetAnaInfoMC(Bool_t b=kTRUE) { fIsMC = b; }

  void AddMultEsti(const TString s) {
    if (fMult.IsNull()) {
      fMult = s;
      if (fMultEstDef.IsNull()) fMultEstDef = s;
    } else {
      fMult += Form(":%s",s.Data());
    }
    return;
  }

  void SetMultRange(const Double_t dMin,
                    const Double_t dMax,
                    const TString sEst="V0M") {
    fCutMinMult = dMin;
    fCutMaxMult = dMax;
    fMultEstDef = sEst;
    return;
  }
//=============================================================================

 private :

  AliAnalysisTaskSEPicoV0Filter(const AliAnalysisTaskSEPicoV0Filter &);
  AliAnalysisTaskSEPicoV0Filter& operator=(const AliAnalysisTaskSEPicoV0Filter &);

  Bool_t fIsMC;  //
  TString fMult; //

  TString  fMultEstDef; //
  Double_t fCutMinMult; //
  Double_t fCutMaxMult; //

  TClonesArray *fV0s; //!
  AliPicoHeaderV0 *fPicoHeader; //!

  TClonesArray *fPicoV0sClArr; //!
  TList *fListUserOutputs;     //!

  ClassDef(AliAnalysisTaskSEPicoV0Filter, 2);
};

#endif
