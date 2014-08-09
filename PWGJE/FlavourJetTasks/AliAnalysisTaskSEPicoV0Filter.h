#ifndef ALIANALYSISTASKSEPICOV0FILTER_H
#define ALIANALYSISTASKSEPICOV0FILTER_H
//=============================================================================

#include "AliAnalysisTaskSE.h"

class TString;
class TClonesArray;

class AliPicoHeaderCJ;
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

  void SetAnaInfoMC(Bool_t b=kTRUE) { fIsAnaInfoMC = b; }
//=============================================================================

 private :

  AliAnalysisTaskSEPicoV0Filter(const AliAnalysisTaskSEPicoV0Filter &);
  AliAnalysisTaskSEPicoV0Filter& operator=(const AliAnalysisTaskSEPicoV0Filter &);

  void CreateUserOutputHistograms();

  Bool_t fIsAnaInfoMC; //

  TClonesArray *fV0s; //!

  AliPicoHeaderCJ *fPicoHeaderCJ; //!

  TClonesArray *fPicoV0sClArr; //!

  TList *fListUserOutputs; //!

  ClassDef(AliAnalysisTaskSEPicoV0Filter, 1);
};

#endif
