#ifndef AliAnalysisSigmaBarChargedGen_H
#define AliAnalysisSigmaBarChargedGen_H

#include "AliAnalysisTaskSE.h"
#include "THashList.h"
#include "TString.h"

class AliAnalysisSigmaBarChargedGen : public AliAnalysisTaskSE {

public:
  AliAnalysisSigmaBarChargedGen();
  AliAnalysisSigmaBarChargedGen(const char *name, TString lExtraOptions = "");
  virtual ~AliAnalysisSigmaBarChargedGen();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

private:
  // outputs
  THashList *fHistos_misc; //!<! Output histos

  AliAnalysisSigmaBarChargedGen(
      const AliAnalysisSigmaBarChargedGen &); // not implemented
  AliAnalysisSigmaBarChargedGen &
  operator=(const AliAnalysisSigmaBarChargedGen &); // not implemented

  ClassDef(AliAnalysisSigmaBarChargedGen, 1);
  // version 1
};

#endif
