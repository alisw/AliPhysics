#ifndef ALIANALYSISTASKLEDCHECK_H
#define ALIANALYSISTASKLEDCHECK_H

#include "AliAnalysisTaskSE.h"

class THistManager;

class AliAnalysisTaskLEDCheck: public AliAnalysisTaskSE {
public:
  AliAnalysisTaskLEDCheck();
  AliAnalysisTaskLEDCheck(const char *name);
  virtual ~AliAnalysisTaskLEDCheck();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);

protected:
  THistManager              *fHistos;

private:
  AliAnalysisTaskLEDCheck(const AliAnalysisTaskLEDCheck &ref);
  AliAnalysisTaskLEDCheck &operator=(const AliAnalysisTaskLEDCheck &);

  ClassDef(AliAnalysisTaskLEDCheck, 1);
};

#endif
