#ifndef ALIANALYSISTASKEMCALTRIGGERPOSITION_H
#define ALIANALYSISTASKEMCALTRIGGERPOSITION_H

#include "AliAnalysisTaskSE.h"
#include "AliInputEventHandler.h"

class THistManager;

class AliAnalysisTaskEmcalTriggerPosition : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskEmcalTriggerPosition();
  AliAnalysisTaskEmcalTriggerPosition(const char *name);
  virtual ~AliAnalysisTaskEmcalTriggerPosition();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *opt);

  void SetEG1Threshold(double threshold) { fThresholdEG1 = threshold; }

protected:

  THistManager                  *fHistos;           //!
  Double_t                       fThresholdEG1;

  ClassDef(AliAnalysisTaskEmcalTriggerPosition, 1);
};

#endif /* ALIANALYSISTASKEMCALTRIGGERPOSITION_H_ */
