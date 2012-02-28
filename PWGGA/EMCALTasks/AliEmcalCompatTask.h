#ifndef ALIEMCALCOMPATTASK_H
#define ALIEMCALCOMPATTASK_H

// $Id$

#include "AliAnalysisTaskSE.h"

class AliEmcalCompatTask : public AliAnalysisTaskSE {
 public:
  AliEmcalCompatTask();
  AliEmcalCompatTask(const char *name);
  virtual ~AliEmcalCompatTask();

  void UserExec(Option_t *option);
  void SetDoCent(Bool_t ce) { fDoCent = ce; }
  void SetDoEp(Bool_t ep) { fDoEp = ep; }

 protected:
  Bool_t fDoCent; //
  Bool_t fDoEp;   // 

 private:
  AliEmcalCompatTask(const AliEmcalCompatTask&);            // not implemented
  AliEmcalCompatTask &operator=(const AliEmcalCompatTask&); // not implemented

  ClassDef(AliEmcalCompatTask, 1); // Class to be able to run on skimmed ESDs
};

#endif
