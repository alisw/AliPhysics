#ifndef ALIEMCALPHYSICSSELECTIONTASK_H
#define ALIEMCALPHYSICSSELECTIONTASK_H

// $Id$

#include "AliPhysicsSelectionTask.h"

class AliPhysicsSelection;
class TH1;

class AliEmcalPhysicsSelectionTask : public AliPhysicsSelectionTask {
 public:
  AliEmcalPhysicsSelectionTask();
  AliEmcalPhysicsSelectionTask(const char* opt);
  virtual ~AliEmcalPhysicsSelectionTask() {};

  virtual void   UserExec(const Option_t *opt);
  virtual void   UserCreateOutputObjects();
  virtual void   Terminate(Option_t*);

  void           SetDoWriteHistos(Bool_t b) { fDoWriteHistos = b; }
  Int_t          GetNCalled() const         { return fNCalled;    }
  Int_t          GetNAccepted() const       { return fNAccepted;  }

 protected:
  Bool_t         fDoWriteHistos; //=true then write output
  Int_t          fNCalled;       //!how often was the PS called
  Int_t          fNAccepted;     //!how often was the event accepted
  TH1           *fHAcc;          //!acceptance histo

 private:
  AliEmcalPhysicsSelectionTask(const AliEmcalPhysicsSelectionTask&);
  AliEmcalPhysicsSelectionTask& operator=(const AliEmcalPhysicsSelectionTask&);

  ClassDef(AliEmcalPhysicsSelectionTask, 1); // Emcal physics selection task
};
#endif
