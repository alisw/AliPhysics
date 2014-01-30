#ifndef ALIEMCALDEBUGTASK_H
#define ALIEMCALDEBUGTASK_H

// $Id$

#include "AliAnalysisTaskSE.h"

class AliEmcalDebugTask : public AliAnalysisTaskSE {
 public:
  AliEmcalDebugTask();
  AliEmcalDebugTask(const char *name);
  virtual ~AliEmcalDebugTask();

  void        SetId(UInt_t id)           { fId       = id; }
  void        SetFileTest(const char *n) { fFileTest =  n; }
  void        SetPrintEnv(Bool_t b)      { fPrintEnv = b;  }

 protected:
  void        UserCreateOutputObjects();
  void        UserExec(Option_t *option);

  UInt_t      fId;         //id to be stored in the output file
  TString     fFileTest;   //path name test 
  Bool_t      fPrintEnv;   //print env if true
  TList      *fOutput;     //!output list
  TString     fFileName;   //!current file name
  UInt_t      fRand;       //!random number

 private:
  AliEmcalDebugTask(const AliEmcalDebugTask&);            // not implemented
  AliEmcalDebugTask &operator=(const AliEmcalDebugTask&); // not implemented

  ClassDef(AliEmcalDebugTask, 1); // Class to be able to run on skimmed ESDs
};

#endif
