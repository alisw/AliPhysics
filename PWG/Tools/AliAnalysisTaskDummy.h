#ifndef ALIANALYSISTASKDUMMY_H
#define ALIANALYSISTASKDUMMY_H
/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskSE.h"

/**
 * \class AliAnalysisTaskDummy
 * \brief Dummy analysis task, only trying to read data - no memory allocation. For file corruption checks.
 */
class AliAnalysisTaskDummy : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskDummy() {}
  AliAnalysisTaskDummy(const char *name, Int_t wait_time) : AliAnalysisTaskSE(name) { fWaitTime = wait_time; }
  virtual ~AliAnalysisTaskDummy() {}

  virtual void UserCreateOutputObjects() {}
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {}

  void SetDebugLevel(Int_t level) { fDebugLevel = level; }

private:
  Int_t         fDebugLevel = 0;       ///< Debug level
  Int_t         fWaitTime = 0;         ///< Wait time in milliseconds

  ClassDef(AliAnalysisTaskDummy, 2);
};

#endif /* ALIANALYSISTASKDUMMY_H */
