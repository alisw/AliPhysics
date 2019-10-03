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
  AliAnalysisTaskDummy();
  virtual ~AliAnalysisTaskDummy();

  virtual void UserCreateOutputObjects() {}
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *) {}

  void SetDebugLevel(Int_t level) { fDebugLevel = level; }

private:
  Int_t         fDebugLevel;          ///< Debug level

  ClassDef(AliAnalysisTaskDummy, 1);
};

#endif /* ALIANALYSISTASKDUMMY_H */
