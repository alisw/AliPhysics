#ifndef ALIANALYSISTASKSTAT_H
#define ALIANALYSISTASKSTAT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Author: Andrei Gheata, 20/12/2010

//==============================================================================
//   AliAnalysisTaskStat - Analysis task that ataches an AliAnalysisStatistics
//      to the analysis manager
//==============================================================================

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliAnalysisStatistics;

class AliAnalysisTaskStat : public AliAnalysisTaskSE
{

protected:
  AliAnalysisStatistics      *fStatistics; // Statistics object
  TList                      *fOutputList; // Output list

private:
  AliAnalysisTaskStat(const AliAnalysisTaskStat& other);
  AliAnalysisTaskStat& operator= (const AliAnalysisTaskStat& other);

public:
  AliAnalysisTaskStat() : AliAnalysisTaskSE(), fStatistics(0), fOutputList(0) {}
  AliAnalysisTaskStat(const char *name);
  virtual ~AliAnalysisTaskStat();

  // Static method to add to the analysis manager
  static AliAnalysisTaskStat *AddToManager(UInt_t offlineMask=0);

  // Implementation of interface methods
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);
  
  // Getters
  AliAnalysisStatistics *GetStatistics() const {return fStatistics;}

  ClassDef(AliAnalysisTaskStat, 1); // Statistics task
};
#endif
