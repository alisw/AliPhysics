#ifndef ALIPHOSCLUSTERSELECTIONTASK_CXX
#define ALIPHOSCLUSTERSELECTIONTASK_CXX

 /* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis Task for selecting PHOS clusters for other analysis
// Authors : Henrik Qvigstad
// Date    : 16.01.2014
/* $Id$ */

class AliVCluster;

#include "TArrayI.h"

#include "AliAnalysisTaskSE.h"

class AliPHOSClusterSelectionTask : AliAnalysisTaskSE {
 public:
  AliPHOSClusterSelectionTask(const char* name = "AliPHOSClusterSelectionTask");
  virtual ~AliPHOSClusterSelectionTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  /* virtual void   Terminate(Option_t *); */

  int GetNClusters() const; // Number of clusters in the current event
  AliVCluster* GetCluster(int index) const; // get cluster of index
  TArrayI GetIndexOfSelected(const AliPHOSClusterSelection* selection) const;

 protected:
  AliPHOSClusterSelectionTask(const AliPHOSClusterSelectionTask&); // not implemented
  AliPHOSClusterSelectionTask& operator=(const AliPHOSClusterSelectionTask&); // not implemented

  
  ClassDef(AliPHOSClusterSelectionTask, 1);
};

#endif
