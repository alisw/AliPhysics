#ifndef ALIPHOSCLUSTERSELECTIONLOGBACKTASK_CXX
#define ALIPHOSCLUSTERSELECTIONLOGBACKTASK_CXX

 /* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis Task for selecting PHOS clusters for other analysis, with Logback.
// Allows to access clusters of past events. The clusters of said events are 
// stored in per-event arrays, and the event-arrays are stored in one 
// first-in-last-out list per event selection scope. Said scope is defined by
// AliPHOSEventSelection or (recommended) whatever subclass used in its place.
//
// Authors : Henrik Qvigstad
// Date    : 16.01.2014
/* $Id$ */

class AliVCluster;

#include "TArrayI.h"

#include "AliAnalysisTaskSE.h"

class AliPHOSClusterSelectionLogbackTask : AliPHOSClusterSelectionTask {
 public:
  AliPHOSClusterSelectionLogbackTask(const char* name = "AliPHOSClusterSelectionLogbackTask");
  virtual ~AliPHOSClusterSelectionLogbackTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  /* virtual void   Terminate(Option_t *); */

  int GetNClustersBacklog(const AliPHOSEventSelection* eventSelection, int eventBacklogIndex) const; // Number of clusters in the current event
  AliVCluster* GetClusterBacklog(const int& index, const AliPHOSEventSelection* eventSelection, int eventBacklogIndex) const; // get cluster of index
  TArrayI GetIndexOfSelectedBacklog(const AliPHOSClusterSelection* selection, const AliPHOSEventSelection& eventSelection, int eventBacklogIndex) const;
  
  void AddEventSelectionLogback(const AliPHOSEventSelection* selection);

 protected:
  AliPHOSClusterSelectionLogbackTask(const AliPHOSClusterSelectionLogbackTask&); // not implemented
  AliPHOSClusterSelectionLogbackTask& operator=(const AliPHOSClusterSelectionLogbackTask&); // not implemented

  
  ClassDef(AliPHOSClusterSelectionLogbackTask, 1);
};

#endif
