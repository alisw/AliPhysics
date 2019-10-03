#ifndef ALIPHOSCLUSTERSELECTIONLOGBACKTASK_CXX
#define ALIPHOSCLUSTERSELECTIONLOGBACKTASK_CXX

class AliPHOSEventSelection;
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
class TObjArray;
class AliPHOSEventSelection;

#include "TArrayI.h"

#include "AliPHOSClusterSelectionTask.h"

class AliPHOSClusterSelectionLogbackTask : public AliPHOSClusterSelectionTask {
 public:
  AliPHOSClusterSelectionLogbackTask(const char* name = "AliPHOSClusterSelectionLogbackTask");
  virtual ~AliPHOSClusterSelectionLogbackTask();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  /* virtual void   Terminate(Option_t *); */

  TObjArray* GetPHOSClustersLogback(const AliPHOSEventSelection* eventSelection, UInt_t eventLogbackIndex, const AliPHOSClusterSelection* clusterSelection=0) const;

  void LogEvent(const AliPHOSEventSelection* selection, UInt_t nEventsToLog);

  static AliPHOSClusterSelectionLogbackTask* GetTask(const char* name = "AliPHOSClusterSelectionLogbackTask");

 protected:
  AliPHOSClusterSelectionLogbackTask(const AliPHOSClusterSelectionLogbackTask&); // not implemented
  AliPHOSClusterSelectionLogbackTask& operator=(const AliPHOSClusterSelectionLogbackTask&); // not implemented

  TMap* fMapOfEventLists; // fMapOfEventLists: EventSelection -> EventList, EventList of EventArray in (CluArray, SelMap), SelMap: CluSelection -> CluArray.
  enum EventArrayIndex { kCluArray, kSelMap, kSize };

  ClassDef(AliPHOSClusterSelectionLogbackTask, 1);
};

#endif
