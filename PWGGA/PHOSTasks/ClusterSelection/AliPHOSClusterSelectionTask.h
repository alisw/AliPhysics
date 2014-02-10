#ifndef ALIPHOSCLUSTERSELECTIONTASK_CXX
#define ALIPHOSCLUSTERSELECTIONTASK_CXX

 /* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Analysis Task for selecting PHOS clusters for other analysis
// Makes some basic selections cuts, universal to all PHOS Analysis.
// Additional selection cuts can be applied via AliPHOSClusterSelection.
// 
// Authors : Henrik Qvigstad
// Date    : 16.01.2014
/* $Id$ */

class TRefArray;

class AliVCluster;

#include "AliAnalysisTaskSE.h"

class AliPHOSClusterSelectionTask : AliAnalysisTaskSE {
 public:
  AliPHOSClusterSelectionTask(const char* name = "AliPHOSClusterSelectionTask");
  virtual ~AliPHOSClusterSelectionTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  /* virtual void   Terminate(Option_t *); */

  TRefArray* GetPHOSClusters() const;
  TRefArray* GetPHOSClustersSelected(const AliPHOSClusterSelection* selection, bool useMap=true, bool addMap=true);

  static AliPHOSClusterSelectionTask* GetTask(const char* name = "AliPHOSClusterSelectionTask");

 protected:
  AliPHOSClusterSelectionTask(const AliPHOSClusterSelectionTask&); // not implemented
  AliPHOSClusterSelectionTask& operator=(const AliPHOSClusterSelectionTask&); // not implemented

  TRefArray* DeterminePHOSClustersSelected(const AliPHOSClusterSelection* selection);

  TRefArray* fClusters;
  TMap* fSelectionMap; // maps: ClusterSelection -> RefArray of Clusters

  // cluster cut variables:
  static const Double_t kMinClusterEnergy;
  static const Double_t kMinBCDistance;  //distance to nearest bad channel
  static const Int_t    kMinNCells;
  static const Double_t kMinM02;
  
  ClassDef(AliPHOSClusterSelectionTask, 1);
};

#endif
