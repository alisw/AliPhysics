#ifndef ALIEMCALTRIGGERCLUSTERANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERCLUSTERANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliCutValueRange.h"
#include "AliEMCalTriggerTracksAnalysisComponent.h"

class TString;
class AliVCluster;
class AliVEvent;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData;

class AliEMCalTriggerClusterAnalysisComponent : public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerClusterAnalysisComponent();
  AliEMCalTriggerClusterAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerClusterAnalysisComponent() {}

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  void SetUsePatches(Bool_t usePatches = kTRUE) { fUsePatches = usePatches; }
  void SetEnergyRange(double min, double max) { fEnergyRange.SetLimits(min, max); }

protected:
  void FillHistogram(const TString &histname, const AliVCluster *clust, AliVEvent *ev, Bool_t inMB);

  AliCutValueRange<double>    fEnergyRange;
  Bool_t                      fUsePatches;

  ClassDef(AliEMCalTriggerClusterAnalysisComponent, 1);       // Analysis component for EMCal cluster
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERCLUSTERANALYSISCOMPONENT_H */
