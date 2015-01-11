#ifndef ALIEMCALTRIGGEREVENTCOUNTERANALYSISCOMPONENT_H
#define ALIEMCALTRIGGEREVENTCOUNTERANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliEMCalTriggerTracksAnalysisComponent.h"

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData;

class AliEMCalTriggerEventCounterAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerEventCounterAnalysisComponent();
  AliEMCalTriggerEventCounterAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerEventCounterAnalysisComponent() {}

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

  void SetUsePatches(Bool_t doUse = kTRUE) { fUsePatches = doUse; }

protected:
  void DefineAxis(TAxis& axis, const char* name,
      const char* title, int nbins, double min, double max,
      const char** labels) const;

  Bool_t          fUsePatches;                                  // Use patches for trigger decision

  ClassDef(AliEMCalTriggerEventCounterAnalysisComponent, 1);    // Analysis component for event counting
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGEREVENTCOUNTERANALYSISCOMPONENT_H */
