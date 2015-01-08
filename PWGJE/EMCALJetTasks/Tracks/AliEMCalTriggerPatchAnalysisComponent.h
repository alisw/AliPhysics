#ifndef ALIEMCALTRIGGERPATCHANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERPATCHANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliEMCalTriggerTracksAnalysisComponent.h"

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData;

class AliEMCalTriggerPatchAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerPatchAnalysisComponent();
  AliEMCalTriggerPatchAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerPatchAnalysisComponent() { }

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

protected:

  ClassDef(AliEMCalTriggerPatchAnalysisComponent, 1);     // Component for trigger patch analysis
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERPATCHANALYSISCOMPONENT_H */
