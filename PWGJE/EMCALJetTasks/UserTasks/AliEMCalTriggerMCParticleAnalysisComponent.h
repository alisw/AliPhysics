#ifndef ALIEMCALTRIGGERMCPARTICLEANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERMCPARTICLEANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliEMCalTriggerTracksAnalysisComponent.h"

class TAxis;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerBinningDimension;
class AliEMCalTriggerEventData;

class AliEMCalTriggerMCParticleAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerMCParticleAnalysisComponent();
  AliEMCalTriggerMCParticleAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerMCParticleAnalysisComponent() {}

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

protected:
  TAxis *DefineAxis(const char *name, const AliEMCalTriggerBinningDimension *binning);

  ClassDef(AliEMCalTriggerMCParticleAnalysisComponent, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERMCPARTICLEANALYSISCOMPONENT_H */
