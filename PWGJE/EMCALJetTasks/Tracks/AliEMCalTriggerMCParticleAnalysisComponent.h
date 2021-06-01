#ifndef ALIEMCALTRIGGERMCPARTICLEANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERMCPARTICLEANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include "AliEMCalTriggerTracksAnalysisComponent.h"

class TAxis;

namespace PWGJE {
  
namespace EMCALJetTasks {

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

  bool IsPhysicalPrimary(const AliVParticle *const part, const AliMCEvent *const ev) const;

  ClassDef(AliEMCalTriggerMCParticleAnalysisComponent, 1);
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGERMCPARTICLEANALYSISCOMPONENT_H */
