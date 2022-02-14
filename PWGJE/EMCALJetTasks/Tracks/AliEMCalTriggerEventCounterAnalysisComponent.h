/**
 * \file AliEMCalTriggerEventCounterAnalysisComponent.h
 * \brief
 *
 * \author Markus Fasel <markus.fasel@cern.ch>, Lawrence Berkeley National Laboratory
 * \date Dec. 12, 2014
 */
#ifndef ALIEMCALTRIGGEREVENTCOUNTERANALYSISCOMPONENT_H
#define ALIEMCALTRIGGEREVENTCOUNTERANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


#include "AliEMCalTriggerTracksAnalysisComponent.h"

namespace PWGJE {
  
namespace EMCALJetTasks {

class AliEMCalTriggerEventData;

/**
 * \class AliEMCalTriggerEventCounterAnalysisComponent
 * \brief Event counter analysis component for the trigger analysis
 *
 * Analysis component counting events for different trigger classes. Task needs
 * to be grouped with a global event selection.
 */
class AliEMCalTriggerEventCounterAnalysisComponent: public AliEMCalTriggerTracksAnalysisComponent {
public:
  AliEMCalTriggerEventCounterAnalysisComponent();
  AliEMCalTriggerEventCounterAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerEventCounterAnalysisComponent() {}

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data);

protected:
  void DefineAxis(TAxis& axis, const char* name,
      const char* title, int nbins, double min, double max,
      const char** labels) const;
  Int_t FindAxis(THnSparse *hist, const char *title) const;

  ClassDef(AliEMCalTriggerEventCounterAnalysisComponent, 1);    // Analysis component for event counting
};

} /* namespace EMCALJetTasks */

} /* namespace PWGJE */

#endif /* ALIEMCALTRIGGEREVENTCOUNTERANALYSISCOMPONENT_H */
