#ifndef ALIEMCALTRIGGERTASKGROUP_H
#define ALIEMCALTRIGGERTASKGROUP_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <TNamed.h>
#include <TObjArray.h>

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventSelection;
class AliEMCalTriggerTracksAnalysisComponent;

class AliEMCalTriggerTaskGroup : public TNamed {
public:
  AliEMCalTriggerTaskGroup();
  AliEMCalTriggerTaskGroup(const char *name);
  virtual ~AliEMCalTriggerTaskGroup();

  void SetEventSelection(const AliEMCalTriggerEventSelection *sel){ fEventSelection = sel; }
  void AddAnalysisComponent(AliEMCalTriggerTracksAnalysisComponent * const analysis);

  TList * InitialiseAnalysisComponents();
  void Process(const AliEMCalTriggerEventData * const event);

protected:
  TObjArray                                   *fAnalysisComponents;
  const AliEMCalTriggerEventSelection         *fEventSelection;

  ClassDef(AliEMCalTriggerTaskGroup, 1);    // Group of analysis components with common event selection
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERTASKGROUP_H */
