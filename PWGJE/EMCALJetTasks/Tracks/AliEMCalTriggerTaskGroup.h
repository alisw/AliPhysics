#ifndef ALIEMCALTRIGGERTASKGROUP_H
#define ALIEMCALTRIGGERTASKGROUP_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <TNamed.h>
#include <TObjArray.h>

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerAnaTriggerDecision;
class AliEMCalTriggerBinningComponent;
class AliEMCalTriggerEventSelection;
class AliEMCalTriggerKineCuts;
class AliEMCalTriggerTracksAnalysisComponent;

class AliEMCalTriggerTaskGroup : public TNamed {
public:
  AliEMCalTriggerTaskGroup();
  AliEMCalTriggerTaskGroup(const char *name);
  virtual ~AliEMCalTriggerTaskGroup();

  void SetEventSelection(const AliEMCalTriggerEventSelection *sel){ fEventSelection = sel; }
  void SetGlobalBinning(const AliEMCalTriggerBinningComponent *const binning) { fBinning = binning; }
  void SetTriggerDecision(const AliEMCalTriggerAnaTriggerDecision *trigger);
  void SetKineCuts(const AliEMCalTriggerKineCuts *cuts) { fKineCuts = cuts; }
  void AddAnalysisComponent(AliEMCalTriggerTracksAnalysisComponent * const analysis);

  TList * InitialiseAnalysisComponents();
  void Process(const AliEMCalTriggerEventData * const event);

protected:
  TObjArray                                   *fAnalysisComponents;
  const AliEMCalTriggerEventSelection         *fEventSelection;
  const AliEMCalTriggerBinningComponent       *fBinning;
  const AliEMCalTriggerKineCuts               *fKineCuts;

  ClassDef(AliEMCalTriggerTaskGroup, 1);    // Group of analysis components with common event selection
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERTASKGROUP_H */
