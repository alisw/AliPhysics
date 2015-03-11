#ifndef ALIEMCALTRIGGERWEIGHTHANDLER_H
#define ALIEMCALTRIGGERWEIGHTHANDLER_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel

#include <TObject.h>

class TF1;

class AliMCEvent;

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerWeightHandler : public TObject {
public:
  AliEMCalTriggerWeightHandler();
  virtual ~AliEMCalTriggerWeightHandler() {}

  void SetUsePtHard(bool usePtHard) { fUsePtHard = usePtHard; }
  void SetWeightModel(TF1 *model) { fWeightModel = model; }
  double GetEventWeight(const AliMCEvent *const event) const;

private:
  TF1               *fWeightModel;    /// Weight model
  bool               fUsePtHard;      /// Calculate weight using pt-hard

  ClassDef(AliEMCalTriggerWeightHandler, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERWEIGHTHANDLER_H */
