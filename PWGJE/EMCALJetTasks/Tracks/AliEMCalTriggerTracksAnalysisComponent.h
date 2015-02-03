#ifndef ALIEMCALTRIGGERTRACKSANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERTRACKSANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <vector>
#include <string>
#include <TNamed.h>
#include "AliEMCalHistoContainer.h"

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerAnaTriggerDecision;
class AliEMCalTriggerBinningComponent;
class AliEMCalTriggerBinningDimension;
class AliEMCalTriggerEventData;
class AliEMCalTriggerKineCuts;

class AliEMCalTriggerTracksAnalysisComponent : public TNamed {
public:
  AliEMCalTriggerTracksAnalysisComponent();
  AliEMCalTriggerTracksAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerTracksAnalysisComponent();

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data) = 0;

  THashList *GetHistList() const { return fHistos->GetListOfHistograms(); }
  void SetBinning(const AliEMCalTriggerBinningComponent * const binning) { fBinning = binning; }
  void SetKineCuts(const AliEMCalTriggerKineCuts * const cuts) { fKineCuts = cuts; }
  void SetTriggerDecision(const AliEMCalTriggerAnaTriggerDecision *trigger) { fTriggerDecision = trigger; }
  void SetComponentDebugLevel(int debuglevel) { fComponentDebugLevel = debuglevel; }

protected:
  TAxis *DefineAxis(const char *name, const AliEMCalTriggerBinningDimension *binning);
  TAxis *DefineAxis(const char *name, int nbins, double min, double max);
  void GetMachingTriggerNames(std::vector<std::string> &triggernames, Bool_t usePatches);
  void PrintTriggerNames(const std::vector<std::string> &, const std::string &componentName) const;

  AliEMCalHistoContainer                      *fHistos;
  const AliEMCalTriggerBinningComponent       *fBinning;
  const AliEMCalTriggerKineCuts               *fKineCuts;
  const AliEMCalTriggerAnaTriggerDecision     *fTriggerDecision;

  Int_t                                         fComponentDebugLevel;

  ClassDef(AliEMCalTriggerTracksAnalysisComponent, 1)
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERTRACKSANALYSISCOMPONENT_H */
