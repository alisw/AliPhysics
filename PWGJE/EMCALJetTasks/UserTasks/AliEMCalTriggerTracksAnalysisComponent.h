#ifndef ALIEMCALTRIGGERTRACKSANALYSISCOMPONENT_H
#define ALIEMCALTRIGGERTRACKSANALYSISCOMPONENT_H
/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#include <TObject.h>
#include "AliEMCalHistoContainer.h"

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerEventData;

class AliEMCalTriggerTracksAnalysisComponent : public TObject {
public:
  AliEMCalTriggerTracksAnalysisComponent();
  AliEMCalTriggerTracksAnalysisComponent(const char *name);
  virtual ~AliEMCalTriggerTracksAnalysisComponent();

  virtual void CreateHistos();
  virtual void Process(const AliEMCalTriggerEventData * const data) = 0;

  THashList *GetHistList() const { return fHistos->GetListOfHistograms(); }

protected:
  AliEMCalHistoContainer *fHistos;

  ClassDef(AliEMCalTriggerTracksAnalysisComponent, 1)
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERTRACKSANALYSISCOMPONENT_H */
