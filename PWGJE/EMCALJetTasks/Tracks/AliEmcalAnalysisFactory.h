/*
 * AliEmcalAnalysisFactory.h
 *
 *  Created on: Feb 23, 2016
 *      Author: markus
 */

#ifndef ALIEMCALANALYSISFACTORY_H
#define ALIEMCALANALYSISFACTORY_H

#include <TObject.h>
#include <TString.h>
#include "AliEmcalTriggerOfflineSelection.h"

class AliEmcalTrackSelection;

namespace EMCalTriggerPtAnalysis {

class AliEmcalAnalysisFactory : public TObject {
public:
  AliEmcalAnalysisFactory() {}
  virtual ~AliEmcalAnalysisFactory(){}

  static TString ClusterContainerNameFactory(Bool_t isAOD);
  static AliEmcalTrackSelection *TrackCutsFactory(TString name, Bool_t isAOD);
  static AliEmcalTriggerOfflineSelection *TriggerSelectionFactory(Double_t el0, Double_t eg1, Double_t eg2, Double_t ej1, Double_t ej2, AliEmcalTriggerOfflineSelection::EmcalEnergyDefinition_t endef = AliEmcalTriggerOfflineSelection::kFEEEnergy);

  ClassDef(AliEmcalAnalysisFactory, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALANALYSISFACTORY_H */
