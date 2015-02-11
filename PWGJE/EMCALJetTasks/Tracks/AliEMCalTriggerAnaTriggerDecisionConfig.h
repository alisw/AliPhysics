/* Copyright(c) 1998-2015, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

// Author: Markus Fasel
#ifndef ALIEMCALTRIGGERANATRIGGERDECISIONCONFIG_H
#define ALIEMCALTRIGGERANATRIGGERDECISIONCONFIG_H

#include "AliEMCalTriggerAnaTriggerDecision.h"

namespace EMCalTriggerPtAnalysis {

class AliEMCalTriggerAnaTriggerDecisionConfig : public TObject {
public:
	AliEMCalTriggerAnaTriggerDecisionConfig();
	virtual ~AliEMCalTriggerAnaTriggerDecisionConfig() {}

	void SetSwapThresholds(Bool_t doSwap = kTRUE) { fSwapThresholds = doSwap;}
	void SetUseOfflinePatches(Bool_t doUse = kTRUE ) { fUseOfflinePatches = kTRUE; }
	void SetEnergyThreshold(AliEMCalTriggerAnaTriggerDecision::ETATriggerType trigger, double threshold){
		fEnergyThresholds[trigger] = threshold;
	}

	Bool_t IsSwapThresholds() const { return fSwapThresholds; }
	Bool_t IsUsingOfflinePatches() const { return fUseOfflinePatches; }
	double GetEnergyThreshold(AliEMCalTriggerAnaTriggerDecision::ETATriggerType trigger) const {
		return fEnergyThresholds[trigger];
	}

	void ConfigureTriggerDecision(AliEMCalTriggerAnaTriggerDecision &dec) const;

private:
	Bool_t 				fSwapThresholds;
	Bool_t				fUseOfflinePatches;
	Double_t			fEnergyThresholds[4];

	ClassDef(AliEMCalTriggerAnaTriggerDecisionConfig, 1);
};

} /* namespace EMCalTriggerPtAnalysis */

#endif /* ALIEMCALTRIGGERANATRIGGERDECISIONCONFIG_H */
