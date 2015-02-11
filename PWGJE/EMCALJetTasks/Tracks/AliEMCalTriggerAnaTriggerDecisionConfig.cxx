/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/*
 * Configuration handler for the trigger decision builder
 *
 *   Author: Markus Fasel
 */

#include "AliEMCalTriggerAnaTriggerDecisionConfig.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerAnaTriggerDecisionConfig)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerAnaTriggerDecisionConfig::AliEMCalTriggerAnaTriggerDecisionConfig():
	fSwapThresholds(kFALSE),
	fUseOfflinePatches(kFALSE)
{
	/*
	 * Default constructor
	 */
	memset(fEnergyThresholds, 0, sizeof(double) * 4);
}

//______________________________________________________________________________
void AliEMCalTriggerAnaTriggerDecisionConfig::ConfigureTriggerDecision(AliEMCalTriggerAnaTriggerDecision& dec) const {
	/*
	 * Configure trigger decision builder
	 *
	 * @param dec: the trigger decision builder
	 */
	dec.SetSwapThresholds(fSwapThresholds);
	dec.SetUseOfflinePatches(fUseOfflinePatches);
	for(int i = 0; i < 4; i++){
		dec.SetOfflineEnergyThreshold(AliEMCalTriggerAnaTriggerDecision::ETATriggerType(i), fEnergyThresholds[i]);
	}
}

} /* namespace EMCalTriggerPtAnalysis */
