#include "AliEMCalTriggerAnaTriggerDecisionConfig.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerAnaTriggerDecisionConfig)

namespace EMCalTriggerPtAnalysis {

/**
 * Default constructor
 */
AliEMCalTriggerAnaTriggerDecisionConfig::AliEMCalTriggerAnaTriggerDecisionConfig():
	fSwapThresholds(kFALSE),
	fUseOfflinePatches(kFALSE),
	fEnergyType(kAmplitudeOnline)
{
	memset(fEnergyThresholds, 0, sizeof(double) * 4);
}

}
