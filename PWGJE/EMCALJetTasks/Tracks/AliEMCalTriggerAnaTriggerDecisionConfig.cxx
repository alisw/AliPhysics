#include "AliEMCalTriggerAnaTriggerDecisionConfig.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerAnaTriggerDecisionConfig)

using namespace PWGJE::EMCALJetTasks;

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
