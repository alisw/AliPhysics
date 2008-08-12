/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Artur Szostak <artursz@iafrica.com>                   *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTTrigger.h"

ClassImp(AliHLTTrigger)


AliHLTTrigger::AliHLTTrigger() :
	AliHLTProcessor(),
	fEventData(NULL),
	fTriggerData(NULL),
	fDecisionMade(false),
	fTriggerEventResult(0)
{
  // Default constructor sets pointers to NULL.
}


AliHLTTrigger::~AliHLTTrigger()
{
  // Default destructor.
}


int AliHLTTrigger::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData)
{
  /// Sets the pointers to the evtData and trigData, then calls the DoTrigger to
  /// execute the actual trigger algorithm.

  fEventData = &evtData;
  fTriggerData = &trigData;
  fDecisionMade = false;
  fTriggerEventResult = 0;
  
  int result = DoTrigger();
  if (result != 0) return result;
  
  // Fill in a default decision of false if none was made.
  if (not fDecisionMade)
  {
    TriggerEvent(false);
  }
  // Cleanup
  fEventData = NULL;
  fTriggerData = NULL;
  return fTriggerEventResult;
}


void AliHLTTrigger::TriggerEvent(bool value)
{
  /// Sets the trigger decision for the current event.
  
  if (fTriggerEventResult != 0) return;  // Do not do anything if a previous call failed.
  TObjString triggerResult(GetTriggerName());
  triggerResult.SetBit(BIT(14), value);  // Use bit 14 for the boolean decision.
  fTriggerEventResult = PushBack(triggerResult, kAliHLTDataTypeTObject|kAliHLTDataOriginOut);
}

