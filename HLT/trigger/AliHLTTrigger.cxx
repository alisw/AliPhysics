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

/// @file   AliHLTTrigger.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   12 Aug 2008
/// @brief  Implementation of the AliHLTTrigger base component class.
///
/// The AliHLTTrigger class is the base class from which all HLT trigger components
/// should be derived.

#include "AliHLTTrigger.h"
#include "AliHLTTriggerDecision.h"

ClassImp(AliHLTTrigger)


AliHLTTrigger::AliHLTTrigger() :
	AliHLTProcessor(),
	fEventData(NULL),
	fTriggerData(NULL),
	fDecisionMade(false),
	fTriggerEventResult(0),
	fDescription(),
	fReadoutList(),
	fTriggerDomain()
{
  // Default constructor sets pointers to NULL.
}


AliHLTTrigger::~AliHLTTrigger()
{
  // Default destructor.
}


void AliHLTTrigger::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // Returns output data size estimate.
  // See header file for more details.

  constBase = sizeof(AliHLTTriggerDecision);
  inputMultiplier = 1;
}


int AliHLTTrigger::DoEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData)
{
  // Sets the pointers to the evtData and trigData, then calls the DoTrigger to
  // execute the actual trigger algorithm.

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
  // Sets the trigger decision for the current event.
  // See header file for more details.
  
  if (fTriggerEventResult != 0) return;  // Do not do anything if a previous call failed.
  AliHLTTriggerDecision triggerResult(value, GetTriggerName(), fReadoutList, fTriggerDomain, fDescription);
  fTriggerEventResult = PushBack(&triggerResult, kAliHLTDataTypeTObject|kAliHLTDataOriginOut);
  if (fTriggerEventResult == 0) fDecisionMade = true;
}


void AliHLTTrigger::TriggerEvent(
    AliHLTTriggerDecision* result, const AliHLTComponentDataType& type,
    AliHLTUInt32_t spec
  )
{
  // Sets a custom trigger decision for the current event.
  // See header file for more details.
  
  if (fTriggerEventResult != 0) return;  // Do not do anything if a previous call failed.
  fTriggerEventResult = PushBack(result, type, spec);
  if (fTriggerEventResult == 0) fDecisionMade = true;
}


void AliHLTTrigger::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // Calls the const version of this method.
  
  // Assign to const temporary variable to make sure we call the constant version
  // of the GetOutputDataTypes method.
  const AliHLTTrigger* t = this;
  t->GetOutputDataTypes(list);
}


int AliHLTTrigger::GetOutputDataTypes(AliHLTComponentDataTypeList& list)
{
  // Calls the const version of this method.
  
  // Assign to const temporary variable to make sure we call the constant version
  // of the GetOutputDataTypes method.
  const AliHLTTrigger* t = this;
  t->GetOutputDataTypes(list);
  list.push_back(kAliHLTDataTypeTObject|kAliHLTDataOriginOut);
  return list.size();
}

