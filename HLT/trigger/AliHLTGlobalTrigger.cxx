// $Id$
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

/// @file   AliHLTGlobalTrigger.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Dec 2008
/// @brief  Implementation of the AliHLTGlobalTrigger base class.
///
/// The AliHLTGlobalTriggerComponent class is an abstract class from which a
/// derived class is constructed by AliHLTTriggerMenu on the fly. The derived
/// class then implements triggering based on the particular trigger menu.

#include "AliHLTGlobalTrigger.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "AliHLTCTPData.h"
#include "TArrayL64.h"
#include "TClonesArray.h"
#include <cstring>
#include <cassert>

ClassImp(AliHLTGlobalTrigger)

// Static factory array.
AliHLTGlobalTrigger::Factory*
AliHLTGlobalTrigger::Factory::fFactory[AliHLTGlobalTrigger::Factory::kMaxFactories]
  = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};


AliHLTGlobalTrigger::AliHLTGlobalTrigger() :
  AliHLTLogging()
  , fCTPDecisions(NULL)
  , fCounters(NULL)
{
  // Default constructor.
}


AliHLTGlobalTrigger::~AliHLTGlobalTrigger()
{
  // Default destructor.
  if (fCounters) {
    delete fCounters;
  }
  if (fCTPDecisions) {
    fCTPDecisions->Delete();
    delete fCTPDecisions;
  }
}


AliHLTGlobalTrigger* AliHLTGlobalTrigger::Factory::CreateNew(const char* name)
{
  // Creates a new instance of the named trigger class.
  
  for (int i = 0; i < kMaxFactories; i++)
  {
    if (fFactory[i] != NULL)
    {
      if (strcmp(fFactory[i]->ClassName(), name) == 0)
      {
        return fFactory[i]->New();
      }
    }
  }
  return NULL;
}


AliHLTGlobalTrigger::Factory::Factory() : AliHLTLogging()
{
  // Default constructor resisters the class factory.
  
  for (int i = 0; i < kMaxFactories; i++)
  {
    if (fFactory[i] == NULL)
    {
      fFactory[i] = this;
      return;
    }
  }
  
  HLTFatal("Trying to register too many global trigger factories.");
}


AliHLTGlobalTrigger::Factory::~Factory()
{
  // The default destructor deregisters the factory.
  
  for (int i = 0; i < kMaxFactories; i++)
  {
    if (fFactory[i] == this)
    {
      fFactory[i] = NULL;
      return;
    }
  }
  
  HLTFatal("Could not find factory to deregister.");
}


void AliHLTGlobalTrigger::ResetCounters(UInt_t number)
{
  // Resets the trigger counters.
  
  if (!fCounters) fCounters = new TArrayL64(number);
  if (!fCounters) return;

  fCounters->Set(number);
  for (UInt_t i = 0; i < number; i++)
  {
    (*fCounters)[i] = 0;
  }
}

void AliHLTGlobalTrigger::IncrementCounter(UInt_t i) 
{
  // increment a specific counter
  if (fCounters && i<(unsigned)fCounters->GetSize()) ++(*fCounters)[i]; 
}

Long64_t AliHLTGlobalTrigger::GetCounter(UInt_t i) const
{
  if (fCounters && i<(unsigned)fCounters->GetSize()) return (*fCounters)[i];
  return 0;
}

int AliHLTGlobalTrigger::AddCTPDecisions(const AliHLTCTPData* pCTPData, const AliHLTComponentTriggerData* trigData)
{
  // add trigger decisions for the valid CTP classes
  if (!pCTPData) return 0;

  AliHLTUInt64_t triggerMask=pCTPData->Mask();
  AliHLTUInt64_t bit0=0x1;
  if (!fCTPDecisions) {
    fCTPDecisions=new TClonesArray(AliHLTTriggerDecision::Class(), gkNCTPTriggerClasses);
    if (!fCTPDecisions) return -ENOMEM;

    fCTPDecisions->ExpandCreate(gkNCTPTriggerClasses);
    for (int i=0; i<gkNCTPTriggerClasses; i++) {
      const char* name=pCTPData->Name(i);
      if (triggerMask&(bit0<<i) && name) {
	AliHLTTriggerDecision* pDecision=dynamic_cast<AliHLTTriggerDecision*>(fCTPDecisions->At(i));
	assert(pDecision);
	if (!pDecision) return -ENOENT;
	pDecision->Name(name);
      }
    }
  }

  for (int i=0; i<gkNCTPTriggerClasses; i++) {
    const char* name=pCTPData->Name(i);
    if ((triggerMask&(bit0<<i))==0 || name==NULL) continue;
    AliHLTTriggerDecision* pDecision=dynamic_cast<AliHLTTriggerDecision*>(fCTPDecisions->At(i));
    HLTDebug("updating CTP trigger decision %d %s (%p casted %p)", i, name, fCTPDecisions->At(i), pDecision);
    if (!pDecision) return -ENOENT;

    bool result=false;
    if (trigData) result=pCTPData->EvaluateCTPTriggerClass(name, *trigData);
    else result=pCTPData->EvaluateCTPTriggerClass(name);
    pDecision->Result(result);
    pDecision->TriggerDomain().Clear();
    pDecision->TriggerDomain().Add(pCTPData->ReadoutList(*trigData));

    Add(fCTPDecisions->At(i), kAliHLTDataTypeTriggerDecision, kAliHLTVoidDataSpec);
  }

  return 0;
}

