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
#include <cstring>

ClassImp(AliHLTGlobalTrigger)

// Static factory array.
AliHLTGlobalTrigger::Factory*
AliHLTGlobalTrigger::Factory::fFactory[AliHLTGlobalTrigger::Factory::kMaxFactories]
  = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL};


AliHLTGlobalTrigger::AliHLTGlobalTrigger() :
  AliHLTLogging(),
  fCounters()
{
  // Default constructor.
}


AliHLTGlobalTrigger::~AliHLTGlobalTrigger()
{
  // Default destructor.
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
  
  fCounters.Set(number);
  for (UInt_t i = 0; i < number; i++)
  {
    fCounters[i] = 0;
  }
}

