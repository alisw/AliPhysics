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
/// The AliHLTGlobalTrigger base class is an abstract class from which a
/// derived class is constructed from AliHLTTriggerMenu on the fly. The derived
/// class then implements triggering based on the particular trigger menu.

#include "AliHLTGlobalTrigger.h"
#include "AliHLTGlobalTriggerWrapper.h"
#include "TClass.h"

ClassImp(AliHLTGlobalTrigger)


AliHLTGlobalTrigger* AliHLTGlobalTrigger::CreateNew(const char* name)
{
  // Creates a new instance of the named global trigger class.
  
  TClass* c = TClass::GetClass(name);
  if (c == NULL) return NULL;
  if (c->GetDeclFileLine() == -1 and c->GetImplFileLine() == -1)
  {
    // Could not find the implementation lines which should be there if the code
    // was compiled. So assuming that this is an interpreted class. In this case
    // we need to use a interface wrapper class to make things work properly.
    AliHLTGlobalTriggerWrapper* trigger = new AliHLTGlobalTriggerWrapper(name);
    if (not trigger->IsValid())
    {
      delete trigger;
      return NULL;
    }
    return trigger;
  }
  else
  {
    return static_cast<AliHLTGlobalTrigger*>(c->New());
  }
}

