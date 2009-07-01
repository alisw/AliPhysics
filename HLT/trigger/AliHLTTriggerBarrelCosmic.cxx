// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/// @file   AliHLTTriggerBarrelCosmic.cxx
/// @author Matthias Richter
/// @date   2009-06-30
/// @brief  HLT cosmics trigger component for the central barrel region.

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTTriggerBarrelCosmic.h"
#include "AliESDEvent.h"
#include "AliHLTTriggerDecision.h"
#include "AliHLTDomainEntry.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTriggerBarrelCosmic)

AliHLTTriggerBarrelCosmic::AliHLTTriggerBarrelCosmic()
  : AliHLTTrigger()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTriggerBarrelCosmic::~AliHLTTriggerBarrelCosmic()
{
  // see header file for class documentation
}

const char* AliHLTTriggerBarrelCosmic::GetTriggerName() const
{
  // see header file for class documentation
  return "BarrelCosmicsTrigger";
}

AliHLTComponent* AliHLTTriggerBarrelCosmic::Spawn()
{
  // see header file for class documentation
  return new AliHLTTriggerBarrelCosmic;
}

int AliHLTTriggerBarrelCosmic::DoTrigger()
{
  // see header file for class documentation
  SetDescription("central barrel cosmics trigger needs to be implemented");
  TriggerEvent(false);
  return 0;
}
