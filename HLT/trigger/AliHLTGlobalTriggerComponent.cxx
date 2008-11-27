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

/// @file   AliHLTGlobalTriggerComponent.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   26 Nov 2008
/// @brief  Implementation of the AliHLTGlobalTriggerComponent component class.
///
/// The AliHLTGlobalTriggerComponentComponent class applies the global HLT trigger to all
/// trigger information produced by components deriving from AliHLTTrigger.

#include "AliHLTGlobalTriggerComponent.h"

ClassImp(AliHLTGlobalTriggerComponent)


AliHLTGlobalTriggerComponent::AliHLTGlobalTriggerComponent() :
	AliHLTTrigger()
{
  // Default constructor.
}


AliHLTGlobalTriggerComponent::~AliHLTGlobalTriggerComponent()
{
  // Default destructor.
}


void AliHLTGlobalTriggerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  // Returns the output data size estimate.

  constBase = strlen(GetTriggerName()) + 1;
  inputMultiplier = 1;
}


int AliHLTGlobalTriggerComponent::DoTrigger()
{
  // This method will apply the global trigger decision.

  //TODO
  return 0;
}
