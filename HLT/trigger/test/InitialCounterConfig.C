// $Id:  $
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

/// @file   InitialCounterConfig.C
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   12 Nov 2009
/// @brief  Implementation of initial counter configuration routine for testing.
///
/// The InitialCounterConfig routine generates the initial counter configuration
/// to test the AliHLTTriggerCounterComponent class.

void InitialCounterConfig()
{
	TMap& counters = AliHLTTriggerCounterComponent::InitialCounterConfig();
	counters.Add(new TObjString("CINTA"), new TObjString("Input interaction trigger"));
	TObjString* name = new TObjString("TRIGGER-A");
	name->SetBit(1<<14); // mark as output trigger
	counters.Add(name, new TObjString("Output trigger A"));
}
