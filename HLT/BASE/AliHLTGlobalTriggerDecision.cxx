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

/// @file   AliHLTGlobalTriggerDecision.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   26 Nov 2008
/// @brief  Implementation of the AliHLTGlobalTriggerDecision class.
/// 
/// The global trigger decision class stores the global HLT decision.

#include "AliHLTGlobalTriggerDecision.h"
#include "Riostream.h"

ClassImp(AliHLTGlobalTriggerDecision)


AliHLTGlobalTriggerDecision::AliHLTGlobalTriggerDecision() :
  AliHLTTriggerDecision(),
  fContributingTriggers(AliHLTTriggerDecision::Class()),
  fInputObjects(),
  fCounters()
{
  // Default constructor.
}


AliHLTGlobalTriggerDecision::AliHLTGlobalTriggerDecision(
    bool result, const AliHLTTriggerDomain& triggerDomain, const char* description
  ) :
  AliHLTTriggerDecision(result, "HLTGlobalTrigger", triggerDomain, description),
  fContributingTriggers(AliHLTTriggerDecision::Class()),
  fInputObjects(),
  fCounters()
{
  // Constructor specifying multiple information fields.
  
  Result(result);
  fInputObjects.SetOwner(kTRUE);
}


AliHLTGlobalTriggerDecision::~AliHLTGlobalTriggerDecision()
{
  // Default destructor.
}


void AliHLTGlobalTriggerDecision::Print(Option_t* option) const
{
  // Prints the contents of the trigger decision.
  
  TString opt(option);
  if (opt.Contains("compact"))
  {
    cout << "Global ";
    AliHLTTriggerDecision::Print("");
  }
  else if (opt.Contains("short"))
  {
    cout << "Global ";
    AliHLTTriggerDecision::Print(option);
    cout << "#################### Input trigger decisions ####################" << endl;
    for (Int_t i = 0; i < NumberOfTriggerInputs(); i++)
    {
      TriggerInput(i)->Print(option);
    }
    if (NumberOfTriggerInputs() == 0)
    {
      cout << "(none)" << endl;
    }
  }
  else if (opt.Contains("counters"))
  {
    cout << "Counter\tValue" << endl;
    for (Int_t i = 0; i < fCounters.GetSize(); i++)
    {
      cout << i << "\t" << fCounters[i] << endl;
    }
    if (fCounters.GetSize() == 0)
    {
      cout << "(none)" << endl;
    }
  }
  else
  {
    cout << "Global ";
    AliHLTTriggerDecision::Print(option);
    cout << "#################### Input trigger decisions ####################" << endl;
    for (Int_t i = 0; i < NumberOfTriggerInputs(); i++)
    {
      cout << "-------------------- Input trigger decision " << i << " --------------------" << endl;
      TriggerInput(i)->Print(option);
    }
    if (NumberOfTriggerInputs() == 0)
    {
      cout << "(none)" << endl;
    }
    cout << "###################### Other input objects ######################" << endl;
    for (Int_t i = 0; i < NumberOfInputObjects(); i++)
    {
      cout << "------------------------ Input object " << i << " ------------------------" << endl;
      InputObject(i)->Print(option);
    }
    if (NumberOfInputObjects() == 0)
    {
      cout << "(none)" << endl;
    }
    cout << "#################### Event class counters ####################" << endl;
    cout << "Counter\tValue" << endl;
    for (Int_t i = 0; i < fCounters.GetSize(); i++)
    {
      cout << i << "\t" << fCounters[i] << endl;
    }
    if (fCounters.GetSize() == 0)
    {
      cout << "(none)" << endl;
    }
  }
}

