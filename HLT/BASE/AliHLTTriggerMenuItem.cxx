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

/// @file   AliHLTTriggerMenuItem.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Dec 2008
/// @brief  Implementation of the AliHLTTriggerMenuItem class.
///
/// The AliHLTTriggerMenuItem contains information about a entry in the global
/// HLT trigger menu.

#include "AliHLTTriggerMenuItem.h"
#include "Riostream.h"

ClassImp(AliHLTTriggerMenuItem)


AliHLTTriggerMenuItem::AliHLTTriggerMenuItem() :
  TObject(),
  fDescription(),
  fConditionExpr(),
  fDomainExpr(),
  fPrescalar(0),
  fPriority(0),
  fScaleDown(1)
{
  // Default constructor.
  
  DefaultResult(true); // The default result for the item is always true.
}


AliHLTTriggerMenuItem::~AliHLTTriggerMenuItem()
{
  // Default destructor.
}


void AliHLTTriggerMenuItem::Print(Option_t* option) const
{
  // Prints the contents of the trigger menu item.
  
  TString opt = option;
  if (opt.Contains("compact"))
  {
    cout << "{fConditionExpr = \"" << fConditionExpr.Data()
         << "\", fDomainExpr = \"" << fDomainExpr.Data()
         << "\", fPrescalar = " << fPrescalar
         << ", fScaleDown = " << fScaleDown
         << ", fPriority = " << fPriority
         << ", default result = " << (DefaultResult() ? "true" : "false")
         << "}" << endl;
  }
  else
  {
    cout << "                    Description = " << fDescription.Data() << endl;
    cout << "   Trigger condition expression = " << fConditionExpr.Data() << endl;
    cout << "Trigger domain merge expression = " << fDomainExpr.Data() << endl;
    cout << "                     Pre-scalar = " << fPrescalar << endl;
    cout << "                       Priority = " << fPriority << endl;
    cout << "                     Scale-down = " << fScaleDown << endl;
    cout << "  Default global trigger result = " << (DefaultResult() ? "true" : "false") << endl;
  }
}

