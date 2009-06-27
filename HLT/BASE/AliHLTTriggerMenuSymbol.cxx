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

/// @file   AliHLTTriggerMenuSymbol.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 Dec 2008
/// @brief  Implementation of the AliHLTTriggerMenuSymbol class.
///
/// The AliHLTTriggerMenuSymbol contains information about a symbol used in the
/// global HLT trigger menu.

#include "AliHLTTriggerMenuSymbol.h"
#include "Riostream.h"

ClassImp(AliHLTTriggerMenuSymbol)


AliHLTTriggerMenuSymbol::AliHLTTriggerMenuSymbol() :
  TObject(),
  fName(),
  fType(),
  fBlockType(kAliHLTAnyDataType),
  fClass(),
  fAssignExpr(),
  fDefaultValue()
{
  // Default constructor.
}


AliHLTTriggerMenuSymbol::~AliHLTTriggerMenuSymbol()
{
  // Default destructor.
}


void AliHLTTriggerMenuSymbol::Print(Option_t* option) const
{
  // Prints the contents of the trigger menu item.
  
  TString opt = option;
  if (opt.Contains("compact"))
  {
    cout << setw(15) << fName << " | "
         << setw(20) << fType << " | ";
    fBlockType.Print("noendl");
    cout << " | " << setw(20) << fClass.Data()
         << " | " << setw(25) << fAssignExpr.Data()
         << " | " << setw(25) << fDefaultValue.Data()
         << setw(0) << endl;
  }
  else
  {
    cout << "                    Name = " << fName.Data() << endl;
    cout << "               Data type = " << fType.Data() << endl;
    cout << "              Block type = "; fBlockType.Print();
    cout << "              Class type = " << fClass.Data() << endl;
    cout << "   Assignment expression = " << fAssignExpr.Data() << endl;
    cout << "Default value expression = " << fDefaultValue.Data() << endl;
  }
}

