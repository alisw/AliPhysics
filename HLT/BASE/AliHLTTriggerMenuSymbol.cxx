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
  fDefaultValue(),
  fRealName()
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
    cout << "{fRealName = \"" << fRealName.Data()
         << "\", fType = \"" << fType.Data()
         << "\", fBlockType = \"";
    fBlockType.Print("noendl");
    cout << "\", fClass = \"" << fClass
         << "\", fAssignExpr = \"" << fAssignExpr
         << "\", fDefaultValue = \"" << fDefaultValue
         << "\"}" << endl;
  }
  else
  {
    cout << "                    Name = " << Name() << endl;
    cout << "                RealName = " << RealName() << endl;
    cout << "               Data type = " << Type() << endl;
    cout << "              Block type = "; fBlockType.Print();
    cout << "              Class type = " << ObjectClass() << endl;
    cout << "   Assignment expression = " << AssignExpression() << endl;
    cout << "Default value expression = " << DefaultValue() << endl;
  }
}


const char* AliHLTTriggerMenuSymbol::RealName() const
{
  // return the real name of the symbol which can contain '-'

  // backward compatibility with old versions of the class
  if (fRealName.IsNull()) return fName;
  return fRealName;
}

void AliHLTTriggerMenuSymbol::Name(const char* value)
{
  // sets the name and real name of the symbol
  // replaces '-' and '.' in the real name with '_' in order to comply with
  // C++ conventions.
  fRealName = value;
  fName = value;
  fName.ReplaceAll("-", "_");
  fName.ReplaceAll(".", "_");
}
