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

/// @file   AliHLTTriggerMenu.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Dec 2008
/// @brief  Implementation of the AliHLTTriggerMenu base class.
///
/// The AliHLTTriggerMenu class implements the HLT global trigger menu,
/// which defines how and on what events the HLT triggers.

#include "AliHLTTriggerMenu.h"
#include "Riostream.h"

ClassImp(AliHLTTriggerMenu)


AliHLTTriggerMenu::AliHLTTriggerMenu() :
  TObject(),
  fName("Unknown"),
  fItems(AliHLTTriggerMenuItem::Class(), 100)
{
  // Default constructor.
}


AliHLTTriggerMenu::~AliHLTTriggerMenu()
{
  // Default destructor.
}


AliHLTTriggerMenu::AliHLTTriggerMenu(const AliHLTTriggerMenu& obj) :
  TObject(obj),
  fName(obj.fName),
  fItems(AliHLTTriggerMenuItem::Class(), obj.fItems.GetEntriesFast())
{
  // Copy constructor performs a deep copy.
  
  for (UInt_t i = 0; i < obj.NumberOfItems(); i++)
  {
    AddItem(*obj.Item(i));
  }
}


AliHLTTriggerMenu& AliHLTTriggerMenu::operator = (const AliHLTTriggerMenu& obj)
{
  // Assignment operator performs a deep copy.
  
  if (this != &obj)
  {
    TObject::operator = (obj);
    fName = obj.fName;
    fItems.Clear();
    for (UInt_t i = 0; i < obj.NumberOfItems(); i++)
    {
      AddItem(*obj.Item(i));
    }
  }
  return *this;
}


void AliHLTTriggerMenu::Print(Option_t* option) const
{
  // Prints the contents of the trigger menu.
  
  cout << "HLT Trigger Menu: " << fName.Data();
  TString opt = option;
  if (opt.Contains("short"))
  {
    cout << ", contains " << NumberOfItems() << " entries." << endl;
    return;
  }
  cout << endl;
  cout << setw(10) << "Prescalar" <<  " | "
       << setw(30) << "Trigger condision" << " | "
       << setw(30) << "Domain merge expression" << setw(0) << endl;
  cout << setfill('-') << setw(10) << "-" <<  "-+-"
       << setw(30) << "-" << "-+-"
       << setw(30) << "-" << setw(0) << endl;
  for (UInt_t i = 0; i < NumberOfItems(); i++)
  {
    Item(i)->Print("compact");
  }
}

