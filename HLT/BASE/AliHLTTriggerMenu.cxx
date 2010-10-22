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

/// @file   AliHLTTriggerMenu.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Dec 2008
/// @brief  Implementation of the AliHLTTriggerMenu base class.
///
/// The AliHLTTriggerMenu class implements the HLT global trigger menu,
/// which defines how and on what events the HLT triggers.

#include "AliHLTTriggerMenu.h"
#include "Riostream.h"
#include <sstream>

ClassImp(AliHLTTriggerMenu)


AliHLTTriggerMenu::AliHLTTriggerMenu() :
  TObject(),
  fName("Unknown"),
  fSymbols(AliHLTTriggerMenuSymbol::Class(), 100),
  fItems(AliHLTTriggerMenuItem::Class(), 100),
  fDefaultDescription(),
  fDefaultDomain(),
  fDefaultConditionOperator("||"),
  fDefaultDomainOperator("|")
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
  fSymbols(AliHLTTriggerMenuSymbol::Class(), obj.fSymbols.GetEntriesFast()),
  fItems(AliHLTTriggerMenuItem::Class(), obj.fItems.GetEntriesFast()),
  fDefaultDescription(obj.fDefaultDescription),
  fDefaultDomain(obj.fDefaultDomain),
  fDefaultConditionOperator(obj.fDefaultConditionOperator),
  fDefaultDomainOperator(obj.fDefaultDomainOperator)
{
  // Copy constructor performs a deep copy.
  
  for (UInt_t i = 0; i < obj.NumberOfSymbols(); i++)
  {
    AddSymbol(*obj.Symbol(i));
  }
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
    fSymbols.Clear();
    for (UInt_t i = 0; i < obj.NumberOfSymbols(); i++)
    {
      AddSymbol(*obj.Symbol(i));
    }
    fItems.Clear();
    for (UInt_t i = 0; i < obj.NumberOfItems(); i++)
    {
      AddItem(*obj.Item(i));
    }
    fDefaultDescription = obj.fDefaultDescription;
    fDefaultDomain = obj.fDefaultDomain;
    fDefaultConditionOperator = obj.fDefaultConditionOperator;
    fDefaultDomainOperator = obj.fDefaultDomainOperator;
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
  
  size_t prescalarColWidth = 9;
  size_t scaledownColWidth = 9;
  size_t priorityColWidth = 8;
  size_t resultColWidth = 6;
  size_t conditionColWidth = 17;
  size_t mergexprColWidth = 23;
  size_t descColWidth = 11;
  // Find out the required widths of the columns:
  for (UInt_t i = 0; i < NumberOfItems(); i++)
  {
    std::stringstream buf;
    buf.copyfmt(cout);
    buf << Item(i)->PreScalar();
    if (buf.str().length() > prescalarColWidth) prescalarColWidth = buf.str().length();
    buf.str("");
    buf << Item(i)->ScaleDown();
    if (buf.str().length()+1 > scaledownColWidth) scaledownColWidth = buf.str().length()+1;
    buf.str("");
    buf << Item(i)->Priority();
    if (buf.str().length() > priorityColWidth) priorityColWidth = buf.str().length();
    size_t len = strlen(Item(i)->TriggerCondition());
    if (len > conditionColWidth and len < 128) conditionColWidth = len;
    len = strlen(Item(i)->MergeExpression());
    if (len > mergexprColWidth and len < 128) mergexprColWidth = len;
    len = strlen(Item(i)->Description());
    if (len > descColWidth and len < 128) descColWidth = len;
  }
  cout << setw(prescalarColWidth) << "Prescalar" << " | "
       << setw(scaledownColWidth) << "Scaledown" << " | "
       << setw(priorityColWidth) << "Priority" << " | "
       << setw(resultColWidth) << "Result" << " | "
       << setw(conditionColWidth) << "Trigger condition" << " | "
       << setw(mergexprColWidth) << "Domain merge expression" << " | "
       << setw(descColWidth) << "Description" << setw(0) << endl;
  cout << setfill('-') << setw(prescalarColWidth) << "-" << "-+-"
       << setw(scaledownColWidth) << "-" << "-+-"
       << setw(priorityColWidth) << "-" << "-+-"
       << setw(resultColWidth) << "-" << "-+-"
       << setw(conditionColWidth) << "-" << "-+-"
       << setw(mergexprColWidth) << "-" << "-+-"
       << setw(descColWidth) << "-" << setw(0) << setfill(' ') << endl;
  for (UInt_t i = 0; i < NumberOfItems(); i++)
  {
    if (Item(i)->PreScalar() == 0)
      cout << setw(prescalarColWidth) << "off" << " | ";
    else
      cout << setw(prescalarColWidth) << Item(i)->PreScalar() << " | ";
    if (Item(i)->ScaleDown() == 1)
      cout << setw(scaledownColWidth) << "none" << " | ";
    else
      cout << setw(scaledownColWidth-1) << Item(i)->ScaleDown()*100 << "% | ";
    cout << setw(priorityColWidth) << Item(i)->Priority() << " | "
         << setw(resultColWidth) << (Item(i)->DefaultResult() ? "true" : "false") << " | "
         << setw(conditionColWidth) << Item(i)->TriggerCondition() << " | "
         << setw(mergexprColWidth) << Item(i)->MergeExpression() << " | "
         << setw(descColWidth) << Item(i)->Description() << setw(0) << endl;
  }
  if (NumberOfItems() == 0) cout << "(none)" << endl;
  
  cout << "Symbol list:" << endl;
  size_t nameColWidth = 4;
  size_t typeColWidth = 9;
  size_t classColWidth = 10;
  size_t assignColWidth = 14;
  size_t valueColWidth = 13;
  // Find out the required widths of the columns:
  for (UInt_t i = 0; i < NumberOfSymbols(); i++)
  {
    size_t len = strlen(Symbol(i)->RealName());
    if (len > nameColWidth and len < 128) nameColWidth = len;
    len = strlen(Symbol(i)->Type());
    if (len > typeColWidth and len < 128) typeColWidth = len;
    len = strlen(Symbol(i)->ObjectClass());
    if (len > classColWidth and len < 128) classColWidth = len;
    len = strlen(Symbol(i)->AssignExpression());
    if (len > assignColWidth and len < 128) assignColWidth = len;
    len = strlen(Symbol(i)->DefaultValue());
    if (len > valueColWidth and len < 128) valueColWidth = len;
  }
  cout << setw(nameColWidth) << "Name"
       << " | " << setw(typeColWidth) << "Data type"
       << " | " << setw(24) << "Block type & spec"
       << " | " << setw(classColWidth) << "Class name"
       << " | " << setw(assignColWidth) << "Assigned value"
       << " | " << setw(valueColWidth) << "Default value"
       << setw(0) << endl;
  cout << setw(nameColWidth) << setfill('-') << "-"
       << "-+-" << setw(typeColWidth) << "-"
       << "-+-" << setw(24) << "-"
       << "-+-" << setw(classColWidth) << "-"
       << "-+-" << setw(assignColWidth) << "-"
       << "-+-" << setw(valueColWidth) << "-"
       << setw(0) << setfill(' ') << endl;
  for (UInt_t i = 0; i < NumberOfSymbols(); i++)
  {
    cout << setw(nameColWidth) << Symbol(i)->RealName() << " | "
         << setw(typeColWidth) << Symbol(i)->Type() << " | ";
    Symbol(i)->BlockType().Print("noendl");
    cout << " | " << setw(classColWidth) << Symbol(i)->ObjectClass()
         << " | " << setw(assignColWidth) << Symbol(i)->AssignExpression()
         << " | " << setw(valueColWidth) << Symbol(i)->DefaultValue()
         << setw(0) << endl;
  }
  if (NumberOfSymbols() == 0) cout << "(none)" << endl;
  
  cout << "Default trigger condition operator: " << fDefaultConditionOperator << endl;
  cout << "   Default trigger domain operator: " << fDefaultDomainOperator << endl;
  cout << "       Default trigger description: \"" << fDefaultDescription << "\"" << endl;
  cout << "     Default global trigger result: " << (DefaultResult() ? "true" : "false") << "" << endl;
  cout << "Default "; fDefaultDomain.Print();
}


void AliHLTTriggerMenu::Clear(Option_t* option)
{
  // Clears the internal symbol and items arrays.
  
  fSymbols.Clear(option);
  fItems.Clear(option);
}


void AliHLTTriggerMenu::AddSymbol(const AliHLTTriggerMenuSymbol& entry)
{
  // Adds a new symbol to the trigger menu.
  
  for (Int_t i = 0; i < fSymbols.GetEntriesFast(); i++)
  {
    const char* symbolname = static_cast<AliHLTTriggerMenuSymbol*>(fSymbols.UncheckedAt(i))->Name();
    if ( strcmp(symbolname, entry.Name()) == 0 )
    {
      Warning("AliHLTTriggerMenu::AddSymbol",
              "The symbol '%s' already exists in the trigger menu. Will not add the new entry.",
              entry.RealName()
             );
      return;
    }
  }
  new (fSymbols[fSymbols.GetEntriesFast()]) AliHLTTriggerMenuSymbol(entry);
}

