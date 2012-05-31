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

/// @file   AliHLTGlobalTriggerConfig.cxx
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 Dec 2008
/// @brief  Implementation of the AliHLTGlobalTriggerConfig class.
///
/// The AliHLTGlobalTriggerConfigComponent class is an interface class used to create
/// global HLT trigger menu configurations.

#include "AliHLTGlobalTriggerConfig.h"
#include "AliHLTTriggerMenu.h"
#include "Riostream.h"

ClassImp(AliHLTGlobalTriggerConfig)

AliHLTTriggerMenu* AliHLTGlobalTriggerConfig::fgMenu = NULL;


AliHLTGlobalTriggerConfig::AliHLTGlobalTriggerConfig(const char* name)
{
  // Default constructor.
  
  NewMenu(name);
}


AliHLTGlobalTriggerConfig::~AliHLTGlobalTriggerConfig()
{
  // Default destructor.
}


void AliHLTGlobalTriggerConfig::NewMenu(const char* name)
{
  // Creates a new trigger menu.
  
  if (fgMenu != NULL)
  {
    delete fgMenu;
    fgMenu = NULL;
  }
  fgMenu = new AliHLTTriggerMenu;
  fgMenu->Name(name);
}


void AliHLTGlobalTriggerConfig::Clear()
{
  // Deletes the current trigger menu.
  
  if (fgMenu != NULL)
  {
    delete fgMenu;
    fgMenu = NULL;
  }
}


void AliHLTGlobalTriggerConfig::AddSymbol(
    const char* name, const char* type, const char* defaultExpr
  )
{
  // Adds a new constant symbol to the trigger menu.
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuSymbol entry;
  entry.Name(name);
  entry.Type(type);
  entry.ObjectClass("");
  entry.AssignExpression("");
  entry.DefaultValue(defaultExpr);
  fgMenu->AddSymbol(entry);
}


void AliHLTGlobalTriggerConfig::AddSymbol(
    const char* name, const char* type, const char* assignExpr,
    const char* defaultExpr, const char* className
  )
{
  // Adds a new symbol to the trigger menu.
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuSymbol entry;
  entry.Name(name);
  entry.Type(type);
  entry.ObjectClass(className);
  entry.AssignExpression(assignExpr);
  entry.DefaultValue(defaultExpr);
  fgMenu->AddSymbol(entry);
}


void AliHLTGlobalTriggerConfig::AddSymbol(
    const char* name, const char* type, const char* assignExpr,
    const char* defaultExpr, const char* className,
    const AliHLTComponentDataType& blockType
  )
{
  // Adds a new symbol to the trigger menu.
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuSymbol entry;
  entry.Name(name);
  entry.Type(type);
  entry.ObjectClass(className);
  entry.AssignExpression(assignExpr);
  entry.DefaultValue(defaultExpr);
  entry.BlockType(blockType);
  fgMenu->AddSymbol(entry);
}


void AliHLTGlobalTriggerConfig::AddSymbol(
    const char* name, const char* type, const char* assignExpr,
    const char* defaultExpr, const char* className,
    const AliHLTComponentDataType& blockType, UInt_t spec
  )
{
  // Adds a new symbol to the trigger menu.
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuSymbol entry;
  entry.Name(name);
  entry.Type(type);
  entry.ObjectClass(className);
  entry.AssignExpression(assignExpr);
  entry.DefaultValue(defaultExpr);
  entry.BlockType(blockType, spec);
  fgMenu->AddSymbol(entry);
}


void AliHLTGlobalTriggerConfig::AddSymbol(
    const char* name, const char* type, const char* assignExpr,
    const char* defaultExpr, const char* className,
    const char* blockType, const char* origin
  )
{
  // Adds a new symbol to the trigger menu.
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuSymbol entry;
  entry.Name(name);
  entry.Type(type);
  entry.ObjectClass(className);
  entry.AssignExpression(assignExpr);
  entry.DefaultValue(defaultExpr);
  entry.BlockType(blockType, origin);
  fgMenu->AddSymbol(entry);
}


void AliHLTGlobalTriggerConfig::AddSymbol(
    const char* name, const char* type, const char* assignExpr,
    const char* defaultExpr, const char* className,
    const char* blockType, const char* origin, UInt_t spec
  )
{
  // Adds a new symbol to the trigger menu.
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuSymbol entry;
  entry.Name(name);
  entry.Type(type);
  entry.ObjectClass(className);
  entry.AssignExpression(assignExpr);
  entry.DefaultValue(defaultExpr);
  entry.BlockType(blockType, origin, spec);
  fgMenu->AddSymbol(entry);
}
    
    
void AliHLTGlobalTriggerConfig::AddItem(
    UInt_t priority, const char* conditionExpr, const char* domainExpr,
    UInt_t prescalar, const char* description, bool defaultResult
  )
{
  // Adds a new entry to the trigger menu with a particular priority.
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuItem entry;
  entry.TriggerCondition(conditionExpr);
  entry.MergeExpression(domainExpr);
  entry.PreScalar(prescalar);
  entry.Priority(priority);
  entry.DefaultResult(defaultResult);
  if (description != NULL) entry.Description(description);
  fgMenu->AddItem(entry);
}


void AliHLTGlobalTriggerConfig::AddItem(
    UInt_t priority, const char* conditionExpr, const char* domainExpr,
    const char* description, Double_t scaledown, bool defaultResult
  )
{
  // Adds a new entry to the trigger menu with a particular priority.
  
  if (scaledown < 0)
  {
    cerr << "ERROR: Cannot have a scale-down value smaller than 0. But a value of "
         << scaledown << " was specified. The valid range is [0..100]." << endl;
    return;
  }
  if (scaledown < 0)
  {
    cerr << "ERROR: Cannot have a scale-down value larger than 100. But a value of "
         << scaledown << " was specified. The valid range is [0..100]." << endl;
    return;
  }
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuItem entry;
  entry.TriggerCondition(conditionExpr);
  entry.MergeExpression(domainExpr);
  entry.Priority(priority);
  entry.ScaleDown(scaledown / 100.);
  entry.DefaultResult(defaultResult);
  if (description != NULL) entry.Description(description);
  fgMenu->AddItem(entry);
}


void AliHLTGlobalTriggerConfig::AddItem(
    const char* conditionExpr, const char* domainExpr, UInt_t prescalar,
    const char* description
  )
{
  // Adds a new entry to the trigger menu.
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuItem entry;
  entry.TriggerCondition(conditionExpr);
  entry.MergeExpression(domainExpr);
  entry.PreScalar(prescalar);
  if (description != NULL) entry.Description(description);
  fgMenu->AddItem(entry);
}


void AliHLTGlobalTriggerConfig::AddItem(
    const char* conditionExpr, const char* domainExpr, const char* description
  )
{
  // Adds a new entry to the trigger menu.
  
  if (fgMenu == NULL) NewMenu("");
  
  AliHLTTriggerMenuItem entry;
  entry.TriggerCondition(conditionExpr);
  entry.MergeExpression(domainExpr);
  if (description != NULL) entry.Description(description);
  fgMenu->AddItem(entry);
}


void AliHLTGlobalTriggerConfig::SetDefaultTriggerDescription(const char* description)
{
  // Sets the default trigger decription.
  
  if (fgMenu == NULL) NewMenu("");
  fgMenu->DefaultDescription(description);
}


void AliHLTGlobalTriggerConfig::SetDefaultTriggerDomain(const AliHLTTriggerDomain& domain)
{
  // Sets the default trigger domain.
  
  if (fgMenu == NULL) NewMenu("");
  fgMenu->DefaultTriggerDomain(domain);
}


AliHLTTriggerDomain& AliHLTGlobalTriggerConfig::DefaultTriggerDomain()
{
  // Returns the default trigger domain for the current trigger menu.
  
  if (fgMenu == NULL) NewMenu("");
  return fgMenu->DefaultTriggerDomain();
}


void AliHLTGlobalTriggerConfig::SetDefaultConditionOperator(const char* op)
{
  // Sets the default operator for trigger condition merging.
  
  if (fgMenu == NULL) NewMenu("");
  fgMenu->DefaultConditionOperator(op);
}


void AliHLTGlobalTriggerConfig::SetDefaultDomainOperator(const char* op)
{
  // Sets the default operator for trigger domain merging.
  
  if (fgMenu == NULL) NewMenu("");
  fgMenu->DefaultDomainOperator(op);
}


void AliHLTGlobalTriggerConfig::SetDefaultResult(bool value)
{
  // Sets the default result when no item is matched.
  
  if (fgMenu == NULL) NewMenu("");
  fgMenu->DefaultResult(value);
}


void AliHLTGlobalTriggerConfig::Print(Option_t* option) const
{
  // Prints the contents of the current trigger menu being manipulated.
  
  if (fgMenu != NULL)
  {
    fgMenu->Print(option);
  }
  else
  {
    cout << "No trigger menu currently being configured, it is empty." << endl;
  }
}

