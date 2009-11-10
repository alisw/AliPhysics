//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERMENUITEM_H
#define ALIHLTTRIGGERMENUITEM_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTTriggerMenuItem.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Dec 2008
/// @brief  Declaration of the AliHLTTriggerMenuItem class.

#include "TObject.h"
#include "TString.h"
#include "TArrayL.h"

/**
 * \class AliHLTTriggerMenuItem
 * A trigger menu item is used to store the information for a single entry in the
 * HLT global trigger menu AliHLTTriggerMenu.
 * It stores information about the trigger condition, trigger domain merging
 * expression, trigger priority and the prescalar to apply.
 * The trigger condition is an expression which indicates what must be true
 * for the trigger menu entry to be fired. A fired item will then use the trigger
 * domain merging expression for the computation of the final global trigger domain.
 * All expressions must be valid C++.
 *
 * The symbols used in the trigger condition expressions are assumed to be AliHLTTrigger
 * names, unless they are predefined in the trigger menu symbols table. All symbols should
 * be valid C++ symbol names. However, the '-' character is allowed as a special extention.
 * The '-' character must not be the first character of the symbol and there cannot be
 * any spaces between it and the alphanumeric characters. If there are any spaces then
 * the '-' character is treated as the normal C++ minus operator. For example,
 * "abc-xyz" is a single whole symbol, while "abc - xyz" are two symbols, abc and xyz,
 * separated by a minus operator.
 *
 * Merging expressions can use all the symbols defined in the trigger menu symbols table
 * including all the implicit symbols used in the trigger conditions which are assumed
 * to be AliHLTTrigger names. If a AliHLTTrigger name is not used in a trigger condition
 * expression, but one wants to use the trigger domain in a merging expression, then a
 * predefined symbol must be added to the trigger menu symbol table. As an example, in
 * the following manner:
 * \code
 * AliHLTGlobalTriggerConfig config("test config");
 * config.AddSymbol("myTriggerName", "bool", "this->Result()", "0", "AliHLTTriggerDecision");
 * \endcode
 * The trigger name "myTriggerName" should be replaced with the actual name of the
 * AliHLTTrigger from which one wants to use the trigger domain result.
 * Symbols with the '-' sign are be handled automatically and will be replaced
 * by their appropriate versions with the minus signs replaced by underscores.
 * This means that a minus sign in any other location is always treated as an operator.
 * If uncertain then just put spaces around the minus operator.
 *
 * \note The following symbol names are reserved and should not be used in either
 * the trigger condition or merging expressions:
 *   _domain_
 *   _description_
 *   _item_result_
 *   _group_result_
 *   _previous_match_
 *   _trigger_matched_
 *   FillFromMenu
 *   NewEvent
 *   Add
 *   CalculateTriggerDecision
 *   GetCounters
 *   SetCounters
 *   CreateNew
 */
class AliHLTTriggerMenuItem : public TObject
{
 public:
  
  /**
   * Default constructor.
   */
  AliHLTTriggerMenuItem();
  
  /**
   * Default destructor.
   */
  virtual ~AliHLTTriggerMenuItem();
  
  /**
   * Inherited from TObject, this prints the contents of the menu item.
   * \param option  Can be "compact", which will print in the compact format.
   */
  virtual void Print(Option_t* option = "") const;
  
  /**
   * Returns the optional comment string.
   */
  const char* Description() const { return fDescription.Data(); }
  
  /**
   * Set the optional comment string.
   */
  void Description(const char* value) { fDescription = value; }
  
  /**
   * Returns the trigger condition expression.
   */
  const char* TriggerCondition() const { return fConditionExpr.Data(); }
  
  /**
   * Set the trigger condition expression.
   */
  void TriggerCondition(const char* value) { fConditionExpr = value; }
  
  /**
   * Returns the trigger domain merging expression.
   */
  const char* MergeExpression() const { return fDomainExpr.Data(); }
  
  /**
   * Set the trigger domain merging expression.
   */
  void MergeExpression(const char* value) { fDomainExpr = value; }
  
  /**
   * Returns the pre-scalar value.
   */
  UInt_t PreScalar() const { return fPrescalar; }
  
  /**
   * Set the pre-scalar value. A value of zero turns off the prescalar.
   */
  void PreScalar(UInt_t value) { fPrescalar = value; }
  
  /**
   * Returns the priority value.
   */
  UInt_t Priority() const { return fPriority; }
  
  /**
   * Set the priority value. Higher values give a higher priority.
   */
  void Priority(UInt_t value) { fPriority = value; }

 private:
  
  TString fDescription;  /// Optional description or comment string.
  TString fConditionExpr;  /// The trigger condition expression.
  TString fDomainExpr;  /// Trigger domain merging expression.
  UInt_t fPrescalar;  /// Pre-scalar value used to optionally reduce the trigger rate. Every modulus n'th event is triggered, where n equals the pre-scalar value.
  UInt_t fPriority;  /// Priority of the trigger menu item. Higher values have higher priority.
  
  ClassDef(AliHLTTriggerMenuItem, 3) // Trigger menu item for global HLT trigger.
};

#endif // ALIHLTTRIGGERMENUITEM_H

