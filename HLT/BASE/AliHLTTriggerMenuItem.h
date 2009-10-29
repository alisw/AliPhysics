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
 * The trigger condition is an expression which indicates that must be true
 * for the trigger menu entry to be fired. A fired item will then use the trigger
 * domain merging expression for the computation of the final global trigger domain.
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

