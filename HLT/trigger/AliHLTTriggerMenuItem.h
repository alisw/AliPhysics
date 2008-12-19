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
 * TODO
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
   * Returns the trigger condition expression.
   */
  const char* TriggerCondision() const { return fConditionExpr.Data(); }
  
  /**
   * Set the trigger condition expression.
   */
  void TriggerCondision(const char* value) { fConditionExpr = value; }
  
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

 private:
  
  TString fConditionExpr;  /// The trigger condition expression.
  TString fDomainExpr;  /// Trigger domain merging expression.
  UInt_t fPrescalar;  /// Pre-scalar value used to optionally reduce the trigger rate. Every modulus n'th event is triggered, where n equals the pre-scalar value.
  
  ClassDef(AliHLTTriggerMenuItem, 1) // Trigger menu item for global HLT trigger.
};

#endif // ALIHLTTRIGGERMENUITEM_H

