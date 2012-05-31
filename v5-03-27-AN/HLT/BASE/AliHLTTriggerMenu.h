//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERMENU_H
#define ALIHLTTRIGGERMENU_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTTriggerMenu.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   19 Dec 2008
/// @brief  Declaration of the AliHLTTriggerMenu base class.

#include "TObject.h"
#include "TString.h"
#include "TClonesArray.h"
#include "AliHLTTriggerMenuSymbol.h"
#include "AliHLTTriggerMenuItem.h"
#include "AliHLTTriggerDomain.h"

/**
 * \class AliHLTTriggerMenu
 * The trigger menu specifies the HLT global trigger component configuration.
 * The global trigger has a list of individual input trigger components deriving
 * from AliHLTTrigger. Each one of these triggers is named. Thus, the trigger menu
 * is a list of trigger condition expressions, where the variables inside the
 * expressions are the names of the input triggers. In this way the global trigger
 * can be configured in a powerful manner using C++ expressions.
 * Attached to each trigger condition expression is a trigger domain merging
 * expression, which indicates how the final global trigger domain should be
 * calculated from the fragments coming from each individual trigger component.
 * Each entry in the trigger menu has a priority. These are set explicitly and
 * set to zero by default. The higher the priority number for a menu item the
 * higher its priority. Multiple items can have the same priority values.
 *
 * An important concept is the trigger priority group. This is a number of menu
 * items that all have the same priority. If all trigger menu items have the
 * same priority value then there is only one priority group with all the items
 * being members of the same group. Otherwise there can be any number of priority
 * groups with any number of one or more trigger menu items in a group.
 * A trigger menu item belongs to priority group N if it has a priority number equal
 * to N. The trigger menu is then evaluated by the global trigger component, such
 * that the highest priority trigger groups (largest N) are evaluated first.
 * Trigger items within a priority group are then evaluated in the order they
 * were added to the trigger menu. Thus, the first item added to the menu in group
 * N is evaluated first, and the last added to group N is evaluated last.
 * Inside a priority group all trigger menu items have their trigger condition
 * expressions evaluated. This is equivalent to evaluating the group's trigger
 * condition expression, where the group's expression is a concatenation of the
 * individual trigger condition expressions of the items in priority group N.
 * This means that the trigger conditions expressions (indeed, also the trigger
 * domain merging expressions) are allowed to have a dangling trailing operator.
 * The trailing operator will then make sense in the full concatenated expression.
 * If no such trailing operator is found then the default trigger conditions operator
 * is used implicitly for the concatenation, as defined in the trigger menu.
 * If the full concatenated condition expression evaluates to true then the priority
 * group's result is also true and the output trigger domain can be calculated.
 * This is done by taking all the merging expressions from only those trigger menu
 * items whose trigger condition expression fragments were true, and concatenating
 * those merging expression fragments together to arrive at the full merging expression.
 * The final trigger domain is calculated by evaluating the merging expression.
 * Note that the concatenation of the merging expression fragments works in the
 * same manner as the trigger condition expressions. So a trailing operator is
 * allowed in each trigger menu item's merging expression, and is implicitly removed
 * if not needed, but used to concatenate with the next triggered expression.
 * The default domain merging operator is used if no trailing operator is present
 * but a concatenation is required.
 * The evaluation of trigger menu items stops at the first priority group whose
 * trigger condition expression evaluated to true. This is important to force
 * mutually exclusive precedence of a higher priority trigger or triggers.
 * The two extremes of this model are:
 * - All trigger menu entries have the same priority so they are all part of the
 *   same priority group, and thus every trigger menu item is evaluated.
 * - Every trigger menu entry has a different priority, so each forms its own priority
 *   group, and the trigger evaluation stops at the first highest priority item
 *   that matches the trigger condition.
 * Another way to look at the model is that priority groups are mutually exclusive.
 * Trigger menu items from two different priority groups cannot be active at the same
 * time. While more than one trigger menu item can be active at the same time if they
 * are from the same priority group.
 * Yet another view at the model is that a priority group forms an explicit trigger
 * condition and trigger domain merging expression, while trigger menu items specify
 * the expression fragments that are concatenated together implicitly. If there is
 * just one trigger menu item in a priority group then the groups expressions are
 * explicit. On the other hand, for multiple items in a group they form implicit
 * expression fragments.
 *
 * \note CTP trigger class names can be used in the trigger menu since the global
 *   trigger will generate and add corresponding trigger decision objects to the
 *   logic on the fly.
 *   In addition, for software triggers, a special SOFTWARE trigger decision is
 *   generated and the SOFTWARE name can be used in the trigger menu for this.
 *   If the software trigger is a calibration event then a trigger decision with
 *   the name CALIBRATION is generated instead. START_OF_DATA and END_OF_DATA
 *   symbols are similarly defined for the start and end of data events respectively.
 */
class AliHLTTriggerMenu : public TObject
{
 public:
  
  /**
   * Default constructor.
   */
  AliHLTTriggerMenu();
  
  /**
   * Default destructor.
   */
  virtual ~AliHLTTriggerMenu();
  
  /**
   * Inherited from TObject, this prints the contents of the trigger menu.
   * \param option  Can be "short" which will print the short format.
   */
  virtual void Print(Option_t* option = "") const;
  
  /**
   * This method removes all items and symbols from the trigger menu.
   * \note The name and default values are not changed. Only the items and symbol
   *    lists are cleared.
   * \param  option  This is passed onto the internal fSymbols and fItems TClonesArrays.
   * The method is inherited from TObject.
   */
  virtual void Clear(Option_t* option = "");

  /**
   * Copy constructor performs a deep copy of the object.
   * \param  obj  Object to copy from.
   */
  AliHLTTriggerMenu(const AliHLTTriggerMenu& obj);
  
  /**
   * Assignment operator performs a deep copy of the object.
   * \param  obj  Object to copy from.
   * \return  This object is returned after being replaced by a copy of <i>obj</i>.
   */
  AliHLTTriggerMenu& operator = (const AliHLTTriggerMenu& obj);
  
  /**
   * Inherited from TObject. Returns the name of the trigger menu.
   */
  virtual const char* GetName() const { return fName.Data(); }
  
  /**
   * Returns the name of the trigger menu.
   */
  const char* Name() const { return fName.Data(); }
  
  /**
   * Sets the name of the trigger menu.
   */
  void Name(const char* name) { fName = name; }
  
  /**
   * Returns the number of symbols in the trigger menu.
   */
  UInt_t NumberOfSymbols() const { return UInt_t(fSymbols.GetEntriesFast()); }
  
  /**
   * Fetches the i'th trigger menu symbol.
   */
  const AliHLTTriggerMenuSymbol* Symbol(UInt_t i) const
  {
    if (i >= UInt_t(fSymbols.GetEntriesFast())) return NULL;
    return static_cast<const AliHLTTriggerMenuSymbol*>( fSymbols.UncheckedAt(Int_t(i)) );
  }
  
  /**
   * Fetches the i'th trigger menu symbol for editing.
   */
  AliHLTTriggerMenuSymbol* Symbol(UInt_t i)
  {
    if (i >= UInt_t(fSymbols.GetEntriesFast())) return NULL;
    return static_cast<AliHLTTriggerMenuSymbol*>( fSymbols.UncheckedAt(Int_t(i)) );
  }
  
  /**
   * Adds a new symbol to the trigger menu. If the symbol being added already
   * exists in the trigger menu then the new symbol will not be added.
   * \param entry  The new trigger menu symbol being added.
   */
  void AddSymbol(const AliHLTTriggerMenuSymbol& entry);
  
  /**
   * Returns the array of symbols.
   */
  const TClonesArray& SymbolArray() const { return fSymbols; }
  
  /**
   * Returns the number of items in the trigger menu.
   */
  UInt_t NumberOfItems() const { return UInt_t(fItems.GetEntriesFast()); }
  
  /**
   * Fetches the i'th trigger menu item.
   */
  const AliHLTTriggerMenuItem* Item(UInt_t i) const
  {
    if (i >= UInt_t(fItems.GetEntriesFast())) return NULL;
    return static_cast<const AliHLTTriggerMenuItem*>( fItems.UncheckedAt(Int_t(i)) );
  }
  
  /**
   * Fetches the i'th trigger menu item for editing.
   */
  AliHLTTriggerMenuItem* Item(UInt_t i)
  {
    if (i >= UInt_t(fItems.GetEntriesFast())) return NULL;
    return static_cast<AliHLTTriggerMenuItem*>( fItems.UncheckedAt(Int_t(i)) );
  }
  
  /**
   * Adds a new entry to the trigger menu.
   */
  void AddItem(const AliHLTTriggerMenuItem& entry)
  {
    new (fItems[fItems.GetEntriesFast()]) AliHLTTriggerMenuItem(entry);
  }
  
  /**
   * Returns the array of menu items.
   */
  const TClonesArray& ItemsArray() const { return fItems; }
  
  /**
   * Sets the default trigger description to use if the global trigger does not
   * fire and returns a negative result.
   */
  void DefaultDescription(const char* value) { fDefaultDescription = value; }
  
  /**
   * Returns the default trigger description to use if the global trigger does not
   * fire and returns a negative result.
   */
  const char* DefaultDescription() const { return fDefaultDescription.Data(); }
  
  /**
   * Sets the default trigger domain to use if the global trigger does not
   * fire and returns a negative result.
   */
  void DefaultTriggerDomain(const AliHLTTriggerDomain& value) { fDefaultDomain = value; }
  
  /**
   * Returns the default trigger domain to use if the global trigger does not
   * fire and returns a negative result.
   */
  const AliHLTTriggerDomain& DefaultTriggerDomain() const { return fDefaultDomain; }
  
  /**
   * Returns the default trigger domain for modification.
   */
  AliHLTTriggerDomain& DefaultTriggerDomain() { return fDefaultDomain; }
  
  /**
   * Sets the default operator used to merge trigger conditions that are matched from
   * the same trigger menu priority group.
   */
  void DefaultConditionOperator(const char* value) { fDefaultConditionOperator = value; }
  
  /**
   * Returns the default operator used to merge trigger conditions that are matched from
   * the same trigger menu priority group.
   */
  const char* DefaultConditionOperator() const { return fDefaultConditionOperator.Data(); }
  
  /**
   * Sets the default operator used to merge trigger domains that are matched from
   * the same trigger menu priority group.
   */
  void DefaultDomainOperator(const char* value) { fDefaultDomainOperator = value; }
  
  /**
   * Returns the default operator used to merge trigger domains that are matched from
   * the same trigger menu priority group.
   */
  const char* DefaultDomainOperator() const { return fDefaultDomainOperator.Data(); }
  
  /**
   * Returns the default result for the global trigger if no item is matched.
   */
  bool DefaultResult() const { return TestBit(BIT(15)) == 1; }
  
  /**
   * Set the default result for the global trigger if no item is matched.
   */
  void DefaultResult(bool value) { SetBit(BIT(15), value); }
  
 private:
  
  TString fName;  /// Name of the trigger menu.
  TClonesArray fSymbols;  /// List of symbols used in trigger expressions.
  TClonesArray fItems;  /// List of trigger menu items.
  TString fDefaultDescription; /// The default trigger description to use for negative global triggers.
  AliHLTTriggerDomain fDefaultDomain;  /// The default trigger domain to use for negative global triggers.
  TString fDefaultConditionOperator;  /// The default operator to use to merge trigger conditions from the same priority group.
  TString fDefaultDomainOperator;  /// The default operator to use to merge trigger domains from the same priority group.
  
  ClassDef(AliHLTTriggerMenu, 4) // Trigger menu for the global HLT trigger.
};

#endif // ALIHLTTRIGGERMENU_H

