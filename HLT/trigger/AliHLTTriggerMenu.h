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

/**
 * \class AliHLTTriggerMenu
 * This class is an abstract class. Classes which derive from this class should
 * implement the logic for a particular trigger menu. The AliHLTTriggerMenu class
 * creates a class deriving from AliHLTTriggerMenu on the fly to implement the
 * trigger logic for that particular trigger menu.
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
   * Adds a new symbol to the trigger menu.
   */
  void AddSymbol(const AliHLTTriggerMenuSymbol& entry)
  {
    new (fSymbols[fSymbols.GetEntriesFast()]) AliHLTTriggerMenuSymbol(entry);
  }
  
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
  
 private:
  
  TString fName;  /// Name of the trigger menu.
  TClonesArray fSymbols;  /// List of symbols used in trigger expressions.
  TClonesArray fItems;  /// List of trigger menu items.
  
  ClassDef(AliHLTTriggerMenu, 2) // Trigger menu for the global HLT trigger.
};

#endif // ALIHLTTRIGGERMENU_H

