//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTTRIGGERMENUSYMBOL_H
#define ALIHLTTRIGGERMENUSYMBOL_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTTriggerMenuSymbol.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 Dec 2008
/// @brief  Declaration of the AliHLTTriggerMenuSymbol class.

#include "TObject.h"
#include "TString.h"
#include "AliHLTDataTypes.h"
#include "AliHLTDomainEntry.h"

/**
 * \class AliHLTTriggerMenuSymbol
 * The trigger menu symbol is used to either label a TObject variable found in
 * one of the input data blocks or create a constant object that can be used in
 * various trigger menu expressions.
 * It essentially represents a C++ variable with a name and type. This is filled
 * with the information found in an input TObject, as specified by the
 * <i>AssignExpression()</i> method. If nothing can be assigned because an
 * appropriate TObject cannot be found, then a default value is used.
 * The correct TObject is found by searching for a TObject with a matching class
 * type as given by <i>ObjectClass()</i> from the input data blocks. The input
 * data block must also have the data block type and specification match those
 * given in the <i>BlockType()</i> method.
 *
 * The symbol name should normally be a valid C++ symbol name. This normally
 * excludes the minus sign from being used in the name. However, since CTP trigger
 * names use the minus sign extensively, the minus sign is also allowed in HLT
 * trigger menu symbol names. This is implicitly converted to an underscore in the
 * name. However, the original name with minus signs is retained and available
 * with the <i>RealName()</i> method. In addition, the <i>RealName</i> is used
 * for display purposes.
 * \note Trigger menu symbols names that use minus signs are synonymous to those
 * that use understores in the same locations. For this reason it is better to use
 * one or the other form, but not both at the same time.
 * \note The following symbol names are reserved and should not be used:
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
class AliHLTTriggerMenuSymbol : public TObject
{
 public:
  
  /**
   * Default constructor.
   */
  AliHLTTriggerMenuSymbol();
  
  /**
   * Default destructor.
   */
  virtual ~AliHLTTriggerMenuSymbol();
  
  /**
   * Inherited from TObject, this prints the contents of the symbol.
   * \param option  Can be "compact", which will print in the compact format.
   */
  virtual void Print(Option_t* option = "") const;
  
  /**
   * Returns the real name which can differ from the symbol name
   * by '-' characters instead of '_'. This approach has been introduced
   * to allow trigger and item names containing '-' characters and thus
   * violating the C++ naming conventions.
   */
  const char* RealName() const;
  
  /**
   * Returns the valid C++ symbol name.
   */
  const char* Name() const { return fName.Data(); }

  /**
   * Set the symbol name. It can contain '-' characters.
   */
  void Name(const char* value);
  
  /**
   * Returns the symbol data type.
   */
  const char* Type() const { return fType.Data(); }
  
  /**
   * Set the symbol data type.
   */
  void Type(const char* value) { fType = value; }
  
  /**
   * Returns the data block type and specification from which the symbol is fetched.
   */
  const AliHLTDomainEntry& BlockType() const { return fBlockType; }
  
  /**
   * Set the data block type and specification from which the symbol is fetched.
   * \param type  The data block type and origin to use.
   */
  void BlockType(const AliHLTComponentDataType& type)
  {
    fBlockType = AliHLTDomainEntry(type);
  }
  
  /**
   * Set the data block type and specification from which the symbol is fetched.
   * \param blocktype  The data block type string of the data block. The value
   *    kAliHLTAnyDataTypeID can be used to specify the any type wild card value.
   * \param origin  The origin of the data block, such as the detector name.
   *    The value kAliHLTDataOriginAny can be used to specify the any origin
   *    wild card value.
   */
  void BlockType(const char* blocktype, const char* origin)
  {
    fBlockType = AliHLTDomainEntry(blocktype, origin);
  }
  
  /**
   * Set the data block type and specification from which the symbol is fetched.
   * \param type  The data block type and origin to use.
   * \param spec  The data block specification to use.
   */
  void BlockType(const AliHLTComponentDataType& type, UInt_t spec)
  {
    fBlockType = AliHLTDomainEntry(type, spec);
  }
  
  /**
   * Set the data block type and specification from which the symbol is fetched.
   * \param blocktype  The data block type string of the data block. The value
   *    kAliHLTAnyDataTypeID can be used to specify the any type wild card value.
   * \param origin  The origin of the data block, such as the detector name.
   *    The value kAliHLTDataOriginAny can be used to specify the any origin
   *    wild card value.
   * \param spec  The data block specification to use.
   */
  void BlockType(const char* blocktype, const char* origin, UInt_t spec)
  {
    fBlockType = AliHLTDomainEntry(blocktype, origin, spec);
  }
  
  /**
   * Returns the class name of the object in the data block.
   */
  const char* ObjectClass() const { return fClass.Data(); }
  
  /**
   * Set the class name of the object in the data block.
   */
  void ObjectClass(const char* value) { fClass = value; }
  
  /**
   * Returns the expression to assign the symbol value.
   */
  const char* AssignExpression() const { return fAssignExpr.Data(); }
  
  /**
   * Set the expression to assign the symbol value.
   * The keyword 'this' is used as the place holder of the object in the matched
   * data block. For example if we want to get a public attribute names xyz from
   * the object to assign to the symbol we would write the assignment expression
   * as "this->xyz".
   */
  void AssignExpression(const char* value) { fAssignExpr = value; }
  
  /**
   * Returns the default value expression.
   */
  const char* DefaultValue() const { return fDefaultValue.Data(); }
  
  /**
   * Set the default value expression.
   */
  void DefaultValue(const char* value) { fDefaultValue = value; }

 private:
  
  TString fName;  /// The name of the symbol (Must be a valid C++ variable name but can contain the '-' character as an extention).
  TString fType;  /// The data type of the symbol (Must be a valid C++ type name).
  AliHLTDomainEntry fBlockType;  /// The data block type and specification this symbol is fetched from.
  TString fClass;  /// The class name of the object to read from (Must be a valid C++ class name).
  TString fAssignExpr;  /// The expression to assign to the symbol (Must be a valid C++ expression).
  TString fDefaultValue;  /// The default value this symbol is set to (this must be a valid C++ expression).
  TString fRealName; /// The name of the symbol, differs from fName by replaced '-' chars.
  
  ClassDef(AliHLTTriggerMenuSymbol, 4) // Trigger menu item for global HLT trigger.
};

#endif // ALIHLTTRIGGERMENUSYMBOL_H

