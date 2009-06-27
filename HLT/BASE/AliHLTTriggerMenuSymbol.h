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
 * TODO
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
   * Returns the symbol name.
   */
  const char* Name() const { return fName.Data(); }
  
  /**
   * Set the symbol name.
   */
  void Name(const char* value) { fName = value; }
  
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
   * data block.
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
  
  TString fName;  /// The name of the symbol (Must be a valid C++ variable name).
  TString fType;  /// The data type of the symbol (Must be a valid C++ type name).
  AliHLTDomainEntry fBlockType;  /// The data block type and specification this symbol is fetched from.
  TString fClass;  /// The class name of the object to read from (Must be a valid C++ class name).
  TString fAssignExpr;  /// The expression to assign to the symbol (Must be a valid C++ expression).
  TString fDefaultValue;  /// The default value this symbol is set to (this must be a valid C++ expression).
  
  ClassDef(AliHLTTriggerMenuSymbol, 2) // Trigger menu item for global HLT trigger.
};

#endif // ALIHLTTRIGGERMENUSYMBOL_H

