#ifndef ALIHLTGLOBALTRIGGERCONFIG_H
#define ALIHLTGLOBALTRIGGERCONFIG_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTGlobalTriggerConfig.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   28 Dec 2008
/// @brief  Declaration of the AliHLTGlobalTriggerConfig class.

#include "TObject.h"

class AliHLTComponentDataType;
class AliHLTTriggerMenu;

/**
 * \class AliHLTGlobalTriggerConfig
 * This class is a user interface class used to make it easy to create HLT global
 * trigger configurations (trigger menus).
 */
class AliHLTGlobalTriggerConfig
{
 public:
  
  /**
   * Default constructor.
   */
  AliHLTGlobalTriggerConfig(const char* name = NULL);
  
  /**
   * Default destructor.
   */
  virtual ~AliHLTGlobalTriggerConfig();
  
  /**
   * Creates a new trigger menu. If a trigger menu is already active then the existing
   * one is replaced with the new menu.
   * \param name  The name of the new trigger menu.
   */
  static void NewMenu(const char* name);
  
  /**
   * Deletes the current trigger menu.
   */
  static void Clear();
  
  /**
   * Returns the current trigger menu.
   */
  static const AliHLTTriggerMenu* Menu() { return fgMenu; }
  
  /**
   * Adds a new symbol to the current trigger menu.
   * \param  name  The name of the symbol. It must be a valid C++ variable name.
   * \param  type  The data type of the symbol. It must be a valid C++ data type.
   * \param  assignExpr  The assignment expression for the symbol. It must be a
   *     valid C++ expression. The 'this' keyword is used as a place holder for
   *     the object found in the data block.
   * \param  defaultExpr  The default value to use for the symbol. It must be a
   *     valid C++ expression.
   * \param className  The class name of the object found in the data block from
   *     which to fetch the symbol's value.
   */
  static void AddSymbol(
      const char* name, const char* type, const char* assignExpr,
      const char* defaultExpr, const char* className = "TObject"
    );
  
  /**
   * Adds a new symbol to the current trigger menu.
   * \param  name  The name of the symbol. It must be a valid C++ variable name.
   * \param  type  The data type of the symbol. It must be a valid C++ data type.
   * \param  assignExpr  The assignment expression for the symbol. It must be a
   *     valid C++ expression. The 'this' keyword is used as a place holder for
   *     the object found in the data block.
   * \param  defaultExpr  The default value to use for the symbol. It must be a
   *     valid C++ expression.
   * \param className  The class name of the object found in the data block from
   *     which to fetch the symbol's value.
   * \param blockType  The data block type and origin from which to fetch the
   *     symbol's value.
   */
  static void AddSymbol(
      const char* name, const char* type, const char* assignExpr,
      const char* defaultExpr, const char* className,
      const AliHLTComponentDataType& blockType
    );
  
  /**
   * Adds a new symbol to the current trigger menu.
   * \param  name  The name of the symbol. It must be a valid C++ variable name.
   * \param  type  The data type of the symbol. It must be a valid C++ data type.
   * \param  assignExpr  The assignment expression for the symbol. It must be a
   *     valid C++ expression. The 'this' keyword is used as a place holder for
   *     the object found in the data block.
   * \param  defaultExpr  The default value to use for the symbol. It must be a
   *     valid C++ expression.
   * \param className  The class name of the object found in the data block from
   *     which to fetch the symbol's value.
   * \param blockType  The data block type and origin from which to fetch the
   *     symbol's value.
   * \param spec  The data block specification to use.
   */
  static void AddSymbol(
      const char* name, const char* type, const char* assignExpr,
      const char* defaultExpr, const char* className,
      const AliHLTComponentDataType& blockType, UInt_t spec
    );
  
  /**
   * Adds a new symbol to the current trigger menu.
   * \param  name  The name of the symbol. It must be a valid C++ variable name.
   * \param  type  The data type of the symbol. It must be a valid C++ data type.
   * \param  assignExpr  The assignment expression for the symbol. It must be a
   *     valid C++ expression. The 'this' keyword is used as a place holder for
   *     the object found in the data block.
   * \param  defaultExpr  The default value to use for the symbol. It must be a
   *     valid C++ expression.
   * \param className  The class name of the object found in the data block from
   *     which to fetch the symbol's value.
   * \param blockType  The data block type string of the data block. The value
   *    kAliHLTAnyDataTypeID can be used to specify the any type wild card value.
   * \param origin  The origin of the data block, such as the detector name.
   *    The value kAliHLTDataOriginAny can be used to specify the any origin
   *    wild card value.
   */
  static void AddSymbol(
      const char* name, const char* type, const char* assignExpr,
      const char* defaultExpr, const char* className,
      const char* blockType, const char* origin
    );
  
  /**
   * Adds a new symbol to the current trigger menu.
   * \param  name  The name of the symbol. It must be a valid C++ variable name.
   * \param  type  The data type of the symbol. It must be a valid C++ data type.
   * \param  assignExpr  The assignment expression for the symbol. It must be a
   *     valid C++ expression. The 'this' keyword is used as a place holder for
   *     the object found in the data block.
   * \param  defaultExpr  The default value to use for the symbol. It must be a
   *     valid C++ expression.
   * \param className  The class name of the object found in the data block from
   *     which to fetch the symbol's value.
   * \param blockType  The data block type string of the data block. The value
   *    kAliHLTAnyDataTypeID can be used to specify the any type wild card value.
   * \param origin  The origin of the data block, such as the detector name.
   *    The value kAliHLTDataOriginAny can be used to specify the any origin
   *    wild card value.
   * \param spec  The data block specification to use.
   */
  static void AddSymbol(
      const char* name, const char* type, const char* assignExpr,
      const char* defaultExpr, const char* className,
      const char* blockType, const char* origin, UInt_t spec
    );
  
  /**
   * Adds a new trigger menu item to the current trigger menu.
   * \param  conditionExpr  The trigger condition expression. It must be a valid
   *     C++ expression, where the symbol names must be either defined in the
   *     menu or the names of AliHLTTrigger components.
   * \param  domainExpr  The trigger domain merging expression. It must be a
   *     valid C++ expression, where the symbol names must be either defined in
   *     the menu or the names of AliHLTTrigger components.
   * \param  prescalar  The prescalar value to use (Zero if not used).
   * \param  description  Optional description string which will be used in the
   *     global result.
   */
  static void AddItem(
      const char* conditionExpr, const char* domainExpr, UInt_t prescalar,
      const char* description = NULL
    );
  
  /**
   * Adds a new trigger menu item to the current trigger menu.
   * \param  conditionExpr  The trigger condition expression. It must be a valid
   *     C++ expression, where the symbol names must be either defined in the
   *     menu or the names of AliHLTTrigger components.
   * \param  domainExpr  The trigger domain merging expression. It must be a
   *     valid C++ expression, where the symbol names must be either defined in
   *     the menu or the names of AliHLTTrigger components.
   * \param  description  Optional description string which will be used in the
   *     global result.
   */
  static void AddItem(
      const char* conditionExpr, const char* domainExpr,
      const char* description = NULL
    );
  
 private:
  
  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerConfig(const AliHLTGlobalTriggerConfig& obj);
  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerConfig& operator = (const AliHLTGlobalTriggerConfig& obj);
  
  static AliHLTTriggerMenu* fgMenu;  /// Trigger menu which is created by this interface class.
  
  ClassDef(AliHLTGlobalTriggerConfig, 0) // Interface class used to construct a global HLT trigger menu configuration.
};

#endif // ALIHLTGLOBALTRIGGERCONFIG_H

