#ifndef ALIHLTGLOBALTRIGGERCOMPONENT_H
#define ALIHLTGLOBALTRIGGERCOMPONENT_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/// @file   AliHLTGlobalTriggerComponent.h
/// @author Artur Szostak <artursz@iafrica.com>
/// @date   26 Nov 2008
/// @brief  Declaration of the AliHLTGlobalTriggerComponent component class.

#include "AliHLTTrigger.h"

class AliHLTTriggerMenu;
class AliHLTGlobalTrigger;
class TClonesArray;

/**
 * \class AliHLTGlobalTriggerComponent
 * This class applies the global HLT trigger to all trigger information produced
 * by components deriving from AliHLTTrigger.
 * Any information delivered by other components in data blocks that contain
 * TObjects can also be used for the trigger algorithm. In such cases a symbol
 * needs to be defined in the global trigger menu which can then be used inside
 * the trigger condition expressions.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b HLTGlobalTrigger <br>
 * Library: \b libAliHLTTrigger.so   <br>
 * Input Data Types: ::kAliHLTAnyDataType <br>
 * Output Data Types: kAliHLTDataTypeTObject|kAliHLTDataOriginOut <br>
 *
 * <h2>Mandatory arguments:</h2>
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -config <i>filename</i> <br>
 *      Indicates the configuration macro file to use for the global HLT trigger menu.
 * \li -includepath <i>path</i> <br>
 *      Indicates the include path to use if the automatically generated code that
 *      implements the global HLT trigger requires non-standard includes.
 * \li -include <i>filename</i> <br>
 *      Indicates a file name that should be included in the automatically generated
 *      trigger implementation code.
 * \li -debug <br>
 *      If specified the automatically generated class will contain extra debugging
 *      code and the ACLiC system will have debugging compilation turned on.
 * \li -usecode <i>filename</i> <i>classname</i> <br>
 *      Used to force the component to use an existing class for the global HLT trigger
 *      class implementation, with the name of <i>classname</i> and found in the file
 *      <i>filename</i>.
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * Configuration by component arguments.
 *
 * <h2>Default CDB entries:</h2>
 * TODO
 *
 * <h2>Performance:</h2>
 * TODO
 *
 * <h2>Memory consumption:</h2>
 * TODO
 *
 * <h2>Output size:</h2>
 * TODO
 *
 * \ingroup alihlt_trigger_components
 */
class AliHLTGlobalTriggerComponent : public AliHLTTrigger
{
 public:
 
  AliHLTGlobalTriggerComponent();
  virtual ~AliHLTGlobalTriggerComponent();
  
  /**
   * Inherited from AliHLTTrigger.
   * @return string containing the global trigger name.
   */
  virtual const char* GetTriggerName() const { return "HLTGlobalTrigger"; };
  
  /**
   * Get a ratio by how much the data volume is shrunk or enhanced.
   * The method returns a size proportional to the trigger name string length
   * for constBase, and 1 for inputMultiplier.
   * @param constBase        <i>[out]</i>: additive part, independent of the
   *                                   input data volume  
   * @param inputMultiplier  <i>[out]</i>: multiplication ratio
   */
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);
  
  /**
   * Initialise the component.
   * \param argc  The number of arguments in argv.
   * \param argv  Array of component argument strings.
   * \returns  Zero on success and negative number on failure.
   */
  virtual Int_t DoInit(int argc, const char** argv);
  
  /**
   * Cleanup the component.
   * \returns  Zero on success and negative number on failure.
   */
  virtual Int_t DoDeinit();
  
  /**
   * Spawn function creates a new object.
   * @return new class instance.
   */
  virtual AliHLTComponent* Spawn();

 protected:

  /**
   * Applies the global HLT trigger.
   * @return Zero is returned on success and a negative error code on failure.
   */
  virtual int DoTrigger();
  
 private:

  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerComponent(const AliHLTGlobalTriggerComponent& obj);
  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerComponent& operator = (const AliHLTGlobalTriggerComponent& obj);
  
  /**
   * Generates the code for the global trigger to apply the given trigger menu.
   * The code will then be compiled on the fly and loaded. The name of the new
   * class is returned so that a new instance of the class can be created via:
   * \code
   *  AliHLTGlobalTrigger::CreateNew(name)
   * \endcode
   * where name is the name of the generated class as returned by this method.
   * \param menu <i>in</i> The trigger menu to create the global trigger class from.
   * \param name <i>out</i> The name of the generated class.
   * \param includePaths <i>in</i> The list of include path strings.
   * \param includeFiles <i>in</i> The list of include file strings.
   * \returns  The error code suitable to return in DoInit. Zero on success.
   */
  int GenerateTrigger(
      const AliHLTTriggerMenu* menu, TString& name,
      const TClonesArray& includePaths, const TClonesArray& includeFiles
    );
  
  /**
   * Loads the code for the generated HLT global trigger class. The code is compiled
   * on the fly if possible, otherwise the CINT interpreter is used to interpret
   * the class.
   * \param filename  The name of the file containing the code for the global trigger class.
   * \param includePaths <i>in</i> The list of include path strings.
   * \returns  The error code suitable to return in DoInit. Zero on success.
   */
  int LoadTriggerClass(const char* filename, const TClonesArray& includePaths);
  
  /**
   * Searches for the specified symbol name in the given list.
   * \param name  The name of the symbol to find.
   * \param list  The list to search for the symbol in.
   * \returns  The position (index) of the symbol found or -1 if it was not found.
   */
  int FindSymbol(const char* name, const TClonesArray& list);
  
  /**
   * Builds the list of symbols to use in the custom global trigger menu
   * implementation class.
   * \param  menu  The trigger menu to create the global trigger class from.
   * \param  list  The list to fill with symbols.
   * \returns  The error code suitable to return in DoInit. Zero on success.
   */
  int BuildSymbolList(const AliHLTTriggerMenu* menu, TClonesArray& list);
  
  AliHLTGlobalTrigger* fTrigger;  //! Trigger object which implements the global trigger menu.
  bool fDebugMode;  //! Indicates if the generated global trigger class should be in debug mode.
  
  ClassDef(AliHLTGlobalTriggerComponent, 0) // Global HLT trigger component class which produces the final trigger decision and readout list.
};

#endif // ALIHLTGLOBALTRIGGERCOMPONENT_H

