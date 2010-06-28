//-*- Mode: C++ -*-
// $Id$
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
#include "TClonesArray.h"

class AliHLTTriggerMenu;
class AliHLTGlobalTrigger;

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
 * Output Data Types: kAliHLTDataTypeGlobalTrigger and kAliHLTDataTypeReadoutList <br>
 *
 * <h2>Mandatory arguments:</h2>
 * None.
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
 * \li -cint <br>
 *      Use CINT to interprete the generated global trigger instead of compiling it.
 *      This will also be the case if no compiler is available.
 * \li -usecode <i>filename</i> <i>classname</i> <br>
 *      Used to force the component to use an existing class for the global HLT trigger
 *      class implementation, with the name of <i>classname</i> and found in the file
 *      <i>filename</i>.
 * \li -skipctp <br>
 *      Indicates that the CTP data should not be added to the global HLT trigger decision.
 * \li -forward-input <br>
 *      Forward the input objects instead of adding them to the global HLT trigger decision.
 *      This will also add a short info on the input objects and decisions, like
 *      -include-input=short, to switch off -include-input=none can be placed after the
 *      parameter
 * \li -include-input[=none,short,objects,both] <br>
 *      Steer adding of input objects to the global HLT trigger decision.
 *      Options: none - include nothing, short - include a short TNames array,
 *               objects - include objects, by default on
 *               both - include both objects and short info
 * \li -process-all-events <br>
 *      Indicates that all events should be processed with the global trigger logic and
 *      not just the data events. The default is not to process just the data events.
 *
 * <h2>Configuration:</h2>
 * Configured from CDB but can be overridden with the -config argument.
 *
 * <h2>Default CDB entries:</h2>
 * HLT/ConfigHLT/HLTGlobalTrigger - Contains the global trigger menu.
 *
 * <h2>Performance:</h2>
 * This is a linear function of the number of input triggers (AliHLTTrigger) that
 * need to be processed.
 * For a modest trigger menu configurations the processing time per event should
 * be on the order of a few milliseconds.
 *
 * <h2>Memory consumption:</h2>
 * Memory consumption is minimal. It should be on the order of 2 or 3 MBytes.
 *
 * <h2>Output size:</h2>
 * This will depend almost linearly on the number of intput triggers and summary
 * data objects used. Thus, for every trigger (AliHLTTrigger object) specified
 * in the trigger menu the output size will require about 1 kBytes.
 * Then for every summary data object (i.e. TObject symbol defined in the trigger
 * menu configuration) one will need an extra few kBytes, depending on the size
 * of the summary objects.
 * In total one would expect no more than a MByte output size for a large trigger
 * configuration and typically only a few kBytes for a small or optimised one.
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
   * Inherited from AliHLTTrigger.
   * This returns kAliHLTDataTypeGlobalTrigger by default.
   * @param list <i>[out]</i>: The list of data types to be filled.
   */
  virtual void GetOutputDataTypes(AliHLTComponentDataTypeList& list) const;
  
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

  enum StatusBits {
    kForwardInput       = BIT(14),  // forward input objects instead of adding them to the decision object
    kIncludeInput       = BIT(15),  // include input objects in the decision object
    kIncludeShort       = BIT(16),  // include short description of input objects: name, title, decision
    kSkipCTP            = BIT(17),  // skip CTP data object in the decision object
  };

  void   SetBit(AliHLTUInt32_t f, bool set) {
    if (set) SetBit(f);
    else ResetBit(f);
  }
  void   SetBit(AliHLTUInt32_t f) { fBits |= f; }
  void   ResetBit(AliHLTUInt32_t f) { fBits &= ~f; }
  bool   TestBit(AliHLTUInt32_t f) const { return (bool) ((fBits & f) != 0); }

 protected:

  /**
   * Applies the global HLT trigger.
   * @return Zero is returned on success and a negative error code on failure.
   */
  virtual int DoTrigger();
  
  /**
   * Reconfigures the component by loading the trigger menu from the given
   * CDB entry.
   * \param cdbEntry  The CDB path to the trigger menu to load.
   * \param chainId   The ID of the component in the chain.
   * \returns  Zero on success and non-zero values otherwise.
   */
  virtual int Reconfigure(const char* cdbEntry, const char* chainId);
  
 private:

  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerComponent(const AliHLTGlobalTriggerComponent& obj);
  /// Not implemented. Do not allow copying of this object.
  AliHLTGlobalTriggerComponent& operator = (const AliHLTGlobalTriggerComponent& obj);
  
  /**
   * Loads a trigger menu object from the CDB.
   * \param cdbPath <i>in</i> The path in the CDB to load the trigger menu object from.
   * \param menu  <i>out</i> A pointer that gets filled with the new trigger menu object.
   * \returns  Zero if the trigger menu object was found and the pointer to it
   *   set in the <i>menu</i> variable. If a non-zero error code is returned then
   *   the <i>menu</i> variable is not changed at all.
   */
  int LoadTriggerMenu(const char* cdbPath, const AliHLTTriggerMenu*& menu);
  
  /**
   * Generates a file name for the generated on the fly code using a UUID.
   * \param name <i>out</i> The name of the class to use.
   * \param filename <i>out</i> The name of the file containing the code.
   */
  void GenerateFileName(TString& name, TString& filename) const;
  
  /**
   * Generates the code for the global trigger to apply the given trigger menu.
   * The code will then be compiled on the fly and loaded. The name of the new
   * class is returned so that a new instance of the class can be created via:
   * \code
   *  AliHLTGlobalTrigger::CreateNew(name)
   * \endcode
   * where name is the name of the generated class as returned by this method.
   *
   * The name of the generated code file is stored in the variable fCodeFileName
   * and the fDeleteCodeFile is set to true.
   *
   * \param menu <i>in</i> The trigger menu to create the global trigger class from.
   * \param name <i>out</i> The name of the generated class.
   * \param filename <i>out</i> The name of the generated file containing the code.
   * \param includePaths <i>in</i> The list of include path strings.
   * \param includeFiles <i>in</i> The list of include file strings.
   * \returns  The error code suitable to return in DoInit. Zero on success.
   */
  int GenerateTrigger(
      const AliHLTTriggerMenu* menu, TString& name, TString& filename,
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
   * Unloads the code that was previously loaded by LoadTriggerClass.
   * \param filename  The name of the file containing the global trigger class logic to be unloaded.
   * \returns  The error code suitable to return in DoInit. Zero on success.
   */
  int UnloadTriggerClass(const char* filename);
  
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
  
  /**
   * Extracts the trailing operator in a C++ expression and returns the found
   * operator in a separate output string.
   * [in/out] \param  expr  The C++ expression to check. The trailing operator
   *      is removed from the expression if found.
   * [out] \param  op   The output variable which will be filled with the
   *      operator found in the expression.
   * \return  true if the trailing operator was found in the expression and
   *      false otherwise.
   */
  bool ExtractedOperator(TString& expr, TString& op);

  /**
   * Add trigger decisions according to the active CTP trigger classes
   * An internal TclonesArray holds the trigger decisions to be added. The trigger
   * decisions are updated according to the active CTP trigger mask.
   * \param pTrigger  The instance of the global trigger
   * \param pCTPData  Instance of the CTP data
   * \param trigData  Current trigger data, if NULL, the active trigger data from the CTP data is used
   */
  int AddCTPDecisions(AliHLTGlobalTrigger* pTrigger, const AliHLTCTPData* pCTPData, const AliHLTComponentTriggerData* trigData);

  /**
   * Print some statistics based on the trigger counters
   */
  int PrintStatistics(const AliHLTGlobalTrigger* pTrigger, AliHLTComponentLogSeverity level=kHLTLogInfo, int offset=1) const;
  
  AliHLTGlobalTrigger* fTrigger;  //! Trigger object which implements the global trigger menu.
  bool fDebugMode;  //! Indicates if the generated global trigger class should be in debug mode.
  bool fRuntimeCompile;  //! Indicates if the generated global trigger class should be compiled
  bool fDeleteCodeFile; //! If true then the code file indicated by fCodeFileName should be deleted during DoDeinit.
  TString fCodeFileName; //! base file name of the generated code for the global trigger
  TString fClassName;  //! The generated/loaded trigger class name.
  TClonesArray* fCTPDecisions; //! AliHLTTriggerDecision objects for the CTP classes
  unsigned long fBufferSizeConst; //! Constant size estimate for GetOutputDataSize.
  double fBufferSizeMultiplier; //! Buffer size multiplier estimate for GetOutputDataSize.
  TClonesArray fIncludePaths; //! Paths specified by the -includepath command line option.
  TClonesArray fIncludeFiles; //! Files specified by the -include command line option.
  TString fLibStateAtLoad; //! This stores the loaded libraries just before we tell CINT to load the interpreted file.
  AliHLTUInt32_t fBits; //! Status bits
  bool fDataEventsOnly; //! Flag indicating if only data events are processed with trigger logic.

  static const char* fgkTriggerMenuCDBPath; //! The path string to read the trigger menu from the CDB.
  
  ClassDef(AliHLTGlobalTriggerComponent, 0) // Global HLT trigger component class which produces the final trigger decision and readout list.
};

#endif // ALIHLTGLOBALTRIGGERCOMPONENT_H

