//-*- Mode: C++ -*-
// @(#) $Id$
#ifndef ALIHLTSAMPLECOMPONENT1_H
#define ALIHLTSAMPLECOMPONENT1_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTSampleComponent1.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  A sample processing component for the HLT.
*/

#include "AliHLTProcessor.h"

/**
 * @class AliHLTSampleComponent1
 * An HLT sample component.
 * This component does not any data processing at all. It just
 * illustrates the existence of several components in ine library and
 * allows to set up a very simple chain with different components.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b Sample-component1 <br>
 * Library: \b libAliHLTSample.so     <br>
 * Input Data Types: @ref kAliHLTAnyDataType <br>
 * Output Data Types: none <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -mandatory1     <i> teststring   </i> <br>
 *      an argument with one parameter
 * \li -mandatory2                           <br>
 *      an argument without parameters
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -optional1      <i> teststring   </i> <br>
 *      an argument with one parameter
 * \li -optional2                            <br>
 *      an argument without parameters
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -config1      <i> teststring   </i> <br>
 *      a configuration argument with one parameter
 * \li -config2                            <br>
 *      a configuration argument without parameters
 *
 * <h2>Default CDB entries:</h2>
 * The component has just one default CDB entry in 
 * <tt>HLT/ConfigSample/SampleComponent1</tt>.
 * It does not load any configuration from the global <tt>ConfigHLT</tt>
 * folder.
 * \li -TObjString object holding a string with the configuration parameters
 *      explained above
 *
 * <h2>Performance:</h2>
 * The component does not any event data processing.
 *
 * <h2>Memory consumption:</h2>
 * The component does not any event data processing.
 *
 * <h2>Output size:</h2>
 * The component has no output data.
 *
 * Furthermore it illustrates the component argument scanning and the
 * component configuration. There are actually two methods to init/
 * configure a component:
 * - via command line arguments. The arguments are specified in the HLT
 *   chain configuration and are passed to the component during
 *   initialization in @ref DoInit()
 * - from a CDB entry. The CDB can contain configuration objects for a
 *   component and the component can implement the handling
 *
 * The component implements the @ref alihltcomponent-low-level-interface.
 * for data processing.
 *
 * Using the latter case, a component can also be reconfigured. Special
 * events are propageted through the chain in order to trigger the re-
 * configuration. The component needs to implement the function
 * @ref Reconfigure(). The simplest version of a configuration object is
 * a string object (TObjString) containing configuration arguments.
 *
 * @ingroup alihlt_tutorial
 */
class AliHLTSampleComponent1 : public AliHLTProcessor {
public:
  AliHLTSampleComponent1();
  virtual ~AliHLTSampleComponent1();

  // AliHLTComponent interface functions
  const char* GetComponentID() { return "Sample-component1";}
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
    list.push_back(kAliHLTAnyDataType);
  }
  AliHLTComponentDataType GetOutputDataType() {return kAliHLTVoidDataType;}
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {constBase = 0;inputMultiplier = 0;};

  // Spawn function, return new class instance
  AliHLTComponent* Spawn() {return new AliHLTSampleComponent1;};

 protected:
  // AliHLTComponent interface functions
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		       AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );
  int Reconfigure(const char* cdbEntry, const char* chainId);
  int ReadPreprocessorValues(const char* modules);

  using AliHLTProcessor::DoEvent;

private:
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   *
   * This function illustrates the scanning of an argument string. The string
   * was presumably fetched from the CDB.
   */
  int Configure(const char* arguments);

  ClassDef(AliHLTSampleComponent1, 1)
};
#endif
