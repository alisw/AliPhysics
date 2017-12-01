//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTOFAGENT_H
#define ALIHLTTOFAGENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTTOFAgent.h
    @author  Jochen Thaeder <jochen@thaeder.de>
    @brief   Agent of the libAliHLTTOF library
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTModuleAgent.h"

/**
 * @class AliHLTTOFAgent
 * This is the agent for the AliHLTTOF library.<br>
 *
 * The sample agent implements all interface function provided by @ref
 * AliHLTModuleAgent :
 * - CreateConfigurations() <br>
 *   The method gets an instance of the AliHLTConfigurationHanler to add
 *   configurations, e.g. 
 *   <pre>
 *   handler->CreateConfiguration("my-puplisher"  , "FilePublisher", NULL , "data.bin");
 *   ...
 *   handler->CreateConfiguration("my-analysis-chain"  , "FileWriter", "my-processor" , "my arguments");
 *   </pre>
 * - GetReconstructionChains() <br>
 *   returns a string of blank separated configurations to be run during
 *   local event reconstruction.
 *   <pre>
 *   return "my-data-sink my-analysis-chain";
 *   </pre>
 * - GetRequiredComponentLibraries() <br>
 *   returns a string of blank separated libraries which have to be loaded
 *   in addition in order to load all required components. <br>
 *   @note Not the right place for library dependencies.
 *   <pre>
 *   return "libAliHLTUtil.so";
 *   </pre>
 * - RegisterComponents() <br>
 *   registers the components: AliHLTDummyComponent, AliHLTSampleComponent1,
 *   AliHLTSampleComponent2, and AliHLTSampleMonitoringComponent<br>
 * - GetHandlerDescription() <br>
 *   Handles HLTOUT data blocks of type {DDL_RAW,SMPL}
 *   <pre>
 *   if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginSample)) {
 *     desc=AliHLTOUTHandlerDesc(kRawReader, dt, GetModuleId());
 *     return 1;
 *   }
 *   </pre>
 * - GetOutputHandler() <br>
 *   Returns handler AliHLTOUTHandlerEquId for HLTOUT data blocks of
 *   type {DDL_RAW,SMPL}
 *   <pre>
 *   if (dt==(kAliHLTDataTypeDDLRaw|kAliHLTDataOriginSample)) {
 *     return new AliHLTOUTHandlerEquId;
 *   }
 *   </pre>
 * - DeleteOutputHandler() <br>
 *   Deletes the output handler. In this case there is no special handling
 *   needed.
 *
 * In order to hook the sample library up to the HLT system on global object
 * @ref gAliHLTTOFAgent of the agent is defined in the source code.
 * 
 * @ingroup alihlt_system
 */
class AliHLTTOFAgent : public AliHLTModuleAgent {
 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTTOFAgent();

  /** destructor */
  virtual ~AliHLTTOFAgent();

  /*
   * ---------------------------------------------------------------------------------
   * Public functions to implement AliHLTModuleAgent's interface.
   * These functions are required for the registration process
   * ---------------------------------------------------------------------------------
   */

  UInt_t GetDetectorMask() const;

  /**
   * Register all configurations belonging to the TOF library with the
   * AliHLTConfigurationHandler. 
   * @param handler   [in] the configuration handler
   * @param rawReader [in] AliRoot RawReader instance 
   * @param runloader [in] AliRoot runloader
   * @return neg. error code if failed
   */
  Int_t CreateConfigurations(AliHLTConfigurationHandler* handler,
			     AliRawReader* rawReader=NULL,
			     AliRunLoader* runloader=NULL) const;

  /**
   * Get the top configurations for local event reconstruction.
   * A top configuration describes a processing chain. It can simply be
   * described by the last configuration(s) in the chain. 
   * @param rawReader [in] AliRoot RawReader instance 
   * @param runloader [in] AliRoot runloader
   * @return string containing the top configurations separated by blanks
   */
  const Char_t* GetReconstructionChains(AliRawReader* rawReader=NULL,
					AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const Char_t* GetRequiredComponentLibraries() const;

  /**
   * Register components for the AliHLTTOF library.
   * @param pHandler  [in] instance of the component handler          
   */
  Int_t RegisterComponents(AliHLTComponentHandler* pHandler) const;
  
  /** interface function, see @ref AliHLTModuleAgent for description */
  Int_t GetHandlerDescription(AliHLTComponentDataType dt, AliHLTUInt32_t spec,
			      AliHLTOUTHandlerDesc& desc) const;

  /** interface function, see @ref AliHLTModuleAgent for description */
  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt, AliHLTUInt32_t spec);
  
  /** interface function, see @ref AliHLTModuleAgent for description */
  Int_t DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  /** interface function, see @ref AliHLTModuleAgent for description */
  AliHLTModulePreprocessor* GetPreprocessor();

 protected:
 
  ///////////////////////////////////////////////////////////////////////////////////
  
private:

  /*
   * ---------------------------------------------------------------------------------
   * Private functions to implement AliHLTModuleAgent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */

  /** copy constructor prohibited */
  AliHLTTOFAgent(const AliHLTTOFAgent&);

  /** assignment operator prohibited */
  AliHLTTOFAgent& operator=(const AliHLTTOFAgent&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
  /** ROOT specific member definition */
  ClassDef(AliHLTTOFAgent, 0);
};

#endif
