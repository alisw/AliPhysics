//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTZDCAGENT_H
#define ALIHLTZDCAGENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTZDCAgent.h
    @author Chiara Oppedisano <Chiara.Oppedisano@to.infn.it>
    @brief   Agent of the libAliHLTZDC library
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTModuleAgent.h"

/**
 * @class AliHLTZDCAgent
 * This is the agent for the AliHLTZDC library.<br>
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
 * @ref gAliHLTZDCAgent of the agent is defined in the source code.
 * 
 * @ingroup alihlt_system
 */
class AliHLTZDCAgent : public AliHLTModuleAgent {
 public:

  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTZDCAgent();
  
/** destructor */
  virtual ~AliHLTZDCAgent();

  /**
   * Register all configurations belonging to the ZDC library with the
   * AliHLTConfigurationHandler. 
   * @param handler   [in] the configuration handler
   * @param rawReader [in] AliRoot RawReader instance 
   * @param runloader [in] AliRoot runloader
   * @return neg. error code if failed
   */
  int CreateConfigurations(AliHLTConfigurationHandler* handler,
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
  const char* GetReconstructionChains(AliRawReader* rawReader=NULL,
				      AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const char* GetRequiredComponentLibraries() const;

  /**
   * Register components for the AliHLTZDC library.
   * @param pHandler  [in] instance of the component handler          
   */
  int RegisterComponents(AliHLTComponentHandler* pHandler) const;
  
  /** interface function, see @ref AliHLTModuleAgent for description */
  int GetHandlerDescription(AliHLTComponentDataType dt,
			    AliHLTUInt32_t spec,
			    AliHLTOUTHandlerDesc& desc) const;

  /** interface function, see @ref AliHLTModuleAgent for description */
  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt, AliHLTUInt32_t spec);
  
  /** interface function, see @ref AliHLTModuleAgent for description */				     
  int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  /** interface function, see @ref AliHLTModuleAgent for description */
  AliHLTModulePreprocessor* GetPreprocessor();

 private:
  /** copy constructor prohibited */
  AliHLTZDCAgent(const AliHLTZDCAgent&);

  /** assignment operator prohibited */
  AliHLTZDCAgent& operator=(const AliHLTZDCAgent&);

  /** ROOT specific member definition */
  ClassDef(AliHLTZDCAgent, 0);
};

#endif
