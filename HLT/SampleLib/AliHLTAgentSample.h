// @(#) $Id$

#ifndef ALIHLTAGENTSAMPLE_H
#define ALIHLTAGENTSAMPLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTAgentSample.h
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTSample library
*/

#include "AliHLTModuleAgent.h"

/**
 * @class AliHLTAgentSample
 * This is the agent for the AliHLTSample library.<br>
 * The AliHLTSample library illustrates usage of the HLT framework. The
 * agent provides information on the features of the sample components
 * and the configuration which should be run during AliRoot reconstruction.
 *
 * The sample agent implements all interface function provided by @ref
 * AliHLTModuleAgent :
 * - @ref CreateConfigurations <br>
 *   The method gets an instance of the AliHLTConfigurationHanler to add
 *   configurations, e.g. 
 *   <pre>
 *   handler->CreateConfiguration("my-puplisher"  , "FilePublisher", NULL , "data.bin");
 *   ...
 *   handler->CreateConfiguration("my-analysis-chain"  , "FileWriter", "my-processor" , "my arguments");
 *   </pre>
 * - @ref GetReconstructionChains <br>
 *   returns a string of blank separated configurations to be run during
 *   local event reconstruction.
 *   <pre>
 *   return "my-data-sink my-analysis-chain";
 *   </pre>
 * - @ref GetRequiredComponentLibraries <br>
 *   returns a string of blank separated libraries which have to be loaded
 *   in addition in order to load all required components. <br>
 *   @note Not the right place for library dependencies.
 *   <pre>
 *   return "libAliHLTUtil.so";
 *   </pre>
 * - not implemented are the in iterface methods
 *   - @ref AliHLTModuleAgent::RegisterComponents
 *
 * In order to hook the sample library up to the HLT system on global object
 * @ref gAliHLTAgentSample of the agent is defined in the source code.
 * 
 * @ingroup alihlt_system
 */
class AliHLTAgentSample : public AliHLTModuleAgent {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTAgentSample();
  /** destructor */
  virtual ~AliHLTAgentSample();

  /**
   * Register all configurations belonging to the sample library with the
   * AliHLTConfigurationHandler. The agent can adapt the configurations
   * to be registered to the current AliRoot setup by checking the
   * runloader.
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
   * The agent can adapt the configurations to be registered to the current
   * AliRoot setup by checking the runloader.
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
   * Register components for the AliHLTSample library.
   * @param pHandler  [in] instance of the component handler          
   */
  int RegisterComponents(AliHLTComponentHandler* pHandler) const;
 protected:

 private:
  /** file name of the generated test data*/
  static const char* fgkAliHLTAgentSampleData;                      //!transient

  /** file name of the output file */
  static const char* fgkAliHLTAgentSampleOut;                       //!transient

  /** ROOT specific member definition */
  ClassDef(AliHLTAgentSample, 1);
};

#endif
