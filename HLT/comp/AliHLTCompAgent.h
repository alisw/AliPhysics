// @(#) $Id$

#ifndef ALIHLTCOMPAGENT_H
#define ALIHLTCOMPAGENT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTCompAgent.h
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTComp library
*/

// see class description below
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTModuleAgent.h"

/**
 * @class AliHLTCompAgent
 * This is the agent for the AliHLTComp library.
 *
 * @ingroup alihlt_system
 */
class AliHLTCompAgent : public AliHLTModuleAgent {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTCompAgent();
  /** destructor */
  virtual ~AliHLTCompAgent();

  /**
   * Register all configurations belonging to this module with the
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
   * Get the top configurations belonging to this module.
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
   * Register components for the AliHLTSim library.
   * @param pHandler  [in] instance of the component handler          
   */
  int RegisterComponents(AliHLTComponentHandler* pHandler) const;
 protected:

 private:
  ClassDef(AliHLTCompAgent, 0);
};

#endif
