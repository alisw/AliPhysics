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
 * This is the agent for the AliHLTSample library.
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
   * Register all configurations belonging to this module with the
   * AliHLTConfigurationHandler. The agent can adapt the configurations
   * to be registered to the current AliRoot setup by checking the
   * runloader.
   * @param handler      the configuration handler
   * @param runloader    AliRoot runloader
   * @return neg. error code if failed
   */
  int CreateConfigurations(AliHLTConfigurationHandler* handler,
			   AliRunLoader* runloader=NULL) const;

  /**
   * Get the top configurations belonging to this module.
   * A top configuration describes a processing chain. It can simply be
   * described by the last configuration(s) in the chain. 
   * The agent can adapt the configurations to be registered to the current
   * AliRoot setup by checking the runloader.
   * @param runloader    AliRoot runloader
   * @return number of configurations, neg. error code if failed
   */
  const char* GetTopConfigurations(AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const char* GetRequiredComponentLibraries() const;

 protected:

 private:
  ClassDef(AliHLTAgentSample, 0);
};

#endif
