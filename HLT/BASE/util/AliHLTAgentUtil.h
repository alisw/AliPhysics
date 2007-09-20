// @(#) $Id$

#ifndef ALIHLTAGENTUTIL_H
#define ALIHLTAGENTUTIL_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTAgentUtil.h
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTUtil library
*/

#include "AliHLTModuleAgent.h"

/**
 * @class AliHLTAgentUtil
 * This is the agent for the AliHLTUtil library.
 *
 * @ingroup alihlt_system
 */
class AliHLTAgentUtil : public AliHLTModuleAgent {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTAgentUtil();
  /** destructor */
  virtual ~AliHLTAgentUtil();

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
  const char* GetLocalRecConfigurations(AliRunLoader* runloader=NULL) const;

  /**
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const char* GetRequiredComponentLibraries() const;

 protected:

 private:
  ClassDef(AliHLTAgentUtil, 0);
};

#endif
