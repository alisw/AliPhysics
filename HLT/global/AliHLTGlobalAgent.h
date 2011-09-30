//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTGLOBALAGENT_H
#define ALIHLTGLOBALAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

//  @file   AliHLTGlobalAgent.h
//  @author Matthias Richter
//  @date   
//  @brief  Agent of the libAliHLTGlobal library
// 

#include "AliHLTModuleAgent.h"

class AliHLTOUTHandler;

/**
 * @class AliHLTGlobalAgent
 * This is the agent for the AliHLTGlobal library.
 *
 * @ingroup alihlt_system
 */
class AliHLTGlobalAgent : public AliHLTModuleAgent {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTGlobalAgent();
  /** destructor */
  virtual ~AliHLTGlobalAgent();

  int RegisterComponents(AliHLTComponentHandler* pHandler) const;

  /**
   * Inherited from AliHLTModuleAgent
   * Register all configurations belonging to this module with the
   * AliHLTConfigurationHandler. The agent can adapt the configurations
   * to be registered to the current AliRoot setup by checking the
   * runloader.
   * @param handler      the configuration handler
   * @param rawReader    AliRawReader instance
   * @param runloader    AliRoot runloader
   * @return neg. error code if failed
   */
  int CreateConfigurations(AliHLTConfigurationHandler* handler,
			   AliRawReader* rawReader=NULL,
			   AliRunLoader* runloader=NULL) const;

  /**
   * Inherited from AliHLTModuleAgent
   * Get the top configurations belonging to this module.
   * A top configuration describes a processing chain. It can simply be
   * described by the last configuration(s) in the chain. 
   * The agent can adapt the configurations to be registered to the current
   * AliRoot setup by checking the runloader.
   * @param rawReader    AliRawReader instance
   * @param runloader    AliRoot runloader
   * @return number of configurations, neg. error code if failed
   */
  const char* GetReconstructionChains(AliRawReader* rawReader=NULL,
				      AliRunLoader* runloader=NULL) const;

  /**
   * Inherited from AliHLTModuleAgent
   * Component libraries which the configurations of this agent depend on.
   * @return list of component libraries as a blank-separated string.
   */
  const char* GetRequiredComponentLibraries() const;

  int GetHandlerDescription(AliHLTComponentDataType dt,
			    AliHLTUInt32_t spec,
			    AliHLTOUTHandlerDesc& desc) const;

  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt, AliHLTUInt32_t spec);

  int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

  AliHLTModulePreprocessor* GetPreprocessor();

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTGlobalAgent(const AliHLTGlobalAgent&);
  /** assignment operator prohibited */
  AliHLTGlobalAgent& operator=(const AliHLTGlobalAgent&);

  ClassDef(AliHLTGlobalAgent, 0);
};

#endif
