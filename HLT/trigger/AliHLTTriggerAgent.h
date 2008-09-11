// $Id$

#ifndef ALIHLTTRIGGERAGENT_H
#define ALIHLTTRIGGERAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTTriggerAgent.h
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTTrigger library
*/

#include "AliHLTModuleAgent.h"

class AliHLTOUTHandler;

/**
 * @class AliHLTTriggerAgent
 * This is the agent for the AliHLTTrigger library.
 *
 * @ingroup alihlt_system
 */
class AliHLTTriggerAgent : public AliHLTModuleAgent {
 public:
  /**
   * standard constructor. The agent is automatically registered in the
   * global agent manager
   */
  AliHLTTriggerAgent();
  /** destructor */
  virtual ~AliHLTTriggerAgent();

  int RegisterComponents(AliHLTComponentHandler* pHandler) const;

  int GetHandlerDescription(AliHLTComponentDataType dt,
			    AliHLTUInt32_t spec,
			    AliHLTOUTHandlerDesc& desc) const;

  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt, AliHLTUInt32_t spec);

  int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTTriggerAgent(const AliHLTTriggerAgent&);
  /** assignment operator prohibited */
  AliHLTTriggerAgent& operator=(const AliHLTTriggerAgent&);

  ClassDef(AliHLTTriggerAgent, 0);
};

#endif
