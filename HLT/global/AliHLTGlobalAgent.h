// $Id$

#ifndef ALIHLTGLOBALAGENT_H
#define ALIHLTGLOBALAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               */

/** @file   AliHLTGlobalAgent.h
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTGlobal library
*/

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

  int GetHandlerDescription(AliHLTComponentDataType dt,
			    AliHLTUInt32_t spec,
			    AliHLTOUTHandlerDesc& desc) const;

  AliHLTOUTHandler* GetOutputHandler(AliHLTComponentDataType dt, AliHLTUInt32_t spec);

  int DeleteOutputHandler(AliHLTOUTHandler* pInstance);

 protected:

 private:
  /** copy constructor prohibited */
  AliHLTGlobalAgent(const AliHLTGlobalAgent&);
  /** assignment operator prohibited */
  AliHLTGlobalAgent& operator=(const AliHLTGlobalAgent&);

  ClassDef(AliHLTGlobalAgent, 0);
};

#endif
