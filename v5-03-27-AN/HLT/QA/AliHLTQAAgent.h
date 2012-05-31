//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTQAAGENT_H
#define ALIHLTQAAGENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTQAAgent.h
    @author Matthias Richter
    @date   2010-03-10
    @brief  Module Agent for HLTqadm library
*/

#include "AliHLTModuleAgent.h"

class AliHLTQAAgent: public AliHLTModuleAgent {
 public:
  AliHLTQAAgent();
  virtual ~AliHLTQAAgent();

  /// inherited from AliHLTModuleAgent
  const char* GetQAPlugins() const;

 protected:
 private:
  /** copy constructor prohibited */
  AliHLTQAAgent(const AliHLTQAAgent&);   
  /** assignment operator prohibited */
  AliHLTQAAgent& operator = (const AliHLTQAAgent&);

  ClassDef(AliHLTQAAgent,0)  // Module Agent for HLTqadm library
};

#endif // ALIHLTQAAGENT_H
