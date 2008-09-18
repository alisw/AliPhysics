// @(#) $Id: AliHLTAgentSample.h 25820 2008-05-16 11:47:09Z richterm $

#ifndef ALIHLTAGENTFMD_H
#define ALIHLTAGENTFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTAgentSample.h
    @author Matthias Richter
    @date   
    @brief  Agent of the libAliHLTSample library
*/

#include "AliHLTModuleAgent.h"


class AliHLTAgentFMD : public AliHLTModuleAgent {
 public:

  AliHLTAgentFMD();
  /** destructor */
  virtual ~AliHLTAgentFMD();

  const char* GetRequiredComponentLibraries() const;

  int RegisterComponents(AliHLTComponentHandler* pHandler) const;

 protected:

 private:

  ClassDef(AliHLTAgentFMD, 1);
};

#endif
