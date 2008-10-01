
#ifndef ALIHLTAGENTFMD_H
#define ALIHLTAGENTFMD_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTAgentSample.h
    @author Hans Hjersin Dalsgaard
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
